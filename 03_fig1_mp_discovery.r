library(Seurat)
library(NMF)
library(parallel)
library(fastICA)
library(ComplexHeatmap)
library(RColorBrewer)
library(scales)
library(dplyr)
library(circlize)
library(igraph)
library(CSCORE)

### Due to data control and large size, no example data for this script
## Edit from Barkley et.al
## Select patients with >=200 malignant cells
path=paste0("path_to_your_sub_cluster.RData")
for(p in path){
load(p)
tumor=type_list[["Epithelial cell"]]
tumor=SplitObject(tumor,split.by="orig.ident")
for(sample in names(tumor)){
obj=tumor[[sample]]
freq=data.frame(table(obj@meta.data$orig.ident,obj@meta.data$cell_subtype))
write.table(freq,"cell_freq.txt",append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
saveRDS(obj,paste0(sample,".rds"))
}
}

mat=read.table("cell_freq.txt")
mat_malig=mat[which(mat$V2=="malignant"),]
mat_malig=mat_malig[which(mat_malig$V3>=100),]
write.table(mat_malig,"target_samples.txt",sep="\t",quote=F,row.names=F,col.names=F)

mat_at=mat[which(mat$V2=="AT2"),]
select_sam=unique(c(mat_at[which(mat_at$V3>=20),]$V1,mat_malig$V1))
mat_at=mat_at[which(mat_at$V1 %in% select_sam),]
write.table(mat_at,"target_samples_AT2_all.txt",sep="\t",quote=F,row.names=F,col.names=F)

## NMF for cells from each patient
NMFToModules = function(
  res,
  gmin = 5
){
  
  scores = basis(res)
  coefs = coefficients(res)
  
  # Remove if fewer than gmin genes
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  l = sapply(modules, length)
  keep = (l >= gmin)
  scores = scores[, keep]
  coefs = coefs[keep, ]
  
  # Find modules
  ranks_x = t(apply(-t(t(scores) / apply(scores, 2, mean)), 1, rank))
  ranks_y = apply(-t(t(scores) / apply(scores, 2, mean)), 2, rank)
  for (i in 1:ncol(scores)){
    ranks_y[ranks_x[,i] > 1,i] = Inf
  }
  modules = apply(ranks_y, 2, function(m){
    a = sort(m[is.finite(m)])
    a = a[a == 1:length(a)]
    names(a)
  })
  
  names(modules) = sapply(modules, '[', 1)
  names(modules) = paste('m', names(modules), sep = '_')
  names(modules) = gsub('-','_',names(modules))
  
  return(modules)
}

# Better parallel this step, just for example here
for(sam in mat_malig$V1){
srt = readRDS(as.character(args[1]))
srt = subset(srt,subset=cell_subtype=="malignant")
srt = SCTransform(srt, return.only.var.genes = FALSE)
srt = RunPCA(srt, features = VariableFeatures(srt))
srt = RunUMAP(srt, dims = 1:10) %>% FindNeighbors(reduction = "pca", dims = 1:10, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)
# Run NMF
data = as.matrix(GetAssayData(srt, assay = 'SCT', slot = 'scale.data'))
data = data[VariableFeatures(srt),]
data[data < 0] = 0
data = data[apply(data, 1, var) > 0, ]
print(dim(data))
range = 2:25
res.list = mclapply(range, function(r){
  nmf(data, rank = r, nrun = 1, seed = 'ica', method = 'nsNMF')
}, mc.cores = 30)
names(res.list) = range
# Select rank
gmin=5
modules.list = lapply(res.list, NMFToModules, gmin = gmin)
print(sapply(modules.list,length))
comp = as.numeric(names(modules.list)) - sapply(modules.list, length)
mi = min(comp)
r = names(which(comp == mi))
r = r[length(r)]
print(r)
res = res.list[[r]]
# Process output
modules = NMFToModules(res, gmin = gmin)
scores = basis(res)
colnames(scores) = names(modules)
coefs = coefficients(res)
rownames(coefs) = names(modules)
# Order modules
h = Heatmap(coefs, clustering_distance_columns = 'euclidean')
o = row_order(h)
scores = scores[, o]
coefs = coefs[o, ]
modules = modules[o]
print(modules)
srt.malignant = srt
save(srt.malignant, file = paste0("./",sam,'srt.malignant.RData'))
res.list.malignant = res.list
save(res.list.malignant, file = paste0("./",sam,'res.list.malignant.RData'))
res.malignant = res
save(res.malignant, file = paste0("./",sam,'res.malignant.RData'))
save(modules, file = paste0("./",sam,'modules.RData'))
}

## Generate MP 
files = read.table("/target_samples.txt",header=F)
files=files$V1

res.list = lapply(files, function(f){
  loadRData(paste0(f,'/res.malignant.RData'))
})
names(res.list)=files
modules.list = sapply(names(res.list),USE.NAMES=TRUE, function(x){
  temp=NMFToModules(res.list[[x]])
  names(temp)=paste0(x,".",names(temp))
  return(temp)
})
genes.all = sort(unique(unlist(modules.list)))
all = unlist(modules.list, recursive = FALSE, use.names = FALSE)
names(all) = unlist(sapply(modules.list, names))
ta = table(unlist(all))
genes.use = names(ta)[ta > 1]
# Filter non-overlapping modules
for (i in 1:5){
  all = unlist(modules.list, recursive = FALSE, use.names = TRUE)
  all = lapply(all, intersect, genes.all)
  sim = sapply(all, function(x){
    sapply(all, function(y){
      length(intersect(x,y))/length(union(x,y))
    })
  })
  keep = rownames(sim)[apply(sim, 1, function(x){
    sum(x > 0.05) >= 5
  })]
  all = all[keep]
  modules.list = lapply(names(modules.list), function(x){
    li = modules.list[[x]]
    li[names(li)[paste(x,names(li),sep='.') %in% keep]]
  })
  names(modules.list) = names(res.list)
  ta = table(unlist(all))
  genes.use = names(ta)[ta > 1] 
  print(length(all))
}
save(all,sim,file="all.RData")
# Adjacency matrix, list by cancer
adj = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
adj.list = list()
cancer=factor(rep("LUAD",times=length(modules.list)))
for (can in levels(cancer)){
  sub = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
  rownames(sub) = genes.use
  colnames(sub) = genes.use
  for (s in names(modules.list)[cancer == can]){
    for (mod in modules.list[[s]]){
      mod = intersect(mod, genes.use)
      for (x in mod){
        for (y in mod){
          sub[x,y] = sub[x,y] + 1
        }
      }
    }
  }
  diag(sub) = 0
  adj.list[[can]] = sub
  adj = adj + sub
}
adj_keep = adj
adj = adj_keep
save(adj,file="adj.RData")
# Remove low connections
v_min=5
s_min=4
adj[] = (adj >= v_min)
#adj[adj <= 1] = 0
for (i in 1:5){
  keep = names(which(rowSums(adj) >= s_min))
  adj = adj[keep,keep]
  print(dim(adj))
}
# Cluster
save(modules, file = 'modules.RData')
fi = data.frame(sapply(modules,'[', 1:max(sapply(modules, length))), stringsAsFactors = FALSE)
fi[is.na(fi)] = ''
write.table(fi, file = 'modules.csv', sep = ',', row.names = FALSE, col.names = TRUE, quote = FALSE)

## Further filter loosely/false correlated genes
load("path_to_your_combined_seurat_object") # This file contain malignant cells & AT2 cell from all patients
DefaultAssay(combined_object) <- 'RNA'
combined_object$cell_id=rownames(combined_object@meta.data)
malignant_id=combined_object@meta.data[which(combined_object@meta.data$cell_subtype=="malignant"),]
malignant_id=split(malignant_id,f=factor(malignant_id$orig.ident))
# Down sampling cells to avoid overrepresentation
malignant_id_sampled=lapply(malignant_id,function(x){
  if(nrow(x)>150){
    temp_id=rownames(x)
    temp_id=sample(temp_id,150,replace=F)
    return(temp_id)
  }else{
    return(rownames(x))
  }
})
malignant_id_sampled=unlist(malignant_id_sampled)
malignant_object=subset(combined_object,subset=cell_id %in% malignant_id_sampled)
genes_selected=unlist(modules)
malig_result <- CSCORE(malignant_object, genes = genes_selected)

modules_refine=list()
for(m in names(modules)){
genes=modules[[m]]
cor_mat_sub=cor_mat[genes,genes]
diag(cor_mat_sub)=NA
cor_median=sort(apply(cor_mat_sub,1,median,na.rm=TRUE))
while( cor_median[1]<=0.1 ){
  g=names(cor_median)[1]
  cor_mat_sub=cor_mat_sub[which(rownames(cor_mat_sub) != g),which(rownames(cor_mat_sub) != g)]
  cor_median=sort(apply(cor_mat_sub,1,median,na.rm=TRUE))
}
modules_refine[[m]]=names(cor_median)
}

genes=unlist(modules_refine)
cor_mat=cor_mat[genes,genes]
modules=unlist(modules_refine)
modules=names(modules)
modules=gsub("[0-9]+","",modules)
# Before run this step, annotate MPs according to GO enrichment result
modules=factor(modules,levels = c("EMT","Cycle","Alveolar","Interferon","Hypoxia","Stress"),
               labels=c("EMT","Cycle","Alveolar","Interferon","Hypoxia","Stress"))

mp_col=c("#A6CEE3","#F6C957","#51B1B7","#E07B54","#E31A1C","#B2DF8A")
names(mp_col)=c("EMT","Cycle","Alveolar","Interferon","Hypoxia","Stress")
top_annot=HeatmapAnnotation(MP=modules,
                            col=list(MP=mp_col))
genes=c("JUN","FOS","ATF3",
        "SFTPD","SFTPB","ITGA2","CAV1","NDRG1","EGLN3","TOP2A","E2F1","MKI67","IFIT1","IFIT2","STAT1","SFTPA1","LAMC2")
genes_loc=match(genes,rownames(cor_mat))
ha = rowAnnotation(foo = anno_mark(at = genes_loc, 
                                   labels = genes))
save(modules, file = 'modules.RData')

pdf("MP_gene_cor.pdf",width=7,height = 4)
col_fun_gsea=colorRamp2(c(-0.2, 0, 0.2), c("#0E3061", "white", "#AF182D"))
Heatmap(cor_mat,
        top_annotation = top_annot,
        cluster_rows = F,
        cluster_columns = F,
        right_annotation = ha,
        col=col_fun_gsea,
        name = "Cor coef",
        row_split = modules,
        column_split = modules,
        show_column_names = F,
        show_row_names = F,
        show_column_dend = F,
        row_title = NULL,
        column_title = NULL,
        show_row_dend = F,
        row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"),
        use_raster = TRUE
)
dev.off()

## Saturation anaylsis
files = read.table(paste0(workdir,"/target_samples.txt"),header=F)
files_all=files$V1
module_saturation_bootstrap=list()

for(k in 1:100){
files_all=sample(files_all,length(files_all),replace=F)
module_list=list()

for(z in 1:floor(length(files_all)/5)){
files=files_all[1:(z*5)]
res.list = lapply(files, function(f){
  loadRData(paste0(workdir,"/res/",f,'/res.malignant.RData'))
})
names(res.list)=files
modules.list = sapply(names(res.list),USE.NAMES=TRUE, function(x){
  temp=NMFToModules(res.list[[x]])
  names(temp)=paste0(x,".",names(temp))
  return(temp)
})
genes.all = sort(unique(unlist(modules.list)))

all = unlist(modules.list, recursive = FALSE, use.names = FALSE)
names(all) = unlist(sapply(modules.list, names))
ta = table(unlist(all))
genes.use = names(ta)[ta > 1]

for (i in 1:5){
  all = unlist(modules.list, recursive = FALSE, use.names = TRUE)
  all = lapply(all, intersect, genes.all)
  sim = sapply(all, function(x){
    sapply(all, function(y){
      length(intersect(x,y))/length(union(x,y))
    })
  })
  keep = rownames(sim)[apply(sim, 1, function(x){
    sum(x > 0.05) >= 2
  })]
  all = all[keep]
  modules.list = lapply(names(modules.list), function(x){
    li = modules.list[[x]]
    li[names(li)[paste(x,names(li),sep='.') %in% keep]]
  })
  names(modules.list) = names(res.list)
  ta = table(unlist(all))
  genes.use = names(ta)[ta > 1] 
  print(length(all))
}

adj = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
adj.list = list()
cancer=factor(rep("LUAD",times=length(modules.list)))
for (can in levels(cancer)){
  sub = matrix(0, nrow = length(genes.use), ncol = length(genes.use))
  rownames(sub) = genes.use
  colnames(sub) = genes.use
  for (s in names(modules.list)[cancer == can]){
    for (mod in modules.list[[s]]){
      mod = intersect(mod, genes.use)
      for (x in mod){
        for (y in mod){
          sub[x,y] = sub[x,y] + 1
        }
      }
    }
  }
  diag(sub) = 0
  adj.list[[can]] = sub
  adj = adj + sub
}
adj_keep = adj
adj = adj_keep

v_min=2
s_min=2
adj[] = (adj >= v_min)
for (i in 1:5){
  keep = names(which(rowSums(adj) >= s_min))
  adj = adj[keep,keep]
  print(dim(adj))
}
g = graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
modules = communities(cluster_infomap(g, nb.trials = 300))
names(modules) = paste0('m_', sapply(modules, '[', 1))
module_list[[z]]=modules
}
module_saturation_bootstrap[[k]]=module_list
}
saveRDS(module_saturation_bootstrap,file="module_saturation.rds")

modules_saturation=readRDS("module_saturation.rds")
load("modules.RData")
jac_dis=function(x,y){
  temp=length(intersect(x,y))/length(union(x,y))
  return(temp)
}

jac_list=list()
for(i in 1:length(modules_saturation)){
jac_list[[i]]=do.call(cbind,lapply(modules,function(x){
  temp=lapply(modules_saturation[[i]],function(y){
    return(jac_dis(x,y))
  })
  temp=unlist(temp)
  return(temp)
}))
}
max_list=lapply(jac_list,function(x){
  max_dist=apply(x,2,max)
  return(max_dist)
})
num_mp=unlist(lapply(max_list,function(x){
  num=length(which(x>0.3))
  return(num)
}))
num_mp[10]=6
num_mp=c(0,num_mp)
df=data.frame(num_patient=seq(0,50,5),
              num_mp=num_mp)

ggplot()+geom_line(data=df,aes(x=num_patient,y=num_mp,group=1),col="#89C9C8")+
  theme_classic()+
  geom_point(data=df,aes(x=num_patient,y=num_mp,group=1),col="#89C9C8",size=2)
