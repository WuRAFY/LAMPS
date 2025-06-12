library(SpaCET)
library(Seurat)
library(RhpcBLASctl)
library(ggpubr)
library(spdep)
library(ggplot2)
library(tidyverse)
library(parallel)
library(ggrastr)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)

Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1", R_MAX_NUM_THREADS = "1")
blas_set_num_threads(1)
omp_set_num_threads(1)

## Run SpaCET
# Large computation resources required, do not run on personal computer
# Due to data control and large file size, example data was not supplied for this step
# Result from SAW and convert to RDS format
load("./input/modules.RData")
path=c("ST1_bin100.rds",
       "ST2_bin100.rds",
       "ST3_bin100.rds",
       "ST4_bin100.rds",
       "ST5_bin100.rds")
samples=c("ST1","ST2","ST3","ST4","ST5")
spacet=list()
i=1
for(p in path){
  seurat=readRDS(p)
  DefaultAssay(seurat)="Spatial"
  SpaCET_obj=convert.Seurat(seurat)
  SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes=200)
  res=SpaCET.deconvolution(
    SpaCET_obj,
    cancerType="LUAD",
    adjacentNormal = FALSE,
    coreNo = 30
  )
  res=SpaCET.GeneSetScore(res, GeneSets = modules)
  spacet[[samples[i]]]=res
  i=i+1
}
saveRDS(spacet,file="spacet_bin100.rds")

## Malignant cell deconvolution & CCI analysis (edit from SpaCET code)
source("deconvote.r")
source("SpaCET_edit.r")
for(i in names(spacet)){
  res=spacet[[i]]
  res=SpaCET.deconvolution.malignant_edit(res, coreNo = 15, malignantCutoff = 0.7)
  tumor_spot <- res@results$deconvolution$propMat["Malignant",]
  tumor_spot <- names(tumor_spot[tumor_spot>=0.7])
  scores=res@results$GeneSetScore
  scores[,-which(colnames(scores) %in% tumor_spot)]=NA
  scores=t(apply(scores,1,med_scale))
  res@results$GeneSetScore_raw=res@results$GeneSetScore
  res@results$GeneSetScore=scores
  res=SpaCET.CCI.colocalization_edit(res)
  res=SpaCET.CCI.LRNetworkScore(res,coreNo=15)
  spacet[[i]]=res
}
saveRDS(spacet,file="spacet_clone.rds")

pdf("spacet_composition.pdf",width=50,height=5)
for(s in names(spacet)){
res=spacet[[s]]
plot_list=list()
for(cell_type in c("Malignant cell state A","Malignant cell state B","Malignant cell state C","Endothelial","Macrophage","T CD4","T CD8","B cell","Plasma","CAF")){
p=SpaCET.visualize.spatialFeature(
  res,
  spatialType = "CellFraction", 
  spatialFeatures=cell_type,
  imageBg = FALSE,
  pointSize = 0.1,
  colors=c("#2B8DC4","cyan","#FFF200","pink","red","red")
  )+theme(aspect.ratio = 1)
plot_list[[cell_type]]=ggrastr::rasterise(p,layers='Point',dpi = 300)}
p_arr=ggarrange(plotlist=plot_list,nrow=1,ncol=10)
plot(p_arr)
}
dev.off()

## Calculate bivariate Moran's I
spacet=readRDS("./input/spacet_clone.rds")
I_all=list()
p_all=list()
j=1
for(r in c(1,2,5,10,20,30)){
k=1
I_aggr=list()
p_aggr=list()
for(i in names(spacet)){
SpaCET_obj=spacet[[i]]
res_deconv <- SpaCET_obj@results$deconvolution$propMat
tme_deconv = res_deconv[-grep("Malignant",rownames(res_deconv)),]
tme_deconv = tme_deconv[-which(rownames(tme_deconv) %in% c("Unidentifiable","cDC","pDC","cDC1 CLEC9A","cDC2 CD1C","cDC3 LAMP3","Mast")),]
tme_deconv_sub = res_deconv[c("CAF","Endothelial","Plasma","B cell","T CD4","T CD8","NK","Macrophage","Neutrophil"),]
tme_deconv_sub = colSums(tme_deconv_sub)
cell_type=c(rownames(tme_deconv),"Lymphocyte")
tme_deconv = t(scale(tme_deconv,center=F,scale=tme_deconv_sub))
tme_deconv = cbind(tme_deconv,Lymphocyte=rowSums(tme_deconv[,c("Plasma","B cell","T CD4","T CD8")]))
tme_deconv[is.na(tme_deconv)]=0
tme_deconv[is.infinite(tme_deconv)]=0
malig_deconv = res_deconv[grep("Malignant cell state",rownames(res_deconv)),]
malig_deconv = data.frame(t(scale(malig_deconv,center=F, scale=colSums(malig_deconv))))
malig_deconv[is.na(malig_deconv)]=0

cord=SpaCET_obj@input$spotCoordinates
cord$pxl_row=(cord$pxl_row-1)/100
cord$pxl_col=(cord$pxl_col-1)/100
coords=as.matrix(cord[,c("pxl_row","pxl_col")])

nb <- dnearneigh(coords,d1=0,d2=r)
empty_neighbors <- which(card(nb) == 0)
if (length(empty_neighbors) > 0) {
  cat("Removing", length(empty_neighbors), "spots without neighbors.\n")
  coords_filter <- coords[-empty_neighbors,]
  nb <- dnearneigh(coords_filter, d1 = 0, d2 = r)  # Recalculate neighbors
}
lw <- nb2listw(nb, style = "W")

I_list=list()
p_list=list()
for(state in colnames(malig_deconv)){
  temp=mclapply(colnames(tme_deconv),function(x){
    res=moran_bv(x=malig_deconv[,state], 
             y=tme_deconv[,x], lw, nsim = 200, scale = TRUE)
    I=res$t0
    p=mean(res$t>=res$t0)
    p=min(p,1-p)
    return(c(I,p))
  },mc.cores=28)
  temp=do.call(rbind,temp)
  temp_I=temp[,1]
  temp_p=temp[,2]
  I_list[[state]]=temp_I
  p_list[[state]]=temp_p
}
  I_list=do.call(cbind,I_list)
  p_list=do.call(cbind,p_list)
  rownames(I_list)=cell_type
  rownames(p_list)=cell_type
  colnames(I_list)=paste0("ST",as.character(k),"_",gsub("Malignant.cell.state.","",colnames(I_list)))
  colnames(p_list)=paste0("ST",as.character(k),"_",gsub("Malignant.cell.state.","",colnames(p_list)))
  I_aggr[[i]]=I_list
  p_aggr[[i]]=p_list
  k=k+1
}
I_all[[j]]=do.call(cbind,I_aggr)
p_all[[j]]=do.call(cbind,p_aggr)
j=j+1
}
scores=list()
composition=list()
cluster=list()
for(i in names(spacet)){
SpaCET_obj=spacet[[i]]
scores[[i]]=SpaCET_obj@results$GeneSetScore_raw
composition[[i]]=SpaCET_obj@results$deconvolution$propMat
cluster[[i]]=SpaCET_obj@results$deconvolution$clustering
}
saveRDS(list(I_all,p_all,scores,composition,cluster),file="moran_multi.rds")

moran=readRDS("./input/moran_multi.rds")
moran_I=moran[[1]]
moran_p=moran[[2]]
scores=moran[[3]]
composition=moran[[4]]
spot=moran[[5]]

coef_list=list()
k=1
# Calculate MP scores using Spearman correlation
for(i in names(scores)){
  temp_scores=scores[[i]]
  temp_comp=composition[[i]]
  temp_comp=temp_comp[grep("Malignant cell state",rownames(temp_comp)),]
  temp_spot=names(spot[[i]])
  temp_scores=temp_scores[,temp_spot]
  temp_comp=temp_comp[,temp_spot]
  res=rcorr(t(temp_comp),t(temp_scores),type="spearman")
  res=res$r
  res=res[grep("Malignant cell state",rownames(res)),c("Alveolar","Interferon","Cycle","Hypoxia","EMT","Stress")]
  rownames(res)=paste0("ST",as.character(k),"_",gsub("Malignant cell state ","",rownames(res)))
  coef_list[[i]]=res
  k=k+1
}
coef_list=t(do.call(rbind,coef_list))
col_score=colorRamp2(c(-0.2, 0, 0.2), c("#2B8DC4", "white", "#FD4C46"))
col_moran=colorRamp2(c(-0.15, 0, 0.15), c("#2B8DC4", "white", "#FD4C46"))
moran_1=moran_I[[1]][c("Endothelial","Macrophage","Neutrophil","NK","B cell","T CD4","T CD8","Plasma","CAF"),]
moran_2=moran_I[[2]][c("Endothelial","Macrophage","Neutrophil","NK","B cell","T CD4","T CD8","Plasma","CAF"),]
moran_5=moran_I[[3]][c("Endothelial","Macrophage","Neutrophil","NK","B cell","T CD4","T CD8","Plasma","CAF"),]
moran_10=moran_I[[4]][c("Endothelial","Macrophage","Neutrophil","NK","B cell","T CD4","T CD8","Plasma","CAF"),]
moran_20=moran_I[[5]][c("Endothelial","Macrophage","Neutrophil","NK","B cell","T CD4","T CD8","Plasma","CAF"),]
moran_30=moran_I[[6]][c("Endothelial","Macrophage","Neutrophil","NK","B cell","T CD4","T CD8","Plasma","CAF"),]

h1=Heatmap(coef_list,col = col_score,
           cluster_rows = F,
           column_split = 3,
           name = "MP score",
           cluster_columns = T,
           column_title = c("Precursor","IFN-responsive","Proliferative"))
h2=Heatmap(moran_1,col = col_moran,
           name = "Moran's I (r=50 um)",
           cluster_rows = F,
           cluster_columns = T)
h3=Heatmap(moran_2,col = col_moran,
           name = "Moran's I (r=100 um)",
           cluster_rows = F,
           cluster_columns = F)
h4=Heatmap(moran_5,col = col_moran,
           name = "Moran's I (r=250 um)",
           cluster_rows = F,
           cluster_columns = F)
h5=Heatmap(moran_10,col = col_moran,
           name = "Moran's I (r=500 um)",
           cluster_rows = F,
           cluster_columns = T)
h6=Heatmap(moran_20,col = col_moran,
           name = "Moran's I (r=1000 um)",
           cluster_rows = F,
           cluster_columns = F)
h7=Heatmap(moran_30,col = col_moran,
           name = "Moran's I (r=1500 um)",
           cluster_rows = F,
           cluster_columns = F)
pdf("moran_multi_scale.pdf",width=8,height = 16)
h1 %v% h2 %v% h3 %v% h4 %v% h5 %v% h6 %v% h7 
dev.off()

tme_cell=c("CAF","Endothelial","Plasma","B cell","T CD4","T CD8","NK","Macrophage","Neutrophil")
k=1
cci=list()
for(i in names(spacet)){
  SpaCET_obj=spacet[[i]]
  cell_type=rownames(SpaCET_obj@results$deconvolution$propMat)
  state=cell_type[grep("Malignant cell state",cell_type)]
  cell_pair=data.frame(expand.grid(state,tme_cell))
  res=apply(cell_pair,1,function(x){
    temp=SpaCET.CCI.cellTypePair_edit(SpaCET_obj,cellTypePair=c(x[1],x[2]))
    return(temp)
  })
  cell_pair=cbind(cell_pair,t(res))
  colnames(cell_pair)=c("Clone","TME","D","p")
  cell_pair$Clone=gsub("Malignant cell state ",paste0("ST",as.character(k),"_"),cell_pair$Clone)
  cci[[i]]=cell_pair
  k=k+1
}
saveRDS(cci,file="spacet_cci.rds")

hr <- hclust(dist(t(coef_list), method = "euclidean"), method="complete") 
hr = cutree(hr,3)
hr=factor(hr,level=c(1,2,3),labels=c("Precursor","Proliferative","IFN-responsive"))
cci=readRDS("./input/spacet_cci.rds")
cci_d=lapply(cci,function(x){
  temp=matrix(x$D,length(unique(x$Clone)),length(unique(x$TME)))
  rownames(temp)=unique(x$Clone)
  colnames(temp)=unique(x$TME)
  return(temp)
})
cci=lapply(cci,function(x){
  temp=matrix(x$p,length(unique(x$Clone)),length(unique(x$TME)))
  rownames(temp)=unique(x$Clone)
  colnames(temp)=unique(x$TME)
  return(temp)
})
cci=do.call(rbind,cci)
cci_d=do.call(rbind,cci_d)
cci=t(cci)
cci=-log10(cci)
cci=cci[c("Endothelial","Macrophage","Neutrophil","NK","B cell","T CD4","T CD8","Plasma","CAF"),]

col_cci=colorRamp2(c(0, 1, 2), c("white", "white", "#FD4C46"))
h_cci=Heatmap(cci,col = col_cci,
              cluster_rows = F,
              name = "CCI significance\n(-log10(p value))")
pdf("moran&cci.pdf",width=8,height = 8)
h1 %v% h3 %v% h_cci
dev.off()