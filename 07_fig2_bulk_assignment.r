library(BayesPrism)
library(GSVA)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggsignif)
library(ggpubr)

## Input files for this script were supplied from the second step.

## Using Bayesprism to deconvolute malignant cell profile
# Do not use all cells as reference, one dataset is enough
# Large computing resources required, do not run on personal computer
load("sub_cluster.RData")
combined_object=merge(type_list[[1]],y=type_list[-1])
meta=data.frame(patient_id=combined_object@meta.data$orig.ident,
                cell_type=combined_object@meta.data$cell_type,
                cell_subtype=combined_object@meta.data$cell_subtype,
                cluster=combined_object@meta.data$seurat_clusters)
rownames(meta)=rownames(combined_object@meta.data)
meta=as.data.frame(meta)
#reference matrix
ensg_id=read.table("features.tsv",sep="\t")
ref.dat=GetAssayData(combined_object,slot="counts",assay="RNA")
rownames(ref.dat)=ensg_id$V1
ref.dat=t(as.matrix(ref.dat))
#drop rare cell type
enough_cell_subtype=names(which(table(meta$cell_subtype)>20))
meta=meta[meta$cell_subtype %in% enough_cell_subtype,]
ref.dat=ref.dat[rownames(meta),]
#filter genes
ref.dat <- cleanup.genes (input=ref.dat,
                           input.type="count.matrix",
                           species="hs", 
                           gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                           exp.cells=5)
#bulk matrix from STARSEM flow
mat=read.table("bulk_expr_count.txt",sep="\t",header=T,row.names=2)
mat=t(mat[,-c(1,2,3)])
ensg_id=as.character(colnames(mat))
ensg_id=unlist(strsplit(ensg_id,split="[.]"))
ensg_id=ensg_id[grep(ensg_id,pattern="ENSG*")]
colnames(mat)=ensg_id
#run BayesPrism
myPrism <- new.prism(
  reference=ref.dat, 
  mixture=mat,
  input.type="count.matrix", 
  cell.type.labels = meta$cell_type, 
  cell.state.labels = meta$cell_subtype,
  key="malignant",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)
res <- run.prism(prism = myPrism, n.cores=as.numeric(args[1]))
save(res,meta,file="bayesprism.rdata")
st=c("type","state")
theta=c("final","first")
for(z in 1:2){
composition=get.fraction (bp=res,
            which.theta=theta[z],
            state.or.type=st[z])
celltype=colnames(composition)
expr_mat=c()
expr_ident=c()
expr_celltype=c()
for(i in celltype){
expr=get.exp(bp=res,
             state.or.type=st[z],
             cell.name=i)
expr=t(expr)
expr_ident=c(expr_ident,colnames(expr))
expr_celltype=c(expr_celltype,rep(i,times=ncol(expr)))
colnames(expr)=paste(colnames(expr),i,sep="_")
expr_mat=cbind(expr_mat,expr)
write.table(expr,paste0("expr_",i,".txt"),sep="\t",quote=F)
}
}

## Patient scoring&assignment
# The output of bayesprism is count matrix, convert to TPM matrix before run this step
load("./input/modules.RData")
# ssGSEA was recommended to run on tpm/rpkm data
med_scale=function(x){
  temp_med=median(x)
  temp_mad=mad(x)
  return((x-temp_med)/temp_mad)
}
mat_tpm=read.table("./input/expr_malignant_tpm.txt",sep="\t",header=T)
res_ssgsea= gsva(
   expr = as.matrix(mat_tpm),
   gset.idx.list = modules,
   method='ssgsea',
   kcdf='Gaussian',
   abs.ranking=TRUE)
res_ssgsea_scale=t(apply(res_ssgsea,1,med_scale))

# Patient assigment
x_test=function(x,y){
  temp=chisq.test(x,y,simulate.p.value = TRUE)
  signif = symnum(temp$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0,0.001,0.01, 0.05, 1),
                  symbols = c("***", "**", "*", " "))
  return(signif)}
kw_test=function(x,g){
  temp=kruskal.test(x,g,na.action=na.exclude)
  signif = symnum(temp$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0,0.001,0.01, 0.05, 1),
                  symbols = c("***", "**", "*", " "))
  return(signif)}

combined_object_sampled=readRDS("./input/trajectory.rds")
cell_metadata=combined_object_sampled@meta.data
#read clinical data and two scores
malignant_metadata=cell_metadata
sc_scores=t(malignant_metadata[,c('Stress','Alveolar','Cycle','EMT','Hypoxia','Interferon')])
bulk_scores=res_ssgsea_scale[c('Stress','Alveolar','Cycle','EMT','Hypoxia','Interferon'),]
#correlation
cor_mat=cor(bulk_scores,sc_scores,method = "pearson")
cor_mat=data.frame(t(cor_mat))
cor_mat_list=split(cor_mat,f=cell_metadata$branch)
cor_median=lapply(cor_mat_list,function(x){
  temp=apply(x,2,median)
  return(temp)
})
cor_median=data.frame(do.call(cbind,cor_median))
cor_median$cluster=colnames(cor_median)[apply(cor_median,1,which.max)]
#read clinical info
info_df=read.table("./input/bulk_info.txt",sep="\t",header=T)
#info_df$Branch=cor_median$cluster

kw_sig=apply(info_df[,c("Purity","TMB")],2,kw_test,g=info_df$Branch)
x_sig=apply(info_df[,c("RNA_subtype","Smoking","Stage","Sex","Ancestry","EGFR_status","KRAS_status","TP53_status")],2,x_test,y=info_df$Branch)
mp_sig=apply(bulk_scores,1,kw_test,g=info_df$Branch)
info_df=info_df[order(info_df$TP53_status,info_df$EGFR_status,info_df$KRAS_status,info_df$TMB),]
info_df$Branch=factor(info_df$Branch,levels = c("Precursor","IFN-responsive","Proliferative"))
bulk_scores= bulk_scores[,info_df$Patient]
annot_top=HeatmapAnnotation(
  Branch=info_df$Branch,
  col=list(Branch = c("Precursor"="#51B1B7","IFN-responsive"="#FF8C00","Proliferative"="#DC143C")))
annot_bot=HeatmapAnnotation(
  TP53_status=info_df$TP53_status,
  EGFR_status=info_df$EGFR_status,
  KRAS_status=info_df$KRAS_status,
  TMB=info_df$TMB,
  RNA_subtype=info_df$RNA_subtype,
  purity=info_df$Purity,
  ethnicity=info_df$Ancestry,
  smoking=info_df$Smoking,
  stage=info_df$Stage,
  sex=info_df$Sex,
  annotation_label = paste0(c("TP53","EGFR","KRAS","TMB","RNA subtype","Purity","Ancestry","Smoking","Stage","Sex"),
                            c(x_sig,kw_sig)[c("TP53_status","EGFR_status","KRAS_status","TMB","RNA_subtype","Purity","Ancestry","Smoking","Stage","Sex")]),
  annotation_legend_param = list(smoking=list(title = "Smoking"),
                                 stage=list(title = "Stage"),
                                 sex=list(title = "Sex"),
                                 RNA_subtype=list(title="RNA subtype"),
                                 ethnicity=list(title="Ancestry"),
                                 EGFR_status=list(title = "EGFR"),
                                 KRAS_status=list(title = "KRAS"),
                                 TP53_status=list(title = "TP53"),
                                 purity=list(title = "Purity"),
                                 TMB= list(title = "TMB")
  ),
  gp = gpar(col = "white"),
  col=list(
    RNA_subtype=c("Normal"="#9EC4BE","TRU-I"="#4169E1","PI"="#DC143C","PP"="#FF9797","TRU"="#FFA500"),
    TMB=colorRamp2(c(0,5,10,15,20), c("#E0F3DB","#95CFBB","#7CBDC1","#5CACC7","#43A2CA")),
    purity=colorRamp2(c(0,0.6,0.7,0.8,0.9,1), c("white","#FEEBE2","#F089AC","#DF639E","#CD3690","#C51B8A")),
    sex=c("female"="#FFE4E1","male"="#6CA6CD"),
    stage=c("I"="#DCE9F4","II"="#ABD0F1","III"="#FF9797","IV"="#E56F5E","Normal"="#9EC4BE"),
    smoking=c("Smoker"="#DC143C","Non-smoker"="#FFF0F5"),
    ethnicity=c("Asian"="#FFA500","European"="#FBE8D5"),
    KRAS_status=c("MT"="#DC143C","WT"="#FFF0F5"),
    EGFR_status=c("MT"="#DC143C","WT"="#FFF0F5"),
    TP53_status=c("MT"="#DC143C","WT"="#FFF0F5")))
col_fun = colorRamp2(c(-1.8, 0, 1.8), c("#2B8DC4", "white", "#FD4C46"))
pdf("bulk_assignment.pdf",width =25 ,height = 4)
Heatmap(bulk_scores,name="MP z-score",top_annotation = annot_top,
        bottom_annotation = annot_bot,
        col=col_fun,
        cluster_rows = T,
        row_labels = paste0(rownames(bulk_scores),mp_sig),
        column_split = info_df$Branch,
        row_title_rot = 0,
        show_column_names = F,
        cluster_columns = F)
dev.off()
chi.test <- function(a, b) {
  return(chisq.test(cbind(a, b)))
}
chi_plot=function(data,x,g,comparisons,label,col){
temp_perc=data[!is.na(data[,x]),] %>% group_by_at(c(x,g)) %>% count() %>% group_by_at(c(g)) %>%
  mutate(countT= sum(n)) %>%
  group_by_at(c(x), add=TRUE) %>%
  mutate(per=round(100*n/countT,2))
p=ggplot(data=temp_perc,aes_string(x=g,y="n",fill=x))+
  geom_bar(data=temp_perc,aes_string(x=g,y="per",fill=x),stat="identity",position="stack",alpha = 0.9)+
  scale_fill_manual(values=col)+
  geom_signif(comparisons = comparisons,
              map_signif_level = TRUE,
              y_position = c(105,110,115),
              textsize = 5,
              test = "chi.test")+
  theme_classic()+
  scale_y_continuous(breaks = seq(0, 100, 20),labels = seq(0,1,0.2))+
  labs(x="",y="",title=label)+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1))
return(p)
}
wilcox_test=function(data,x,g,comparisons,label,col){
  p=ggplot(data=data,aes_string(x=g,y=x))+
    geom_boxplot(alpha = 0.9)+
    geom_signif(comparisons = comparisons,
                map_signif_level = TRUE,
                step_increase=0.05,
                textsize = 5,
                test = "wilcox.test")+
    theme_classic()+
    labs(x="",y="",title=label)+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45,hjust = 1))
  return(p)
}
p_gender=chi_plot(info_df,
         "Sex",
         "Branch",
         list(c("Precursor","IFN-responsive"),
              c("IFN-responsive","Proliferative"),
              c("Precursor","Proliferative")),
         label = "Gender",
         col=c("female"="#FFE4E1","male"="#6CA6CD"))
p_smoking=chi_plot(info_df,
                  "Smoking",
                  "Branch",
                  list(c("Precursor","IFN-responsive"),
                       c("IFN-responsive","Proliferative"),
                       c("Precursor","Proliferative")),
                  label = "Smoking",
                  col=c("Smoker"="#DC143C","Non-smoker"="#FFF0F5"))
p_stage=chi_plot(info_df,
                 "Stage",
                 "Branch",
                 list(c("Precursor","IFN-responsive"),
                      c("IFN-responsive","Proliferative"),
                      c("Precursor","Proliferative")),
                 label = "Stage",
                 col=c("I"="#DCE9F4","II"="#ABD0F1","III"="#FF9797","IV"="#E56F5E","Normal"="#9EC4BE"))
p_ancestry=chi_plot(info_df,
                     "Ancestry",
                     "Branch",
                     list(c("Precursor","IFN-responsive"),
                          c("IFN-responsive","Proliferative"),
                          c("Precursor","Proliferative")),
                     label = "Ancestry",
                     col=c("Asian"="#FFA500","European"="#FBE8D5"))
p_ethnicity_2=chi_plot(info_df,
                     "Branch",
                     "Ancestry",
                     list(c("Asian","European")),
                     label = "Ancestry",
                     col=c("Precursor"="#51B1B7","IFN-responsive"="#FF8C00","Proliferative"="#DC143C"))
p_tp53=chi_plot(info_df,
                "TP53_status",
                "Branch",
                list(c("Precursor","IFN-responsive"),
                     c("IFN-responsive","Proliferative"),
                     c("Precursor","Proliferative")),
                label = "TP53",
                col=c("MT"="#DC143C","WT"="#FFF0F5"))
p_egfr=chi_plot(info_df,
                "EGFR_status",
                "Branch",
                list(c("Precursor","IFN-responsive"),
                     c("IFN-responsive","Proliferative"),
                     c("Precursor","Proliferative")),
                label = "EGFR",
                col=c("MT"="#DC143C","WT"="#FFF0F5"))
p_kras=chi_plot(info_df,
                "KRAS_status",
                "Branch",
                list(c("Precursor","IFN-responsive"),
                     c("IFN-responsive","Proliferative"),
                     c("Precursor","Proliferative")),
                label = "KRAS",
                col=c("MT"="#DC143C","WT"="#FFF0F5"))
p_rna=chi_plot(info_df,
               "RNA_subtype",
               "Branch",
               list(c("Precursor","IFN-responsive"),
                    c("IFN-responsive","Proliferative"),
                    c("Precursor","Proliferative")),
               label = "Trancriptomic subtype",
               col=c("Normal"="#9EC4BE","TRU-I"="#4169E1","PI"="#DC143C","PP"="#FF9797","TRU"="#FFA500"))
p_rna2=chi_plot(info_df,
               "Branch",
               "RNA_subtype",
               list(c("TRU-I","TRU"),
                    c("PI","PP")),
               label = "Trancriptomic subtype",
               col=c("Precursor"="#51B1B7","IFN-responsive"="#FF8C00","Proliferative"="#DC143C"))

p_tmb=ggplot(data=info_df,aes_string(x="Branch",y="TMB"))+
  geom_violin(aes(fill=Branch),alpha=0.9)+
  geom_boxplot(fill="white",outliers = FALSE,width=0.15)+
  scale_fill_manual(values=c("Precursor"="#51B1B7","IFN-responsive"="#FF8C00","Proliferative"="#DC143C"))+
  geom_signif(comparisons =list(c("Precursor","IFN-responsive"),
                                c("IFN-responsive","Proliferative"),
                                c("Precursor","Proliferative")) ,
              map_signif_level = TRUE,
              step_increase=0.05,
              textsize = 5,
              test = "wilcox.test")+
  theme_classic()+
  labs(x="",y="",title="TMB")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1))
p_purity=ggplot(data=info_df,aes_string(x="Branch",y="Purity"))+
  geom_violin(aes(fill=Branch),alpha=0.9)+
  geom_boxplot(fill="white",outliers = FALSE,width=0.15)+
  scale_fill_manual(values=c("Precursor"="#51B1B7","IFN-responsive"="#FF8C00","Proliferative"="#DC143C"))+
  geom_signif(comparisons =list(c("Precursor","IFN-responsive"),
                                c("IFN-responsive","Proliferative"),
                                c("Precursor","Proliferative")) ,
              textsize = 5,
              map_signif_level = TRUE,
              step_increase=0.05,
              test = "wilcox.test")+
  theme_classic()+
  labs(x="",y="",title="Purity")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1))
plot_list=list(p_gender,p_smoking,p_stage,p_ethnicity,p_tp53,p_egfr,p_kras,p_rna,p_tmb,p_purity)
pdf("bulk_feature.pdf",width=16,height = 8)
ggarrange(plotlist = plot_list,nrow=2,ncol = 5)
dev.off()
