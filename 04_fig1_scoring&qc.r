library(AUCell)
library(clusterCrit)
library(dplyr)
library(tidyverse)
library(parallel)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(ggpubr)
library(ggsignif)

### Due to data control and large size, no example data for this script
## Scoring
combined_object=readRDS("malignant_cell(sampled).RData") # This file should contain all malignant cells from all patients, but a downasampled version was supplied
cells_AUC = AUCell_run(GetAssayData(combined_object,slot="counts",assay="RNA"), modules,BPPARAM=BiocParallel::MulticoreParam(10))
scores=t(as.matrix(getAUC(cells_AUC)))
combined_object@meta.data=cbind(combined_object@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","cell_type","cell_subtype")],scores)

## Optimization and clustering
med_scale=function(x){
  temp_med=median(x)
  temp_mad=mad(x)
  return((x-temp_med)/temp_mad)
}
scores=apply(scores,2,med_scale)
s_vals=c()
ch_vals=c()
db_vals=c()
crit_list=list()

for(k in 2:10){
km=kmeans(scores_sampled,centers=k,iter.max = 100)
res=km$cluster
crit=intCriteria(scores_sampled,as.integer(res),c("Silhouette","Calinski_Harabasz","Davies_Bouldin"))
crit_list[[as.character(k)]]=crit
}
df=do.call(rbind,crit_list)
colnames(df)=c("Silhouette","Calinski_Harabasz","Davies_Bouldin")
rownames(df)=as.character(2:10)
df= df %>% as_tibble(rownames = "k") %>% pivot_longer(-c(k),names_to = c("idx"),values_to = "val")
df$k=as.numeric(df$k)
df=data.frame(df)
df$val=unlist(df$val)
pdf("km_index.pdf")
p=ggplot(data=df)+geom_line(aes(x=k,y=val,col=idx,group=idx))+
  theme_classic()+
  facet_grid(idx~.,scales="free")
plot(p)
dev.off()
# k=8 was selected
km=kmeans(scores,centers=8,iter.max = 1000)
res=km$cluster
combined_object$cluster=res

## Heatmap
cell_metadata=combined_object@meta.data
scores=cell_metadata[,c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")]
scores=t(apply(scores,2,med_scale))

# Prepare patient metadata file first
patient_metadata=read.table("patient_metadata.txt",sep="\t",header=T)
cell_metadata$cell_id=rownames(cell_metadata)
combined_metadata=merge(cell_metadata,patient_metadata,by.x="orig.ident",by.y="sample",all.x = T)
combined_metadata[which(combined_metadata$tissue %in% c("normal","nLung","adjacent_normal","distant_normal","intermediate_normal")),]$stage="Normal"
rownames(combined_metadata)=combined_metadata$cell_id
combined_metadata=combined_metadata[!is.na(combined_metadata$smoking),]
combined_metadata$smoking=factor(combined_metadata$smoking,levels = c("0","1"),labels = c("non_smoker","smoker"))
combined_metadata$cluster=factor(combined_metadata$cluster,levels=c("1","2","3","4","5","6","7","8"),
                                labels = c("C1_Low RNA count","C2_Interferon high","C3_Cycle","C4_Stress","C5_Interferon low","C6_Normal like","C7_Intermediate","C8_Hypoxia"))
# Downsample cells to avoid overrepresentation and RAM leak
cell_sampled=split(combined_metadata,f=factor(combined_metadata$orig.ident))
cell_sampled=lapply(cell_sampled,function(x){
  if(nrow(x)>50){
    temp_id=rownames(x)
    temp_id=sample(temp_id,50,replace=F)
    return(temp_id)
  }else{
    return(rownames(x))
  }
})
cell_sampled=unlist(cell_sampled)

cell_sampled=unlist(lapply(unique(combined_metadata$cluster),function(x){
  temp_df=combined_metadata[which(combined_metadata$cluster==x),]
  temp_df_sampled=temp_df[sample(1:nrow(temp_df),200,replace = F),]
  return(temp_df_sampled$cell_id)
}))

cell_metadata_sampled=cell_metadata[cell_sampled,]
scores_sampled=scores[,cell_sampled]
res_sampled=factor(res[colnames(scores_sampled)])
metadata_sampled=combined_metadata[colnames(scores_sampled),]

scores_sampled=scores_sampled[,rownames(metadata_sampled)]
res_sampled=res_sampled[rownames(metadata_sampled)]

annot=HeatmapAnnotation(
  smoking=metadata_sampled$smoking,
  stage=metadata_sampled$stage,
  sex=metadata_sampled$sex,
  ethnicity=metadata_sampled$ethnicity,
  cell_type=metadata_sampled$cell_subtype,
  `log2(count_RNA)`=log2(metadata_sampled$nCount_RNA),
  annotation_legend_param = list( smoking=list(title = "Smoking"),
                                  stage=list(title = "Stage"),
                                  sex=list(title = "Sex"),
                                  ancestry=list(title="Ancestry"),
                                  cell_type=list(title="Cell type")
  ),
  col=list(
    sex=c("female"="#DC143C","male"="#4169E1"),
    cell_type=c("AT2 cell"="#54bbd5","Malignant cell"="#f17f73"),
    stage=c("I"="#DCE9F4","II"="#ABD0F1","III"="#FF9797","IV"="#E56F5E","Normal"="#9EC4BE"),
    smoking=c("Smoker"="#DC143C","Non-smoker"="#FFF0F5"),
    ethnicity=c("Asian"="#FFA500","European"="#FBE8D5"),
    `log2(count_RNA)`=colorRamp2(c(8, 12, 16), c("#377eb8", "white", "#e4211c"))
    ))
col_fun_gsea=colorRamp2(c(-1.5, 0, 1.5), c("#377eb8", "white", "#e4211c"))

pdf("cluster_heatmap_k=8.pdf",width=12.5,height=3.7)
Heatmap(scores_sampled,name="MP z-score",top_annotation = annot,col=col_fun_gsea,
             cluster_rows = T,
             column_split = metadata_sampled$cluster,
             column_title_gp = grid::gpar(fontsize = 8),
             row_title_rot = 0,
             show_column_names = F,
             cluster_columns = F)
dev.off()

combined_metadata$`log2(count_RNA)`=log2(combined_metadata$nCount_RNA)
ggplot()+geom_violin(data=combined_metadata,aes(x=cluster,y=`log2(count_RNA)`),fill="#DC143C")+
  geom_boxplot(data=combined_metadata,aes(x=cluster,y=`log2(count_RNA)`),width=0.3)+
  theme_classic()

pdf("cluster_RNA_count.pdf",width=8,height=4)
my_comparisons <- list( c("C1_Low RNA count", "C2_Interferon high"), 
                        c("C1_Low RNA count", "C3_Cycle"), 
                        c("C1_Low RNA count", "C4_Stress"),
                        c("C1_Low RNA count", "C5_Interferon low"),
                        c("C1_Low RNA count", "C6_Normal like"),
                        c("C1_Low RNA count", "C7_Intermediate"),
                        c("C1_Low RNA count", "C8_Hypoxia")
                        )
ggplot(data=combined_metadata,aes(x=cluster,y=`log2(count_RNA)`))+
  geom_violin(fill="#ABD0F1")+
  geom_boxplot(outlier.size = 0.2,width=0.3)+
  stat_compare_means(comparisons = my_comparisons,test = "wilcox.test",p.adjust.methods="fdr")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))+
  labs(x=NULL)
dev.off()

## Cell cluster with low RNA count should be removed here
combined_object=subset(combined_object,subset=cluster != "C1_Low RNA count")
save(combined_object,file-"malignant_cell.RData")
