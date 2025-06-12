library(ggplot2)
library(NMF)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

## TME clustering
mat=read.table("./input/TME_composition.txt",sep="\t",header=T,row.names = 1)
mat=as.matrix(mat)
mat=t(mat/rowSums(mat))
#run NMF
W_list=list()
H_list=list()
cla_list=list()
for(rank in 2:12){
res=nmf(mat,rank,nrun=10)
W_list[[as.character(rank)]]=basis(res)
H_list[[as.character(rank)]]=coef(res)
cla_list[[as.character(rank)]]=predict(res)
}
cla=do.call(cbind,cla_list)
write.table(cla,file=paste0("./class_multi.txt"),sep="\t",quote=FALSE)

## Plot
x_test=function(x,y){
  temp=chisq.test(x,y,simulate.p.value = TRUE)
  signif = symnum(temp$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0,0.001,0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
  return(signif)}
kw_test=function(x,g){
  temp=kruskal.test(x,g,na.action=na.exclude)
  signif = symnum(temp$p.value, corr = FALSE, na = FALSE,
                  cutpoints = c(0,0.001,0.01, 0.05, 0.1, 1),
                  symbols = c("***", "**", "*", ".", " "))
  return(signif)}

#load metadata
info_df=read.table("./input/bulk_info.txt",sep="\t",header = T)
rownames(info_df)=info_df$Patient
cell_type=c("Lymphloid","Lymphloid","Endothelial","Stromal","Myeloid","Myeloid","Myeloid","Lymphloid","Myeloid","Lymphloid","Lymphloid","Lymphloid","Lymphloid")
cell_type=factor(cell_type,levels=c("Lymphloid","Myeloid","Stromal","Epithelial","Endothelial"))
#load TME subtype from A. Bagaev
load("./input/tme_score_bageav.rdata")
fge_func=c("Anti-tumor Immune infiltrate","Anti-tumor Immune infiltrate",
           "Anti-tumor Immune infiltrate","Anti-tumor Immune infiltrate",
           "Anti-tumor Immune infiltrate","Anti-tumor Immune infiltrate",
           "Anti-tumor Immune infiltrate","Anti-tumor Immune infiltrate",
           "Pro-tumor Immune infiltrate","Anti-tumor Immune infiltrate",
           "Anti-tumor Immune infiltrate","Anti-tumor Immune infiltrate",
           "Pro-tumor Immune infiltrate","Pro-tumor Immune infiltrate",
           "Angiogenesis&Fibroblasts","Angiogenesis&Fibroblasts",
           "Pro-tumor Immune infiltrate","Pro-tumor Immune infiltrate",
           "Pro-tumor Immune infiltrate","Pro-tumor Immune infiltrate",
           "Pro-tumor Immune infiltrate","Pro-tumor Immune infiltrate",
           "Angiogenesis&Fibroblasts","Angiogenesis&Fibroblasts",
           "Angiogenesis&Fibroblasts","Angiogenesis&Fibroblasts",
           "Angiogenesis&Fibroblasts","Tumor Proliferation",
           "EMT")
fge_func=factor(fge_func,levels=c("Angiogenesis&Fibroblasts","Pro-tumor Immune infiltrate","Anti-tumor Immune infiltrate","EMT","Tumor Proliferation"))
#load Ecotyper subtype
load("./input/Ecotpyer.RData")
info_df$ce_subtype=factor(ce_subtype$Carcinoma.Ecotype,levels=c(paste0("CE",1:10)))

#load article immune subtype from  V. Thorsson
load("./input/tme_score_thorsson.rdata")


#reorder
info_df$TME.subtype..k.4.=factor(info_df$TME.subtype..k.4.,
                                 levels=c("Initiating","Fibrotic","Immune recruiting","Immune activated"))
info_df$RNA_subtype=factor(info_df$RNA_subtype,levels = c("TRU-I","TRU","PI","PP"))
info_df=info_df[order(info_df$TME.subtype..k.4.,info_df$Ancestry,info_df$RNA_subtype,info_df$Purity),]
mat=mat[,rownames(info_df)]
fge_score=fge_score[,rownames(info_df)]
abundance_mat=abundance_mat[,rownames(info_df)]
immune_score=immune_score[,rownames(info_df)]

col_fun = colorRamp2(c(0, 0.075, 0.15), c("#377eb8", "white", "#e4211c"))
col_fun_gsea=colorRamp2(c(-1.2, 0, 1.2), c("#377eb8", "white", "#e4211c"))
col_fun_eco=colorRamp2(c(0, 0.1, 0.2), c("#377eb8", "white", "#e4211c"))
ce_col=c(brewer.pal(10,"Paired"))
names(ce_col)=levels(info_df$ce_subtype)
pdf(paste0("TME_heatmap.pdf"),width =25 ,height = 20)
kw_sig=apply(info_df[,c("Purity","TMB")],2,kw_test,g=factor(info_df$TME.subtype..k.4.))
x_sig=apply(info_df[,c("RNA_subtype","ce_subtype","Smoking","Stage","Sex","Ancestry","EGFR_status","KRAS_status","TP53_status")],2,x_test,y=factor(info_df$TME.subtype..k.4.))
annot_tme=HeatmapAnnotation(
  TME_subtype=info_df$TME.subtype..k.4.,
  RNA_subtype=info_df$RNA_subtype,
  CE_subtype=info_df$ce_subtype,
  smoking=info_df$Smoking,
  stage=info_df$Stage,
  sex=info_df$Sex,
  ancestry=info_df$Ancestry,
  EGFR_status=info_df$EGFR_status,
  KRAS_status=info_df$KRAS_status,
  TP53_status=info_df$TP53_status,
  purity=info_df$Purity,
  TMB=info_df$TMB,
  annotation_label = paste0(c("TME subtype","RNA subtype","CE subtype"
                              ,"Smoking","Stage","Sex","Ancestry","EGFR_status","KRAS_status","TP53_status","Purity","TMB"),
                            c("",x_sig,kw_sig)),
  annotation_legend_param = list(RNA_subtype=list(title = "RNA subtype"),
                                 CE_subtype=list(title = "Ecotyper"),
                                 smoking=list(title = "Smoking"),
                                 stage=list(title = "Stage"),
                                 sex=list(title = "Sex"),
                                 ancestry=list(title = "Ancestry"),
                                 EGFR_status=list(title = "EGFR"),
                                 KRAS_status=list(title = "KRAS"),
                                 TP53_status=list(title = "TP53"),
                                 purity=list(title = "Purity"),
                                 TMB= list(title = "TMB")
  ),
  col=list(TME_subtype=c("Initiating"="#80BA8A","Fibrotic"="#ED9F9B","Immune recruiting"="#6BB7CA","Immune activated"="#F4CEB4"),
           purity=colorRamp2(c(0,0.6,0.7,0.8,0.9,1), c("white","#FEEBE2","#F089AC","#DF639E","#CD3690","#C51B8A")),
           CE_subtype=ce_col,
           TMB=colorRamp2(c(0,5,10,15,20), c("#E0F3DB","#95CFBB","#7CBDC1","#5CACC7","#43A2CA")),
           sex=c("female"="#DC143C","male"="#4169E1"),
           stage=c("I"="#DCE9F4","II"="#ABD0F1","III"="#FF9797","IV"="#E56F5E","Normal"="#9EC4BE"),
           smoking=c("Smoker"="#DC143C","Non-smoker"="#FFF0F5"),
           ancestry=c("Asian"="#FFA500","European"="#FBE8D5"),
           KRAS_status=c("MT"="#DC143C","WT"="#FFF0F5"),
           EGFR_status=c("MT"="#DC143C","WT"="#FFF0F5"),
           TP53_status=c("MT"="#DC143C","WT"="#FFF0F5"),
           RNA_subtype=c("TRU-I"="#4169E1","PI"="#DC143C","PP"="#FF9797","TRU"="#FFA500")))
kanssandra_sig=apply(mat,1,kw_test,g=factor(info_df$TME.subtype..k.4.))
p_kanssandra=Heatmap(mat,name="TME cell percentage", top_annotation = annot_tme,col=col_fun,
           row_split=cell_type,
           column_split = info_df$TME.subtype..k.4.,
           row_labels = paste0(row.names(mat),kanssandra_sig),
           cluster_rows = F,
           row_title_rot = 0,
           cluster_columns = F,
           show_column_names = F)
ie_sig=apply(fge_score,1,kw_test,g=factor(info_df$TME.subtype..k.4.))
p_ie=Heatmap(fge_score,name="FGE z-score",col=col_fun_gsea,
             row_split = fge_func,
             row_labels = paste0(row.names(fge_score),ie_sig),
             cluster_rows = F,
             row_title_rot = 0,
             show_column_names = F,
             cluster_columns = F)
ce_sig=apply(abundance_mat,1,kw_test,g=factor(info_df$TME.subtype..k.4.))
p_ce=Heatmap(abundance_mat,name="Ecotyper cell abundance",col=col_fun_eco,
        cluster_rows = F,
        row_labels = paste0(row.names(abundance_mat),ce_sig),
        row_split = cell_type_mapping$Ecotype,
        row_title_rot = 0,
        show_column_names = F,
        cluster_columns = F)
imm_sig=apply(immune_score,1,kw_test,g=factor(info_df$TME.subtype..k.4.))
p_imm=Heatmap(immune_score,name="Immune score",col=colorRamp2(c(-1.2, 0, 1.2), c("#377eb8", "white", "#e4211c")),
              row_labels = paste0(row.names(immune_score),imm_sig),
              cluster_rows = T,
              row_title_rot = 0,
              show_column_names = F,
              cluster_columns = F)
p=p_kanssandra %v% p_ie %v% p_ce %v% p_imm
draw(p)
dev.off()
