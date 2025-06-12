library(igraph)
library(clusterCrit)
library(rms)
library(ComplexHeatmap)
library(nnet)
library(circlize)
library(tidyverse)
library(circlize)
library(ggradar)
library(DDRTree)
library(scales)
library(ggpubr)

med_scale=function(x){
  temp_med=median(x)
  temp_mad=mad(x)
  return((x-temp_med)/temp_mad)
}

bulk_info=read.table("input/bulk_info.txt",sep="\t",stringsAsFactors = T,header=T)
TME_mat=read.table("./input/TME_composition.txt",sep="\t",header=T,row.names = 1)
mp_mat=read.table("./input/MP_score.txt",sep="\t",header=T,row.names = 1)
combined_mat=cbind(TME_mat,mp_mat)
combined_mat=log(combined_mat+1)
combined_mat=scale(combined_mat)

## Feature correlation
cor_mat=cor(combined_mat,method = "pearson")
col_fun=colorRamp2(c(-1,-0.3, 0, 0.3,1), c("#2B8DC4","#2B8DC4","white","#FD4C46","#FD4C46"))
pdf("EM_heatmap.pdf",width=7,height = 6)
dev.off()

#Feature clustering
cor_mat=cor_mat+1
g=graph_from_adjacency_matrix(cor_mat, mode = "undirected", weighted=TRUE, diag = F)
g_cluster=cluster_leiden(g,resolution_parameter=1.2,n_iterations=500)
res=g_cluster$membership
colnames(cor_mat)[which(res==5)]
pdf("EM_clustering.pdf",width=6,height = 6)
par(mfrow=c(2,2))
for(k in c(1:3,5)){
cor_partial=cor_mat[which(res==k),which(res==k)]
temp_g=graph_from_adjacency_matrix(cor_partial, mode = "undirected", weighted=TRUE, diag = F)
V(temp_g)$size=30
V(temp_g)$color<-"grey"
V(temp_g)$label.color <- "black"
V(temp_g)$label.cex<-0.8
E(temp_g)$width<- 5*abs(E(temp_g)$weight)
p=plot.igraph(temp_g,
            edge.curved=0,
            vertex.shape='circle',
            layout=layout.circle)
}
dev.off()

## RSR scoring
combined_rank=data.frame(apply(combined_mat,2,rank))
rownames(combined_rank)=rownames(combined_mat)
colnames(combined_rank)=colnames(combined_mat)
combined_rank$patient=rownames(combined_rank)
em=res
em=paste0("EM",as.character(em))
names(em)=colnames(cor_mat)
combined_rank=combined_rank %>% pivot_longer(-c("patient"),values_to = "rank",names_to = "cell_type")
combined_rank$EM=factor(em[combined_rank$cell_type])
combined_rank=combined_rank %>% group_by(EM,patient) %>% summarise(sum(rank))
em_total_rank=as.vector(table(em)*nrow(combined_mat))
names(em_total_rank)=paste0("EM",1:6)
combined_rank$total_rank=em_total_rank[combined_rank$EM]
combined_rank$RSR=combined_rank$`sum(rank)`/combined_rank$total_rank
combined_rank=merge(combined_rank,bulk_info,all.x = TRUE,by.x = "patient",by.y = "Patient")
combined_rank=combined_rank[which(combined_rank$EM %in% c("EM1","EM2","EM3","EM5")),]
combined_rank$EM=factor(combined_rank$EM,levels=c("EM2","EM5","EM1","EM3"),
                        labels = c("EM1","EM2","EM3","EM4"))
pdf("EM_RSR.pdf",width=10,height = 3)
ggplot(data=combined_rank,aes(x=Ancestry,y=RSR))+geom_violin(aes(x=Ancestry,y=RSR,fill=Ancestry))+
  geom_boxplot(width=0.3,outlier.size = 0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5),
        panel.border = element_rect(fill=NA,colour = "black"))+
  labs(x=NULL)+
  scale_fill_manual(values = c("Asian"="#FFA500","European"="#FBE8D5"))+
  stat_compare_means(comparisons = list(c("Asian","European")),map_signif_level=TRUE,
                     test = "wilcox.test",p.adjust.methods="fdr")+
  facet_wrap(vars(EM),nrow=1)
dev.off()

## Purity regression
rsr_cor=combined_rank %>% group_by(EM) %>% summarise(correlation = cor(RSR,Purity))
rsr_df=combined_rank %>% pivot_wider(id_cols = patient,names_from = EM,values_from = c(RSR,Purity))
rsr_df=rsr_df[,1:6]
colnames(rsr_df)=c("patient","EM3","EM1","EM4","EM2","purity")
model <- lm(purity ~ EM1 + EM2 + EM3 + EM4, data = rsr_df)
coefficients <- coef(summary(model))[-1, "Estimate"]  # Remove intercept
p_values <- coef(summary(model))[-1, "Pr(>|t|)"]
plot_data <- data.frame(
  Variable = names(coefficients),
  Coefficient = coefficients,
  P_Value = p_values
)
pdf("EM_purity_regression.pdf",width=5,height = 3)
ggplot(plot_data, aes(x = Variable, y = Coefficient, fill = Coefficient > 0)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("p=%.3f", P_Value)), 
            vjust = ifelse(plot_data$Coefficient > 0, -1.2, 1.5), size = 4) +
  scale_fill_manual(
    values = c("TRUE" = "lightpink1", "FALSE" = "lightblue"),
    labels = c("Negative Coefficient", "Positive Coefficient")
  ) +
  labs(
    title = "Purity regression",
    x = "Variable",
    y = "Coefficient",
    fill = "Coefficient Sign"
  ) +
  theme_classic2()
dev.off()
