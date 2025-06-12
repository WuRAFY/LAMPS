library(survminer)
library(survival)
library(compositions)
library(rms)
library(nnet)
library(igraph)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(glmnet)
library(tidyverse)

## Feature_selection
info_cox_full_model=readRDS("./input/cox_input_unfilter.rds")
info_cox_selection=info_cox_full_model[which(info_cox_full_model$overall_survival>0),]
info_eas=info_cox_selection[which(info_cox_selection$Ancestry=="Asian"),]
info_eur=info_cox_selection[which(info_cox_selection$Ancestry=="European"),]

X = model.matrix(~ Sex + Age + Stage + Smoker+
                    TMB+Purity+Ploidy+GII+
                    SI+Num.Clones+MATH+Num.driver+
                    EGFR+TP53+RBM10+KRAS+APC+SETD2+NF1+CTNND2+FAT4+CTNNA2+ANK1+
                    strata(Ancestry), data = info_cox_selection)
X = X[,-1]
coef_selected=lapply(1:500,function(x){
  cv_fit <- cv.glmnet(
    x = X,
    y = Surv(info_cox_selection[rownames(X),]$overall_survival,info_cox_selection[rownames(X),]$deceased),
    family = "cox",
    alpha = 1,
    nfolds = 10
  )
  final_coefs <- as.numeric(coef(cv_fit, s = cv_fit$lambda.min))
  final_coefs=ifelse(final_coefs==0,0,1)
  return(final_coefs)
})

cv_fit <- cv.glmnet(
  x = X,
  y = Surv(info_cox_selection[rownames(X),]$overall_survival,info_cox_selection[rownames(X),]$deceased),
  family = "cox",
  alpha = 1,
  nfolds = 10
)
final_coefs <- as.matrix(coef(cv_fit, s = cv_fit$lambda.min))
coef_selected=do.call(cbind,coef_selected)
rownames(coef_selected)=rownames(final_coefs)
coef_selected=rowSums(coef_selected)/500
coef_selected=coef_selected[coef_selected>=0.3]

## COX model regressed out ancestry(batch effect)
fit.coxph_selected=coxph(Surv(overall_survival, deceased) ~ Sex + Age + Stage + Smoker+
                           Purity+MATH+
                           TP53+APC+SETD2+FAT4+CTNNA2+
                           EM1+EM2+EM3+EM4+strata(Ancestry), 
                         data = info_cox_selection)
fit.coxph_combined=coxph(Surv(overall_survival, deceased) ~ Sex + Age + Stage + Smoker+
                           Purity+MATH+
                           TP53+APC+SETD2+FAT4+CTNNA2+
                         EM1+EM2+EM3+EM4+Ancestry, 
                         data = info_cox_selection)
cox_importance_combined = car::Anova(fit.coxph_combined,type=2)
cox_importance_combined = cox_importance_combined[-16,]
cox_importance_combined <- data.frame(importance=cox_importance_combined[,1]/sum(cox_importance_combined[,1])*100) %>% 
  `rownames<-`(row.names(cox_importance_combined))
cox_importance_combined$feature=rownames(cox_importance_combined)

pdf("COX.pdf",width=8,height = 8)
ggforest(fit.coxph_combined, data = info_cox_selection, main="Combined (Ancestry regressed)")
dev.off()

cor_p_val=lapply(info_cox_full_model[,3:ncol(info_cox_full_model)],function(x){
  p_val=lapply(info_cox_full_model[,3:ncol(info_cox_full_model)],function(y){
    class_x=ifelse(is.factor(x),1,0)
    class_y=ifelse(is.factor(y),1,0)
    if(class_x+class_y==0){
      test=cor.test(x,y,method="spearman",exact = FALSE)
      p=test$p.value
    }else if(class_x==1&class_y==0){
      test=kruskal.test(y,x)
      p=test$p.value
    }else if(class_x==0&class_y==1){
      test=kruskal.test(x,y)
      p=test$p.value}
    else{
      test=chisq.test(x,y)
      p=test$p.value
    }
    return(p)
  })
  p_val=unlist(p_val)
  return(p_val)
})
cor_p_val=do.call(rbind,cor_p_val)
diag(cor_p_val)=NA
cor_p_val_adj=matrix(p.adjust(as.vector(as.matrix(cor_p_val)), method='fdr'),ncol=30)
colnames(cor_p_val_adj)=colnames(cor_p_val)
rownames(cor_p_val_adj)=rownames(cor_p_val)

pdf("COX_correlation.pdf",width = 7.5,height = 6)
FDR_log=-log10(cor_p_val_adj)
col_fun = colorRamp2(c(0, 5, 10), c("#2B8DC4", "white", "#FD4C46"))
Heatmap(FDR_log,col=col_fun,name = "-log10(FDR)",
        cluster_rows = T,
        cluster_columns = T)
dev.off()

info_ord=c("Sex","Smoker","Age","Stage",
           "CTNNA2","TP53","SETD2","FAT4","APC",
           "Purity","MATH",
           "EM1","EM2","EM3","EM4")
hazard_ratio=exp(fit.coxph_selected$coefficients)
hazard_ratio=hazard_ratio[-c(3,4)]
names(hazard_ratio)=gsub("^Smoker","",names(hazard_ratio))
names(hazard_ratio)=gsub("European","",names(hazard_ratio))
names(hazard_ratio)=gsub("male","",names(hazard_ratio))
names(hazard_ratio)=gsub("MT","",names(hazard_ratio))
names(hazard_ratio)=gsub("IV","",names(hazard_ratio))
hazard_ratio=hazard_ratio[info_ord]
importance=cox_importance_combined$importance
names(importance)=cox_importance_combined$feature
importance=importance[info_ord]
importance[which(hazard_ratio<1)]=-importance[which(hazard_ratio<1)]
FDR_log=FDR_log[info_ord,info_ord]
diag(FDR_log)=0
FDR_log[FDR_log<=2]=0
p_val=summary(fit.coxph_selected)$coefficients[,5]
names(p_val)=gsub("^Smoker","",names(p_val))
names(p_val)=gsub("European","",names(p_val))
names(p_val)=gsub("male","",names(p_val))
names(p_val)=gsub("MT","",names(p_val))
names(p_val)=gsub("IV","",names(p_val))
p_val=p_val[info_ord]
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}
sizeTransform <- function(x) return(log(x+1)+1)
myHRplot <- function(coords, v=NULL, params) {
  #plot up for size larger than 1 and down otherwise.
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.border <-  params("vertex", "border")
  if (length(vertex.border) != 1 && !is.null(v)) {
    vertex.border <- vertex.border[v]
  }
  vertex.size <-  params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  #Set point size here (cex)
  for (i in 1:length(vertex.size)){
    if (vertex.size[i]>=0){
      points(x=coords[i,1], y=coords[i,2], col=vertex.border[i], bg=vertex.color[i], 
             pch=21, cex=sizeTransform(vertex.size[i]), lwd=4)
    } else points(x=coords[i,1], y=coords[i,2],col=vertex.border[i], bg=vertex.color[i], 
                  pch=23, cex=sizeTransform(abs(vertex.size[i])),lwd=4)
  }
}
add_shape("myHRplot", clip=shapes("circle")$clip,
          plot=myHRplot)
g <- graph_from_adjacency_matrix(FDR_log, mode = "undirected", weighted = T)
color.use = ifelse(hazard_ratio[colnames(FDR_log)]>1,"#EB9184","#80ACF9")
V(g)$size<-importance
V(g)$shape="myHRplot"
V(g)$border=ifelse(p_val<0.05,"grey30",NA)
E(g)$color <- scales::rescale(log(E(g)$weight),c(1,100)) %>% round %>% 
  colorRampPalette(c("grey90","black"))(100)[.] %>% adjustcolor(alpha.f=0.6)
label.locs <- radian.rescale(x=1:length(V(g)), direction=-1, start=0)

pdf("COX_circle.pdf",width = 12,height = 6)
coords <- layout_in_circle(g)
plot(g,layout=coords,
     #vertex.size=log(+1),
     vertex.label.dist=3,
     vertex.label.degree=label.locs,
     vertex.label.color="black",
     vertex.label.family="sans",
     vertex.label.font=2,
     vertex.color=color.use,
     #vertex.shape="myHRplot",
     #vertex.shape="circle",
     edge.width=scales::rescale(E(g)$weight,c(0.1,5)),
     edge.curved=.15)
leg.scale <- c(1,5,10,20,30)
leg <- legend('bottomright',legend=leg.scale,
              pt.cex=sizeTransform(leg.scale),col='white',
              pch=21, pt.bg='white',ncol=5,bty = "n",
              title="Univariate Cox model\nhazard ratio")
leg.x <- leg$text$x
leg.y <- leg$text$y+0.1
points(leg.x[1:5],leg.y[1:5],pch=21,lwd=2,bg="grey80",
       cex=sizeTransform(leg.scale[1:5]))
legend('right',legend=c("<0.05",">=0.05"),
       pt.cex=sizeTransform(2.2),col='black',pt.lwd=2,
       pch=21, pt.bg='grey80',ncol=2,bty = "n",
       title="Univariate Cox model\np-value")
max.q <- max(FDR_log)
min.q <- min(FDR_log)
q.legend <- -log(c(0.01,10^-10,10^-20,10^-40))
col.line <- scales::rescale(log(c(min.q,q.legend,max.q)),c(1,100)) %>% round %>% 
  colorRampPalette(c("grey94","grey40"))(100)[.] %>% .[2:5] 
lwd.line <- scales::rescale(c(min.q,q.legend,max.q),c(0.1,5)) %>% .[2:5]
legend('topright',legend=c("0.01",expression('10'^-10),expression('10'^-20),expression('10'^-40)),
       col=col.line,lty=1,lwd=lwd.line,box.lty=0,bty="n",
       title="Correlation\nq-value")
dev.off()

## Feature importance
library(car)
cox_importance=car::Anova(fit.coxph_combined,type=2)
cox_importance=cox_importance[1:15,]
cox_importance <- data.frame(importance=cox_importance[,1]/sum(cox_importance[,1])*100) %>% 
  `rownames<-`(row.names(cox_importance))
feature_factor=factor(c("Clinical","Clinical","Clinical","Clinical",
                        "Molecular","Molecular",
                        "Driver","Driver","Driver","Driver","Driver",
                        "Covel","Covel","Covel","Covel"))
names(feature_factor)=rownames(cox_importance)
col=c("Clinical"="#80BA8A",
      "Molecular"="#6BB7CA",
      "Driver"="#9CD1C8",
      "ITH"="#C6B3D3",
      "Covel"="#ED9F9B")
feature_col=col[as.character(feature_factor)]
cox_importance$class=feature_factor
cox_importance=cox_importance[order(cox_importance$importance),]
cox_importance$feature=rownames(cox_importance)
cox_importance$feature=factor(cox_importance$feature,levels = cox_importance$feature)
pdf("COX_feature_importance.pdf.pdf",width=4,height = 6)
ggplot(data=cox_importance)+
geom_col(aes(x=feature, y = importance, fill = class),width=0.9) +
  theme_classic()+
  scale_fill_manual(values=col)+
  ylim(c(0,40))+
  coord_flip()
dev.off()

cox_importance_collapse=split(cox_importance,f=cox_importance$class)
cox_importance_collapse=lapply(cox_importance_collapse,function(x){
  return(sum(x$importance))
})
cox_importance_collapse=unlist(cox_importance_collapse)
cox_importance_collapse=data.frame(feature=names(cox_importance_collapse),importance=cox_importance_collapse)
cox_importance_collapse=cox_importance_collapse[order(cox_importance_collapse$importance),]
cox_importance_collapse$feature=factor(cox_importance_collapse$feature,level=cox_importance_collapse$feature)
pdf("COX_feature_importance_collapse.pdf",width = 4,height = 3)
ggplot(data=cox_importance_collapse)+
geom_col(aes(x=feature, y = importance, fill = feature),width=0.9) +
  theme_classic()+
  scale_fill_manual(values=col)+
  ylim(c(0,50))+
  coord_flip()
dev.off()

## Harrel's C index
samplecindexWrap <- function(i,input.factors,data,p=0.7,...){
  outp <- sapply(
    i,function(x){
      #cat(i," ")
      set.seed(x)
      n <- nrow(data)
      if(p<=1) s <- sample(1:n,round(n*p))
      else s <- sample(1:n,p)
      train <- data[s,]
      test <- data[-s,]
      cindexWrap(input.factors,test,train,...)
      #cat("Done! ")
    }
  )
  return(outp)
}
cindexWrap <- function(input.factors,testdf,traindf,
                       cindex="Harrell",k=NULL,...){
  tr <- coxWrap(input.factors,traindf,...)
  pred <- predict(tr,testdf)
  #if (!is.null(k)) pred <- c(1:k)[discretize(pred,sort(pred)[seq(1,length(pred),length.out=k+1)[-c(1,k+1)]])]
  if (!is.null(k)) pred <- c(1:k)[discretize(order(pred),seq(1,length(pred),length.out=k+1)[-c(1,k+1)])]
  if (cindex=="Harrell") {
    library(survcomp)
    return(concordance.index(pred,testdf$overall_surviva,testdf$deceased,na.rm=T)$c.index)
  } else if (cindex=="Uno") {
    library(survAUC)
    train.rsp <- Surv(traindf$overall_surviva,traindf$deceased)
    test.rsp <- Surv(testdf$overall_surviva,testdf$deceased)
    return(UnoC(train.rsp,test.rsp,pred))
  }
}
coxWrap <- function(input.factors,survdf,summ=F,...){
  library(survival)
  f <- as.formula(paste0("Surv(overall_survival, deceased)~",paste(input.factors,collapse="+")))
  cox <- coxph(f,data=survdf)
  if(summ)print(summary(cox))
  return(cox)
}
c_index_factors=list(`Only_EM`=c("EM1","EM2","EM3","EM4","strata(Ancestry)"),
                     `Only_Clinical`=c("Age","Sex","Smoker","Stage","strata(Ancestry)"),
                     `EM&Clinical`=c("EM1","EM2","EM3","EM4","Age","Sex","Smoker","Stage","strata(Ancestry)"),
                     `Only_Molecular`=c("Purity","MATH","strata(Ancestry)"),
                     `Only_Driver`=c("TP53","APC","SETD2","FAT4","CTNNA2","strata(Ancestry)"),
                  `Basic_model`=c("Age","Sex","Smoker","Stage","strata(Ancestry)",
                                "TP53","APC","SETD2","FAT4","CTNNA2",
                                "Purity","MATH"),
                  `Basic_model&EM`=c("Age","Sex","Smoker","Stage","strata(Ancestry)",
                                  "TP53","APC","SETD2","FAT4","CTNNA2",
                                  "Purity","MATH",
                                  "EM1","EM2","EM3","EM4"))
c_index_list=lapply(c_index_factors,function(x){
  temp=samplecindexWrap(1:200,
                        x,
                        info_cox_selection)
  return(temp)
})
c_index=data.frame(model=factor(rep(names(c_index_list),each=200),level=names(c_index_list)),
                   index=unlist(c_index_list))

pdf("COX_C_index.pdf",width=5,height = 5)
ggplot(data=c_index,aes(x=model,y=index))+geom_boxplot(fill="#F9BEBB",lwd=0.6,outlier.size=0.5,
                                                       notch=T,width=0.5,position=position_dodge(0.7))+
  stat_compare_means(comparisons = list(c("Only_EM","Only_Clinical"),
                                        c("Only_Clinical","EM&Clinical"),
                                        c("Basic_model&EM","Basic_model")),
                     test = "wilcox.test",p.adjust.methods="fdr")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.5),
        panel.border = element_rect(fill=NA,colour = "black"))+
  geom_hline(yintercept=0.5,alpha=0.8,col="black",linetype=2,size=0.8)
dev.off()

## Pancer EM COX model
# Score patients with 6 MPs first as previously described
mp_score_list=readRDS("./input/pancancer_mp_score.rds")
# Run Kassandra first or directly download from Kassandra web
celltype=c("CD8_T_cells_PD1_high","CD8_T_cells_PD1_low","Endothelium",
                     "Fibroblasts","Macrophages_M1","Macrophages_M2","Monocytes","NK_cells",
                     "Neutrophils","Non_plasma_B_cells","Plasma_B_cells",
                     "T_helpers","Tregs")
EM=c("EM1","EM1","EM1","EM1","EM1",
     "EM2","EM2",
     "EM3","EM3","EM3","EM3","EM3",
     "EM4","EM4","EM4","EM4")
names(EM)=c("Macrophages_M1","Macrophages_M2","Neutrophils","NK_cells","Endothelium",
            "T_helpers","Non_plasma_B_cells",
            "CD8_T_cells_PD1_high","CD8_T_cells_PD1_low","Plasma_B_cells","Tregs","Interferon",
            "Cycle","EMT","Hypoxia","Fibroblasts")
pancancer_type=names(mp_score_list)
# Filter cancer type with over 95% survival rate (5 yr) and blood cancer 
pancancer_type=pancancer_type[-which(pancancer_type %in% c("DLBC","PRAD","TGCT","THCA","THYM","PCPG","CHOL","KICH","UCS"))]
pancancer_list=list()
cox_list=list()
plot_list=list()
cox_list_sep=list()
for(cancer_type in pancancer_type){
  mp_score=mp_score_list[[cancer_type]]
  tme=read.table(paste0("./input/Kassandra_TCGA/",cancer_type,"_BG_deconvolution.tsv"),sep="\t",header=T,row.names = 1)
  comm_sample=intersect(colnames(mp_score),colnames(tme))
  clinical=read.table(paste0("./input/Clinical_TCGA/clinical/",cancer_type,"_clinicalMatrix"),sep="\t",header=T,row.names = 1)
  rownames(clinical)=gsub("-",".",rownames(clinical))
  prognosis=read.table(paste0("./input/Clinical_TCGA/survival/",cancer_type,"_survival.txt"),sep="\t",header=T,row.names = 1)
  rownames(prognosis)=gsub("-",".",rownames(prognosis))
  clinical=clinical[comm_sample,]
  clinical=clinical[which(clinical$sample_type=="Primary Tumor"),]
  comm_sample=rownames(clinical)
  prognosis=prognosis[comm_sample,]
  mp_score=mp_score[,comm_sample]
  tme=tme[celltype,comm_sample]
combined_mat=t(rbind(mp_score,tme))
combined_rank=data.frame(apply(combined_mat,2,rank))
combined_rank$patient=rownames(combined_rank)
combined_rank=combined_rank %>% pivot_longer(-c("patient"),values_to = "rank",names_to = "cell_type")
combined_rank$EM=factor(EM[combined_rank$cell_type])
combined_rank=combined_rank[!is.na(combined_rank$EM),]
combined_rank=combined_rank %>% group_by(EM,patient) %>% summarise(sum(rank))
EM_total_rank=as.vector(table(EM)*nrow(combined_mat))
names(EM_total_rank)=paste0("EM",1:4)
combined_rank$total_rank=EM_total_rank[combined_rank$EM]
combined_rank$RSR=combined_rank$`sum(rank)`/combined_rank$total_rank
rsr_df=combined_rank %>% pivot_wider(id_cols = patient,names_from = EM,values_from = c(RSR))
rsr_df=data.frame(rsr_df)
rownames(rsr_df)=rsr_df$patient
rsr_df=scale(rsr_df[,c(2:5)],center = T)
if(cancer_type %in% c("CESC","OV","UCEC","UCS")){
  stage=gsub("Stage ","",clinical$clinical_stage)
  stage=gsub("A[0-9]*","",stage)
  stage=gsub("B[0-9]*","",stage)
}else if(cancer_type %in% c("GBM","LGG","PCPG","PRAD","SARC")){
  stage=NA 
}else if(cancer_type %in% c("THYM")){
  stage=gsub("Stage ","",clinical$masaoka_stage)
  stage=gsub("A[0-9]*","",stage)
  stage=gsub("B[0-9]*","",stage)
}else{
  stage=gsub("Stage ","",clinical$pathologic_stage)
  stage=gsub("A[0-9]*","",stage)
  stage=gsub("B[0-9]*","",stage)
}
df=data.frame(sample=clinical$X_INTEGRATION,
              tissue=clinical$sample_type,
              individual=clinical$X_PATIENT,
              age=clinical$age_at_initial_pathologic_diagnosis,
              sex=tolower(clinical$gender),
              stage=factor(stage,levels = c("I","II","III","IV")),
              overall_survival=prognosis$OS.time,
              deceased=prognosis$OS,
              combined_mat,
              rsr_df[comm_sample,c("EM1","EM2","EM3","EM4")])
pancancer_list[[cancer_type]]=df

quanti=summary(df$EM4)
df_surv=df[which(df$EM4>=quanti[5]|df$EM4<=quanti[2]),]
df_surv$EM4_content=ifelse(df_surv$EM4>=median(df$EM4),"High","Low")
fit <- survfit(Surv(overall_survival, deceased) ~ EM4_content, data = df_surv)
p=ggsurvplot(fit,
           data = df_surv,
           pval = T,title = cancer_type,
           risk.table = T)
plot_list[[cancer_type]]=p
cox_model=coxph(Surv(overall_survival, deceased) ~ EM1+EM2+EM3+EM4,
                data = df)
cox_list[[cancer_type]]=cox_model
cox_model_sep=coxph(Surv(overall_survival, deceased) ~ Cycle+EMT+Hypoxia+Fibroblasts,
                data = df)
cox_list_sep[[cancer_type]]=cox_model_sep
}

cox_coefficients=lapply(cox_list,function(x){
  temp=summary(x)
  temp=temp$coefficients[,2]
  return(temp)
})
cox_coefficients=do.call(rbind,cox_coefficients)
cox_coefficients=data.frame(cox_coefficients)
cox_coefficients$cancer_type=rownames(cox_coefficients)
cox_pval=lapply(cox_list,function(x){
  temp=summary(x)
  temp=temp$coefficients[,5]
  return(temp)
})
cox_pval=do.call(rbind,cox_pval)
cox_pval=data.frame(cox_pval)
cox_pval$cancer_type=rownames(cox_pval)
pdf("EM_prognosis_filtered.pdf",width = 30,height = 20)
arrange_ggsurvplots(plot_list,nrow = 4,ncol = 8)
dev.off()

cox_pval=cox_pval %>% pivot_longer(c("EM1","EM2","EM3","EM4"),values_to = "p_val",names_to = "EM")
cox_pval$FDR=p.adjust(cox_pval$p_val,method = "fdr")
cox_coefficients=cox_coefficients %>% pivot_longer(c("EM1","EM2","EM3","EM4"),values_to = "hazard",names_to = "EM")
df=data.frame(cancer_type=cox_pval$cancer_type,
              EM=cox_pval$EM,
              q=cox_pval$p_val,
              hazard=cox_coefficients$hazard)
df$q_cate=ifelse(df$q>0.05,"nosig","sig")
df$hazard_cate=ifelse(df$hazard>1,"pos","neg")
df$hazard_q=df$hazard_cate
df[df$q_cate=="nosig",]$hazard_q="nosig"
df$hazard=abs(df$hazard-1)

pdf("pancancer_EM.pdf",width = 4,height = 6)
ggplot(df, aes(x =EM , y = cancer_type)) +
  geom_point(aes(col=hazard_q,shape=hazard_cate),size=3.5) +
  scale_color_manual(values = c("pos"="#EB9184","neg"="#80ACF9","nosig" = "grey90")) +
  scale_shape_manual(values=c("pos"=16,"neg"=18))+
  labs(
    title = "Pancancer EM prognosis",
    x = "EM",
    y = "Cancer Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
df_sum=data.frame(table(df$hazard_q,df$EM))
df_sum=df_sum[which(df_sum$Var1 != "nosig"),]
ggplot(df_sum)+geom_col(aes(x=Var2,y=Freq,fill=Var1))+
  scale_fill_manual(values = c("pos"="#EB9184","neg"="#80ACF9","nosig" = "grey90"))+
  theme_classic()
dev.off()