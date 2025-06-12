library(Seurat)
library(ggplot2)
library(DDRTree)
library(parallel)
library(scales)
library(slingshot)
library(SingleCellExperiment)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

### Due to data control and large size, no example data for this script
med_scale=function(x){
  temp_med=median(x)
  temp_mad=mad(x)
  return((x-temp_med)/temp_mad)
}

load("malignant_cell.RData")
load("modules.RData")

modules=modules[c('Stress','Alveolar','Cycle','EMT','Hypoxia','Interferon')]
colors.module = c(hue_pal()(length(modules)))
names(colors.module) = names(modules)
modules=lapply(modules,function(x){return(x)})

combined_object@meta.data[which(combined_object@meta.data$cell_subtype=="AT2"),]$state="AT2"
combined_object$cell_id=rownames(combined_object@meta.data)

patient_list=data.frame(table(combined_object@meta.data$cell_subtype,combined_object@meta.data$orig.ident))
patient_list=patient_list[which(patient_list$Var1=="malignant"),]
patient_list=patient_list[which(patient_list$Freq>=100),]

srt_at=sample(colnames(subset(combined_object,subset=cell_subtype=="AT2")),50,replace=F)

obj_list=mclapply( as.character(patient_list$Var2),function(patient){
combined_object_sampled=subset(combined_object,subset=(orig.ident==patient & cell_subtype =="malignant")|(cell_id %in% srt_at))
cell_metadata=combined_object_sampled@meta.data

my.pca = combined_object_sampled@meta.data[,c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")]
my.pca = t(apply(my.pca,2,med_scale))

ncells_limit=100
ncells=ncol(combined_object_sampled)
dm = DDRTree(my.pca,dimensions=2)

tmp = t(data.matrix(data.frame(dm$Z)))
colnames(tmp)=paste0("Component_",as.character(1:2))
rownames(tmp)=colnames(my.pca)

tmp=tmp[rownames(combined_object_sampled@meta.data),]
combined_object_sampled[["DDRTree"]] <- CreateDimReducObject(embeddings = tmp, key="Component_", assay=DefaultAssay(combined_object_sampled))
combined_object_sampled@meta.data[,c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")]=t(my.pca)
return(combined_object_sampled)
},mc.cores=20,mc.set.seed = TRUE)
save(obj_list,file="trajectory_individual.RData")

## Lineage identification using Slingshot
# Filter patient with not enough malignant cells
obj_list_raw=obj_list
obj_list=list()
z=1
for(i in 1:length(obj_list_raw)){
    if(ncol(obj_list_raw[[i]])>=150){
      obj_list[[z]]=obj_list_raw[[i]]
      z=z+1
    }
  }
rm(obj_list_raw)
#Slingshot require clustering first, you can set number of clusters here for each patient
clu=rep(6,length(obj_list))
p_time_list=list()
pdf("lineage.pdf")
for(i in 1:length(obj_list)){
plot_list=list()
combined_object=obj_list[[i]]
sce=as.SingleCellExperiment(combined_object)
rd1=sce@int_colData$reducedDims$DDRTREE
cl2 <- kmeans(rd1, centers = clu[i],iter.max = 3000,nstart = 150)$cluster
colData(sce)$kmeans <- cl2
colData(sce)$GMM <- cl2

sta_clu=data.frame(table(sce$state,sce$GMM))
sta_clu=sta_clu[which(sta_clu$Var1=="AT2"),]
sta_clu=as.character(sta_clu[which.max(sta_clu$Freq),]$Var2)

sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'DDRTREE',start.clus=sta_clu,extend='n')

df=data.frame(reducedDims(sce)$DDRTREE,cluster=factor(sce$GMM))
curves <- slingCurves(sce, as.df = TRUE)
lineages = slingMST(sce, as.df = TRUE)
p=ggplot(data=df) +
  geom_point(aes(x = Component_1, y = Component_2,col = cluster),size=1) + 
  geom_point(data = lineages,aes(x=Component_1,y=Component_2), size = 4)+
  geom_path(data = lineages %>% arrange(Order),aes(x=Component_1,y=Component_2,group = Lineage))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )
p=set_panel_size(p,width  = unit(4, "inch"),
                 height = unit(4, "inch"))
plot_list[[1]]=p

w=data.frame(slingCurveWeights(sce))
w=scale(t(w),center=FALSE,scale=colSums(t(w)))
w=as.matrix(t(w))
p_time=as.matrix(sce@colData[,grep("slingPseudotime",colnames(sce@colData))])
p_time[is.na(p_time)]=0
p_aggr=w*p_time
p_aggr=rowSums(p_aggr)/ncol(p_aggr)
df$Pseudotime=p_aggr
df=cbind(df,sce@colData[,c("cell_subtype","Alveolar","Cycle","Interferon")])
df_text=curves[which(curves$Order==150),]
patient=gsub("_[A-Z]+-1","",rownames(df)[nrow(df)])

colors <- colorRampPalette(c("#1a2a6c","#b21f1f","#fdbb2d"))(100)
p= ggplot(data=df) +
  geom_point(aes(x = Component_1, y = Component_2,col = Pseudotime),size=1) + 
  scale_colour_gradientn(colours=colors)+
  geom_text(data = df_text,aes(x=Component_1,y=Component_2,label=Lineage),size=8)+
  geom_path(data = curves %>% arrange(Order),aes(x=Component_1,y=Component_2,group = Lineage))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )
p=set_panel_size(p,width  = unit(4, "inch"),
                 height = unit(4, "inch"))
plot_list[[2]]=p

df=df[rev(rownames(df)),]
p=ggplot()+geom_point(data=df,aes(x=Component_1,y=Component_2,col=cell_subtype),size=1)+labs(title = patient)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )
p=set_panel_size(p,width  = unit(4, "inch"),
                 height = unit(4, "inch"))
plot_list[[3]]=p

p_time=data.frame(as.matrix(sce@colData[,grep("slingPseudotime",colnames(sce@colData)),drop=FALSE]))
all_lineage=length(grep("slingPseudotime",colnames(sce@colData)))
p_time$patient=patient
p_time=cbind(p_time,sce@colData[,c("Alveolar","Cycle","Interferon","Stress","Hypoxia","EMT")])
p_time_list[[i]]=p_time

z=4
for(pse_time in 1:all_lineage){
  temp_df=p_time
  pse_time=temp_df[,pse_time]
  temp_df$pse_time=(pse_time-min(pse_time,na.rm=TRUE))/(max(pse_time,na.rm=TRUE)-min(pse_time,na.rm=TRUE))
  temp_df=pivot_longer(temp_df,-c(colnames(temp_df)[grep("slingPseudotime",colnames(temp_df))],"patient","pse_time"),names_to = "mp",values_to = "score")
  temp_df$mp=factor(temp_df$mp)
  p=ggplot()+
    geom_smooth(data=temp_df,aes(x=pse_time,y=score,color=mp),method = "lm",alpha=0.4)+
    scale_color_manual(values=c("Interferon"="#FF8C00",
                                "Cycle"="#DC143C",
                                "Alveolar"="#51B1B7",
                                "EMT"="#FFDAB9",
                                "Hypoxia"="#FFB6C1",
                                "Stress"="#B2DF8A"))+
    labs(x="Pseudo-time scaled",y="Score",title=paste0(patient,"_lineage_",as.character(z-3)))+
    xlim(c(0,1.25))+
    theme_classic()+
    theme(legend.position="none")
  p.smoothedmaxes <- 
    ggplot_build(p)$data[[1]] %>% 
    group_by(group) %>% 
    mutate( xmean = mean(x)) %>% 
    filter( x == max(x))
  p.smoothedmaxes$mp=sort(c("Alveolar","EMT","Hypoxia","Interferon","Cycle","Stress"))
  p=p +
    geom_text( data = p.smoothedmaxes, 
               mapping = aes(x = x, y = y, label = mp), 
               nudge_x = 0.1,
               col = p.smoothedmaxes$colour,
               inherit.aes = FALSE)
  p=set_panel_size(p,width  = unit(4, "inch"),
                   height = unit(4, "inch"))
  plot_list[[z]]=p
  z=z+1
}
p=ggpubr::ggarrange(plotlist = plot_list,ncol=3,nrow=4)
print(p)
}
dev.off()
saveRDS(p_time_list,file="p_time_list.rds")

## Lineage clustering
obj_list=readRDS("p_time_list.rds")

obj_list_coef=lapply(obj_list,function(x){
  p_time=x[,grep("sling",colnames(x)),drop=FALSE]
  scores=x[,c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")]
  cor_val=apply(p_time,2,function(y){
    cor(scores,y,method = "spearman",use="complete.obs")
  })
  rownames(cor_val)=c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")
  patient=unique(x$patient)
  colnames(cor_val)=paste(patient,colnames(cor_val),sep="_")
  return(cor_val)
})
cor_val=do.call(cbind,obj_list_coef)

dist_mat <- as.dist(1 - cor(cor_val,method = "pearson"))
hc <- hclust(dist_mat, method = "ward.D2")
branch=cutree(hc,2)
col_fun = colorRamp2(c(-0.7, 0, 0.7), c("#2F6294", "white", "#BF3A38"))
ht=Heatmap(cor_val,
           column_split = branch,
           clustering_distance_columns = "pearson",
           clustering_method_columns = "ward.D2",
           border = "white",
           col = col_fun,
           show_column_names = F
)
pdf("lineage_heatmap.pdf")
plot(ht)
dev.off()

cor_df_1=data.frame(t(apply(cor_val[,which(branch==1)],1,quantile,probs=c(0.25,0.5,0.75))))
cor_df_1$MP=rownames(cor_df_1)
cor_df_1$lineage="LG2"
colnames(cor_df_1)[1:3]=c("p25","p50","p75")
cor_df_1$MP=factor(cor_df_1$MP,levels = rev(c("Interferon","EMT","Hypoxia","Cycle","Stress","Alveolar")))

cor_df_2=data.frame(t(apply(cor_val[,which(branch==2)],1,quantile,probs=c(0.25,0.5,0.75))))
cor_df_2$MP=rownames(cor_df_2)
cor_df_2$lineage="LG1"
colnames(cor_df_2)[1:3]=c("p25","p50","p75")
cor_df_2$MP=factor(cor_df_2$MP,levels = rev(c("Interferon","EMT","Hypoxia","Cycle","Stress","Alveolar")))

p1=ggplot()+
  geom_point(data=cor_df_1,aes(y=as.numeric(MP)+0.1,x=p50,col=lineage), shape=15, size=3) +
  geom_linerange(data=cor_df_1,aes(y=as.numeric(MP)+0.1,xmin=p25, xmax=p75,col=lineage)) +
  geom_point(data=cor_df_2,aes(y=as.numeric(MP)-0.1,x=p50,col=lineage), shape=15, size=3) +
  geom_linerange(data=cor_df_2,aes(y=as.numeric(MP)-0.1,xmin=p25, xmax=p75,col=lineage)) +
  scale_color_manual(values=c("LG1"="#E07B54","LG2"="#F6C957"),guide="none")+
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_classic()+
  theme(
        axis.line.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.text.y= element_blank(),
        axis.title.x= element_blank(),
        axis.title.y= element_blank())
p2=ggplot(data=data.frame(MP=c("Interferon","EMT","Hypoxia","Cycle","Stress","Alveolar"),y=rev(c(1:6))))+
  geom_text(aes(x = 0,y=y, label = MP), hjust = 0)+
  theme_void()
p3=ggplot(data=data.frame(q_val=paste0(c(signif(p.adjust(wil[c("Interferon","EMT","Hypoxia","Cycle","Stress","Alveolar")],method ="fdr"),digits = 3)),
                                    c(" ***"," "," **"," ***"," "," *")),y=rev(c(1:6))))+
  geom_text(aes(x = 0,y=y, label = q_val), hjust = 0)+
  theme_void()

layout <- c(
  area(t = 0, l = 0, b = 30, r = 3), # left plot, starts at the top of the page (0) and goes 30 units down and 3 units to the right
  area(t = 1, l = 4, b = 30, r = 9), # middle plot starts a little lower (t=1) because there's no title. starts 1 unit right of the left plot (l=4, whereas left plot is r=3), goes to the bottom of the page (30 units), and 6 units further over from the left plot (r=9 whereas left plot is r=3)
  area(t = 0, l = 9, b = 30, r = 11))

pdf("lineage_forest.pdf")
p2+p1+p3+plot_layout(design = layout)
dev.off()
