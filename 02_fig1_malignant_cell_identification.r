library(infercnv)

### Due to data control and large size, no example data for this script

## Cell distribution
load("sub_cluster.RData")

sub_object=type_list[["Epithelial cell"]]
sub_object= NormalizeData(sub_object) %>% FindVariableFeatures(selection.method="vst",nfeatures=2000) %>% ScaleData() %>%
   RunPCA(verbose=FALSE,npcs=30%>%
   RunUMAP(reduction="pca",dims=1:30)%>%
   FindNeighbors(reduction="pca",dims=1:30) %>% 
   FindClusters(resolution=0.5)
img=DimPlot(sub_object,group.by="orig.ident",pt.size=0.1)+labs(title="Epithelial cell")
png(paste0("Epithelial cell_UMAP.png"),width=800,height=800)
plot(img2)
dev.off()

## Cell markers
cell_markers=list(
`Epithelial cell`=list(
      `Epithelial`=c("EPCAM"),
      `AT1`=c("AGER","CAV1","PDPN"),
      `AT2`=c("SFTPB","SFTPC","SFTPD"),
      `Ciliated`=c("FOXJ1","TPPP3","CAPS"),
      `Club`=c("SCGB3A2","SCGB1A1","SCGB3A1")
)
)
for(m in names(cell_markers[["Epithelial cell"]])){
    png("Epithelial cell_canonical_markers.png",width=1500,height=1500)
    img1=FeaturePlot(sub_object, features = cell_markers[["Epithelial cell"]][[m]],pt.size=0.1)+labs(caption=m)
    plot(img1)
    dev.off()
    png("Epithelial cell_canonical_violin.png",width=1500,height=1500)
    img2=VlnPlot(sub_object,features = cell_markers[["Epithelial cell"]][[m]])+labs(caption=m)
    plot(img2)
    dev.off()
}
png("Epithelial cell_dot.png",width=2000,height=1000)
dot=DotPlot(sub_object,cols=c("RdYlBu"),group.by="seurat_clusters",cluster.idents=TRUE,features=cell_markers[["Epithelial cell"]])+RotatedAxis()+scale_size(range = c(5, 20))
plot(dot)
dev.off()
png("Epithelial cell_feature_count.png",width=2000,height=1000)
feature=FeaturePlot(sub_object, features = c("nCount_RNA","nFeature_RNA"),pt.size=0.1)
plot(feature)
dev.off()
png("Epithelial cell_feature_violin.png",width=2000,height=1000)
feature=VlnPlot(sub_object, features = c("nCount_RNA","nFeature_RNA"),pt.size=0.1)
plot(feature)
dev.off()

## GII calculation
#Before run this step, run infercnv first
seg=read.table("17_HMM_predHMMi6.leiden.hmm_mode-subclusters.pred_cnv_regions.dat",header=T,sep="\t")
infercnv_obj = readRDS("run.final.infercnv_obj")
seg=split(seg,f=seg$cell_group_name)
GII=lapply(seg,function(segment){
           gii=sum((segment$end-segment$start)*abs(segment$state-3))/3099734149
           return(gii)
           })
cell_indice=unlist(infercnv_obj@tumor_subclusters$subclusters,recursive=F)
cell_indice=lapply(names(cell_indice),function(x){
  cell=names(cell_indice[[x]])
  if(x %in% names(GII)){
  s=GII[[x]]}else{
  s=0}
  df=data.frame(cell,s)
  return(df)
})
cell_indice=do.call(rbind,cell_indice)
write.table(cell_indice,"genome_score_seg.txt",append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
