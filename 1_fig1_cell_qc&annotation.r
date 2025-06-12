library(Seurat)
library(DoubletFinder)
library(dplyr)
library(ggplot2)
library(harmony)
library(future)

### Due to data control and large size, no example data for this script

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 10000 * 1024^2)
##1.Read cellranger output and perferom quality control
args=commandArgs(T)
path=args[1]
sample=args[2]
doublet_rate=as.numeric(args[3])
# Read cell ranger outputs
mat=Read10X(path)
# Filter cell debris
seu_object=CreateSeuratObject(counts=mat,min.genes=200,project=sample)
# UMI-gene plot and mito-RNA percentage
seu_object[["percent.mt"]]=PercentageFeatureSet(seu_object,pattern="^MT-")
pdf(paste0(sample,"_qc.pdf"))
VlnPlot(seu_object,features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol=3)
FeatureScatter(seu_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seu_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
seu_object=subset(seu_object,subset=nFeature_RNA>200&percent.mt<20)
# DF does not support SCTransform v2 for now
seu_object=NormalizeData(seu_object) %>%
    FindVariableFeatures(selection.method="vst",nfeatures=2000)%>%
    ScaleData()%>%
    RunPCA(npcs = 30,verbose=FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30)
# DoubletFinder
# When PC=30 neuroendocrine cells will be filtered but get a clear map.
sweep.res.list_LUAD=paramSweep_v3(seu_object,PCs=1:30,num.cores=12)
sweep.stats_LUAD=summarizeSweep(sweep.res.list_LUAD,GT=FALSE)
bcmvn_LUAD=find.pK(sweep.stats_LUAD)
pK=as.character(bcmvn_LUAD$pK[which.max(bcmvn_LUAD$BCmetric)])
pK=as.numeric(pK)
nExp_poi=round(doublet_rate*nrow(seu_object@meta.data))
seu_object=doubletFinder_v3(seu_object,PCs=1:30,pN=0.25,pK=pK,nExp=nExp_poi,reuse.pANN=FALSE)
names(seu_object@meta.data)[grep("DF*",names(seu_object@meta.data))]="DF.classification"
seu_object=subset(seu_object,subset=DF.classification=="Singlet")
dev.off()
saveRDS(seu_object,file=paste0(sample,".rds"))

##2.Merge samples from same dataset and perform dimension reduction
args=commandArgs(T)
# build a sample list text file for each dataset before this step
sample=read.table(args[1],sep="\t")
sample=sample[,1]
obj_list=c()
for(sam in sample){
    temp_obj=(paste0(sam,".rds"))
    obj_list=c(obj_list,temp_obj)
}
combined_object=merge(obj_list[[1]],y=obj_list[-1],add.cell.ids=sample)
combined_object=NormalizeData(combined_object) %>% FindVariableFeatures() %>% ScaleData(verbose=FALSE)
combined_object=RunPCA(combined_object,npcs=30,verbose = FALSE)
combined_object = RunHarmony(combined_object,assay.use="RNA",dims.use=1:30,group.by.vars="orig.ident")
combined_object=RunUMAP(combined_object,reduction="harmony",dims=1:30)%>%
   FindNeighbors(reduction="harmony",dims=1:30)%>%
   FindClusters(resolution=0.8)
clu=levels(combined_object@meta.data$seurat_clusters)
save(combined_object,file="cluster.RData")

## Markers for cell type annotation
DimPlot(combined_object,reduction="umap",label=TRUE,repel=TRUE,pt.size=0.1)
DefaultAssay(combined_object) = "RNA"
# DEG
for(i in clu){
tryCatch({
markers = FindConservedMarkers(combined_object,ident.1=i,assay="RNA",grouping.var="orig.ident",verbose=FALSE,min.cells.group=3)
top_markers = rownames(markers)[1:9]
if(length(which(!is.na(top_markers)))>0){
img1=FeaturePlot(combined_object, features = top_markers,pt.size=0.1)+labs(caption=paste0("Cluster ",i))
png(paste0("markers_cluster_",as.character(i),".png"),width=1000,height=1000)
plot(img1)
dev.off()
png(paste0("violin_cluster_",as.character(i),".png"),width=1000,height=1000)
img2=VlnPlot(combined_object,features=top_markers,pt.size=0)+labs(caption=paste0("Cluster ",i))
plot(img2)
dev.off()}},error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
# Gene markers
png(paste0("umap_individual.png"),width=1000,height=1000)
features=c("EPCAM","CLDN5","COL1A2","PTPRC","CD3E","CD4","CD8A","CD19","MZB1","MS4A2","G0S2","APOE","FCN1","LYZ","AIF1","MKI67")
DimPlot(combined_object,reduction="umap",group.by="orig.ident",pt.size=0.1)
dev.off()
png(paste0("umap.png"),width=1000,height=1000)
DimPlot(combined_object,reduction="umap",label=TRUE,repel=TRUE,pt.size=0.1)
dev.off()
png(paste0("canonical_markers.png"),width=1500,height=1500)
FeaturePlot(combined_object, features = features,pt.size=0.1)
dev.off()
png(paste0("canonical_markers_violin.png"),width=1500,height=1500)
VlnPlot(combined_object,features = features,pt.size=0.1)
dev.off()
png(paste0("canonical_markers_dot.png"),width=1000,height=1000)
DotPlot(combined_object,cols=c("RdYlBu"),group.by="seurat_clusters",cluster.idents=TRUE,features=features)+RotatedAxis()
dev.off()

## Cell subtype annotation
# You should assign cell type to each cell before this step
sub_type=c("Epithelial cell","Lymphoid","Myeloid","Endothelial cell","Stromal cell","Proliferating")
dims=c(30,30,30,30,30,30)
res=c(0.7,0.9,0.7,0.5,0.5,0.5)
names(dims)=sub_type
names(res)=sub_type

type_list=list()
combined_object=DietSeurat(combined_object,counts=TRUE,data=TRUE,scale.data=FALSE)
for(s in sub_type){
sub_object=subset(combined_object,subset=cell_type==s)
sub_object=NormalizeData(sub_object) %>% FindVariableFeatures(selection.method="vst",nfeatures=2000) %>% ScaleData() %>%
   RunPCA(verbose=FALSE,npcs=dims[s])%>%
   RunHarmony(assay.use="RNA",dims.use=1:dims[s],group.by.vars="orig.ident")%>%
   RunUMAP(reduction="harmony",dims=1:dims[s])%>%
   FindNeighbors(reduction="harmony",dims=1:dims[s]) %>% 
   FindClusters(resolution=res[s])
img1=DimPlot(sub_object,group.by="orig.ident",pt.size=0.1)+labs(title=s)
img2=DimPlot(sub_object,group.by="ident",label=TRUE,pt.size=0.1)+labs(title=s)
png(paste0(s,"_umap_individual_res",as.character(res[s]),".png"),width=800,height=800)
plot(img1)
dev.off()
png(paste0(s,"_umap_res",as.character(res[s]),".png"),width=800,height=800)
plot(img2)

DefaultAssay(sub_object) = "RNA"
for(i in levels(sub_object@meta.data$seurat_clusters)){
markers= FindConservedMarkers(sub_object,ident.1=i,assay="RNA",grouping.var="orig.ident",verbose=FALSE)
top_markers = rownames(markers)[1:9]
img1=FeaturePlot(sub_object, features = top_markers,pt.size=0.1)+labs(caption=paste0("Cluster ",i))
png(paste0(s,"_markers_res_",as.character(res[s]),"_cluster_",as.character(i),".png"),width=1500,height=1500)
plot(img1)
dev.off()
png(paste0(s,"_violin_res_",as.character(res[s]),"_cluster_",as.character(i),".png"),width=1500,height=1500)
img2=VlnPlot(sub_object, features= top_markers,pt.size=0.1)+labs(caption=paste0("Cluster ",i))
plot(img2)
dev.off()
}
dev.off()
type_list[[s]]=sub_object
}
save(type_list,file="sub_cluster.RData")
