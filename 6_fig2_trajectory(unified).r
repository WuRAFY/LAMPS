library(Seurat)
library(ggplot2)
library(DDRTree)
library(parallel)
library(scales)
library(monocle3)
library(igraph)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## Construct the unified trajectory
med_scale=function(x){
  temp_med=median(x)
  temp_mad=mad(x)
  return((x-temp_med)/temp_mad)
}

combined_object=readRDS("./input/trajectory.rds")
# The output file trajectory.rds can be used as example data
#combined_object=readRDS("trajectory.rds")
load("./input/modules.RData")
colors.module = c(hue_pal()(length(modules)))
names(colors.module) = names(modules)
modules=lapply(modules,function(x){return(x)})
combined_object$cell_id=rownames(combined_object@meta.data)
# Down sample cells to avoid over representation, this step should be skipped when using example data
AT2_id=combined_object@meta.data[which(combined_object@meta.data$cell_subtype=="AT2 cell"),]
AT2_id=split(AT2_id,f=factor(AT2_id$orig.ident))
malignant_id=combined_object@meta.data[which(combined_object@meta.data$cell_subtype=="Malignant cell"),]
malignant_id=split(malignant_id,f=factor(malignant_id$orig.ident))
AT2_id_sampled=lapply(AT2_id,function(x){
  if(nrow(x)>30){
    temp_id=rownames(x)
    temp_id=sample(temp_id,30,replace=F)
    return(temp_id)
  }else{
    return(rownames(x))
  }
})
AT2_id_sampled=unlist(AT2_id_sampled)
malignant_id_sampled=lapply(malignant_id,function(x){
  if(nrow(x)>500){
    temp_id=rownames(x)
    temp_id=sample(temp_id,500,replace=F)
    return(temp_id)
  }else{
    return(rownames(x))
  }
})
malignant_id_sampled=unlist(malignant_id_sampled)
combined_id=c(AT2_id_sampled,malignant_id_sampled)
combined_object_sampled=subset(combined_object,subset=cell_id %in% combined_id)
my.pca = combined_object_sampled@meta.data[,c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")]
my.pca = t(apply(my.pca,2,med_scale))
dm = DDRTree(my.pca,dimensions=2,ncenter=700)
tmp = t(data.matrix(data.frame(dm$Z)))
colnames(tmp)=paste0("Componenet_",as.character(1:2))
rownames(tmp)=colnames(my.pca)
tmp=tmp[rownames(combined_object_sampled@meta.data),]
combined_object_sampled[["DDRTree"]] <- CreateDimReducObject(embeddings = tmp, key="Component_", assay=DefaultAssay(combined_object_sampled))
combined_object_sampled@meta.data[,c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")]=t(my.pca)
#save(combined_object_sampled,file="trajectory.rds")

## Pseudotime ordering using Monocle3
get_earliest_principal_node <- function(cds){
  cell_ids <- which(colData(cds)[, "state"] == "AT2")
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

#Pseudtotime ordering with AT2 cell as root
expression_matrix=GetAssayData(combined_object_sampled,slot="counts",assay="RNA")
emb=combined_object_sampled[["DDRTree"]]@cell.embeddings
cell_metadata=combined_object_sampled@meta.data[colnames(expression_matrix),c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon","orig.ident","nCount_RNA","nFeature_RNA","percent.mt","cell_subtype","state")]
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata)
cds <- preprocess_cds(cds)
cds@int_colData@listData$reducedDims$UMAP <- as.matrix(emb[,c(1,2)])
cds@clusters$UMAP_so$clusters <- cell_metadata[,"state"]
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds,use_partition=FALSE)
cds <- order_cells(cds,root_pr_nodes=get_earliest_principal_node(cds))

# Cut trajectory into 3 branches (edit from Monocle3 code)
branch_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  branch_points <- which(igraph::degree(g) > 2)
  branch_points = branch_points[branch_points %in% root_nodes(cds, reduction_method) == FALSE]
  return(branch_points)
}
leaf_nodes <- function(cds,reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  leaves <- which(igraph::degree(g) == 1)
  leaves = leaves[leaves %in% root_nodes(cds, reduction_method) == FALSE]
  return(leaves)
}
root_nodes <- function(cds, reduction_method="UMAP"){
  g = principal_graph(cds)[[reduction_method]]
  root_pr_nodes <- which(names(igraph::V(g)) %in%
                           cds@principal_graph_aux[[reduction_method]]$root_pr_nodes)
  names(root_pr_nodes) <-
    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes
  return(root_pr_nodes)
}
buildBranchCellDataSet <- function(cds,
                                   progenitor_method = c('sequential_split', 'duplicate'), 
                                   branch_states = NULL, 
                                   branch_point = 1,
                                   branch_labels = NULL, 
                                   stretch = TRUE)
{
  if (!is.null(branch_labels) & !is.null(branch_states)) {
    if(length(branch_labels) != length(branch_states))
      stop("length of branch_labels doesn't match with that of branch_states")
    branch_map <- setNames(branch_labels, as.character(branch_states))
  }
  pr_graph_cell_proj_mst <- principal_graph(cds)[['UMAP']]
  root_cell_point_in_Y <- cds@principal_graph_aux$UMAP$root_pr_nodes
  root_cell <- names(which(degree(pr_graph_cell_proj_mst, v = root_cell_point_in_Y, mode = "all")==1, useNames = T))[1]
  paths_to_root <- list()
  mst_branch_nodes <- branch_nodes(cds, "UMAP")
  branch_cell <- names(mst_branch_nodes[branch_point])
  mst_no_branch_point <- pr_graph_cell_proj_mst - V(pr_graph_cell_proj_mst)[branch_cell]
  path_to_ancestor <- shortest_paths(pr_graph_cell_proj_mst, branch_cell, root_cell)
  path_to_ancestor <- names(unlist(path_to_ancestor$vpath))
  for (backbone_nei in V(pr_graph_cell_proj_mst)[suppressWarnings(nei(branch_cell))]$name){
    descendents <- bfs(mst_no_branch_point, V(mst_no_branch_point)[backbone_nei], unreachable=FALSE)
    descendents <- descendents$order[!is.na(descendents$order)]
    descendents <- V(mst_no_branch_point)[descendents]$name
    path_to_root <- unique(c(path_to_ancestor, branch_cell, descendents))
    closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
    path_to_root <- row.names(closest_vertex)[which(V(pr_graph_cell_proj_mst)[closest_vertex]$name %in% path_to_root)]
    closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
    path_to_root <- intersect(path_to_root, colnames(cds))
    paths_to_root[[backbone_nei]] <- path_to_root
  }
  all_cells_in_subset <- c()
  names(paths_to_root) <- branch_labels
  for (path_to_ancestor in paths_to_root){
    if (length(path_to_ancestor) == 0){
      stop("Error: common ancestors between selected State values on path to root State")
    }
    all_cells_in_subset <- c(all_cells_in_subset, path_to_ancestor)
  }
  all_cells_in_subset <- unique(all_cells_in_subset)
  common_ancestor_cells <- intersect(paths_to_root[[1]], paths_to_root[[2]])
  cds <- cds[, row.names(pData(cds[,all_cells_in_subset]))] #or just union(ancestor_cells, branch_cells)
  Pseudotime <- pseudotime(cds)
  pData <- pData(cds)
  pData$Pseudotime=Pseudotime
  if(stretch) {
    max_pseudotime <- -1
    for (path_to_ancestor in paths_to_root){
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime)  
      if (max_pseudotime < max_pseudotime_on_path){
        max_pseudotime <- max_pseudotime_on_path
      }
    }
    branch_pseudotime <- max(pData[common_ancestor_cells,]$Pseudotime)
    for (path_to_ancestor in paths_to_root){
      max_pseudotime_on_path <- max(pData[path_to_ancestor,]$Pseudotime) 
      path_scaling_factor <-(max_pseudotime - branch_pseudotime) / (max_pseudotime_on_path - branch_pseudotime)
      if (is.finite(path_scaling_factor)){
        branch_cells <- setdiff(path_to_ancestor, common_ancestor_cells)
        pData[branch_cells,]$Pseudotime <- ((pData[branch_cells,]$Pseudotime - branch_pseudotime) * path_scaling_factor + branch_pseudotime)
      }
    }
    pData$Pseudotime <- 100 * pData$Pseudotime / max_pseudotime
  }
  pData$original_cell_id <- row.names(pData)
  branch=lapply(names(paths_to_root),function(x){
    df=data.frame(cell_id=paths_to_root[[x]],Branch=x)
    rownames(df)=df$cell_id
    return(df)
  })
  branch=do.call(rbind,branch)
  pData(cds)$Branch=branch[rownames(pData(cds)),]$Branch
  pData(cds)[common_ancestor_cells, "Branch"] <- names(paths_to_root)[1]
  pData(cds)$Pseudotime=pseudotime(cds)
  return (cds)
}
cds_branch=buildBranchCellDataSet(cds,
                                  progenitor_method = c('duplicate'), 
                                  branch_states = NULL, 
                                  branch_point = 9, # Edit the branching point according to your trajectory
                                  branch_labels = c("Precursor","IFN-responsive", "Proliferative"), 
                                  stretch = TRUE)
combined_object_sampled$pseudotime=pData(cds_branch)[combined_object_sampled$cell_id,]$Pseudotime
combined_object_sampled$branch=pData(cds_branch)[combined_object_sampled$cell_id,]$Branch
saveRDS(combined_object_sampled,file="trajectory.rds")

## Plot
cell_metadata=combined_object_sampled@meta.data
cell_metadata=cbind(cell_metadata,combined_object_sampled[["DDRTree"]]@cell.embeddings)
cell_metadata=cell_metadata[rev(rownames(cell_metadata)),]
dp_mst <- cds@principal_graph[["UMAP"]]
ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>%
  as.data.frame() %>%
  dplyr::mutate(sample_name = rownames(.),
                sample_state = rownames(.))
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  dplyr::select(source = "from", target = "to") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(source="sample_name",
                                   source_prin_graph_dim_1="Component_1",
                                   source_prin_graph_dim_2="Component_2"),
                   by = "source") %>%
  dplyr::left_join(ica_space_df %>%
                     dplyr::select(target="sample_name",
                                   target_prin_graph_dim_1="Component_1",
                                   target_prin_graph_dim_2="Component_2"),
                   by = "target")

pdf("trajectory_split.pdf",width=30,height = 30)
ggplot()+geom_segment(aes_string(x="source_prin_graph_dim_1",
                                 y="source_prin_graph_dim_2",
                                 xend="target_prin_graph_dim_1",
                                 yend="target_prin_graph_dim_2"),
                      size=0.2,
                      color="black",
                      linetype="solid",
                      na.rm=TRUE,
                      data=edge_df)+
  geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=cell_subtype),size=0.3)+
  scale_color_manual(values = c("AT2 cell"="#54bbd5","Malignant cell"="#f17f73"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )+facet_wrap(vars(orig.ident_m),nrow=13)+theme(legend.position="none")
dev.off()

cell_metadata_order=cell_metadata[sample(rownames(cell_metadata)),]
h1=ggplot()+geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=Alveolar),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Alveolar")+
  scale_color_gradientn(colours = c("grey", "grey","#51AAAE", "#51AAAE"), 
                        values = rescale(x = c(min(cell_metadata$Alveolar), 0, 0.5, max(cell_metadata$Alveolar)),
                                         from = c(min(cell_metadata$Alveolar), max(cell_metadata$Alveolar))))
h2=ggplot()+geom_point(data=cell_metadata[order(cell_metadata$Stress),],aes(x=Component_1,y=Component_2,col=Stress),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Stress")+
  scale_color_gradientn(colours = c("grey", "grey", "#ADD487","#ADD487"), 
                        values = rescale(x = c(min(cell_metadata$Stress), 0, 1, max(cell_metadata$Stress)),
                                         from = c(min(cell_metadata$Stress), max(cell_metadata$Stress))))
h3=ggplot()+geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=Cycle),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Cycle")+
  scale_color_gradientn(colours = c("grey", "grey", "#DC143C","#DC143C"), 
                        values = rescale(x = c(min(cell_metadata$Cycle), 0, 2, max(cell_metadata$Cycle)),
                                         from = c(min(cell_metadata$Cycle), max(cell_metadata$Cycle))))
h4=ggplot()+geom_point(data=cell_metadata[order(cell_metadata$EMT),],aes(x=Component_1,y=Component_2,col=EMT),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "EMT")+
  scale_color_gradientn(colours = c("grey", "grey","#FFDAB9","#FFDAB9"), 
                        values = rescale(x = c(min(cell_metadata$EMT), 0, 1, max(cell_metadata$EMT)),
                                         from = c(min(cell_metadata$EMT), max(cell_metadata$EMT))))
h5=ggplot()+geom_point(data=cell_metadata[order(cell_metadata$Hypoxia),],aes(x=Component_1,y=Component_2,col=Hypoxia),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Hypoxia")+
  scale_color_gradientn(colours = c("grey", "grey", "#FFB6C1", "#FFB6C1"), 
                        values = rescale(x = c(min(cell_metadata$Hypoxia), 0, 2.5, max(cell_metadata$Hypoxia)),
                                         from = c(min(cell_metadata$Hypoxia), max(cell_metadata$Hypoxia))))
h6=ggplot()+geom_point(data=cell_metadata[order(cell_metadata$Interferon),],aes(x=Component_1,y=Component_2,col=Interferon),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Interferon")+
  scale_color_gradientn(colours = c("grey", "grey", "#FF8C00", "#FF8C00"), 
                        values = rescale(x = c(min(cell_metadata$Interferon), 0, 5, max(cell_metadata$Interferon)),
                                         from = c(min(cell_metadata$Interferon), max(cell_metadata$Interferon))))

plot_list=list(h1,h2,h3,h4,h5,h6)
pdf("trajectory_score.pdf",width=15,height=7)
ggarrange(plotlist = plot_list,nrow=2,ncol=3)
dev.off()

pdf("trajectory_metadata.pdf",width = 5,height = 3.5)
ggplot()+geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=pseudotime),size=0.3)+
  geom_segment(aes_string(x="source_prin_graph_dim_1",
                          y="source_prin_graph_dim_2",
                          xend="target_prin_graph_dim_1",
                          yend="target_prin_graph_dim_2"),
               size=0.5,
               color="black",
               linetype="solid",
               na.rm=TRUE,
               data=edge_df)+
  scale_color_gradientn(colors=c("#1a2a6c","#b21f1f","#fdbb2d"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )

ggplot()+geom_point(data=cell_metadata[which(cell_metadata$cell_subtype != "AT2 cell"),],aes(x=Component_1,y=Component_2,fill = state,col=state),shape=21,size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )
cell_metadata_at=cell_metadata[sample(rownames(cell_metadata),nrow(cell_metadata),replace=F),]
ggplot()+geom_point(data=cell_metadata_at,aes(x=Component_1,y=Component_2,col=cell_subtype),size=0.3)+
  scale_color_manual(values = c("AT2 cell"="#54bbd5","Malignant cell"="#f17f73"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )+
  labs(title = "Cell type")
cell_metadata_stage=cell_metadata[which(cell_metadata$cell_subtype=="Malignant cell"),]
cell_metadata_stage=cell_metadata_stage[!is.na(cell_metadata_stage$stage),]
cell_metadata_stage=cell_metadata_stage[order(cell_metadata_stage$stage),]
ggplot()+geom_point(data=cell_metadata_stage,aes(x=Component_1,y=Component_2,col=stage),size=0.3)+
  scale_color_manual(values = c("I"="#DCE9F4","II"="#ABD0F1","III"="#FF9797","IV"="#E56F5E","Normal"="#9EC4BE"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )+
  labs(title = "Stage(Only malignant cells)")
ggplot()+geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=sex),size=0.3)+
  scale_color_manual(values = c("female"="#DC143C","male"="#4169E1"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )+
  labs(title = "Sex")

ggplot()+geom_point(data=cell_metadata[!is.na(cell_metadata$smoking),],aes(x=Component_1,y=Component_2,col=smoking),size=0.3)+
  scale_color_manual(values = c("Smoker"="#DC143C","Non-smoker"="#FFF0F5"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )+
  labs(title = "Smoking")

ggplot()+geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=ethnicity),size=0.3)+
  scale_color_manual(values = c("Asian"="#FFA500","European"="#FBE8D5"))+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )+
  labs(title = "Ethnicity")
cell_metadata_gii=cell_metadata
cell_metadata_gii$branch=as.character(cell_metadata_gii$branch)
cell_metadata_gii$branch[which(cell_metadata_gii$cell_subtype == "AT2 cell")]="AT2 cell"
ggplot()+geom_point(data=cell_metadata_gii[which(cell_metadata_gii$cell_subtype != "AT2 cell"),],aes(x=Component_1,y=Component_2,col=gii_score),size=0.3)+
  scale_color_gradientn(colours = c("white","white","#DC143C","#DC143C"), 
                        values = rescale(x = c(0,0.3,0.6,1),
                                         from = c(0, 1)))+
  labs(title = "GII score (Only malignant cells)")+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  )
cell_metadata_gii$branch=factor(cell_metadata_gii$branch,level=c("AT2 cell","Precursor","IFN-responsive","Proliferative"))
ggplot(data=cell_metadata_gii,aes(x=branch,y=gii_score))+geom_violin(aes(fill=branch))+
  geom_boxplot(width=0.3)+
  scale_fill_manual(values =c("AT2 cell"="#97C8AF","Precursor"="#51B1B7","IFN-responsive"="#FF8C00","Proliferative"="#DC143C"))+
  stat_compare_means(comparisons = list(c("AT2 cell","Precursor"),
                                        c("Precursor","IFN-responsive"),
                                        c("IFN-responsive","Proliferative"),
                                        c("Precursor","Proliferative")
  ),
  test = "wilcox.test",p.adjust.methods="fdr")+
  theme_classic()
ggplot()+geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=Alveolar),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Alveolar")+
  scale_color_gradientn(colours = c("grey", "grey","#51AAAE", "#51AAAE"), 
                        values = rescale(x = c(min(cell_metadata$Alveolar), 0, 0.5, max(cell_metadata$Alveolar)),
                                         from = c(min(cell_metadata$Alveolar), max(cell_metadata$Alveolar))))
ggplot()+geom_point(data=cell_metadata[order(cell_metadata$Stress),],aes(x=Component_1,y=Component_2,col=Stress),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Stress")+
  scale_color_gradientn(colours = c("grey", "grey", "#ADD487","#ADD487"), 
                        values = rescale(x = c(min(cell_metadata$Stress), 0, 1, max(cell_metadata$Stress)),
                                         from = c(min(cell_metadata$Stress), max(cell_metadata$Stress))))
ggplot()+geom_point(data=cell_metadata,aes(x=Component_1,y=Component_2,col=Cycle),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Cycle")+
  scale_color_gradientn(colours = c("grey", "grey", "#DC143C","#DC143C"), 
                        values = rescale(x = c(min(cell_metadata$Cycle), 0, 2, max(cell_metadata$Cycle)),
                                         from = c(min(cell_metadata$Cycle), max(cell_metadata$Cycle))))
ggplot()+geom_point(data=cell_metadata[order(cell_metadata$EMT),],aes(x=Component_1,y=Component_2,col=EMT),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "EMT")+
  scale_color_gradientn(colours = c("grey", "grey","#FFDAB9","#FFDAB9"), 
                        values = rescale(x = c(min(cell_metadata$EMT), 0, 1, max(cell_metadata$EMT)),
                                         from = c(min(cell_metadata$EMT), max(cell_metadata$EMT))))
ggplot()+geom_point(data=cell_metadata[order(cell_metadata$Hypoxia),],aes(x=Component_1,y=Component_2,col=Hypoxia),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Hypoxia")+
  scale_color_gradientn(colours = c("grey", "grey", "#FFB6C1", "#FFB6C1"), 
                        values = rescale(x = c(min(cell_metadata$Hypoxia), 0, 2.5, max(cell_metadata$Hypoxia)),
                                         from = c(min(cell_metadata$Hypoxia), max(cell_metadata$Hypoxia))))
ggplot()+geom_point(data=cell_metadata[order(cell_metadata$Interferon),],aes(x=Component_1,y=Component_2,col=Interferon),size=0.3)+
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
  ) +
  labs(title = "Interferon")+
  scale_color_gradientn(colours = c("grey", "grey", "#FF8C00", "#FF8C00"), 
                        values = rescale(x = c(min(cell_metadata$Interferon), 0, 5, max(cell_metadata$Interferon)),
                                         from = c(min(cell_metadata$Interferon), max(cell_metadata$Interferon))))
dev.off()

# Branch summary
mp_mat=cell_metadata[,c("branch",'Stress','Alveolar','Cycle','EMT','Hypoxia','Interferon')]
mp_mat=mp_mat %>%
  group_by(branch ) %>% 
  summarise_at(vars("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon"), mean)
mp_mat=data.frame(mp_mat)
rownames(mp_mat)=mp_mat$branch
mp_mat=mp_mat[,-1]
mp_mat=data.frame(mp_mat)
col_fun_gsea=colorRamp2(c(-1, 0, 1), c("#0E3061", "white", "#AF182D"))
pdf("trajectory_branch_summary.pdf",width=6,height = 3)
Heatmap(mp_mat,col=col_fun_gsea,name = "MP z-score")
dev.off()
