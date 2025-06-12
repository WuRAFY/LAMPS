library(tidyverse)
library(phangorn)
library(readxl)
library(ggtree)
library(ggplot2)
library(patchwork)
library(cowplot)
library(phytools)
library(igraph)
library(ggraph)
library(tidygraph)
library(ape)
library(treeio)
library(extrafont)
library(vcfR)
library(ggnewscale)
library(patchwork)
loadfonts()

######################################################## function #######################################################
#######################
## build the nj tree ##
#######################
rtrw_build <- function(mutation_code_sample){
  # phyDat1 <- phyDat(t(mutation_code_sample), type = "AA")
  phyDat <- phyDat(t(mutation_code_sample), type = "USER", levels = c("R", "M"))
  # hamming1 <- dist.hamming(phyDat1, ratio = F)
  hamming <- dist.hamming(phyDat, ratio = F)
  
  trw <- NJ(hamming)
  rtrw=root(trw, outgroup = "N", resolve.root = T)
  plot(rtrw, main = patient, )
  return(rtrw)
}

#########################################################################################################
#####################
## Fitch Algorithm ##
#####################
{
  # get the preliminary state set
  fitch_bottom_up <- function(tree, tip_states) {
    # 初始化节点状态和替换次数
    tree <- reorder(tree, "postorder")
    node_states <- list()
    total_changes <- 0
    all_nodes <- unique(c(tree$edge[, 2], tree$edge[, 1]))
    for (i in all_nodes) {
      if (i <= length(tree$tip.label)) {
        node_states[[i]] <- tip_states[[tree$tip.label[i]]]
      } else {
        children <- tree$edge[tree$edge[, 1] == i, 2]
        child_states <- list(node_states[[children[1]]], node_states[[children[2]]])
        intersection <- intersect(child_states[[1]], child_states[[2]])
        if (length(intersection) > 0) {
          node_states[[i]] <- intersection
        } else {
          node_states[[i]] <- union(child_states[[1]], child_states[[2]])
          total_changes <- total_changes + 1
        }
      }
    }
    return(list(node_states = node_states, total_changes = total_changes))
  }
  
  ## get final node state set
  check_internal_node <- function(best_result){
    tree <- best_result$tree
    all_nodes <- unique(c(tree$edge[, 1]))
    for(i in all_nodes[-1]){
      children <- tree$edge[tree$edge[, 1] == i, 2]
      parent <- tree$edge[tree$edge[, 2] == i, 1]
      
      parent_stateSet <- best_result$node_states[[parent]]
      preliminary_stateSet <- best_result$node_states[[i]]
      
      children_union <- union(best_result$node_states[[children[1]]], best_result$node_states[[children[2]]])
      children_inter <- intersect(best_result$node_states[[children[1]]], best_result$node_states[[children[2]]])
      
      parent_more <- setdiff(parent_stateSet, preliminary_stateSet)
      
      if(length(parent_more) == 0){
        best_result$node_states[[i]] <- parent_stateSet
      }else if(length(setdiff(preliminary_stateSet, children_inter))!=0){
        stateSet_diff <- setdiff(parent_stateSet, preliminary_stateSet)
        best_result$node_states[[i]] <- c(best_result$node_states[[i]], stateSet_diff)
        
      }else if(sum(parent_more %in% best_result$node_states[[children[1]]], parent_more %in% best_result$node_states[[children[2]]]) > 0){
        children1 <- intersect(best_result$node_states[[children[1]]], parent_stateSet)
        children2 <- intersect(best_result$node_states[[children[2]]], parent_stateSet)
        union_child <- union(union(children1, children2), best_result$node_states[[i]])
        
        best_result$node_states[[i]] <- union_child
        
      }else{
        best_result$node_states[[i]] <- best_result$node_states[[i]]
      }
    }
    return(best_result)
  }
  
  ## get transition score
  get_transition_score <- function(best_result, adjust_result, transition_score){
    preliminary_node <- best_result$node_states
    final_node <- adjust_result$node_states
    tree <- adjust_result$tree
    p_list <- list()
    for (internal_node in unique(tree$edge)) {
      ancestral_states <- final_node[[internal_node]]
      ancestral_score <- rep(0, length(ancestral_states))
      names(ancestral_score) <- ancestral_states
      p_list[[internal_node]] <- ancestral_score
    }
    ancestral_states <- final_node[[unique(tree$edge[, 1])[1]]]
    ancestral_score <- rep(1/length(ancestral_states), length(ancestral_states))
    names(ancestral_score) <- ancestral_states
    p_list[[unique(tree$edge[, 1])[1]]] <- ancestral_score
    
    for(internal_node in unique(tree$edge[, 1])){
      ancestral_states <- final_node[[internal_node]]
      children <- tree$edge[tree$edge[,1]==internal_node, 2]
      
      for(temp_state in ancestral_states){
        for (child in children){
          preliminary_node_child <- preliminary_node[[child]]
          final_node_child <- final_node[[child]]
          if(temp_state %in% preliminary_node_child){
            child_score <- p_list[[internal_node]][temp_state]
            child_state <- temp_state
            p_list[[child]][child_state] <- p_list[[child]][child_state] + child_score
          }else{
            
            if(temp_state %in% final_node_child){
              child_score <- p_list[[internal_node]][temp_state]/length(final_node_child)
              child_state <- final_node_child
              p_list[[child]][child_state] <- p_list[[child]][child_state] + child_score
              final_node_child_more <- setdiff(final_node_child, preliminary_node_child)
              change_node <- setdiff(child_state, temp_state)
              change_node <- setdiff(change_node, final_node_child_more)
              
              transition_score[temp_state, change_node] <- transition_score[temp_state, change_node] + child_score
              
            }else{
              child_score <- p_list[[internal_node]][temp_state]/length(preliminary_node_child)
              child_state <- preliminary_node_child
              p_list[[child]][child_state] <- p_list[[child]][child_state] + child_score
              
              change_node <- child_state
              transition_score[temp_state, change_node] <- transition_score[temp_state, change_node] + child_score
              
            }
          }
        }
      }
    }
    transition_score <- round(transition_score, 2)
    return(c(transition_score = list(transition_score), p_list = list(p_list)))
  }
}

###############################
## plot nj tree With subtype ##
###############################
rtrw_plotWithsubtype_LUAD <- function(rtrw, tree_annotation, patient = substr(rtrw$tip.label[1], 1, 9)){
  tree_annotation$Subtype_Name <- as.character(tree_annotation$Subtype_Name)
  # patient <- substr(rtrw$tip.label[1], 1, 9)
  
  # get the branch.length of "N"
  rtrw_tibble <- as_tibble(rtrw)
  root.edge  <- as.numeric(rtrw_tibble[which(rtrw_tibble$label == "N"), "branch.length"])
  Nnode <- rtrw$Nnode
  
  # delete the root node
  rtrw <- drop.tip(rtrw, "N")
  rtrw_tibble <- as_tibble(rtrw)
  rtrw_tibble$root_label <- NA
  rtrw_tibble$root_label[which(rtrw_tibble$parent == rtrw_tibble$node)] <- "N"
  rtrw_tibble$root_subtype <- NA
  
  max_depth <- max(node.depth.edgelength(rtrw))
  max_depth_root <- max_depth + root.edge
  
  # give the normal sample subtype
  
  rtrw_tibble$root_label[which(rtrw_tibble$parent == rtrw_tibble$node)] <- "N"
  rtrw_tibble$root_subtype[which(rtrw_tibble$parent == rtrw_tibble$node)] <- "Normal" 
  
  # change the data type of tree
  rtrw_annotation <- as.treedata(rtrw_tibble)
  
  if(length(rtrw$tip.label)>20){
    tiplab_size = 2.5
  }else if(length(rtrw$tip.label)>15){
    tiplab_size = 4
  }else{
    tiplab_size = 5
  }
  
  # plot the tree by ggtree
  p1 <- ggtree(rtrw_annotation) %<+% 
    tree_annotation+
    geom_tiplab(aes(color = Subtype_Name), hjust = -0.1, vjust = 0.5, show.legend = F, size = tiplab_size)+
    geom_tippoint(aes(color=Subtype_Name), size = tiplab_size, show.legend = F)+
    geom_rootedge(root.edge)+
    geom_label(aes(x=-(root.edge+2), label=root_label, colour = root_subtype), label.size = NA, show.legend = F)+
    geom_treescale(0, 0)+
    scale_color_manual(values = c("Normal"="#9EC4BE","Initiating"="#80BA8A","Immune recruiting"="#6BB7CA","Immune killing"="#F4CEB4","Fibrosis"="#ED9F9B", "Precursor"="#52A8AC", "IFN-responsive"="#EE8827", "Proliferative"="#CA213E"))+
    ggtitle(patient)
  
  p2 <- gheatmap(p1,
                 tree_annotation[, "branch", drop=F],
                 offset = max_depth_root*1.6,
                 width = max_depth_root/max_depth*0.3,
                 colnames_angle = -45,
                 hjust = 0,
                 # colnames_offset_y = .4,
                 colnames=F) +
    labs(fill = "Branch") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
    scale_fill_manual(values = c("Normal"="#9EC4BE","Initiating"="#80BA8A","Immune recruiting"="#6BB7CA","Immune killing"="#F4CEB4","Fibrosis"="#ED9F9B", "Precursor"="#52A8AC", "IFN-responsive"="#EE8827", "Proliferative"="#CA213E", "Immune cold"="#7AADCC", "Immune hot"="#E5909F"))
  return(p2)
}

############################################
## plot the network of subtype transition ##
############################################
get_network <- function(transition_score){
  edges <- transition_score %>%
    as.data.frame.table(responseName = "weight") %>%
    filter(weight > 0) %>%
    rename(from = Var1, to = Var2) %>%
    mutate(
      from = as.character(from),
      to = as.character(to)
    )
  graph <- as_tbl_graph(
    graph_from_data_frame(edges, directed = TRUE)
  )
  p1 <- ggraph(graph, layout = "circle") + 
    geom_edge_fan(aes(edge_width = weight, label = weight),
                  show.legend = FALSE,
                  arrow = arrow(length = unit(3, 'mm')), 
                  angle_calc = "along",
                  label_dodge = unit(2.5, "mm"),
                  label_size = 4,
                  check_overlap = T,
                  strength = 2,
                  end_cap = circle(8, 'mm'),
                  family = "Times New Roman")+
    geom_node_point(size = 13, shape=21, aes(fill = name, color = name))+ 
    scale_fill_manual(values = c("Normal"="#9EC4BE","Initiating"="#80BA8A","Immune recruiting"="#6BB7CA","Immune killing"="#F4CEB4","Fibrosis"="#ED9F9B", "Precursor"="#52A8AC", "IFN-responsive"="#EE8827", "Proliferative"="#CA213E"))+
    scale_color_manual(values = c("Normal"="#9EC4BE","Initiating"="#80BA8A","Immune recruiting"="#6BB7CA","Immune killing"="#F4CEB4","Fibrosis"="#ED9F9B", "Precursor"="#52A8AC", "IFN-responsive"="#EE8827", "Proliferative"="#CA213E"))+
    scale_edge_width(range=c(0,3))+
    guides(size=F,fill=F)+
    theme_graph()+
    xlim(-1.1,1.1)+
    ylim(-1.1,1.1)+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    guides(fill=FALSE, color = F) 
  return(p1)
}
######################################################## function #######################################################

###########################################################################################
## 一个病人例子
###########################################################################################
input_path <- "./input/"
out_path <- "./result/"
dir.create(out_path)

###################
## load the data ##
###################
patient <- "LUAD_RM03"
DNA_RNA_merge <- readRDS(paste0(input_path, "DNA_RNA_merge_", patient, ".rds"))
# Due to data control, mutations were masked
mutation_code_sample <- readRDS(paste0(input_path, "mutation_code_sample_", patient, ".rds"))

#######################
## build the nj tree ##
#######################
rtrw <- rtrw_build(mutation_code_sample)
write.tree(rtrw, paste0(out_path, patient, ".tree"))

###################################################
##### get the transition path based on Fitch ######
###################################################
rtrw <- read.tree(paste0(out_path, patient, ".tree"))
trw <- drop.tip(rtrw, "N")

fitch_subtype <- matrix(as.character(DNA_RNA_merge$Subtype_Name), ncol = 1)
rownames(fitch_subtype) <- DNA_RNA_merge$mask_sample_id
tip_states <- lapply(fitch_subtype[trw$tip.label,], function(x){x})

best_result <- fitch_bottom_up(trw, tip_states)
best_result <- c(tree = list(trw), best_result)
adjust_result <- check_internal_node(best_result)

## 初始化状态转变矩阵
transition_score <- matrix(0, nrow = 3, ncol = 3)
rownames(transition_score) <- c("IFN-responsive", "Precursor", "Proliferative")
colnames(transition_score) <- c("IFN-responsive", "Precursor", "Proliferative")

## 得到状态转变矩阵
patient_transition <- get_transition_score(best_result, adjust_result, transition_score)
patient_transition <- c(rtrw = list(rtrw), patient_transition)

######################
##### save plot ######
######################
p1 <- rtrw_plotWithsubtype_LUAD(rtrw, DNA_RNA_merge, patient = substr(rtrw$tip.label[1], 1, 9))
print(p1)
p2 <- get_network(patient_transition$transition_score)
print(p2)

cairo_pdf(paste0(out_path, paste0("TreeWithSubtypeTransition_", patient, ".pdf")), width = 10, height = 6)
print(p1 + p2)
dev.off()


