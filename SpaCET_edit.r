SpaCET.deconvolution.malignant_edit <- function(SpaCET_obj, Malignant="Malignant", malignantCutoff=0.7, coreNo=6)
{
  coreNoDect <- parallel::detectCores(logical = FALSE)
  if(coreNoDect<coreNo)
  {
    print(paste0("Since the number of your physical cores is ",coreNoDect,", coreNo=",coreNoDect," is used automatically."))
    coreNo <- coreNoDect
  }
  if(Sys.info()[['sysname']] == "Windows")
  {
    print("Since Windows does not support > 1 core, coreNo=1 is used automatically.")
    coreNo <- 1
  }

  if(length(SpaCET_obj@results$deconvolution$Ref$lineageTree[[Malignant]])>1)
  {
    stop("Your deconvolution results have included multiple malignant cell states. We do not recommend deconvolve malignant cell fraction further.")
  }

  if(is.null(SpaCET_obj@results$deconvolution$propMat))
  {
    stop("Please do the complete deconvolution first by using SpaCET.deconvolution.")
  }else{
    res_deconv <- SpaCET_obj@results$deconvolution$propMat
  }

  if(length(Malignant)>1 | length(Malignant)==0)
  {
    stop("Please input the only one major malignant cell type.")
  }else{
    if(!Malignant%in%rownames(res_deconv))
    {
      stop("The input malignant cell type does not exist in the deconvolution results. Please check whether you input correct the name of malignant cell type. Of note, R language is case sensitive generally.")
    }
  }

  if(malignantCutoff>1 | malignantCutoff<0)
  {
    stop("Please input a value within 0~1 for the cutoff of malignant spots.")
  }
  st.matrix.data <- as.matrix(SpaCET_obj@input$counts)
  st.matrix.data <- st.matrix.data[rowSums(st.matrix.data)>0,]
  st.matrix.data.mal <- st.matrix.data[,res_deconv[Malignant,]>=malignantCutoff]
  st.matrix.data.mal.CPM <- sweep(st.matrix.data.mal, 2, Matrix::colSums(st.matrix.data.mal), "/") *1e5
  st.matrix.data.mal.log <- log2(st.matrix.data.mal.CPM+1)
  # clustering
  set.seed(123)
  suppressPackageStartupMessages(
    library(MUDAN)
  )
  med_scale=function(x){
  temp_med=median(x)
  temp_mad=mad(x)
  return((x-temp_med)/temp_mad)
  }
  #hajack transcriptome clustering, replace PC with scores
  pcs <- t(SpaCET_obj@results$GeneSetScore)[colnames(st.matrix.data.mal),c("Stress","Alveolar","Cycle","EMT","Hypoxia","Interferon")]
  pcs <- apply(pcs,2,med_scale)
  d <- as.dist(1-cor(t(pcs)))
  hc <- hclust(d, method='ward.D')
  cluster_numbers <- 2:9
  clustering <- cutree(hc,k=cluster_numbers)
  clustering <- t(clustering)
  rownames(clustering) <- paste0("c",rownames(clustering))
  # silhouette
  suppressPackageStartupMessages({
    library(factoextra)
    library(NbClust)
    library(cluster)
  })
  v <- c()
  for(i in cluster_numbers)
  {
    clustering0 <- cutree(hc,k=i)
    sil <- silhouette(clustering0, d, Fun=mean)
    v <- c(v, mean(sil[,3]))
    #db = intCriteria(pcs,clustering0,c("Davies_Bouldin"))
    #db = db$davies_bouldin
    #v=c(v,db)
  }
  maxN =3
  silMat <- cbind(cluster=cluster_numbers,silhouette=v)
  silMat <- cbind(silMat,maxN=cluster_numbers%in%maxN)
  clustering <- apply(clustering,1:2,function(x) LETTERS[x])
  Content <- as.character(clustering[maxN-1,])
  names(Content) <- colnames(clustering)
  states <- sort(unique(Content))
  refProfiles <- data.frame()
  sigGenes <- list()
  lineageTree <- list()
  refProfiles[rownames(st.matrix.data.mal.CPM),"Malignant"] <- rowMeans(st.matrix.data.mal.CPM)
  for(i in states)
  {
    refProfiles[rownames(st.matrix.data.mal.CPM),paste0("Malignant cell state ",i)] <- rowMeans(st.matrix.data.mal.CPM[,Content==i])
    tempMarkers <- c()
    for(j in setdiff(states,i))
    {
      library(limma)
      TT <- as.numeric(Content%in%c(i))
      WT <- as.numeric(Content%in%c(j))
      design <- cbind(TT,WT)
      fit <- lmFit(st.matrix.data.mal.log,design)
      cont.matrix <- makeContrasts(TTvsWT=TT-WT,levels=design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      res <- topTable(fit2,coef=1,number=nrow(st.matrix.data.mal.log))
      res <- res[order(res[,"t"],decreasing=T),]
      tempMarkers <- c(tempMarkers, rownames(res)[1:500])
    }
    tempMarkers <- table(tempMarkers)
    sigGenes[[paste0("Malignant cell state ",i)]] <- names(tempMarkers)[tempMarkers==1]
  }
  lineageTree[[Malignant]] <- paste0("Malignant cell state ",states)
  Refnew <- list(refProfiles=refProfiles, sigGenes=sigGenes, lineageTree=lineageTree)
  knownCellTypes <- names(SpaCET_obj@results$deconvolution$Ref$lineageTree)
  knownCellTypes <- setdiff(knownCellTypes,Malignant)
  if("Unidentifiable"%in%rownames(SpaCET_obj@results$deconvolution$propMat))
  {
    knownCellFractions <- c(knownCellTypes,"Unidentifiable")
  }else{
    knownCellFractions <- knownCellTypes
  }
  #res_deconv=res_deconv[-grep("Malignant cell state",rownames(res_deconv)),]
  propMat <- SpatialDeconv(
    ST=st.matrix.data,
    Ref=Refnew,
    malProp=res_deconv[knownCellFractions,],
    malRef=SpaCET_obj@results$deconvolution$Ref$refProfiles[,knownCellTypes],
    mode="deconvMal",
    coreNo=coreNo
  )
  propMat<- rbind(res_deconv,propMat[!rownames(propMat)%in%"Malignant",])
  SpaCET_obj@results$deconvolution$propMat <- propMat
  SpaCET_obj@results$deconvolution$clustering <- Content
  SpaCET_obj
}

SpaCET.CCI.colocalization_edit <- function(SpaCET_obj)
{
  res_deconv <- SpaCET_obj@results$deconvolution$propMat
  res_malig = res_deconv[grep("Malignant",rownames(res_deconv)),]
  res_malig = scale(res_malig,center=F,scale=res_deconv["Malignant",])
  res_tme = res_deconv[-grep("Malignant",rownames(res_deconv)),]
  res_tme = scale(res_tme,center=F,scale=colSums(res_deconv[c("CAF","Endothelial","Plasma","B cell","T CD4","T CD8","NK","cDC","pDC","Macrophage","Mast","Neutrophil","Unidentifiable"),]))
  res_deconv = rbind(res_malig,res_tme)
  res_deconv[is.na(res_deconv)]=0
  res_deconv[is.infinite(res_deconv)]=1
  res_deconv <- res_deconv[!rownames(res_deconv)%in%c("Unidentifiable","Macrophage other"),]
  res_deconv <- round(res_deconv,2)
  overallFraction <- rowMeans(res_deconv)

  cc_corr <- psych::corr.test(t(res_deconv),t(res_deconv),method="spearman",adjust="none",ci=FALSE)

  cc_corr_r <- cc_corr$r
  cc_corr_p <- cc_corr$p

  cc_corr_r.m <- reshape2::melt(cc_corr_r)
  cc_corr_p.m <- reshape2::melt(cc_corr_p)

  summary_df <- data.frame(
    cell_type_1 = as.character(cc_corr_r.m[,1]),
    cell_type_2 = as.character(cc_corr_r.m[,2]),
    fraction_product = overallFraction[as.character(cc_corr_r.m[,1])]*overallFraction[as.character(cc_corr_r.m[,2])],
    fraction_rho = round(cc_corr_r.m[,3],3),
    fraction_pv = cc_corr_p.m[,3]
    )
  rownames(summary_df) <- paste0(summary_df[,1],"_",summary_df[,2])
  summary_df[is.na(summary_df[,"fraction_rho"] ),"fraction_rho"] <- 0
  summary_df[is.na(summary_df[,"fraction_pv"] ),"fraction_pv"] <- 1

  Ref <- SpaCET_obj@results$deconvolution$Ref
  reff <- Ref$refProfiles[unique(unlist(Ref$sigGenes[names(Ref$sigGenes)%in%c(names(Ref$lineageTree),"T cell")])),]
  reff <- reff-rowMeans(reff)

  cc_corr <- psych::corr.test(reff,reff,method="spearman",adjust="none",ci=FALSE)

  cc_corr_r <- cc_corr$r
  cc_corr_p <- cc_corr$p

  cc_corr_r.m <- reshape2::melt(cc_corr_r)
  cc_corr_p.m <- reshape2::melt(cc_corr_p)

  summary_df2 <- data.frame(
    cell_type_1 = as.character(cc_corr_r.m[,1]),
    cell_type_2 = as.character(cc_corr_r.m[,2]),
    reference_rho = round(cc_corr_r.m[,3],3),
    reference_pv = cc_corr_p.m[,3]
  )
  rownames(summary_df2) <- paste0(summary_df2[,1],"_",summary_df2[,2])

  summary_df[rownames(summary_df2),"reference_rho"] <- summary_df2[,"reference_rho"]
  summary_df[rownames(summary_df2),"reference_pv"] <- summary_df2[,"reference_pv"]

  summary_df <- summary_df[ summary_df[,1] != summary_df[,2], ] #remove same cell type

  SpaCET_obj@results$CCI$colocalization <- summary_df

  SpaCET_obj
}

SpaCET.CCI.cellTypePair_edit <- function(SpaCET_obj, cellTypePair){
  groupMat <- data.frame()
  testRes <- data.frame()
  cellTypePair <- sort(cellTypePair)
  res_deconv <- SpaCET_obj@results$deconvolution$propMat
  LRNetworkScoreMat <- SpaCET_obj@results$CCI$LRNetworkScore

  cutoff1 <- summary(res_deconv[cellTypePair[1],])[5]
  cutoff2 <- summary(res_deconv[cellTypePair[2],])[5]

  cutoff11 <- quantile(res_deconv[cellTypePair[1],],0.85)
  cutoff22 <- quantile(res_deconv[cellTypePair[2],],0.85)

  Content <- res_deconv[1,]
  for(i in 1:length(Content))
  {
    if(res_deconv[cellTypePair[1],i]>cutoff11 & res_deconv[cellTypePair[2],i]>cutoff22)
    {
      Content[i] <- "Both"
    }else if(res_deconv[cellTypePair[1],i]>cutoff11 & res_deconv[cellTypePair[2],i]<cutoff2){
      Content[i] <- cellTypePair[1]
    }else if(res_deconv[cellTypePair[2],i]>cutoff22 & res_deconv[cellTypePair[1],i]<cutoff1){
      Content[i] <- cellTypePair[2]
    }else{
      Content[i] <- NA
    }
  }

  groupMat[paste0(cellTypePair[1],"_",cellTypePair[2]),colnames(res_deconv)] <- Content

  if(sum(Content%in%"Both")>5)
  {
    fg.df <- data.frame(group=Content,value=LRNetworkScoreMat[2,],stringsAsFactors=FALSE)
    fg.df <- fg.df[!fg.df[,"group"]%in%NA,]
    fg.df[fg.df[,1]%in%cellTypePair,1] <- "Single"

    cohend_res <- psych::cohen.d(fg.df, group="group", alpha=.05, std=TRUE)
    cd1 <- signif(cohend_res$cohen.d["value","effect"],2)

    n1 <- sum(fg.df[,1]%in%"Both")
    n2 <- sum(!fg.df[,1]%in%"Both")
    pv2 <- signif(wilcox.test(fg.df[fg.df[,1]%in%"Both",2],fg.df[!fg.df[,1]%in%"Both",2])$p.value,2)

    testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"groupCompare_cohen.d"] <- cd1
    testRes[paste0(cellTypePair[1],"_",cellTypePair[2]),"groupCompare_pv"] <- pv2
    return(c(cd1,pv2))
   }else{
    return(c(NA,NA))
   }
}