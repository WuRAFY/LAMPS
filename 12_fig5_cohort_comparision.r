library(scales)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggpubr)
library(tidyverse)
library(survminer)
library(broom)

## Cohort comparision, control for smoking status
chi.test <- function(a, b) {
  return(chisq.test(cbind(a, b)))
}
chi_plot=function(data,x,g,comparisons,label,col){
  temp_perc=data[!is.na(data[,x]),] %>% group_by_at(c(x,g)) %>% count() %>% group_by_at(c(g)) %>%
    mutate(countT= sum(n)) %>%
    group_by_at(c(x), add=TRUE) %>%
    mutate(per=round(100*n/countT,2))
  p=ggplot(data=temp_perc,aes_string(x=g,y="n",fill=x))+
    geom_bar(data=temp_perc,aes_string(x=g,y="per",fill=x),stat="identity",position="stack",alpha = 0.9)+
    scale_fill_manual(values=col)+
    geom_signif(comparisons = comparisons,
                map_signif_level = TRUE,
                y_position = c(105,110,115),
                textsize = 5,
                test = "chi.test")+
    theme_classic()+
    scale_y_continuous(breaks = seq(0, 100, 20),labels = seq(0,1,0.2))+
    labs(x="",y="",title=label)+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45,hjust = 1))
  return(p)
}
wilcox_test=function(data,x,g,comparisons,label,col){
  p=ggplot(data=data,aes_string(x=g,y=x))+
    geom_boxplot(alpha = 0.9)+
    geom_signif(comparisons = comparisons,
                map_signif_level = TRUE,
                step_increase=0.05,
                textsize = 5,
                test = "wilcox.test")+
    theme_classic()+
    labs(x="",y="",title=label)+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45,hjust = 1))
  return(p)
}

info=read.table("./input/bulk_info.txt",sep="\t",header = T, stringsAsFactors = T)
chi_plot(info,x="Smoking",g="Ancestry",comparisons = list(c("Asian","European")),label = "Smoking",
         col=c("Smoker"="#DC143C","Non-smoker"="#FFF0F5"))

info_s=info[which(info$Smoking=="Smoker"),]
info_ns=info[which(info$Smoking=="Non-smoker"),]

p1=chi_plot(info_s,x="Stage",g="Ancestry",comparisons = list(c("Asian","European")),label = "Stage",
         col=c("I"="#DCE9F4","II"="#ABD0F1","III"="#FF9797","IV"="#E56F5E","Normal"="#9EC4BE"))
p2=chi_plot(info_s,x="Sex",g="Ancestry",comparisons = list(c("Asian","European")),label = "Sex",
         col=c("female"="#FFE4E1","male"="#6CA6CD"))
p3=ggplot(data=info_s,aes_string(x="Ancestry",y="TMB"))+
  geom_violin(aes(fill=Ancestry),alpha=0.9)+
  geom_boxplot(fill="white",outliers = FALSE,width=0.15)+
  scale_fill_manual(values=c("Asian"="#FFA500","European"="#FBE8D5"))+
  geom_signif(comparisons =list(c("Asian","European")) ,
              map_signif_level = TRUE,
              step_increase=0.05,
              textsize = 5,
              test = "wilcox.test")+
  theme_classic()+
  labs(x="",y="",title="TMB")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1))
p4=ggplot(data=info_s,aes_string(x="Ancestry",y="Age"))+
  geom_violin(aes(fill=Ancestry),alpha=0.9)+
  geom_boxplot(fill="white",outliers = FALSE,width=0.15)+
  scale_fill_manual(values=c("Asian"="#FFA500","European"="#FBE8D5"))+
  geom_signif(comparisons =list(c("Asian","European")) ,
              map_signif_level = TRUE,
              step_increase=0.05,
              textsize = 5,
              test = "wilcox.test")+
  theme_classic()+
  labs(x="",y="",title="Age")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1))

p5=chi_plot(info_ns,x="Stage",g="Ancestry",comparisons = list(c("Asian","European")),label = "Stage",
            col=c("I"="#DCE9F4","II"="#ABD0F1","III"="#FF9797","IV"="#E56F5E","Normal"="#9EC4BE"))
p6=chi_plot(info_ns,x="Sex",g="Ancestry",comparisons = list(c("Asian","European")),label = "Sex",
            col=c("female"="#FFE4E1","male"="#6CA6CD"))
p7=ggplot(data=info_ns,aes_string(x="Ancestry",y="TMB"))+
  geom_violin(aes(fill=Ancestry),alpha=0.9)+
  geom_boxplot(fill="white",outliers = FALSE,width=0.15)+
  scale_fill_manual(values=c("Asian"="#FFA500","European"="#FBE8D5"))+
  geom_signif(comparisons =list(c("Asian","European")) ,
              map_signif_level = TRUE,
              step_increase=0.05,
              textsize = 5,
              test = "wilcox.test")+
  theme_classic()+
  labs(x="",y="",title="TMB")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1))
p8=ggplot(data=info_ns,aes_string(x="Ancestry",y="age"))+
  geom_violin(aes(fill=Ancestry),alpha=0.9)+
  geom_boxplot(fill="white",outliers = FALSE,width=0.15)+
  scale_fill_manual(values=c("Asian"="#FFA500","European"="#FBE8D5"))+
  geom_signif(comparisons =list(c("Asian","European")) ,
              map_signif_level = TRUE,
              step_increase=0.05,
              textsize = 5,
              test = "wilcox.test")+
  theme_classic()+
  labs(x="",y="",title="Age")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust = 1))

pdf("cohort_comparision.pdf",width=10,height = 10)
ggarrange(plotlist = list(p1,p2,p3,p4,p5,p6,p7,p8),ncol = 4,nrow = 2)
dev.off()

## Branch/TME regression model
source("utilities.r")
info$Sex=factor(info$Sex,levels=c("male","female"),labels=c("Male","Female"))
info$TME=factor(info$TME.subtype..k.4.,levels = c("Initiating","Immune recruiting","Immune activated","Fibrotic"))
ggforest_multi <- function(model, data = NULL,
                     main = "Estimates", cpositions=c(0.02, 0.22, 0.4),
                     fontsize = 0.7, refLabel = "reference", noDigits=2) {
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  df_list=list()
  for(i in unique(coef$y.level)){
  conf.high <- conf.low <- estimate <- NULL
  data  <- .get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  coef = coef[which(coef$y.level == i),]
  gmodel <- glance(model)
  allTerms <- lapply(seq_along(terms), function(i){
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data),
                 pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data),
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[,1:2], 1, paste0, collapse="")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds,])[,c("var", "level", "N", "p.value", "estimate", "conf.low", "conf.high", "pos")]
  toShowExp <- toShow[,5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(toShowExp, digits=noDigits)
  toShowExpClean <- data.frame(toShow,
                               pvalue = signif(toShow[,4],noDigits+1),
                               toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, noDigits+1), " ",
                                 ifelse(toShowExpClean$p.value < 0.05, "*",""),
                                 ifelse(toShowExpClean$p.value < 0.01, "*",""),
                                 ifelse(toShowExpClean$p.value < 0.001, "*",""))
  toShowExpClean$ci <- paste0("(",toShowExpClean[,"conf.low.1"]," - ",toShowExpClean[,"conf.high.1"],")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var = paste0(toShowExpClean$var," (",i,")")
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=",toShowExpClean$N,")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, ]
  toShowExpClean = toShowExpClean[which(toShowExpClean$estimate > -14),]
  df_list[[i]]=toShowExpClean}
  toShowExpClean=do.call(rbind,df_list)
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, na.rm = TRUE)
  breaks <- axisTicks(rangeb, log = FALSE, nint = 7)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + .15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] +  cpositions[1] * width
  y_nlevel <- rangeplot[1]  +  cpositions[2] * width
  y_cistring <- rangeplot[1]  +  cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize * 3
  
  p <- ggplot(toShowExpClean, aes(seq_along(var), estimate)) +
    geom_rect(aes(xmin = seq_along(var) - .5, xmax = seq_along(var) + .5,
                  ymin = rangeplot[1], ymax = rangeplot[2],
                  fill = ordered(seq_along(var) %% 2 + 1))) +
    scale_fill_manual(values = c("#FFFFFF33", "#00000033"), guide = "none") +
    geom_point(pch = 15, size = 4) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
    geom_hline(yintercept = 0, linetype = 3) +
    coord_flip(ylim = rangeplot) +
    ggtitle(main) +
    scale_y_continuous(
      name = "",
      labels = sprintf("%g", breaks),
      expand = c(0.02, 0.02),
      breaks = breaks) +
    theme_light() +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          legend.position = "none",
          panel.border=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    xlab("") +
    annotate(geom = "text", x = x_annotate, y = y_variable,
             label = toShowExpClean$var, fontface = "bold", hjust = 0,
             size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = y_nlevel, hjust = 0,
             label = toShowExpClean$level, vjust = -0.1, size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = y_nlevel,
             label = toShowExpClean$N, fontface = "italic", hjust = 0,
             vjust = ifelse(toShowExpClean$level == "", .5, 1.1),
             size = annot_size_mm) +
    annotate(geom = "text", x = x_annotate, y = y_cistring,
             label = toShowExpClean$estimate.1, size = annot_size_mm,
             vjust = ifelse(toShowExpClean$estimate.1 == "reference", .5, -0.1)) +
    annotate(geom = "text", x = x_annotate, y = y_cistring,
             label = toShowExpClean$ci, size = annot_size_mm,
             vjust = 1.1,  fontface = "italic") +
    annotate(geom = "text", x = x_annotate, y = y_stars,
             label = toShowExpClean$stars, size = annot_size_mm,
             hjust = -0.2,  fontface = "italic") 
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
}
model_branch=nnet::multinom(Branch ~ Ancestry + Smoking + 
              Sex + Stage + Age, data = info,family = "binomial")
model_tme=nnet::multinom(TME ~ Ancestry + Smoking + 
                              Sex + Stage + Age, data = info,family = "binomial")
pdf("branch.pdf",width=6,height = 6)
ggforest_multi(model_branch)
dev.off()
pdf("tme.pdf",width=6,height = 9)
ggforest_multi(model_tme)
dev.off()