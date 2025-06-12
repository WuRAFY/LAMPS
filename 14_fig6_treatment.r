library(ggsankey)
library(oncokbR)
library(dplyr)
library(ggsankey)
library(circlize)
library(ggplot2)
library(cbioportalR)
library(future)
library(future.apply)

## Using OncoKB to annotate mutations
# Due to data control, this data will be supplied upon reasonable request
load("./input/oncokb.RData")
## ICB response prediction
# Due to data control, this data will be supplied upon reasonable request
tide=read.csv("./input/tide.csv")
rownames(tide)=tide$Patient
info$tide=tide[info$sample,]$Responder
info$tide_score=tide[info$sample,]$TIDE
info$r=ifelse(info$tide_score>-1,"N","Y")
info$target=factor(info$target,level=c("LEVEL_1","LEVEL_2","LEVEL_3A","LEVEL_3B","LEVEL_4","No_target"),
                   labels = c("LEVEL_1","LEVEL_2","No_target","No_target","No_target","No_target"))
info$target=paste0(as.character(info$target),"_",info$mut)
df=data.frame(table(info$tme,info$branch))
TME=c("Initiating"="#80BA8A","Immune recruiting"="#F4CEB4","Immune killing"="#ED9F9B","Fibrosis"="#6BB7CA")
branch = c("Precursor"="#51B1B7","IFN-responsive"="#FF8C00","Proliferative"="#DC143C")

pdf("therapy_freq.pdf",width = 6,height = 2.5)
df=info %>% group_by(tme,ethnicity) %>% count(r) %>% group_by(ethnicity) %>% mutate(perc=n/sum(n))
df=df[which(df$r=="Y"),]
ggplot()+geom_bar(data=df,aes(x=tme,y=perc,fill=ethnicity),stat = "identity",position="dodge")+
  scale_fill_manual(values = c("Asian"="#A32A31","European"="#DFB6BC"))+
  ylim(c(0,0.32))+
  theme_classic()
df_t=info %>% group_by(branch,ethnicity) %>% count(target) %>% group_by(ethnicity) %>% mutate(perc=n/sum(n))
df_t=df_t[-grep("No_target",df_t$target),]
ggplot()+geom_bar(data=df_t,aes(x=branch,y=perc,fill=target),stat = "identity",position="stack")+
  #scale_fill_manual(values = c("Asian"="#407BD0","European"="#B7D0EA"))+
  scale_fill_manual(values=c("LEVEL_1_BRAF"="#92B4C8",
                    "LEVEL_1_EGFR"="#ABD3E1",
                    "LEVEL_1_ERBB2"="#E3EDE0",
                    "LEVEL_1_KRAS"="#FFE9BE",
                    "LEVEL_1_MET"="#FAC795",
                    "LEVEL_2_MET.AMPLIFICATION"="#EEA599"))+
  ylim(c(0,0.32))+
  facet_wrap(vars(ethnicity))+
  theme_classic()
dev.off()