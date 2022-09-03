## TCGA fixed script

library(survival)
library(survminer)
library(tidyverse)
library(GSVA)
library(customLayout)
library(mixtools)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

load('//gpfs/output/20201016_Loom_files_for_Bladder_Cancer/analyse/cluster_marker_for_epi_0214_cleaned5_marker_20210227.RData')

inputparam <- commandArgs()

expression_file <- inputparam[6]
clinical_survival_file <- inputparam[7]
probe_file <- inputparam[8]
names <- inputparam[9]


id.annota <- read.delim(probe_file, check.names=F)
id.annota <- id.annota[, 1:2]
TCGAexpr <- readr::read_tsv(expression_file)
TCGAexpr <- as.data.frame(TCGAexpr)
TCGAexpr <- merge(id.annota, TCGAexpr, by.x="id", by.y="Ensembl_ID")
TCGAexpr <- aggregate(TCGAexpr[,-c(1,2)], list(TCGAexpr$gene), FUN=sum)
TCGAexpr <- column_to_rownames(TCGAexpr,"Group.1")
TCGAexpr <- TCGAexpr[rowSums(TCGAexpr)>10, ] %>% as.matrix()
TCGAclin <- read.delim(gzfile(clinical_survival_file), check.names = F)
TCGAclin <- TCGAclin[TCGAclin$sample %in% colnames(TCGAexpr),]
TCGAclin <- TCGAclin[,-3]






genelistall <- list()

genelistall <- lapply(
  unique(epi_0214_cleaned5_marker$cluster),function(x){
    epi_0214_cleaned5_marker$gene[epi_0214_cleaned5_marker$cluster %in% x & epi_0214_cleaned5_marker$p_val_adj < 0.01 & epi_0214_cleaned5_marker$avg_logFC>1 & epi_0214_cleaned5_marker$pct.1 > 0.3]
  }
)
names(genelistall) <- unique(epi_0214_cleaned5_marker$cluster)


customLayout::lay_new(matrix(c(1,2),ncol=1),heights = c(1,1)) -> pleft
customLayout::lay_new(matrix(c(1,2),ncol=2),widths = c(1,2)) -> newfig
customLayout::lay_split_field(newfig,pleft,1) -> newfig

TCGAclin_merged_list <- list()
res.cox.list <- list()
pdf(paste0('Survival.',names,'.pdf'),height=4,width=8)
sink(paste0('Survival.',names,'.output.txt'))
for (genesetid in names(genelistall) ){
  cat ("Working at",genesetid,'\n')
  tryCatch({
    candidateset <- genelistall[[genesetid]]
    TFHset = list(candidateset)
    names(TFHset) <- 'TFHset'
    ES <- gsva(TCGAexpr, TFHset)
    ES <- data.frame(t(ES), check.names=F)
    TCGAclin_merged <- merge(TCGAclin, ES, by.x="sample", by.y=0)
    TCGAclin_merged$TFH.set = ifelse(TCGAclin_merged$TFHset>median(TCGAclin_merged$TFHset),'high','low')
    ggplot(TCGAclin_merged,aes(x=TFHset)) + geom_density() + theme_bw() -> p1
    ggplot(TCGAclin_merged,aes(x=OS.time,y=TFHset)) + geom_density2d(color='black',linetype='dashed') + geom_point(alpha=0.05) + theme_bw() -> p2
    ggsurvplot(survfit(Surv(OS.time, OS)~TFH.set, data=TCGAclin_merged), conf.int=F, pval=TRUE,palette='jco') + ggtitle(genesetid) -> p3
    print(lay_grid(list(p1,p2,p3$plot),lay=newfig))
    res.cox <- coxph(Surv(OS.time, OS)~TFHset, data=TCGAclin_merged)
    print(res.cox)
    TCGAclin_merged -> TCGAclin_merged_list[[genesetid]]
    res.cox -> res.cox.list[[genesetid]]
  },error=function(e) {cat ('failed to work on ',genesetid,'\n')})
  
}
dev.off()

save(list=c('TCGAclin_merged_list','res.cox.list'),file=paste0('Survival.',names,'.results.Rdata'))
