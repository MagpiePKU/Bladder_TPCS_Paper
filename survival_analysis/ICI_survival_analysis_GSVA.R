library(IMvigor210CoreBiologies)
library(DESeq2)
library(limma)
library(edgeR)
library(DT)
library(ComplexHeatmap)
library(corrplot)
library(survival)
library(GSVA)
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
library(glmnet)
library(pheatmap)
library(ArchR)

rm(list=ls())

dir.create('/Users/wing/Desktop/Eulerian/Results/20210211_BLCA_Single_Cell_Paper_Writting_Preparation/Figures/Figure6/Figure6A')
setwd('/Users/wing/Desktop/Eulerian/Results/20210211_BLCA_Single_Cell_Paper_Writting_Preparation/Figures/Figure6/Figure6A')

data(human_gene_signatures)
ind_genes <- human_gene_signatures

# colors
data(color_palettes)

# variables
irf <- "Best Confirmed Overall Response"
ml <- "FMOne mutation burden per MB"
goi <- names(ind_genes)
tablesDir <- system.file("tables", 
                         package="IMvigor210CoreBiologies")

# font sizes
labCex <- 0.9
namesCex <- 0.9
legendCex <- 0.9
titleCex <- 1
axisCex <- 0.9
titleF <- 1

data(cds)  
cds2 <- cds
data(fmone)  
fmi <- fmone

# normalize 
geneNames <- setNames(fData(cds2)$Symbol, 
                      as.character(rownames(fData(cds2))))
voomD <- filterNvoom(counts(cds2),
                     minSamples=ncol(counts(cds2))/10,
                     minCpm=0.25)
m <- voomD$E
m <- t(scale( t( m ),
              center=TRUE, 
              scale=TRUE)
)
# add signature scores to pData()
m2 <- m
rownames(m2) <- geneNames[rownames(m2)]

## Gene expression is in `m2`. ZY

load('/Users/wing/Desktop/Eulerian/Results/20210211_BLCA_Single_Cell_Paper_Writting_Preparation/Data/cluster_marker_for_epi_0214_cleaned5_marker_20210227.RData')

## Now pull the OS 

genelistall <- lapply(
  unique(epi_0214_cleaned5_marker$cluster),function(x){
    epi_0214_cleaned5_marker$gene[epi_0214_cleaned5_marker$cluster %in% x & epi_0214_cleaned5_marker$p_val_adj < 0.01 & epi_0214_cleaned5_marker$avg_logFC>1 & epi_0214_cleaned5_marker$pct.1 > 0.3]
  }
)
names(genelistall) <- unique(epi_0214_cleaned5_marker$cluster)



## Summarise clinical data

sampleNames(cds@phenoData) -> name_vector
cds$`Best Confirmed Overall Response` -> response_vector
cds$`Immune phenotype` -> immune_phenotype_vector
cds$os -> os_vector
cds$censOS -> os_cens_vector
clinical_data <- as.data.frame(cbind(name_vector,as.character(response_vector),as.character(immune_phenotype_vector),as.numeric(as.character(os_vector)),as.numeric(as.character(os_cens_vector))))
colnames(clinical_data) <- c('sample','best_response','immune_phenotype','OS.time','OS')

clinical_data$OS.time <- as.numeric(as.character(clinical_data$OS.time))
clinical_data$OS <- as.numeric(as.character(clinical_data$OS))

## Premake figure

customLayout::lay_new(matrix(c(1,2),ncol=1),heights = c(1,1)) -> pleft
customLayout::lay_new(matrix(c(1,2),ncol=2),widths = c(1,2)) -> newfig
customLayout::lay_split_field(newfig,pleft,1) -> newfig

## Do GSVA

res.cox.list <- list() 
resultlist <- list() 
for (genesetid in unique(names(genelistall))){
  print(genesetid)
  candidateset <- genelistall[[genesetid]]
  TFHset = list(candidateset)
  names(TFHset) <- 'TFHset'
  ES <- gsva(m2, TFHset)
  ES <- data.frame(t(ES), check.names=F)
  TCGAclin_merged <- merge(clinical_data, ES, by.x="sample", by.y=0)
  TCGAclin_merged$TFH.set = ifelse(TCGAclin_merged$TFHset>median(TCGAclin_merged$TFHset),'high','low')
  ggplot(TCGAclin_merged,aes(x=TFHset)) + geom_density() + theme_bw() + facet_wrap(~best_response) -> p1
  ggplot(TCGAclin_merged,aes(x=TFHset)) + geom_density() + theme_bw() + facet_wrap(~immune_phenotype) -> p2
  ggsurvplot(survfit(Surv(OS.time, OS)~TFH.set, data=TCGAclin_merged), conf.int=F, pval=TRUE,palette='jco') + ggtitle(genesetid) -> p3
  pdf(paste0('Survival.',genesetid,'.pdf'),height=4,width=8)
  print(lay_grid(list(p1,p2,p3$plot),lay=newfig))
  dev.off() 
  res.cox <- coxph(Surv(OS.time, OS)~TFHset, data=TCGAclin_merged)
  print(res.cox)
  res.cox.list[[genesetid]] <- res.cox
  TCGAclin_merged -> resultlist[[genesetid]]
}

resultlist[["11"]] -> TCGAclin_merged

TCGAclin_merged$immune_phenotype <- as.character(TCGAclin_merged$immune_phenotype)
TCGAclin_merged$immune_phenotype[is.na(TCGAclin_merged$immune_phenotype)] <- 'unknown'
plist <- list() 
for ( x in as.character(unique(TCGAclin_merged$immune_phenotype))){
  ggsurvplot(survfit(Surv(OS.time, OS)~TFH.set, data=TCGAclin_merged[TCGAclin_merged$immune_phenotype %in% x,]), conf.int=F, pval=TRUE,palette='jco') + ggtitle(x) -> plist[[x]]
}

load('TCGA/Survival.TCGA-BLCA.results.Rdata') 
TCGAclin_merged_list[['11']] -> real_TCGA_data

ggsurvplot(survfit(Surv(OS.time, OS)~TFH.set, data=TCGAclin_merged[TCGAclin_merged$immune_phenotype %ni% 'desert',]), conf.int=F, pval=TRUE,palette='jco') + ggtitle('Not Immune Desert') -> p1
ggsurvplot(survfit(Surv(OS.time, OS)~TFH.set, data=TCGAclin_merged[TCGAclin_merged$immune_phenotype %in% 'desert',]), conf.int=F, pval=TRUE,palette='jco') + ggtitle('Immune Desert') -> p2
ggsurvplot(survfit(Surv(OS.time/30, OS)~TFH.set, data=real_TCGA_data), conf.int=F, pval=TRUE,palette='jco') + ggtitle('Real TCGA no IO') -> p3

pdf('Figure6A.Cluster.11.TPCS.Survival.in.IO.and.TCGA.for.paper.pdf',height=6,width=8)
p1
p2
p3
dev.off() 

pdf('Figure6A.Cluster.11.TPCS.Activity.Distribution.per.immune.infiltration.class.in.IO.pdf',height=3,width=3.5)
table(TCGAclin_merged$immune_phenotype,TCGAclin_merged$TFH.set) -> tb1
tb1 <- tb1/rowSums(tb1)
tb1df <- melt(tb1)
ggplot(tb1df[tb1df$Var1 %ni% 'unknown',],aes(x=Var1,y=value,fill=Var2)) + geom_col() + theme_classic() + scale_fill_manual(values = c('red','cornflowerblue')) + xlab('Immune Phenotype') + ylab('Cluster 11 Activity Class')
dev.off() 

table(TCGAclin_merged[TCGAclin_merged$immune_phenotype %ni% 'unknown',]$immune_phenotype,TCGAclin_merged[TCGAclin_merged$immune_phenotype %ni% 'unknown',]$TFH.set) -> tb2


# 
# lapply(names(resultlist),function(x){
#   y <- resultlist[[x]]
#   colnames(y)[6:7] <- c('expression','classification')
#   y$set <- x
#   y
# }) %>% bind_rows -> resultdf
# resultmat <- dcast(resultdf,sample+best_response+immune_phenotype+OS.time+OS ~ set, value.var='expression',fun.aggregate = sum)
# rownames(resultmat) <- resultmat$sample
# resultmat$DC <- ifelse(resultmat$best_response %in% c("PD",'NE'),'0','1')
# resultmat$OR <- ifelse(resultmat$best_response %in% c("PD",'NE','SD'),'0','1')
# drawmap <- resultmat[,colnames(resultmat) %ni% c('sample','best_response','immune_phenotype','OS.time','OS','DC','OR')]
# 
# pheatmap::pheatmap(drawmap[resultmat$immune_phenotype%in%'inflamed',],annotation_row = resultmat[resultmat$immune_phenotype%in%'inflamed',] %>% dplyr::select(best_response,immune_phenotype,OS.time,OS,DC,OR),clustering_method = 'ward.D2',show_rownames = F,scale='row')
# 
# selected_set <- c('Basal','mNK_negative','NK3','Treg_eff','TPCS','FPS')
# 
# 
# 
# resultmat2 <- resultmat
# # resultmat2$mNK <- -1 * resultmat2$mNK
# resultmat2$death = ifelse(resultmat2$OS==1,0,1)
# 
# coxph(Surv(OS.time, OS)~mNK_negative+TRTm+Treg_eff+Basal, data=resultmat2)
# # 
# # ggsurvplot(survfit(Surv(OS.time, OS)~mNK_negative+TRTm+Treg_eff+Basal, data=TCGAclin_merged), conf.int=F, pval=TRUE,palette='jco') + ggtitle(genesetid) -> p3
# 
# resultmat2$PRCR <- ifelse(resultmat2$best_response %in% c('PR','CR'),1,0)
# 
# annotation_colors <- list(
#   'DC' = c('0'='white','1'='black'),
#   'OR' = c('0'='white','1'='black'),
#   'death' = c('0'='white','1'='black'),
#   'OS.time' = c('0'='white','1'='black'),
#   'immune_phenotype' = c('desert'='white','excluded'='black','inflamed'='red'),
#   'best_response' = c('PD'='black','NE'='black','SD'='gold1','PR'='pink','CR'='red'),
#   'PRCR' = c('0'='white','1'='red'),
#   'Intravesical_BCG_administered' = c('N'='white','Y'='black'),
#   'Received_platinum' = c('N'='white','Y'='black'),
#   "TC_Level" = c('NE'='white','TC0'='gray','TC1'='pink','TC2'='red'),
#   "IC_Level" = c('NE'='white','IC0'='gray','IC1'='pink','IC2'='red'),
#   'TCGA_Subtype' = c('I'='gold1','II'='orange','III'='red','IV'='blueviolet')
# )
# pheatmap::pheatmap(resultmat2[,selected_set],annotation_row = resultmat2 %>% dplyr::select(best_response,immune_phenotype,OS.time,OR,PRCR),clustering_method = 'ward.D2',show_rownames = F,scale='row',annotation_colors = annotation_colors)
# 
# openxlsx::read.xlsx('~/Downloads/thnov10p7002s2.xlsx',sheet='S1use') -> t1
# head(resultmat)
# resultmat$ID <- resultmat$sample
# resultmat3 <- left_join(resultmat,t1)
# rownames(resultmat3) <- resultmat3$ID
# rownames(resultmat2) <- rownames(resultmat)
# 
# resultmat3$scaled_TMB <- scale(resultmat3$FMOne_mutation_burden_per_MB)
# resultmat3$scaled_TNB <- scale(resultmat3$Neoantigen_burden_per_MB)
# 
# selected_set <- c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS')
# 
# resultmat3[,selected_set] -> mat1
# mat1 <- scale(mat1)
# mat1[mat1>1] <- 1
# mat1[mat1< -1] <- -1
# pdf('Heatmap.correlate.pdf',height=40,width=20)
# pheatmap::pheatmap(t(mat1),annotation_col = resultmat3 %>% dplyr::select(best_response,immune_phenotype,OS.time,OR,TCGA_Subtype,Lund,TC_Level,IC_Level,scaled_TMB,scaled_TNB,TumorPurity_estimate,Intravesical_BCG_administered,Received_platinum),clustering_method = 'ward.D2',show_rownames = T,show_colnames = F,scale='row',annotation_colors = annotation_colors,cellwidth = 2, cellheight = 15)
# pheatmap::pheatmap(t(mat1),annotation_col = resultmat3 %>% dplyr::select(best_response,immune_phenotype,OS.time,OR,TCGA_Subtype,Lund,TC_Level,IC_Level,scaled_TMB,scaled_TNB,TumorPurity_estimate,Intravesical_BCG_administered,Received_platinum),clustering_method = 'ward.D2',show_rownames = T,show_colnames = F,scale='row',annotation_colors = annotation_colors,cellwidth = 1, cellheight = 15)
# pheatmap::pheatmap(t(mat1[resultmat3$immune_phenotype %in% c('excluded','inflamed'),]),annotation_col = resultmat3[resultmat3$immune_phenotype %in% c('excluded','inflamed'),] %>% dplyr::select(best_response,immune_phenotype,OS.time,OR,TCGA_Subtype,Lund,TC_Level,IC_Level,scaled_TMB,scaled_TNB,TumorPurity_estimate,Intravesical_BCG_administered,Received_platinum),clustering_method = 'ward.D2',show_rownames = T,show_colnames = F,scale='row',annotation_colors = annotation_colors,cellwidth = 1, cellheight = 15)
# dev.off()
# 
# cor_test_list <- list()
# for (x in c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS','FMOne_mutation_burden_per_MB','Neoantigen_burden_per_MB')){
#   cor.test(resultmat3[,x],as.numeric(resultmat3$best_response %in% c("PR",'CR','SD'))) -> cor_test_list[[x]]
# } 
# cor_df <- data.frame(
#   'pvalue' = lapply(names(cor_test_list),function(x) cor_test_list[[x]]$p.value) %>% unlist(),
#   'estimate' = lapply(names(cor_test_list),function(x) cor_test_list[[x]]$estimate) %>% unlist(),
#   'feature' = names(cor_test_list)
# )
# 
# cor_df <- arrange(cor_df,-pvalue)
# cor_df$feature <- as.character(cor_df$feature)
# cor_df$feature[grepl("M1",cor_df$feature)] <- 'M1_related'
# cor_df$feature[grepl("mutation_burden",cor_df$feature)] <- 'TMB'
# cor_df$feature[grepl("Neoantigen_burden",cor_df$feature)] <- 'TNB'
# cor_df$feature <- factor(cor_df$feature,levels=cor_df$feature,ordered=T)
# pdf('cor.test.pval.pdf',height=4,width=4)
# ggplot(cor_df,aes(x=feature,y=-log(pvalue))) + geom_col() + coord_flip() + theme_classic()
# dev.off() 
# 
# cor_test_list1 <- list()
# for (x in c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS','FMOne_mutation_burden_per_MB','Neoantigen_burden_per_MB')){
#   cor.test(resultmat3[resultmat3$immune_phenotype %in% c('inflamed'),x],as.numeric(resultmat3[resultmat3$immune_phenotype %in% c('inflamed'),]$best_response %in% c("PR",'CR','SD'))) -> cor_test_list1[[x]]
# } 
# cor_df1 <- data.frame(
#   'pvalue' = lapply(names(cor_test_list1),function(x) cor_test_list1[[x]]$p.value) %>% unlist(),
#   'estimate' = lapply(names(cor_test_list1),function(x) cor_test_list1[[x]]$estimate) %>% unlist(),
#   'feature' = names(cor_test_list1)
# )
# 
# cor_df1 <- arrange(cor_df1,-pvalue)
# cor_df1$feature <- as.character(cor_df1$feature)
# cor_df1$feature[grepl("M1",cor_df1$feature)] <- 'M1_related'
# cor_df1$feature[grepl("mutation_burden",cor_df1$feature)] <- 'TMB'
# cor_df1$feature[grepl("Neoantigen_burden",cor_df1$feature)] <- 'TNB'
# cor_df1$feature <- factor(cor_df1$feature,levels=cor_df1$feature,ordered=T)
# pdf('cor.test.pval.Inflamed.pdf',height=4,width=4)
# ggplot(cor_df1,aes(x=feature,y=-log(pvalue))) + geom_col() + coord_flip() + theme_classic()
# dev.off() 
# 
# cor_test_list2 <- list()
# for (x in c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS','FMOne_mutation_burden_per_MB','Neoantigen_burden_per_MB')){
#   cor.test(resultmat3[resultmat3$immune_phenotype %in% c('excluded'),x],as.numeric(resultmat3[resultmat3$immune_phenotype %in% c('excluded'),]$best_response %in% c("PR",'CR','SD'))) -> cor_test_list2[[x]]
# } 
# cor_df2 <- data.frame(
#   'pvalue' = lapply(names(cor_test_list2),function(x) cor_test_list2[[x]]$p.value) %>% unlist(),
#   'estimate' = lapply(names(cor_test_list2),function(x) cor_test_list2[[x]]$estimate) %>% unlist(),
#   'feature' = names(cor_test_list2)
# )
# 
# cor_df2 <- arrange(cor_df2,-pvalue)
# cor_df2$feature <- as.character(cor_df2$feature)
# cor_df2$feature[grepl("M1",cor_df2$feature)] <- 'M1_related'
# cor_df2$feature[grepl("mutation_burden",cor_df2$feature)] <- 'TMB'
# cor_df2$feature[grepl("Neoantigen_burden",cor_df2$feature)] <- 'TNB'
# cor_df2$feature <- factor(cor_df2$feature,levels=cor_df2$feature,ordered=T)
# pdf('cor.test.pval.Excluded.pdf',height=4,width=4)
# ggplot(cor_df2,aes(x=feature,y=-log(pvalue))) + geom_col() + coord_flip() + theme_classic()
# dev.off() 
# 
# 
# cor_test_list3 <- list()
# for (x in c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS','FMOne_mutation_burden_per_MB','Neoantigen_burden_per_MB')){
#   cor.test(resultmat3[resultmat3$immune_phenotype %in% c('desert'),x],as.numeric(resultmat3[resultmat3$immune_phenotype %in% c('desert'),]$best_response %in% c("PR",'CR','SD'))) -> cor_test_list3[[x]]
# } 
# cor_df3 <- data.frame(
#   'pvalue' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$p.value) %>% unlist(),
#   'estimate' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$estimate) %>% unlist(),
#   'feature' = names(cor_test_list3)
# )
# 
# cor_df3 <- arrange(cor_df3,-pvalue)
# cor_df3$feature <- as.character(cor_df3$feature)
# cor_df3$feature[grepl("M1",cor_df3$feature)] <- 'M1_related'
# cor_df3$feature[grepl("mutation_burden",cor_df3$feature)] <- 'TMB'
# cor_df3$feature[grepl("Neoantigen_burden",cor_df3$feature)] <- 'TNB'
# cor_df3$feature <- factor(cor_df3$feature,levels=cor_df3$feature,ordered=T)
# pdf('cor.test.pval.Desert.pdf',height=4,width=4)
# ggplot(cor_df3,aes(x=feature,y=-log(pvalue))) + geom_col() + coord_flip() + theme_classic()
# dev.off() 
# 
# cor_df1$origin <- 'inflamed'
# cor_df2$origin <- 'excluded'
# cor_df3$origin <- 'desert'
# 
# cor_df4 <- rbind(cor_df1,cor_df2,cor_df3)
# 
# pdf('cor.test.pval.all.pdf',height=12,width=9)
# ggplot(cor_df4,aes(x=feature,y=-log(pvalue),fill=origin)) + geom_col(position='dodge') + coord_flip() + theme_classic() + ggthemes::scale_fill_tableau() + theme(legend.position = 'bottom',axis.text.y = element_blank(),axis.title.y = element_blank()) -> p1
# ggplot(cor_df4,aes(x=feature,y=estimate,fill=origin)) + geom_col(position='dodge') + coord_flip() + theme_classic() + ggthemes::scale_fill_tableau()  + theme(legend.position = 'bottom') -> p2
# library(patchwork)
# p2+p1
# dev.off()
# 
# pdf("Survival.M1.scale.pdf",height=4,width=6)
# ggsurvplot(survfit(Surv(OS.time, OS)~Macrophage_M1_binary, data=resultmat3), conf.int=F, pval=TRUE,palette='jco') + ggtitle('M1') 
# dev.off() 
# 
# pheatmap::pheatmap(drawmap[resultmat3$immune_phenotype%in%'excluded',],annotation_row = resultmat3[resultmat3$immune_phenotype%in%'excluded',] %>% dplyr::select(best_response,immune_phenotype,OS.time,OS,IC_Level,TC_Level),clustering_method = 'ward.D2',show_rownames = F,scale='none',annotation_colors = annotation_colors)
# 
# pheatmap::pheatmap(drawmap,annotation_row = resultmat3%>% dplyr::select(best_response,immune_phenotype,OS.time,OS,IC_Level,TC_Level),clustering_method = 'ward.D2',show_rownames = F,scale='none',annotation_colors = annotation_colors)
# 
# resultmat3$IC_score <- 0 
# resultmat3$IC_score[resultmat3$IC_Level%in%'IC1'] <- 1
# resultmat3$IC_score[resultmat3$IC_Level%in%'IC2'] <- 2
# cor_test_list3 <- list()
# for (x in c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS','FMOne_mutation_burden_per_MB','Neoantigen_burden_per_MB','B_n','B_switch','B_tumor_specific','TRTm','TRTl','TRTm_specific','rNK','Treg10','Treg_5','tCD8_1','tCD8_2','Treg_cm','Treg_naive','B_preswitch')){
#   cor.test(resultmat3[,x],as.numeric(resultmat3$IC_score)) -> cor_test_list3[[x]]
# } 
# cor_df3 <- data.frame(
#   'pvalue' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$p.value) %>% unlist(),
#   'estimate' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$estimate) %>% unlist(),
#   'feature' = names(cor_test_list3)
# )
# 
# resultmat3$IC_score <- 0 
# resultmat3$IC_score[resultmat3$IC_Level%in%'IC1'] <- 1
# resultmat3$IC_score[resultmat3$IC_Level%in%'IC2'] <- 2
# cor_test_list3 <- list()
# for (x in c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS','FMOne_mutation_burden_per_MB','Neoantigen_burden_per_MB','B_n','B_switch','B_tumor_specific','TRTm','TRTl','TRTm_specific','rNK','Treg10','Treg_5','tCD8_1','tCD8_2','Treg_cm','Treg_naive','B_preswitch')){
#   cor.test(resultmat3[,x],as.numeric(resultmat3$IC_score)) -> cor_test_list3[[x]]
# } 
# cor_df3 <- data.frame(
#   'pvalue' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$p.value) %>% unlist(),
#   'estimate' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$estimate) %>% unlist(),
#   'feature' = names(cor_test_list3)
# )
# cor_df3 <- arrange(cor_df3,estimate,-pvalue)
# cor_df3$feature <- as.character(cor_df3$feature)
# cor_df3$feature[grepl("M1",cor_df3$feature)] <- 'M1_related'
# cor_df3$feature[grepl("mutation_burden",cor_df3$feature)] <- 'TMB'
# cor_df3$feature[grepl("Neoantigen_burden",cor_df3$feature)] <- 'TNB'
# cor_df3$feature <- factor(cor_df3$feature,levels=cor_df3$feature,ordered=T)
# pdf('cor.test.pval.IC_score.pdf',height=4,width=4)
# ggplot(cor_df3,aes(x=feature,y=-log(pvalue))) + geom_col(position='dodge') + coord_flip() + theme_classic() + ggthemes::scale_fill_tableau() + theme(legend.position = 'bottom',axis.text.y = element_blank(),axis.title.y = element_blank()) -> p1
# ggplot(cor_df3,aes(x=feature,y=estimate)) + geom_col(position='dodge') + coord_flip() + theme_classic() + ggthemes::scale_fill_tableau()  + theme(legend.position = 'bottom') -> p2
# p2+p1
# dev.off()
# 
# 
# 
# resultmat3$TC_score <- 0 
# resultmat3$TC_score[resultmat3$TC_Level%in%'TC1'] <- 1
# resultmat3$TC_score[resultmat3$TC_Level%in%'TC2'] <- 2
# cor_test_list3 <- list()
# for (x in c('mNK_negative','Treg_eff','TRTl_specific','NK3','Macrophage_M1_scale','Basal','TPCS','FPS','FMOne_mutation_burden_per_MB','Neoantigen_burden_per_MB','B_n','B_switch','B_tumor_specific','TRTm','TRTl','TRTm_specific','rNK','Treg10','Treg_5','tCD8_1','tCD8_2','Treg_cm','Treg_naive','B_preswitch')){
#   cor.test(resultmat3[,x],as.numeric(resultmat3$TC_score)) -> cor_test_list3[[x]]
# } 
# cor_df3 <- data.frame(
#   'pvalue' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$p.value) %>% unlist(),
#   'estimate' = lapply(names(cor_test_list3),function(x) cor_test_list3[[x]]$estimate) %>% unlist(),
#   'feature' = names(cor_test_list3)
# )
# cor_df3 <- arrange(cor_df3,estimate,-pvalue)
# cor_df3$feature <- as.character(cor_df3$feature)
# cor_df3$feature[grepl("M1",cor_df3$feature)] <- 'M1_related'
# cor_df3$feature[grepl("mutation_burden",cor_df3$feature)] <- 'TMB'
# cor_df3$feature[grepl("Neoantigen_burden",cor_df3$feature)] <- 'TNB'
# cor_df3$feature <- factor(cor_df3$feature,levels=cor_df3$feature,ordered=T)
# pdf('cor.test.pval.TC_score.pdf',height=4,width=4)
# ggplot(cor_df3,aes(x=feature,y=-log(pvalue))) + geom_col(position='dodge') + coord_flip() + theme_classic() + ggthemes::scale_fill_tableau() + theme(legend.position = 'bottom',axis.text.y = element_blank(),axis.title.y = element_blank()) -> p1
# ggplot(cor_df3,aes(x=feature,y=estimate)) + geom_col(position='dodge') + coord_flip() + theme_classic() + ggthemes::scale_fill_tableau()  + theme(legend.position = 'bottom') -> p2
# p2+p1
# dev.off() 