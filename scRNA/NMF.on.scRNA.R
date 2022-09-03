## Figure 4
## Liger NMF
suppressPackageStartupMessages(library(liger))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(getopt))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
options(future.globals.maxSize=2000000000)
library(reshape2)
library(ArchR)

# Use Liger to perform iterative NMF on samples, to get NMF-decomposed metagenes

runOriginalLiger <- function (object, k, assay = NA, split.by = "orig.ident", lambda = 5,
                              thresh = 1e-06, max.iters = 30, reduction.name = "iNMF_raw",
                              reduction.key = "riNMF_", nrep = 1, H.init = NULL, W.init = NULL,
                              V.init = NULL, rand.seed = 1, print.obj = FALSE, ...)
{
  if(is.na(assay)){
    assay <- DefaultAssay(object = object)
  }
  if (is.character(x = split.by) && length(x = split.by) ==
      1) {
    split.by <- object[[split.by]]
  }
  split.cells <- split(x = colnames(x = object), f = split.by)
  scale.data <- lapply(X = split.cells, FUN = function(x) {
    return(t(x = GetAssayData(object = object, slot = "scale.data",
                              assay = assay)[, x]))
  })
  out <- liger::optimizeALS(object = scale.data, k = k, lambda = lambda,
                            thresh = thresh, max.iters = max.iters, nrep = nrep,
                            H.init = H.init, W.init = W.init, V.init = V.init, rand.seed = rand.seed,
                            print.obj = print.obj)
  colnames(x = out$W) <- VariableFeatures(object = object)
  return(out)
}

## Get Args

dir.create(paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/'))
dir.create(paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/'))
setwd(paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/'))
load('/gpfs/output/20201016_Loom_files_for_Bladder_Cancer/analyse/results/plot_for_paper/saved_RData/total_cell_for_paper_with_annotation_20210419.RData')

load('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/Figure7/Figure7D/Figure7D.new_combined_obj_sub_tissue.seurat.obj.for.cellchat.Rdata')

load('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/Figure4/Figure4F/Figure4F.early_clusters_epi.Rdata')
load('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/Figure1/Figure1B/Figure1B.DCD.normal.bladder.epithelial.cell.data.and.marker.Rdata')


table(total_cell_for_paper$annotation) %>% as.data.frame -> df

total_cell_for_paper@meta.data %>% as.data.frame -> df_target
new_combined_obj_sub_tissue@meta.data %>% as.data.frame -> df_anno
# df_target <- left_join(df_target,df_anno %>% dplyr::select(name_for_cellchat,cell))
# table(df_target$annotation,df_target$name_for_cellchat) -> tb1
# tb1%>% as.data.frame -> df
# table(total_cell_for_paper$annotation) %>% as.data.frame -> df
# write.table(df,file='prepare.celltype.tsv',quote = F,sep = '\t',row.names = F,col.names = T)

read.delim('prepare.celltype.zy.tsv') -> newanno
colnames(newanno) <- c('annotation','freq','zy_annotation')
df_target <- left_join(df_target,newanno %>% dplyr::select(-freq))
df_target$zy_annotation <- as.character(df_target$zy_annotation)
df_target$zy_annotation[df_target$cell %in% df_anno$cell[grepl('arterial',df_anno$name_for_cellchat)] & df_target$zy_annotation %in% 'Endo'] <- 'Endo_Arterial'
df_target$zy_annotation[df_target$cell %in% df_anno$cell[grepl('venous',df_anno$name_for_cellchat)] & df_target$zy_annotation %in% 'Endo'] <- 'Endo_Venous'
df_target$zy_annotation[df_target$cell %in% df_anno$cell[grepl('lymphatic',df_anno$name_for_cellchat)] & df_target$zy_annotation %in% 'Endo'] <- 'Endo_Lymph'
df_target$zy_annotation[df_target$cell %in% df_anno$cell[grepl('capillary',df_anno$name_for_cellchat)] & df_target$zy_annotation %in% 'Endo'] <- 'Endo_Capillary'
df_target$zy_annotation[df_target$cell %in% df_anno$cell[grepl('^Progenitor$',df_anno$name_for_cellchat)]] <- 'Urothelial_Progenitor'
df_target$zy_annotation[df_target$cell %in% dat_dcd_epi@meta.data$cell[grepl('Basal',dat_dcd_epi$cell_type_annotation_details)]] <- 'Urothelial_Basal'
df_target$zy_annotation[df_target$cell %in% dat_dcd_epi@meta.data$cell[grepl('Intermediate|Umbrella',dat_dcd_epi$cell_type_annotation_details)]] <- 'Urothelial_Im/Um'
df_target$zy_annotation[df_target$cell %in% dat_dcd_epi@meta.data$cell[grepl('progenitor',dat_dcd_epi$cell_type_annotation_details)]] <- 'Urothelial_Progenitor'


total_cell_for_paper@meta.data$final_annotation_zy <- df_target$zy_annotation


total_cell_for_paper@meta.data$general_annotation_zy <- gsub('_.+','',total_cell_for_paper@meta.data$final_annotation_zy)
total_cell_for_paper@meta.data$general_annotation_zy[grepl('Mphi|cDC|Mono',total_cell_for_paper@meta.data$final_annotation_zy)] <- 'Myeloid'
total_cell_for_paper@meta.data$general_annotation_zy[grepl('^Urothelial|^Ca',total_cell_for_paper@meta.data$final_annotation_zy)] <- 'Epithelial'
total_cell_for_paper@meta.data$final_annotation_zy[grepl('^Ca_Basal',total_cell_for_paper@meta.data$final_annotation_zy) &total_cell_for_paper@meta.data$tissue_type %in% 'normal_bladder_tissue' ] <- 'Urothelial_Basal_like_remove'
total_cell_for_paper@meta.data$final_annotation_zy[grepl('^Ca_TPCS',total_cell_for_paper@meta.data$final_annotation_zy) &total_cell_for_paper@meta.data$tissue_type %in% 'normal_bladder_tissue' ] <- 'Urothelial_TPCS_like_remove'

total_cell_for_paper <- subset(total_cell_for_paper,subset = final_annotation_zy %ni% c('Urothelial_Basal_like_remove','Urothelial_TPCS_like_remove'))

dcd_obj <- subset(total_cell_for_paper,subset=tissue_type %in% 'normal_bladder_tissue')
dcd_obj <- NormalizeData(dcd_obj) %>%  FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(reduction="pca",dims = 1:10,n.neighbors = 10,min.dist = 0.1) %>% RunTSNE(reduction = "pca",dims = 1:10)

pdf('Figure4.Supp.DCD.UMAP.from.total.obj.pdf',width=10,height=5)
DimPlot(dcd_obj,group.by = 'general_annotation_zy',label = T,reduction = 'umap',pt.size = 0.1) 
DimPlot(dcd_obj,group.by = 'general_annotation_zy',label = T,reduction = 'tsne',pt.size = 0.1)
DimPlot(subset(dcd_obj,subset=general_annotation_zy %in% 'Epithelial'),group.by = 'final_annotation_zy',label = T,reduction = 'umap',pt.size = 0.1) 
DimPlot(subset(dcd_obj,subset=general_annotation_zy %in% 'Epithelial'),group.by = 'final_annotation_zy',label = F,reduction = 'umap',pt.size = 0.1) + facet_wrap(~final_annotation_zy)
DimPlot(subset(dcd_obj,subset=general_annotation_zy %in% 'Epithelial'),group.by = 'final_annotation_zy',label = T,reduction = 'tsne',pt.size = 0.1) 
DimPlot(subset(dcd_obj,subset=general_annotation_zy %in% 'Epithelial'),group.by = 'final_annotation_zy',label = F,reduction = 'tsne',pt.size = 0.1) + facet_wrap(~final_annotation_zy)
dev.off()

save(dcd_obj,file='Figure2.dcd_obj.full.cells.Rdata')




save(total_cell_for_paper,file='/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/RNA.full.Data.20210421.ZY.Rdata')

cancer_tissue_cell_for_paper <- subset(total_cell_for_paper,subset = tissue_type %in% c('cancer_tissue','cancer_adjacent_tissue','metastatic_prostate_tissue','DCD','normal_bladder_tissue'))

table(cancer_tissue_cell_for_paper$final_annotation_zy,cancer_tissue_cell_for_paper$general_annotation_zy)

table(cancer_tissue_cell_for_paper@meta.data[cancer_tissue_cell_for_paper$tissue_type %in% 'normal_bladder_tissue',]$final_annotation_zy)

typelist <- as.character(unique(cancer_tissue_cell_for_paper@meta.data$general_annotation_zy))

for (celltypeid in typelist){
  cat (celltypeid, '\n')
  tryCatch({
    merged_dataset <- subset(cancer_tissue_cell_for_paper,subset = general_annotation_zy %in% celltypeid)
    merged_dataset <- FindVariableFeatures(merged_dataset)
    merged_dataset <- ScaleData(merged_dataset, split.by = "orig.ident", do.center = FALSE)
    k_experience = 20
    table(merged_dataset$orig.ident) -> tb1
    names(tb1[tb1>k_experience]) -> names_to_retain # use only samples > k_experience cells otherwise would hinder the analysis. 
    # names_to_retain
    merged_dataset_samples_good <- subset(merged_dataset,subset=orig.ident %in% names_to_retain)
    # perform liger iNMF
    liger_output_for_merged_dataset <- runOriginalLiger(merged_dataset_samples_good, k = k_experience, lambda = 5, split.by = "orig.ident")
    colnames(liger_output_for_merged_dataset$W) <- rownames(GetAssayData(object = merged_dataset_samples_good,assay='RNA',slot = "scale.data"))
    save(liger_output_for_merged_dataset,file=paste0('Figure4.',celltypeid,'.NMF.original.output.Rdata'))
#   },error=function(e){
#     cat ('Failed\n')
#   })
# }
# 
# k_experience = 20
# for (celltypeid in typelist){
#   cat (celltypeid, '\n')
#   tryCatch({
#     merged_dataset <- subset(cancer_tissue_cell_for_paper,subset = general_annotation_zy %in% celltypeid)
#     load(file=paste0('Figure4.',celltypeid,'.NMF.original.output.Rdata'))
    # k_experience = length(liger_output_for_merged_dataset$H)
    lapply(c(1:k_experience),function(x){
      as.data.frame(liger_output_for_merged_dataset$H[[x]])
    })  %>%  bind_rows() -> liger_output_for_merged_dataset_binded_H_loading_mtx
    colnames(liger_output_for_merged_dataset_binded_H_loading_mtx) <- paste0('NMF_',c(1:ncol(liger_output_for_merged_dataset_binded_H_loading_mtx)))
    liger_output_for_merged_dataset_binded_H_loading_mtx$cell <- rownames(liger_output_for_merged_dataset_binded_H_loading_mtx)
    merged_dataset@meta.data %>% as.data.frame->temp
    temp$cell <- rownames(temp)
    left_join(temp,liger_output_for_merged_dataset_binded_H_loading_mtx,by='cell') -> all_liger
    rownames(all_liger) <- c(1:nrow(all_liger))
    save(all_liger,file=paste0('Figure4.',celltypeid,'.NMF.combined.Rdata'))
  },error=function(e){
    cat ('Failed\n')
  })
}


typelist2 <- c('B','CD8','Myeloid','NK')

for (celltypeid in typelist){
  cat (celltypeid, '\n')
  tryCatch({
    merged_dataset <- subset(cancer_tissue_cell_for_paper,subset = general_annotation_zy %in% celltypeid)
    k_experience = 20
    table(merged_dataset$orig.ident) -> tb1
    names(tb1[tb1>k_experience]) -> names_to_retain # use only samples > k_experience cells otherwise would hinder the analysis. 
    # names_to_retain
    merged_dataset_samples_good <- subset(merged_dataset,subset=orig.ident %in% names_to_retain)
    # perform liger iNMF
    load(file=paste0('Figure4.',celltypeid,'.NMF.original.output.Rdata'))
    lapply(c(1:length(liger_output_for_merged_dataset$H)),function(x){
      as.data.frame(liger_output_for_merged_dataset$H[[x]])
    })  %>%  bind_rows() -> liger_output_for_merged_dataset_binded_H_loading_mtx
    colnames(liger_output_for_merged_dataset_binded_H_loading_mtx) <- paste0('NMF_',c(1:ncol(liger_output_for_merged_dataset_binded_H_loading_mtx)))
    liger_output_for_merged_dataset_binded_H_loading_mtx$cell <- rownames(liger_output_for_merged_dataset_binded_H_loading_mtx)
    merged_dataset_samples_good@meta.data %>% as.data.frame->temp
    temp$cell <- rownames(temp)
    left_join(temp,liger_output_for_merged_dataset_binded_H_loading_mtx,by='cell') -> all_liger
    rownames(all_liger) <- c(1:nrow(all_liger))
    save(all_liger,file=paste0('Figure4.',celltypeid,'.NMF.combined.Rdata'))
  },error=function(e){
    cat ('Failed\n')
  })
}


for (celltypeid in typelist){
  load(file=paste0('Figure4.',celltypeid,'.NMF.original.output.Rdata'))
  load(file=paste0('Figure4.',celltypeid,'.NMF.combined.Rdata'))
  k_experience = nrow(liger_output_for_merged_dataset$W)
  liger_output_for_merged_dataset$W %>% as.data.frame -> df
  df$NMF <- paste0('NMF_',c(1:k_experience))
  df <- melt(df,id.vars='NMF')
  if(celltypeid %in% 'Epithelial'){
    all_liger <- all_liger[grepl("^Ca|^Uro",all_liger$final_annotation_zy),]
  }
  
  pdf(file = paste0('Figure4.',celltypeid,'.NMF.essential.loading.pdf'),height=12,width=(24/5)*k_experience)
  resdf <- all_liger[,c('final_annotation_zy','cell',paste0('NMF_',c(1:k_experience)))] %>% melt(id.vars=c("final_annotation_zy",'cell'))

  lapply(c(paste0('NMF_',c(1:k_experience))),function(x){
    ggplot(resdf[resdf$variable %in% x & !is.na(resdf$value),],aes(x=final_annotation_zy,y=value,fill=final_annotation_zy)) + geom_violin() + geom_boxplot(fill='black',outlier.alpha = 0,width=0.06)  + theme_classic() + theme(legend.position = 'none',axis.title=element_blank(),axis.text.x = element_text(angle=90,hjust=0,vjust=1)) + ggthemes::scale_fill_tableau() 
    # + ggpubr::stat_compare_means(formula = value~final_annotation_zy,method = 'wilcox.test',p.adjust.method = 'none',label='p.signif',hide.ns=T)
    }) %>% patchwork::wrap_plots(nrow=1) -> p1
  lapply(c(paste0('NMF_',c(1:k_experience))),function(x){
    df1 <- df[df$NMF %in% x,]
    df1 <- arrange(df1,value)
    df1$gene <- factor(df1$variable,levels=df1$variable,ordered=T)
    ggplot(df1,aes(x=gene,y=(value))) + geom_point() + ggrepel::geom_label_repel(data=df1[df1$value > quantile(df1$value,seq(0,1,0.001)[993]),],mapping=aes(label=gene,y=value,x=gene),cex=2) + theme_classic() + geom_point(data=df1[df1$value > quantile(df1$value,seq(0,1,0.001)[993]),],color='red') + 
      theme(axis.text = element_blank(),axis.title.y = element_blank()) + xlab(x) 
  }) %>% patchwork::wrap_plots(nrow=1) -> p2
  print(p1/p2 + plot_layout(ncol=1,heights = c(1,3)))
  dev.off()
  cat('Done',celltypeid,'\n')
}

##
## runGsea

runGSEA_mod <- function (object, gene_sets = c(), mat_w = T, mat_v = 0, custom_gene_sets = c())
{
  if (mat_w) {
    gene_loadings <- object$W
    if (mat_v) {
      gene_loadings <- gene_loadings + Reduce("+", lapply(mat_v,
                                                          function(v) {
                                                            object$V[[v]]
                                                          }))
    }
  }
  else {
    gene_loadings <- Reduce("+", lapply(mat_v, function(v) {
      object$V[[v]]
    }))
  }
  gene_ranks <- t(apply(gene_loadings, MARGIN = 1, function(x) {
    rank(x)
  }))
  colnames(gene_ranks) <- sapply(colnames(gene_ranks), toupper)
  gene_id <- as.character(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                                colnames(gene_ranks), "ENTREZID", "SYMBOL"))
  colnames(gene_ranks) <- gene_id
  gene_ranks <- gene_ranks[, !is.na(colnames(gene_ranks))]
  if (inherits((custom_gene_sets)[1], "tbl_df")) {
    pathways <- split(custom_gene_sets, x = custom_gene_sets$entrez_gene,
                      f = custom_gene_sets$gs_name)
    pathways <- lapply(pathways, function(x) {
      as.character(x)
    })
  }
  else if (length(custom_gene_sets)) {
    pathways <- custom_gene_sets
  }
  else {
    pathways <- fgsea::reactomePathways(colnames(gene_ranks))
    if (length(gene_sets)) {
      pathways <- pathways[intersect(gene_sets, names(pathways))]
    }
  }
  gsea <- apply(gene_ranks, MARGIN = 1, function(x) {
    fgsea::fgsea(pathways, x, minSize = 15, maxSize = 500,
                 nperm = 10000)
  })
  gsea <- lapply(gsea, function(x) {
    as.matrix(x[order(x$pval), ])
  })
  return(gsea)
}

lapply(typelist, function(celltypeid){
  load(file=paste0('Figure1S.',celltypeid,'.NMF.original.output.Rdata'))
  load(file=paste0('Figure1S.',celltypeid,'.NMF.combined.Rdata'))
  runGSEA_mod(liger_output_for_merged_dataset) -> liger_output_for_merged_dataset_gsea
  lapply(c(1:length(liger_output_for_merged_dataset_gsea)),function(y){
    liger_output_for_merged_dataset_gsea[[y]] %>% as.data.frame -> x
    x$origin <- celltypeid
    x$nmfid <- paste0("NMF_",y)
    x
  }) 
}) -> GSEA_result 
names(GSEA_result) <- typelist 

## 

library(RcisTarget)
data(motifAnnotations_hgnc)
motifRankings <- importRankings("/gpfs/genomedb/cisTopic/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather") # note not using the long one
# Motif enrichment analysis:
lapply(typelist, function(celltypeid){
  load(file=paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/Figure4.',celltypeid,'.NMF.original.output.Rdata'))
  load(file=paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/Figure4.',celltypeid,'.NMF.combined.Rdata'))
  k_experience <- nrow(liger_output_for_merged_dataset$W)
  lapply(c(1:k_experience),function(x){
    sort(liger_output_for_merged_dataset$W[x,],decreasing=T) -> test1
    test1 <- test1[test1>0]
    head(test1,100) %>% names 
  }) -> ca_only_feature_list_from_NMF_for_cistarget
names(ca_only_feature_list_from_NMF_for_cistarget) <- paste0('NMF_',c(1:k_experience))
lapply(names(ca_only_feature_list_from_NMF_for_cistarget),function(x){
  cat (x,'\n')
  tryCatch({
    t1 <- cisTarget(ca_only_feature_list_from_NMF_for_cistarget[[x]], motifRankings,
                    motifAnnot=motifAnnotations_hgnc,verbose = T)
    return(t1)
  },error=function(e){cat ('Cannot work',x,'\n')})
}) -> motifEnrichmentTable_wGenes_list
names(motifEnrichmentTable_wGenes_list) <- names(ca_only_feature_list_from_NMF_for_cistarget)
return(motifEnrichmentTable_wGenes_list)
}) -> cisTarget_result_list 
 
names(cisTarget_result_list) <- typelist 
names(GSEA_result) <- typelist  

save(cisTarget_result_list,file='Figure4.CisTarget.Result.List.Rdata')

## Draw TF-gene network underlying NMF-features metagene 



lapply(typelist, function(celltypeid){
  load(file=paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/Figure4.',celltypeid,'.NMF.original.output.Rdata'))
  load(file=paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/Figure4.',celltypeid,'.NMF.combined.Rdata'))
  k_experience <- nrow(liger_output_for_merged_dataset$W)
  lapply(c(1:k_experience),function(x){
    sort(liger_output_for_merged_dataset$W[x,],decreasing=T) -> test1
    test1 <- test1[test1>0]
    head(test1,100) %>% names 
  }) -> ca_only_feature_list_from_NMF_for_cistarget
  names(ca_only_feature_list_from_NMF_for_cistarget) <- paste0('NMF_',c(1:k_experience))
  lapply(names(ca_only_feature_list_from_NMF_for_cistarget),function(x){
    cat (x,'\n')
    tryCatch({
      t1 <- cisTarget(ca_only_feature_list_from_NMF_for_cistarget[[x]], motifRankings,
                      motifAnnot=motifAnnotations_hgnc,verbose = T)
      return(t1)
    },error=function(e){cat ('Cannot work',x,'\n')})
  }) -> motifEnrichmentTable_wGenes_list
  names(motifEnrichmentTable_wGenes_list) <- names(ca_only_feature_list_from_NMF_for_cistarget)
  return(motifEnrichmentTable_wGenes_list)
}) -> cisTarget_result_list 


## VisGene

lapply(names(cisTarget_result_list),function(nameid){
  cat (nameid, '\n')
  motifEnrichmentTable_wGenes_list <- cisTarget_result_list[[nameid]]
  celltypeid <- nameid
  load(file=paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/Figure4.',celltypeid,'.NMF.original.output.Rdata'))
  load(file=paste0('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/New/Figure4/NMF/Figure4.',celltypeid,'.NMF.combined.Rdata'))
  k_experience <- nrow(liger_output_for_merged_dataset$W)
  lapply(c(1:k_experience),function(x){
    sort(liger_output_for_merged_dataset$W[x,],decreasing=T) -> test1
    test1 <- test1[test1>0]
    head(test1,100) %>% names 
  }) -> ca_only_feature_list_from_NMF_for_cistarget
  names(ca_only_feature_list_from_NMF_for_cistarget) <- paste0('NMF_',c(1:k_experience))
  
  names(motifEnrichmentTable_wGenes_list) <- paste0('NMF_',c(1:20))
  # motifEnrichmentTable_wGenes_list <- motifEnrichmentTable_wGenes_list[motifEnrichmentTable_wGenes_list!=NULL]
  lapply(names(motifEnrichmentTable_wGenes_list),function(x){
    motifEnrichmentTable_wGenes_list[[x]] -> y
    if(!is.null(y)){
      y$nmfid <- x   
      return(y)
    }
  }) %>% bind_rows -> motifEnrichmentTable_wGenes_df
  motifEnrichmentTable_wGenes_df <- motifEnrichmentTable_wGenes_df[motifEnrichmentTable_wGenes_df$NES>3,]
  motifEnrichmentTable_wGenes_df <- left_join(motifEnrichmentTable_wGenes_df,motifEnrichmentTable_wGenes_df %>% dplyr::group_by(TF_highConf) %>% dplyr::summarise(TF_observed_n=n()))
  motifEnrichmentTable_wGenes_df_specific <- motifEnrichmentTable_wGenes_df[motifEnrichmentTable_wGenes_df$TF_observed_n<=2,]
  for (nmfid in unique(motifEnrichmentTable_wGenes_df_specific$nmfid)){
    tryCatch({
    cat (nameid, nmfid, '\n')
      maxNMFnum <- min(10,length(motifEnrichmentTable_wGenes_df_specific$motif))
      signifMotifNames <- motifEnrichmentTable_wGenes_df_specific$motif[motifEnrichmentTable_wGenes_df_specific$nmfid %in% nmfid][1:maxNMFnum]
      signifMotifMatchTF <- gsub(' .+','',motifEnrichmentTable_wGenes_df_specific$TF_highConf[motifEnrichmentTable_wGenes_df_specific$nmfid %in% nmfid][1:maxNMFnum])
      names(signifMotifMatchTF) <- signifMotifNames
      incidenceMatrix <- getSignificantGenes(ca_only_feature_list_from_NMF_for_cistarget[[nmfid]], 
                                             motifRankings,
                                             signifRankingNames=signifMotifNames,
                                             plotCurve=TRUE, maxRank=5000, 
                                             genesFormat="incidMatrix",
                                             method="aprox")$incidMatrix
      library(igraph)
      library(reshape2)
      edges <- melt(incidenceMatrix)
      edges <- edges[which(edges[,3]==1),1:2]
      colnames(edges) <- c("from","to")
      library(visNetwork)
      motifs <- unique(as.character(edges[,1]))
      genes <- unique(as.character(edges[,2]))
      nodes <- data.frame(id=c(motifs, genes),   
                          label=c(motifs, genes),    
                          title=c(motifs, genes), # tooltip 
                          shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                          color=c(rep("red", length(motifs)), rep("cornflowerblue", length(genes))))
      data.frame(label=as.character(signifMotifNames),new_name=as.character(signifMotifMatchTF)) -> temp
      temp$new_name <- as.character(temp$new_name)
      left_join(nodes,temp,by='label') -> nodes
      nodes$label <- as.character(nodes$label)
      nodes$new_name <- as.character(nodes$new_name)
      nodes$label[!is.na(nodes$new_name)] <- nodes$new_name[!is.na(nodes$new_name)]
      visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                              nodesIdSelection = TRUE,width = "300%",height = "300%") -> vn_obj
      visSave(vn_obj, file = paste0("Figure1M.VisNetwork.NMF.Features.after.Rcistopic.TF.gene.network.",nameid,".",nmfid,".html"), background = "white",selfcontained = T)
    },error=function(e){cat ('No Enrichment for',nameid, nmfid, '\n')})
  }
})

# use TF/target to fetch pathway


lapply(names(cisTarget_result_list),function(nameid){
  cat (nameid, '\n')
  motifEnrichmentTable_wGenes_list <- cisTarget_result_list[[nameid]]
  celltypeid <- nameid
  names(motifEnrichmentTable_wGenes_list) <- paste0('NMF_',c(1:20))
  lapply(names(motifEnrichmentTable_wGenes_list),function(x){
    motifEnrichmentTable_wGenes_list[[x]] -> y
    if(!is.null(y)){
      y$nmfid <- x   
      return(y)
    }
  }) %>% bind_rows %>% as.data.frame() -> motifEnrichmentTable_wGenes_df
  lapply(unique(motifEnrichmentTable_wGenes_df$nmfid),function(nmfid){
    motifEnrichmentTable_wGenes_df_sub <- motifEnrichmentTable_wGenes_df[motifEnrichmentTable_wGenes_df$nmfid %in% nmfid,]
    motifEnrichmentTable_wGenes_df_sub <- motifEnrichmentTable_wGenes_df_sub[motifEnrichmentTable_wGenes_df_sub$NES>3,]
    motifEnrichmentTable_wGenes_df_sub$enrichedGenes %>% stringr::str_split(pattern=";") %>% unlist %>% unique() -> targets
    motifEnrichmentTable_wGenes_df_sub$TF_highConf %>% stringr::str_split(pattern=";") %>% unlist %>% unique() -> TFs
    TFs <- gsub(' .+','',gsub('^ ','',TFs))
    library(spagi)
    data(pathway.path)
    lapply(c(targets,TFs),function(tfname){
      data.frame(pathname=names(pathway.path),Score=
                   lapply(names(pathway.path),function(x){
                     sum(unlist(pathway.path[[x]]) %in% tfname)
                   }) %>% unlist ) -> df
      df$tfname = tfname
      df
    }) %>% bind_rows -> test
    test <- test[test$Score>0,]
    test$nmfid <- nmfid
    test$celltypeid <- celltypeid
    table(test$pathname) %>% sort %>% tail %>% names -> enriched_receptor
    test[test$pathname %in% enriched_receptor,]
  }) %>% bind_rows -> enriched_pathway_for_tf_in_this_celltype
  enriched_pathway_for_tf_in_this_celltype
}) %>% bind_rows -> enriched_pathway_for_all

  
save.image('Figure4.NMF.working.image.Rdata')

