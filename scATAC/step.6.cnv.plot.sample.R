library(ComplexHeatmap)
read.delim('/gpfs/output/Bladder_Cancer_Landscape_Project/CT/CNV_amplicon_analysis/all_putative_amplicon.sorted.bed',col.names = c('chr','start','end','gene','log2FC','depth','probes','weight')) -> putative_amplicon
putative_amplicon$chr <- paste0('chr',putative_amplicon$chr)
putative_amplicon_gr <- reduce(easyLift::easyLiftOver(flank(reduce(makeGRangesFromDataFrame(putative_amplicon,keep.extra.columns = F)),width = 100000),map = 'hg19_hg38') )
putative_amplicon_gr <- reduce(flank(putative_amplicon_gr,width=100000))
read.delim('/gpfs/output/Bladder_Cancer_Landscape_Project/CT/CNV_amplicon_analysis/all_putative_focal_del.sorted.bed',col.names = c('chr','start','end','gene','log2FC','depth','probes','weight')) -> putative_focal_del
putative_focal_del$chr <- paste0('chr',putative_focal_del$chr)
putative_focal_del_gr <- reduce(easyLift::easyLiftOver(flank(reduce(makeGRangesFromDataFrame(putative_focal_del,keep.extra.columns = F)),width = 100000),map = 'hg19_hg38') )
putative_focal_del_gr <- reduce(flank(putative_focal_del_gr,width=100000))

cnv_ranges <- CNVMat@rowRanges
# get_amplicon_rows <- findOverlaps(cnv_ranges,putative_amplicon_gr)@from 
# matrix_of_possible_amplicon <- assays(CNVMat)[['CNVMatrix']][get_amplicon_rows,]
# 
# ampcnv_value_per_cluster_list <- list() 
# for (seurat_id in unique(CNVMat@colData$seurat_clusters)) {
#   matrix_of_possible_amplicon[,CNVMat@colData$seurat_clusters == seurat_id] -> tempvector
#   rowMeans(tempvector) -> mean_cnv_value
#   mean_cnv_value -> ampcnv_value_per_cluster_list[[seurat_id]]
# }

# t(as.matrix(bind_rows(ampcnv_value_per_cluster_list))) -> ampcnv_per_cluster_data
# 
# pdf('testcnvamp.pdf',height=10,width=10)
# pheatmap::pheatmap(ampcnv_per_cluster_data,cluster_cols = F,scale='none',cellheight = 10)
# dev.off() 

colnames(CNVMat@colData)[ncol(CNVMat@colData)] <- 'seurat_clusters'
cnv_value_per_cluster_list <- list() 
cnvtestdata <- list() 
cnv_sds_value_per_cluster_list <- list()
for (seurat_id in unique(CNVMat@colData$seurat_clusters)) {
  assays(CNVMat)[['CNVMatrix']][,CNVMat@colData$seurat_clusters == seurat_id] -> tempvector
  rowMeans(tempvector) -> mean_cnv_value
  colSds(as.matrix(tempvector)) -> sds_cnv_value
  mean_cnv_value -> cnv_value_per_cluster_list[[seurat_id]]
  sds_cnv_value -> cnv_sds_value_per_cluster_list[[seurat_id]]
  if(ncol(tempvector)<200){
    tempvector[,sample(ncol(tempvector),200,replace = T)] -> cnvtestdata[[seurat_id]]
  }else{
    tempvector[,sample(ncol(tempvector),200)] -> cnvtestdata[[seurat_id]]
  }
}

t(as.matrix(bind_rows(cnv_value_per_cluster_list))) -> cnv_per_cluster_data

# pdf('testcnvall.pdf',height=10,width=10)
# pheatmap::pheatmap(cnv_per_cluster_data,cluster_cols = F,scale='none',cellheight = 10,clustering_method = 'ward.D2')
# dev.off() 
# 
# library(mixtools)
# normalmixEM(as.vector(unlist(unlist(cnv_value_per_cluster_list))),k=3) -> cnv_mdl
# lower_bound <- cnv_mdl$mu[2] - 3*cnv_mdl$sigma[2]
# upper_bound <- cnv_mdl$mu[2] + 3*cnv_mdl$sigma[2]

lower_bound <- -0.2
upper_bound <- 0.2 

lapply(cnvtestdata,function(x) as.matrix(x)) -> cnvtestdata2
as.matrix(bind_rows(cnvtestdata2)) -> tdf
as.data.frame(rep(names(cnvtestdata),each=200)) -> group_name_list
rownames(group_name_list) <- c(1:nrow(group_name_list))
colnames(group_name_list) <- 'cluster'
group_name_list$cluster <- factor(group_name_list$cluster,levels=paste0('C',c(1:length(unique(group_name_list$cluster)))),ordered=T)


as.data.frame(CNVMat@rowRanges) -> cnv_ranges_df
rownames(cnv_ranges_df) <- c(1:nrow(cnv_ranges_df))
cnv_ranges_df$seqnames <- factor(cnv_ranges_df$seqnames,levels=paste0('chr',c(1:22,"X")),ordered=T)
cnv_ranges_df$oldid <- 1:nrow(cnv_ranges_df)
cnv_ranges_df <- arrange(cnv_ranges_df,seqnames,start)
cnv_ranges <- makeGRangesFromDataFrame(cnv_ranges_df,keep.extra.columns = T)

tdf[tdf<lower_bound] <- -1
tdf[tdf>upper_bound] <- 1
tdf[tdf<upper_bound & tdf>lower_bound] <- 0 
tdf <- t(tdf)
rownames(tdf) <- c(1:nrow(tdf))
colnames(tdf) <- c(1:ncol(tdf))
tdf <- tdf[order(group_name_list$cluster),order(cnv_ranges_df$oldid)]
group_name_list$cluster <- group_name_list$cluster[order(group_name_list$cluster)]

# pdf('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/Save/Bladder_epithelial_cell_20201212_test_CNV_adjustment/CNV_subsampled_200cells_per_cluster.pdf',height=10,width=40)
# # pheatmap::pheatmap(tdf[,order(cnv_ranges_df$seqnames)],cluster_cols = F,cluster_rows=F,scale='none',annotation_row = group_name_list,annotation_col=cnv_ranges_df %>% dplyr::select(seqnames),show_rownames = F,show_colnames = F)
# Heatmap(tdf[,order(cnv_ranges_df$seqnames)], 
#     top_annotation = HeatmapAnnotation( foo = anno_block(gp = gpar(fill = rep(c('black','gray'),100)[1:23]),
#         # labels = as.vector(unique(cnv_ranges_df$seqnames)[order(unique(cnv_ranges_df$seqnames))]),
#         labels_gp = gpar(col = "white", fontsize = 10))),column_split = cnv_ranges_df[order(cnv_ranges_df$seqnames),] %>% dplyr::select(seqnames),
#     left_annotation = rowAnnotation( foo = anno_block(gp = gpar(fill = 2:(length(as.vector(unique(group_name_list$cluster)+1)))),
#         # labels = as.vector(unique(group_name_list$cluster)),
#         labels_gp = gpar(col = "white", fontsize = 10))),row_split=group_name_list,
#     cluster_rows=F,cluster_columns=F,
# show_column_dend=F,show_column_names=F,
# show_row_dend=F,show_row_names=F,
# gap=unit(0, "npc"),row_gap=unit(0, "npc"))
# dev.off() 


pdf(paste0(workingpath,'/Plots/CNV_subsampled_200cells_per_cluster.pdf'),height=10,width=40)
# pheatmap::pheatmap(tdf[,order(cnv_ranges_df$seqnames)],cluster_cols = F,cluster_rows=F,scale='none',annotation_row = group_name_list,annotation_col=cnv_ranges_df %>% dplyr::select(seqnames),show_rownames = F,show_colnames = F)
Heatmap(tdf[,order(cnv_ranges_df$seqnames)],
    top_annotation = HeatmapAnnotation( foo = anno_block(gp = gpar(fill = rep(c('black','gray'),100)[1:23]),
        # labels = as.vector(unique(cnv_ranges_df$seqnames)[order(unique(cnv_ranges_df$seqnames))]),
        labels_gp = gpar(col = "white", fontsize = 10))),column_split = cnv_ranges_df[order(cnv_ranges_df$seqnames),] %>% dplyr::select(seqnames),
    left_annotation = rowAnnotation( foo = anno_block(gp = gpar(fill = 2:(length(as.vector(unique(group_name_list$cluster)+1)))),
        # labels = as.vector(unique(group_name_list$cluster)),
        labels_gp = gpar(col = "white", fontsize = 10))),row_split=group_name_list,
    cluster_rows=F,cluster_columns=F,
show_column_dend=F,show_column_names=F,
show_row_dend=F,show_row_names=F,
gap=unit(0, "npc"),row_gap=unit(0, "npc"))
dev.off()
