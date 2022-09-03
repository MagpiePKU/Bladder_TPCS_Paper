suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk)) # use the 'loom' branch
suppressPackageStartupMessages(library(velocyto.R))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(reticulate))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(liger))
suppressPackageStartupMessages(library(SeuratWrappers))
suppressPackageStartupMessages(library(metR))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GenomeInfoDb))
suppressPackageStartupMessages(library(ChIPseeker))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(easyLift))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ArchR))
here()

scv <- reticulate::import("scvelo") # Note: add the reticulate pointer to avoid error
ad <- reticulate::import("anndata", convert = FALSE) # Note: add the reticulate pointer to avoid error
# load('/gpfs/genomedb/singleR/ref_Human_all.RData')

options(bitmapType='cairo') # solve the PNG/X11 issue





plot_regulon_function <- function(
  target_gene = 'SOX2',
  target_chromosome='chr3',
  peaks_to_plot = NA,
  mycl,
  myclt,
   GTF,
  use_Groups = 'ALL',
  target_interaction_Groups = 'ALL',
  project_input=proj_arrow,
  subsetproj=subsetproj,
  peakMatrix=peakMatrix,
  cA_good=cA_good,
  cA=cA,
  peakSinCA_gr=peakSinCA_gr,
  plotSummary = c("bulkTrack", "scTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 10, 2, 1.5, 1),
  group_by_keyword = 'seurat_clusters',
  extend_width=0,
  plot_width = 16,
  plot_height = 10,
  prefix=prefix,
  corCutOff = corCutOff,
  scCellsMax=300
) {
        # define groups
        if (use_Groups == 'ALL'){
          use_Groups = unique(project_input$seurat_clusters)
          use_Groups <- as.character(use_Groups)
          use_Groups <- use_Groups[order(use_Groups)]
        }
        # Build annotation
        GTF[GTF$symbol %in% target_gene,] -> sub_gtf
        sub_gtf <- sub_gtf[!grepl('-AS',sub_gtf$X9) & !grepl('-OT',sub_gtf$X9),]
        sub_gtf[sub_gtf$X3 %in% 'gene',] -> sub_gtf_for_genes
        sub_gtf[sub_gtf$X3 %in% 'exon',] -> sub_gtf_for_exons
        sub_gtf[sub_gtf$X3 %in% 'transcript',] -> sub_gtf_for_transcript
        sub_gtf_for_transcript$length <- sub_gtf_for_transcript$X5-sub_gtf_for_transcript$X4
        sub_gtf_for_transcript <- sub_gtf_for_transcript[sub_gtf_for_transcript$length == max(sub_gtf_for_transcript$length),][1,]
        sub_gtf_for_TSS <- sub_gtf_for_transcript
        if(unique(sub_gtf_for_TSS$X7) %in% '-'){
          sub_gtf_for_TSS$X4 <- sub_gtf_for_TSS$X5 - 1
        }else{
          sub_gtf_for_TSS$X5 <- sub_gtf_for_TSS$X4 + 1
        }
        sub_gtf_for_exons <- sub_gtf_for_exons[sub_gtf_for_exons$X7 %in% unique(sub_gtf_for_TSS$X7),]
        genes_gr <- makeGRangesFromDataFrame(sub_gtf_for_genes %>% dplyr::select(X1,X4,X5,X7,symbol),seqnames.field = 'X1',start.field = 'X4',end.field = 'X5',strand.field = 'X7',keep.extra.columns = T)
        exon_gr <- makeGRangesFromDataFrame(sub_gtf_for_exons%>% dplyr::select(X1,X4,X5,X7,symbol),seqnames.field = 'X1',start.field = 'X4',end.field = 'X5',strand.field = 'X7',keep.extra.columns = T)
        tss_gr <- makeGRangesFromDataFrame(sub_gtf_for_TSS%>% dplyr::select(X1,X4,X5,X7,symbol),seqnames.field = 'X1',start.field = 'X4',end.field = 'X5',strand.field = 'X7',keep.extra.columns = T)
        customizedAnnotation <- createGeneAnnotation(genes=genes_gr,exons=exon_gr,TSS=tss_gr)
        genes_gr_extended <- GenomicRanges::flank(genes_gr,width=extend_width)

  
      tryCatch({
        peakMatrix[peakMatrix$nearestGene %in% target_gene,] -> target_peaks
      if(sum(target_peaks$seqnames %in% target_chromosome) == 0){
        cat ('No peaks with target gene\n')
        plot_region_gr <- genes_gr_extended
        loopsobj_to_plot <- getCoAccessibility(subsetproj,resolution = 20,corCutOff = 0.2,returnLoops=T)
      }else{
        target_cluster <- mycl[names(mycl) %in% rownames(target_peaks)]
        target_cluster_uniq <- unique(as.vector(target_cluster))
        target_cluster_with_peak_number <- myclt[myclt$mycl %in% target_cluster_uniq,]
        if(nrow(target_cluster_with_peak_number)>=5){
          topcluster <- target_cluster_with_peak_number[1:2,]
        }else{
          topcluster <- target_cluster_with_peak_number[1,]
        }
        peaks_in_topcluster <- names(mycl[mycl %in% topcluster$mycl])
        initial_peaks <- peaks_in_topcluster
        increase_in_peaks = 1
        step_outside = 0
        while(increase_in_peaks > 0 ){
          cA_good_test[cA_good_test$queryHits %in% initial_peaks | cA_good_test$subjectHits %in% initial_peaks,] -> cluster_extended
          extended_peaks <- unique(c(cluster_extended$queryHits,cluster_extended$subjectHits))
          increase_in_peaks <- length(extended_peaks) - length(initial_peaks)
          initial_peaks <- extended_peaks
          step_outside = step_outside + 1
        }
        cat ('The regulon extension stops at step',step_outside,'with peaks number',length(extended_peaks),'\n')
        
        extended_peaks_matrix <- peakMatrix[extended_peaks,]
        plot_region_dataframe <- as.data.frame(list('chr'=unique(extended_peaks_matrix$seqnames),'start'=min(extended_peaks_matrix$start),'end'=max(extended_peaks_matrix$end)))
        plot_region_dataframe$width <- plot_region_dataframe$end - plot_region_dataframe$start
        if(extend_width > 0 ){
          extension_bp = extend_width
        }else{
          if(plot_region_dataframe$width < 200000){
            extension_bp <- 100000
          }else{
            extension_bp <- round(plot_region_dataframe$width/5)
          }
        }
        plot_region_dataframe$end <- plot_region_dataframe$end + extension_bp
        plot_region_dataframe$start <- plot_region_dataframe$start - extension_bp
        plot_region_dataframe$width <- plot_region_dataframe$end - plot_region_dataframe$start
        plot_region_gr <- makeGRangesFromDataFrame(plot_region_dataframe)
        # Build Loop
        extended_peaks_matrix_gr <- makeGRangesFromDataFrame(extended_peaks_matrix)
        
        # 1. Calculate loop number
        
        loopsobj <- getCoAccessibility(subsetproj,resolution = 20,corCutOff = 0.2,returnLoops=T)
        loopsobj_df <- GenomicRanges::as.data.frame(loopsobj)
        loopsobj_df_start <- loopsobj_df
        loopsobj_df_start$end <- loopsobj_df_start$start + 1
        loopsobj_df_start_gr <- makeGRangesFromDataFrame(loopsobj_df_start)
        findOverlaps(loopsobj_df_start_gr,extended_peaks_matrix_gr) -> start_overlap
        loops_number <- length(unique(start_overlap@from))
        
        # 2. Calculate graph width
        width_for_graph <- plot_region_dataframe$width
        
        # 3. Calculate loop density
        loop_density <- width_for_graph/loops_number
        if(loop_density < 1000 & loops_number > 100){
          resolution_new <- max(50,10^((round(log(loop_density,10)))))
        }else{
          resolution_new <- 20
        }
        
        # 4. Scale co-accessibility resolution based on loop density
        loopsobj <- getCoAccessibility(subsetproj,resolution = resolution_new,corCutOff = 0.2,returnLoops=T)
        loopsobj_df <- GenomicRanges::as.data.frame(loopsobj)
        loops_number <- nrow(loopsobj_df)
        loopsobj_df_start <- loopsobj_df
        loopsobj_df_end <- loopsobj_df
        loopsobj_df_start$end <- loopsobj_df_start$start + 1
        loopsobj_df_end$start <- loopsobj_df_end$end - 1
        loopsobj_df_start_gr <- makeGRangesFromDataFrame(loopsobj_df_start)
        loopsobj_df_end_gr <- makeGRangesFromDataFrame(loopsobj_df_end)
        
        findOverlaps(loopsobj_df_start_gr,extended_peaks_matrix_gr) -> start_overlap
        findOverlaps(loopsobj_df_end_gr,extended_peaks_matrix_gr) -> end_overlap
        loops_df_of_extended_peaks <- loopsobj_df[unique(c(start_overlap@from,end_overlap@from)),]
        
        loopsobj_to_plot <- list('CoAccessibility' = makeGRangesFromDataFrame(loops_df_of_extended_peaks %>% dplyr::select(seqnames,start,end,value),keep.extra.columns = T))
        loops_number <- length(unique(start_overlap@from))
        width_for_graph <- plot_region_dataframe$width
        loop_density_new <- width_for_graph/loops_number
        
        cat ('Before trimming, loop density is 1 loop per',loop_density,'bp, after trimming, loop density is 1 loop per',loop_density_new,'bp with new resolution',resolution_new,'bp for coaccessibility \n')
        
        loopsobj_to_plot <- list('CoAccessibility' = makeGRangesFromDataFrame(loops_df_of_extended_peaks %>% dplyr::select(seqnames,start,end,value),keep.extra.columns = T)) 
      }
      },error=function(e){
        plot_region_gr <- genes_gr_extended
        loopsobj_to_plot <- getCoAccessibility(subsetproj,resolution = 20,corCutOff = 0.2,returnLoops=T)
      })
      
      # Build annotation
      GTF[GTF$symbol %in% target_gene,] -> sub_gtf
      sub_gtf <- sub_gtf[!grepl('-AS',sub_gtf$X9) & !grepl('-OT',sub_gtf$X9),]
      sub_gtf[sub_gtf$X3 %in% 'gene',] -> sub_gtf_for_genes
      sub_gtf[sub_gtf$X3 %in% 'exon',] -> sub_gtf_for_exons
      sub_gtf[sub_gtf$X3 %in% 'transcript',] -> sub_gtf_for_transcript
      sub_gtf_for_transcript$length <- sub_gtf_for_transcript$X5-sub_gtf_for_transcript$X4
      sub_gtf_for_transcript <- sub_gtf_for_transcript[sub_gtf_for_transcript$length == max(sub_gtf_for_transcript$length),][1,]
      sub_gtf_for_TSS <- sub_gtf_for_transcript
      if(unique(sub_gtf_for_TSS$X7) %in% '-'){
        sub_gtf_for_TSS$X4 <- sub_gtf_for_TSS$X5 - 1 
      }else{
        sub_gtf_for_TSS$X5 <- sub_gtf_for_TSS$X4 + 1 
      }
      sub_gtf_for_exons <- sub_gtf_for_exons[sub_gtf_for_exons$X7 %in% unique(sub_gtf_for_TSS$X7),]
      genes_gr <- makeGRangesFromDataFrame(sub_gtf_for_genes %>% dplyr::select(X1,X4,X5,X7,symbol),seqnames.field = 'X1',start.field = 'X4',end.field = 'X5',strand.field = 'X7',keep.extra.columns = T)
      exon_gr <- makeGRangesFromDataFrame(sub_gtf_for_exons%>% dplyr::select(X1,X4,X5,X7,symbol),seqnames.field = 'X1',start.field = 'X4',end.field = 'X5',strand.field = 'X7',keep.extra.columns = T)
      tss_gr <- makeGRangesFromDataFrame(sub_gtf_for_TSS%>% dplyr::select(X1,X4,X5,X7,symbol),seqnames.field = 'X1',start.field = 'X4',end.field = 'X5',strand.field = 'X7',keep.extra.columns = T)
      customizedAnnotation <- createGeneAnnotation(genes=genes_gr,exons=exon_gr,TSS=tss_gr)
      
      # Plot
        p <- plotBrowserTrack(
          ArchRProj = subsetproj,
          region = plot_region_gr,
          plotSummary = plotSummary,
          sizes = sizes,
          useGroups = use_Groups,
          groupBy = group_by_keyword,
          geneAnnotation = customizedAnnotation,
          geneSymbol = target_gene,
          scCellsMax = scCellsMax,
          threads = 10,
          features = peaks_to_plot,
          loops = loopsobj_to_plot
        )

      # print pdf
      plotPDF(p,
              name = paste0("Plot-Tracks-Regulon-with-",prefix,'-',target_gene,".pdf"),
              ArchRProj = subsetproj,
              addDOC = FALSE, width = plot_width, height = plot_height)
}

## plot_regulon
## zy 20201127
## Function to plot a regulon network by co-accesibility, for a group of genes
## the function does:
## 1. filter the gene list (genes without peaks-associated are removed)
## 2. calculate co-clustering of enhancer per chromosome, by hclust
## 3. plot the regulon using plot_regulon_function
## inputs
## target_gene_list: a vector of gene names (character)
## use_Groups: default to 'ALL', or a character vector denoting specific groups of cells to plot
## target_interaction_Groups: default to 'ALL', or a character vector denoting the co-accessibility network to be calculated from. For example, if you only wish to analyze regulatory network in T cells but not cancer cells, use 'T cells' or 'C1' (depending on the cell cluster name you gave)
## extend_width: default to 0, a numeric variable defining how many bp to extend from the gene. If 0, automatically extends the width as defined by the regulon 
## plot_width: size of plot, see plotBrowswerTrack
## plot_height: size of plot, see plotBrowswerTrack
## prefix: a string prefix for saving the pdf
## corCutOff: a numeric variable for correlation cutoff, usually set to 0.5 (50% correlation, above chance)
## scCellsMax: a numeric variable, how many single cells to plot (if plotted)
## plotSummary: a character vector defining which tracks to plot, see plotBrowswerTrack
## sizes: a numeric vector defining height of the tracks to plot, see plotBrowswerTrack
## group_by_keyword: the name of cell group, such as 'seurat_clusters' or 'RNA_integration_clusters', etc
## peaks_to_plot: annotation sets for additional annotation, such as DMR or ATAC peaks. A GRangeList object
## GTF: a GTF for extracting gene information
## 

plot_regulon <- function(
  target_gene_list='SOX2',
  use_Groups = 'ALL',
  target_interaction_Groups = 'ALL',  
  prefix='cancer',
  project_input=proj_arrow,
  peaks_to_plot=DMR_region_to_plot_gr,
  plotSummary = c("bulkTrack", "scTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 10, 2, 1.5, 1),
  group_by_keyword = 'seurat_clusters',
  extend_width=0,
  plot_width = 16,
  plot_height = 10,
  corCutOff = 0.2,
  scCellsMax=300,
  GTF=GTF
){
  # define groups
  if (use_Groups == 'ALL'){
    use_Groups = unique(project_input$seurat_clusters)
  }
  # define list
  final_gene_list <- list() 
  for (geneid in target_gene_list){
    if(sum(peakMatrix$nearestGene %in% geneid) > 0){
      final_gene_list[[geneid]] <- geneid
    }else{
      cat ('no peaks associated with',geneid,'\n')
    }
  }
  GTF[GTF$symbol %in% target_gene_list & GTF$X3 %in% 'gene',] -> subGTF_for_genelist
  target_chr_list <- as.character(unique(subGTF_for_genelist$X1))
  if(length(target_interaction_Groups) > 1){
    target_interaction_Groups <- paste0(target_interaction_Groups,collapse = "|")
  }
  if(target_interaction_Groups == 'ALL'){
    subsetproj = project_input
  }else{
    idxSample <- (grepl(target_interaction_Groups,project_input$seurat_clusters))
    subset_cell_names <- project_input$cellNames[idxSample]
    subsetproj <- project_input[subset_cell_names,]
    subsetproj <- addCoAccessibility(
      ArchRProj = subsetproj,
      reducedDims = "IterativeLSI"
    )
  }
  cA <- getCoAccessibility(
    ArchRProj = subsetproj,
    corCutOff = corCutOff,
    resolution = 50,
    returnLoops = FALSE
  )
  peakMatrix <- as.data.frame(subsetproj@peakSet,row.names=NULL)
  target_peak_matrix <- peakMatrix[peakMatrix$nearestGene %in% final_gene_list,]
  cA_good <- cA[cA$FDR<(1/nrow(as.data.frame(cA))),] # Select FDR unlikely (more than one in the dataset)
  cA_good_test <- cA_good
  unique(c(cA_good$queryHits,cA_good$subjectHits)) -> hits_row
  peakSinCA <- peakMatrix[hits_row,]
  peakSinCA_gr <- makeGRangesFromDataFrame(peakSinCA,keep.extra.columns = T)
  target_peak_matrix <- peakMatrix[peakMatrix$nearestGene %in% names(final_gene_list),]
  target_chr_list <- as.character(unique(target_peak_matrix$seqnames))
  for(target_chr in target_chr_list){
    cat ('Working on chr',target_chr,'\n')
    target_genes_on_target_chr <- unique(target_peak_matrix$nearestGene[target_peak_matrix$seqnames %in% target_chr])
    chr_selected_peaks <- rownames(peakMatrix)[peakMatrix$seqnames %in% target_chr]
    cA_good_test_chr <- cA_good_test[cA_good_test$queryHits %in% chr_selected_peaks | cA_good_test$subjectHits %in% chr_selected_peaks,]
    as.data.table(cA_good_test_chr) -> cA_good_test_chr
    dcast(as.data.table(cA_good_test_chr),queryHits ~ subjectHits,value.var='correlation',fun.aggregate=sum) -> test
    test <- 1-test[,2:ncol(test)]
    hc <- hclust(as.dist(test), method="complete") 
    mycl <- cutree(hc, h=max(hc$height)/1.1);
    myclt <- arrange(as.data.frame(table(mycl)),-Freq)
    for (target_gene in target_genes_on_target_chr){
      cat ('Working on gene',target_gene,'\n')
      tryCatch({
        plot_regulon_function(target_gene = target_gene,target_chromosome = target_chr,peaks_to_plot = peaks_to_plot,mycl=mycl,myclt=myclt,GTF=GTF,use_Groups = use_Groups,target_interaction_Groups = target_interaction_Groups, prefix=prefix,project_input = project_input,peakMatrix=peakMatrix,cA_good=cA_good,cA=cA,peakSinCA_gr=peakSinCA_gr,subsetproj=subsetproj,  plotSummary = plotSummary,sizes = sizes,group_by_keyword = group_by_keyword,  extend_width=extend_width,corCutOff = corCutOff,  scCellsMax=scCellsMax)
      },error=function(e){cat ('No regulon for gene',target_gene,'\n')})
    }
  }
}



## Test code 1

# target_gene_list='ALK'
# use_Groups = 'ALL'
# target_interaction_Groups = 'ALL'
# prefix='cancer'
# project_input=proj_arrow
# peaks_to_plot=AnnotationSet
# plotSummary = c("bulkTrack",  "featureTrack", "loopTrack", "geneTrack")
# sizes = c(6, 2, 1.5, 1)
# group_by_keyword = 'seurat_clusters'
# extend_width=100000
# plot_width = 16
# plot_height = 10
# 
# if (use_Groups == 'ALL'){
#   use_Groups = unique(project_input$seurat_clusters)
# }
# # define list
# final_gene_list <- list() 
# for (geneid in target_gene_list){
#   if(sum(peakMatrix$nearestGene %in% geneid) > 0){
#     final_gene_list[[geneid]] <- geneid
#   }else{
#     cat ('no peaks associated with',geneid,'\n')
#   }
# }
# GTF[GTF$symbol %in% target_gene_list & GTF$X3 %in% 'gene',] -> subGTF_for_genelist
# target_chr_list <- as.character(unique(subGTF_for_genelist$X1))
# if(length(target_interaction_Groups) > 1){
#   target_interaction_Groups <- paste0(target_interaction_Groups,collapse = "|")
# }
# if(target_interaction_Groups == 'ALL'){
#   subsetproj = project_input
# }else{
#   idxSample <- (grepl(target_interaction_Groups,project_input$seurat_clusters))
#   subset_cell_names <- project_input$cellNames[idxSample]
#   subsetproj <- project_input[subset_cell_names,]
#   subsetproj <- addCoAccessibility(
#     ArchRProj = subsetproj,
#     reducedDims = "IterativeLSI"
#   )
# }
# cA <- getCoAccessibility(
#   ArchRProj = subsetproj,
#   corCutOff = 0.1,
#   resolution = 50,
#   returnLoops = FALSE
# )
# peakMatrix <- as.data.frame(subsetproj@peakSet,row.names=NULL)
# target_peak_matrix <- peakMatrix[peakMatrix$nearestGene %in% final_gene_list,]
# cA_good <- cA[cA$FDR<(1/nrow(as.data.frame(cA))),] # Select FDR unlikely (more than one in the dataset)
# cA_good_test <- cA_good
# unique(c(cA_good$queryHits,cA_good$subjectHits)) -> hits_row
# peakSinCA <- peakMatrix[hits_row,]
# peakSinCA_gr <- makeGRangesFromDataFrame(peakSinCA,keep.extra.columns = T)
# target_peak_matrix <- peakMatrix[peakMatrix$nearestGene %in% names(final_gene_list),]
# target_chr_list <- as.character(unique(target_peak_matrix$seqnames))
# target_chr = target_chr_list[1]
# cat ('Working on chr',target_chr,'\n')
# target_genes_on_target_chr <- unique(target_peak_matrix$nearestGene[target_peak_matrix$seqnames %in% target_chr])
# chr_selected_peaks <- rownames(peakMatrix)[peakMatrix$seqnames %in% target_chr]
# cA_good_test_chr <- cA_good_test[cA_good_test$queryHits %in% chr_selected_peaks | cA_good_test$subjectHits %in% chr_selected_peaks,]
# as.data.table(cA_good_test_chr) -> cA_good_test_chr
# dcast(as.data.table(cA_good_test_chr),queryHits ~ subjectHits,value.var='correlation',fun.aggregate=sum) -> test
# test <- 1-test[,2:ncol(test)]
# hc <- hclust(as.dist(test), method="complete") 
# mycl <- cutree(hc, h=max(hc$height)/1.1);
# myclt <- arrange(as.data.frame(table(mycl)),-Freq)
# 
# target_gene = target_gene_list[1]
# plot_regulon_function(target_gene = target_gene,target_chromosome = target_chr,peaks_to_plot = peaks_to_plot,mycl=mycl,myclt=myclt,GTF=GTF,use_Groups = use_Groups,target_interaction_Groups = target_interaction_Groups, prefix=prefix,project_input = project_input,peakMatrix=peakMatrix,cA_good=cA_good,cA=cA,peakSinCA_gr=peakSinCA_gr,subsetproj=subsetproj,  plotSummary = plotSummary,sizes = sizes,group_by_keyword = group_by_keyword,  extend_width=extend_width)

## Test code 2 

# plot_regulon(
#   target_gene_list='SOX2',
#   use_Groups = 'ALL',
#   target_interaction_Groups = 'ALL',  
#   prefix='cancer',
#   project_input=proj_arrow,
#   peaks_to_plot=AnnotationSet,
#   plotSummary = c("bulkTrack", 'scTrack', "featureTrack", "loopTrack", "geneTrack"),
#   sizes = c(6, 6, 2, 1.5, 1),
#   group_by_keyword = 'seurat_clusters',
#   extend_width=100000,
#   plot_width = 16,
#   plot_height = 14
# )


markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
  "CD14", "CEBPB", "MPO", #Monocytes
  "IRF8", 'PTPRB', # Immune
  "FOXP3",'AIRE', # Treg
  'PECAM1','PLVAP','COL6A2','COL1A1','CD10','PTPRB','EPCAM','CD45','PTPRC','MME',
  'CCR7','LEF1','SELL','TCF7','CD27','CD28','S1PR1','KLRG1','NCR3','KLRF1','CD27','CD28','ID2','PRF1','HOPX','GZMB','GZMM','GZMK','CX3CR1','FGGR3A','CCL5','GPR183','XCL1','XCL2','CXCR3','CXCR4','CD44','CD6','S1PR5','CX3CR1','KLF2','ADGRG1','LAIR1','LAIR2','CD69','MYADM','CAPG','RORA','NR4A1','NR4A2','NR4A3','ITGA3','CD160','KIR2DL4','TMIGD2','IKZF2','ICOS','ENTPD1','CXCL13','GZMB','MIR155HG','TNFRSF9','SLC4A10','ZBTB16','RORC','RORA','NCR3','CCR7','ANXA1','ANXA2','ICAM1','ICAM2','PTGER2','GNLY','CTSW','RGS1','CXCR6','MYADM','CXCR5','BCL6','ICA1','IL6ST','MAGEH1','BTLA','CD200','RUNX3','EOMES','IL23R','IL17A','IL17F','FURIN','CTSH','CCR6','CAPG','BHLHE40','IGFLR1','IL2RA','IL10RA','IKZF2','RTKN2','CDC25B','S1PR4','CCR4','CCR8','TNFRSF18','CTLA4','BATF','IL21R','TBX21','CD5','CD85','CD56','CD127','CD160','CD161', # ZZM
  "CD3D", "CD8A", "TBX21", "IL7R",'CD4','PRF1', 'GZMB','GZMA','KLRG1', #TCells
  'IFNA','IFNB','IFNG','TNFA','TGFB','IL1A','IL1B','IL2','IL3','IL4','IL5','IL6','CXCL8','IL10','IL12A','IL12B', # Inflammation
  'NCAM1','NKG7','KLRD1','KLRG1','KIR3DL1','LYAR','CD38','AAK1','BHLHE40','NR2E3','TICAM1','JPT1','ELF1','CD226','RASGRP1','FASLG','KLRC1','KLRD1','KLRK1','KRT74',#  NK
  'TNFRSF25','SELL','GDF10',
  'FCRL4','FCRL3','FCRL2','FCRL1','LY9',
  'CCL1','FLT3',
  'FZD4','ROBO4',
  'GJA4','VWA1','TIE1','TAL1','S1PR1',
  'SYT13','SLC37A4','SLC6A13',
  'TGFB2-AS1','TAGLN','TGFB2',
  'HOXD3','HOXB9',
  'MS4A8',
  'CDC42','KIF1B',
  'NLRP9','SGK2','MYBL2',
  'HES4',
  'FCGR2C','FCGR3B','CD160',
  'ATOH7','RASAL1',
  'TBX2',
  'GJB3','TTC22','FOXG1','TGM5','KRT19','KRT15',
  'EPCAM', # Epithelial cells
  'UPK1A','UPK1B','UPK3A', # Urothelial cells
  'TJP1','DSP','CDH1', # Epithelial - not yet mesenchymal
  'VIM','CDH2','FOXC2','SNAI1','SNAI2','TWIST1','GSC','FN1','ITBG6','MMP2','MMP3','MMP9','SOX10', # Mesenchymal
  'SOX2', 'IFITM10','KLF6','FGFR3','PTGES','FOXA1','AR', # MIBC vs NMIBC
  'CD274','CTLA4','HAVCR2','ICOS','ICOSLG','TNFRSF4','TNFRSF9','FASLG','LAG3','PDCD1','TOX','EOMES','TCF7', # Immunosuppress
  'HRAS','FGF4','CCND1' # Patient specific markers
)

track_markers_to_plot_cancer <- unique(c('MDM4','ONECUT','NRG1','ALK','MET','TBX5',
                                  'TESC','TNIK','TERC','PPARG','E2F3','ID4','RIPOR2','CDKAL1','PRL','RICTOR','YWHAZ','ERLIN2','HKR1','MDM2','CDK4','IL4R','AKT2','IQSEC1','CDK12','SOX2','FGFR3','FOXA1','ESRRA','ERBB2','AR','VAV2','CDKN2A','TBX3','SOX17','HRAS','CCND1','FGF4','FGF19','FGF3','SOX4','SOS1','ONECUT1','MYC','MYCL','PIM1','PPIE','VEGFA','NDRG1','FOS','CCNA1','EGFR','MYO5A','ALPK2','TXNL1','NFIB','PRKCB','PPARG','AKAP13','CD24'))

MARKERS <- unique(c(markerGenes,track_markers_to_plot_cancer))

library(igraph)

setwd('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/')
addArchRThreads(threads = 64)
addArchRGenome("hg38")
proj_arrow <- loadArchRProject(path = "Save/Bladder_20201127")

load('peakMatrix.Rdata')
load('peakSinCA.Rdata')
load('peakSinCA_temp.with.distToGene.Rdata')
load('tss_peakSinCA.with.enhancer.link.info.Rdata')
load('testdflist.gene.associated.enhancers.list.Rdata')
load('AnnotationSet.DMR.ATAC.bedfiles.all.list.Rdata')
load('cA_good_with_range.500kb.coaccessible.pairs.Rdata')
load('markers_arrow.Rdata')
GTF <- read_tsv('/gpfs/genomedb/b38/gencode.v34.basic.annotation.gtf.gz',col_names = F,comment = '#')
GTF$symbol <- gsub('.+gene_name \\\"','',GTF$X9)
GTF <- separate(GTF,col='symbol',into=c('symbol','useless'))
GTF$useless <- NULL 




