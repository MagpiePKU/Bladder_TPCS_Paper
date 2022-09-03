suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratDisk)) # use the 'loom' branch
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

load('AnnotationSet.DMR.ATAC.bedfiles.all.list.Rdata')

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

## Test regulon with adjacency matrix


plot_regulon_defined_by_matrix <- function(testdflist=testdflist,geneid='PDCD1',useGroups=c('C1','C2','C3','C4'),markerList=markerList,plotwidth=12,plotheight=6){
        # Check the extracted regulon for target gene
            # checkpoint
      # geneid = 'PDCD1'
      # useGroups=c('C1','C2','C3','C4')
      # head(testdflist[[geneid]])
      # tail(testdflist[[geneid]])
      
      # build target gene specific GTF annotation
      GTF[GTF$symbol %in% geneid,] -> sub_gtf
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
      
      # build draw region
      test_target_df <- testdflist[[geneid]]
      test_target_peak_gr <- makeGRangesFromDataFrame(test_target_df)
      
      # subset the loop start/end with target cell type (define 'peaks of celltype' at certain threshold, say as FDR<0.01 & LogFC>1 peaks in the specific cell type. Otherwise not.)
      markerList2 <- lapply(markerList,function(x){
          x <- x[findOverlaps(x,test_target_peak_gr)@from,]
          xdf <- as.data.frame(x)
          xdf <- left_join(xdf,test_target_df %>% dplyr::select(seqnames,start,end,peakId),by = c("seqnames", "start", "end"))
          xdfgr <- makeGRangesFromDataFrame(xdf,keep.extra.columns = T)
          return(xdfgr)
      })
      markerList2 <- markerList2[names(markerList2) %in% useGroups]
      markerdf <- bind_rows(lapply(markerList2,function(x){x=as.data.frame(x);return(x)}))
      
      # polish draw region by defining peaks within the marker range. Otherwise, remove them
      test_target_peak_gr <- test_target_peak_gr[findOverlaps(test_target_peak_gr,makeGRangesFromDataFrame(markerdf))@from,]
      # extend the plot region to include all peaks-to-plot 
      test_target_gr <- makeGRangesFromDataFrame(data.frame('chr'=unique(test_target_df$seqnames),'start'=min(test_target_df$start - 10000),'end'=max(test_target_df$end)+10000))
            
      # build loops
      target_loop <- cA_good_within_distance[
        cA_good_within_distance$subjectHits %in% markerdf$peakId & 
        cA_good_within_distance$queryHits %in% markerdf$peakId , ] # select loops landing with only those markers (in other words, we exclude common links and only plot the differential-peak-associated links, which should be differential links by nature. mainly for visualization purposes. this could be removed)
      target_loop$value <- target_loop$correlation
      target_loop <- as.data.frame(target_loop)
      target_loop$start_1 <- target_loop$start
      target_loop$end_1 <- target_loop$end
      target_loop$start <- rowMins(as.matrix(target_loop %>% dplyr::select(start_1,end_1)))
      target_loop$end <- rowMaxs(as.matrix(target_loop %>% dplyr::select(start_1,end_1)))
      loopsobj_to_plot <- list('CoAccessibility' = makeGRangesFromDataFrame(as.data.frame(target_loop) %>% dplyr::select(seqnames,start,end,value),keep.extra.columns = T)) 
      
      # estimate plot size 
      test_target_gr_size <- max(test_target_df$end)-min(test_target_df$start) + 20000
      customized_region_size <- round(test_target_gr_size/2000)
      
      # intersect the annotation regions to differential marker peaks, in order to ** only plot ** the DMR/ATAC diff region on differential peaks. In other words remove all annotations not associated with DE peaks. 
      CustomizedAnnotationSet <- AnnotationSet
      for (nameid in names(CustomizedAnnotationSet)){
        if(sum(findOverlaps(CustomizedAnnotationSet[[nameid]],test_target_peak_gr)@from)>0){
          CustomizedAnnotationSet[[nameid]] <- CustomizedAnnotationSet[[nameid]][findOverlaps(CustomizedAnnotationSet[[nameid]],test_target_peak_gr)@from,] # only include intersecting annotation region with DE peaks. 
        }else{
          if(sum(findOverlaps(CustomizedAnnotationSet[[nameid]],test_target_gr)@from)>0){
            CustomizedAnnotationSet[[nameid]] <- CustomizedAnnotationSet[[nameid]][- findOverlaps(CustomizedAnnotationSet[[nameid]],test_target_gr)@from,] # if no overlapping annotation region, remove the ones within visualization range. if nothing is overlapping in the range, do nothing (the if parameter)
          }
        }
        CustomizedAnnotationSet[[nameid]] <- flank(CustomizedAnnotationSet[[nameid]],width = customized_region_size * 5)
      }
      # checkpoint
      # unlist(lapply(CustomizedAnnotationSet,function(x) length(x)))
      # unlist(lapply(AnnotationSet,function(x) length(x)))
      
      # draw figure
      p <- plotBrowserTrack(
          ArchRProj = proj_arrow,
          useGroups = useGroups,
          sizes = c(4, 3, 3, 1),
          groupBy = "seurat_clusters",  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
          minCells = 100, region = test_target_gr,
          geneAnnotation = customizedAnnotation ,
          features = CustomizedAnnotationSet,threads = 10,
          loops = loopsobj_to_plot
      )
      # plot pdf
      plotPDF( p,
              name = paste0('Plot-Tracks-',geneid,'-defined_matrix_regulon',".pdf"),
              ArchRProj = proj_arrow,
              addDOC = FALSE, width = plotwidth, height = plotheight)
}

## Find louvain cluster, by defined cA object and a peakSinCA object
# this function returns a peakSinCA object which contains the louvain cluster the peak belongs to in the end. However, NOTE that the louvain cluster is chromosome specific. 

find_associated_ccan <- function(peakSinCA=peakSinCA,cA_good_within_distance=cA_good_within_distance){
  tempdflist <- list() 
  for (chrid in unique(peakSinCA$seqnames)){
    peakSinCA[peakSinCA$seqnames%in%chrid,] -> tempdf
    chr_all_peak <- (as.data.frame(peakSinCA[peakSinCA$seqnames%in%chrid,])$peakId)
    chr_tss_peak <- (as.data.frame(peakSinCA[peakSinCA$seqnames%in%chrid & peakSinCA$peakType %in% 'Promoter',])$peakId)
    dt <- as.data.frame(cA_good_within_distance[
      (cA_good_within_distance$queryHits %in% chr_all_peak & 
         cA_good_within_distance$subjectHits %in% chr_all_peak & 
         cA_good_within_distance$correlation>0.3) | 
        cA_good_within_distance$queryHits %in% chr_tss_peak | 
        cA_good_within_distance$subjectHits %in% chr_tss_peak,])[,c(1,2,4)]
    # build adjacency matrix (sparse) with igraph package
    g <- graph_from_data_frame(dt,directed = F)
    louvain_ccan <- cluster_louvain(g,dt$correlation)
    membership(louvain_ccan)-> membership_from_louvain_ccan
    as.data.frame(as.matrix(membership_from_louvain_ccan)) -> membership_from_louvain_ccan
    membership_from_louvain_ccan$peakId <- as.integer(as.character(rownames(membership_from_louvain_ccan)))
    colnames(membership_from_louvain_ccan)[1] <- 'louvain_cluster'
    tempdf <- left_join(tempdf,membership_from_louvain_ccan,by='peakId')
    tempdf -> tempdflist[[chrid]]
  }
  peakSinCA_return <- bind_rows(tempdflist)
  return(peakSinCA_return)
}


setwd('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/')
addArchRThreads(threads = 64)
addArchRGenome("hg38")

GTF <- read_tsv('/gpfs/genomedb/b38/gencode.v34.basic.annotation.gtf.gz',col_names = F,comment = '#')
GTF$symbol <- gsub('.+gene_name \\\"','',GTF$X9)
GTF <- separate(GTF,col='symbol',into=c('symbol','useless'))
GTF$useless <- NULL 

load('AnnotationSet.DMR.ATAC.bedfiles.all.list.Rdata')

library(igraph)


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

MARKERS <- unique(c(markerGenes,'TOX','MAF','TNFSF9','TNFRSF9','NRP1','CD8B','ID3','ID2','MAX','UPK1A','UPK1B','UPK3','NCAM1','SOX2','AR','ESRRA','KLRD1','KLRG1','KIR3DL1','FOXA1','FOXA2','NKG7','NCR3','AIRE','KLRB1','KLRK1','GZMM','GZMB','GZMA','VAV2','CD19','CD79A','CD8A','KIR2DL4','KIR2DL3','LYAR','LAYN','LY6A'),MARKERS)

meta_jw <- read.delim('/gpfs/output/20201016_Loom_files_for_Bladder_Cancer/analyse/modified.sample.meta.20201203.tsv')
meta_jw$Sample <- toupper(gsub('-.+','',meta_jw$sampleid))
meta_jw <- meta_jw %>% dplyr::select(donor,tissue_type,cancer_type,invasive,`T`,N,M,plot_T,sex,age,operation,tissue_plot,simplified_tissue,Sample) %>% unique
rownames(meta_jw) <- meta_jw$Sample


