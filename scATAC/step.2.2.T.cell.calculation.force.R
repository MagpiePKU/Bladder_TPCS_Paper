
proj_subtype <- addIterativeLSI( ArchRProj = proj_subtype,useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 4,clusterParams = list( resolution = c(0.1, 0.2, 0.4,1), sampleCells = 10000, n.start = 10), varFeatures = 50000, dimsToUse = 1:30,force = T)
# Add cluster based on LSI
proj_subtype <- addClusters(input = proj_subtype,reducedDims = "IterativeLSI",method = "Seurat",name = "seurat_clusters",resolution = 1.5,force=T)
proj_subtype <- addImputeWeights(proj_subtype,reducedDims = "IterativeLSI",sampleCells=30000,nRep=3)
# Add embeding
proj_subtype <- addUMAP(ArchRProj = proj_subtype, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 10, minDist = 0.1,metric = "cosine",force = T)
proj_subtype <- addTSNE(ArchRProj = proj_subtype, reducedDims = "IterativeLSI", name = "TSNE",   perplexity = 30,force = T)
# Add markers
markers_subtype <- getMarkerFeatures(ArchRProj = proj_subtype, useMatrix = "GeneScoreMatrix", groupBy = "seurat_clusters",  bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")
markerList <- getMarkers(markers_subtype, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
# Add grouped coverage
proj_subtype <- addGroupCoverages(ArchRProj = proj_subtype, groupBy = "seurat_clusters",force = T)
# Add peaks
proj_subtype <- addReproduciblePeakSet(
  ArchRProj = proj_subtype,
  groupBy = "seurat_clusters",force = T,
  pathToMacs2 = '/gpfs/bin/anaconda3/bin/macs2' # Note this is a new installation, can run with any account now.
)
proj_subtype <- addPeakMatrix(proj_subtype,force = T)
# Add co-accessibility
proj_subtype <- addCoAccessibility(
  ArchRProj = proj_subtype,corCutOff = 0.3,maxDist = 500000,
  reducedDims = "IterativeLSI"
)
cA_T_cell <- getCoAccessibility(
  ArchRProj = proj_subtype,
  corCutOff = 0.3,
  resolution = 50,
  returnLoops = FALSE
)
cA_good_within_distance <- cA_T_cell[cA_T_cell$FDR<(1/nrow(as.data.frame(cA_T_cell))),] # Select FDR unlikely (more than one in the dataset)

peakMatrix <- as.data.frame(proj_subtype@peakSet,row.names=NULL)
unique(c(cA_good_within_distance$queryHits,cA_good_within_distance$subjectHits)) -> hits_row
peakSinCA <- peakMatrix[hits_row,]
peakSinCA_gr <- makeGRangesFromDataFrame(peakSinCA,keep.extra.columns = T)
peakSinCA_return <- find_associated_ccan(peakSinCA,cA_good_within_distance)
peakSinCA_return_gr <- makeGRangesFromDataFrame(peakSinCA_return,keep.extra.columns = T)

AnnotationSet_registered_with_CANN <- lapply(AnnotationSet,function(x){
  x[findOverlaps(x,peakSinCA_return_gr)@from,] -> xgr
  peakSinCA_return_gr[findOverlaps(x,peakSinCA_return_gr)@to,] -> ygr
  xgr$louvain_cluster <- ygr$louvain_cluster
  xgr$potential_gene <- ygr$nearestGene
  xgr$potential_type <- ygr$peakType
  xgr$distToTSS_potentialgene <- ygr$distToTSS
  xdf <- as.data.frame(xgr)
  xdf <- left_join(xdf,xdf %>% dplyr::group_by(potential_gene) %>% dplyr::summarise(mindist=min(distToTSS_potentialgene)))
  xgr$mindist <- xdf$mindist
  xgr$is_promoter <- ifelse(xgr$potential_type %in% 'Promoter' | xgr$distToTSS_potentialgene <= xgr$mindist + 400 ,1,0)
  xgr_good <- xgr[xgr$is_promoter==1,]
  return(xgr_good)
})

peakSinCA_temp <- left_join(peakSinCA,peakSinCA %>% dplyr::group_by(nearestGene) %>% dplyr::summarise(mindist=min(distToTSS)))

NewAnnotationSet <- AnnotationSet_registered_with_CANN

cell_type_DE_peaks_with_DMR <- lapply(markerList,function(x){
  list_return = lapply(NewAnnotationSet, function(y){
    x <- x[findOverlaps(x,y)@from,]
    xdf <- as.data.frame(x)
    xdf <- left_join(xdf,peakSinCA_temp %>% dplyr::select(seqnames,start,end,peakId),by = c("seqnames", "start", "end"))
    xdfgr <- makeGRangesFromDataFrame(xdf,keep.extra.columns = T)
    return(xdfgr)
  })
  return(list_return)
})

# 2. overlap the DE-peaks with annotation to interaction to find essential (non-false-discovery) peaks
non_redundant_coaccessibility_network_peaks <- unique(c(cA_good_within_distance$queryHits,cA_good_within_distance$subjectHits))
cell_type_DE_peaks_with_DRPS <- lapply(cell_type_DE_peaks_with_DMR,function(x){
  list_return = lapply(x,function(y){
    y <- y[y$peakId %in% non_redundant_coaccessibility_network_peaks,]
    return(y)
  })
  return(list_return)
})

cell_type_DE_peaks_with_DRPS[[4]][[8]] -> CD4DRPS
cell_type_DE_peaks_with_DRPS[[3]][[8]] -> CD8DRPS

AnnotationSet_registered_with_CANN[[8]][findOverlaps(AnnotationSet_registered_with_CANN[[8]],CD4DRPS)@from,] -> CD4DRPSanno
AnnotationSet_registered_with_CANN[[8]][findOverlaps(AnnotationSet_registered_with_CANN[[8]],CD8DRPS)@from,] -> CD8DRPSanno

tb <- table(CD4DRPSanno$potential_gene)
tb[tb>3]
tb <- table(CD8DRPSanno$potential_gene)
tb[tb>3]


