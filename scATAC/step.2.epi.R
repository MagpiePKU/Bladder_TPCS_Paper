proj_subtype <- addIterativeLSI( ArchRProj = proj_subtype,useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 4,clusterParams = list( resolution = c(0.1, 0.2, 0.4,1), sampleCells = 30000, n.start = 10), varFeatures = 50000, dimsToUse = 1:30,force = T,threads=5) # epithelial
# Add cluster based on LSI
proj_subtype <- addClusters(input = proj_subtype,reducedDims = "IterativeLSI",method = "Seurat",name = "seurat_clusters",resolution = 1.5,maxClusters = NULL,force=T)
proj_subtype <- addImputeWeights(proj_subtype,reducedDims = "IterativeLSI",sampleCells=30000,nRep=3)
# Add embeding
proj_subtype <- addUMAP(ArchRProj = proj_subtype, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 10, minDist = 0.1,metric = "cosine",force = T)
proj_subtype <- addTSNE(ArchRProj = proj_subtype, reducedDims = "IterativeLSI", name = "TSNE",   perplexity = 30,force = T)

dir.create(paste0(workingpath,'QualityControl'))
pdf(paste0(workingpath,'QualityControl/UMAP.pdf'),height=12,width=12)
p1 <- plotEmbedding(ArchRProj = proj_subtype, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_subtype, colorBy = "colData", name = "seurat_clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj_subtype, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p4 <- plotEmbedding(ArchRProj = proj_subtype, colorBy = "colData", name = "seurat_clusters", embedding = "TSNE")
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

plotgenes <- unique(c('CD3D','CD8A','CD4','FOXP3','NCAM1','CD19','CD79A','MS4A1','CD163','NCAM1','PDCD1','TOX','EOMES','CD274','AIRE','SELL','TCF7','MAX','ID2','ID3','MAF','TBX21','CD81','CD160','IFNG','TNFRSF9','TNFRSF18','KLF14','ADAM33','CD80','CR1','CREBBP','CSF2','ICAM3','IL1RL2','ITPKA','MAP4K2','MST1R','TAL1','TANK','TNFRSF8','TNFSF14','FASLG','GNLY','GZMM','RUNX3','SELPLG','STAT4','FYN','APP','BTLA','C3AR1','CCR7','CXCR5','DPP4','FAS','HDAC1','ITGA6','LTB','MIR4323','TIRAP','CD244','CCL5','APOBR','FLT3LG','IFI16','IFITM1','MAFF','IL2RB','IFNL1','PRF1','PIK3CD','SIGIRR','TAP1','ZBTB7A','YPEL1','APOBR','APP','BCL2','BTLA','C2','CCR4','CCR7','CD27','CD4','CD55','CPB1','CTLA4','CXCR5','DPP4','FAS','FOS','ICOS','IKBKE','IL17RA','IL23A','IL2RA','IL2RB','IL6R','ITGA1','JAML','LCP1','MAP3K1','MAPK13','NFATC1','PSMB10','SLAMF1','TLR5','STAT1','TNFSF4','TNFSF8','SMAD3','TXK'))

p <- plotEmbedding(ArchRProj = proj_subtype,colorBy = "GeneScoreMatrix",name = plotgenes,embedding = "TSNE",imputeWeights = getImputeWeights(proj_subtype))
p2 <- plotEmbedding(ArchRProj = proj_subtype,colorBy = "GeneScoreMatrix",name = plotgenes,embedding = "UMAP",imputeWeights = getImputeWeights(proj_subtype))

pdf(paste0(workingpath,'QualityControl/Diagnostic.MAGIC.gene.expression.marker.pdf'),width=5,height=5)
for (id in names(p)){print(p[[id]])}
for (id in names(p2)){print(p2[[id]])}
dev.off()


# Add markers
markers_subtype <- getMarkerFeatures(ArchRProj = proj_subtype, useMatrix = "GeneScoreMatrix", groupBy = "seurat_clusters",  bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",threads = 5)
markerList <- getMarkers(markers_subtype, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markers_peak_subtype <- getMarkerFeatures(ArchRProj = proj_subtype, useMatrix = "PeakMatrix", groupBy = "seurat_clusters",  bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon",threads = 5)
markerList_peak <- getMarkers(markers_peak_subtype, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")


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
cA <- getCoAccessibility(
  ArchRProj = proj_subtype,
  corCutOff = 0.3,
  resolution = 50,
  returnLoops = FALSE
)
cA_good_within_distance <- cA[cA$FDR<(1/nrow(as.data.frame(cA))),] # Select FDR unlikely (more than one in the dataset)

peakMatrix <- as.data.frame(proj_subtype@peakSet,row.names=NULL)
peakMatrix$peakId <- c(1:nrow(peakMatrix))
rownames(peakMatrix) <- peakMatrix$peakId
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

cell_type_DE_peaks_with_DMR <- lapply(markerList_peak,function(x){
  list_return = lapply(NewAnnotationSet, function(y){
    tryCatch({
      # x$startOld <- x$start
    # x$endOld <- x$end
    # x <- as.data.frame(x,row.names = NULL) %>% dplyr::select(-strand) %>% unique
    # x$start <- rowMins(as.matrix(x%>%dplyr::select(startOld,endOld)))
    # x$end <- rowMaxs(as.matrix(x%>%dplyr::select(startOld,endOld)))
    x <- makeGRangesFromDataFrame(as.data.frame(x,row.names = NULL) %>% dplyr::select(seqnames,start,end) %>% unique,keep.extra.columns = F)
    x <- x[unique(findOverlaps(x,y)@from),]
    xdf <- as.data.frame(x)
    xdf <- left_join(xdf,peakSinCA_temp %>% dplyr::select(seqnames,start,end,peakId),by = c("seqnames", "start", "end"))
    xdfgr <- makeGRangesFromDataFrame(xdf,keep.extra.columns = T)
    return(xdfgr)
    },error=function(e){cat ('')})
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



