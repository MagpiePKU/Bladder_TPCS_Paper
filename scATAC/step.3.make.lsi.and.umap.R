
proj_arrow <- addIterativeLSI( ArchRProj = proj_arrow,useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 4,clusterParams = list( resolution = c(0.1, 0.2, 0.4,1), sampleCells = 10000, n.start = 10), varFeatures = 50000, dimsToUse = 1:30)
# Add cluster based on LSI
proj_arrow <- addClusters(input = proj_arrow,reducedDims = "IterativeLSI",method = "Seurat",name = "seurat_clusters",resolution = 1,force=T)
proj_arrow <- addImputeWeights(proj_arrow,reducedDims = "IterativeLSI",sampleCells=30000,nRep=3)
# Add embeding
proj_arrow <- addUMAP(ArchRProj = proj_arrow, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5,metric = "cosine")
proj_arrow <- addTSNE(ArchRProj = proj_arrow, reducedDims = "IterativeLSI", name = "TSNE",   perplexity = 30)

p1 <- plotEmbedding(ArchRProj = proj_subtype, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p1 <- plotEmbedding(ArchRProj = proj_subtype, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = proj_subtype, colorBy = "colData", name = c("TSSEnrichment","PromoterRatio","DoubletScore",'seurat_clusters'),   embedding = "UMAP")
pdf(paste0(workingpath,'/QualityControl/UMAP.pdf'),height=8,width=8)
print(p1)
print(p2)
dev.off()

