workingpath='/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/Save/Bladder_20201225/'


p1 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "colData", name = c("TSSEnrichment","PromoterRatio","DoubletScore",'seurat_clusters'),   embedding = "UMAP")
pdf(paste0(workingpath,'/QualityControl/UMAP.pdf'),height=8,width=8)
print(p1)
print(p3)
print(p2)
dev.off()

