
proj_arrow <- addIterativeLSI( ArchRProj = proj_arrow,useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 3,clusterParams = list( resolution = c(0.1, 0.2, 0.4,1), sampleCells = 10000, n.start = 10), varFeatures = 50000, dimsToUse = 1:30,force = T)
# Add cluster based on LSI
proj_arrow <- addClusters(input = proj_arrow,reducedDims = "IterativeLSI",method = "Seurat",name = "seurat_clusters",resolution = 1.2,force=T)
proj_arrow <- addImputeWeights(proj_arrow,reducedDims = "IterativeLSI",sampleCells=50000,nRep=3)
# Add embeding
proj_arrow <- addUMAP(ArchRProj = proj_arrow, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5,metric = "cosine",force = T)
proj_arrow <- addTSNE(ArchRProj = proj_arrow, reducedDims = "IterativeLSI", name = "TSNE",   perplexity = 30,force = T)

p1 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
p2 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "colData", name = c("TSSEnrichment","PromoterRatio","DoubletScore",'seurat_clusters'),   embedding = "UMAP")
pdf(paste0(workingpath,'/QualityControl/UMAP.pdf'),height=8,width=8)
print(p1)
print(p2)
print(p3)
dev.off()
saveArchRProject(ArchRProj = proj_arrow, load = FALSE)


p1 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggplot(data=p1$data,mapping=aes(x=x,y=y)) + geom_density2d(linetype='dashed',size=0.5) + geom_point(alpha=0.2) + facet_wrap(~color,ncol=8) -> p_per_sample_distr
pdf(paste0(workingpath,'/QualityControl/UMAP.density.per.sample.pdf'),height=50,width=56)
print(p_per_sample_distr)
dev.off() 

pmagic <- plotEmbedding(ArchRProj = proj_arrow,colorBy = "GeneScoreMatrix",name = unique(c('ZEB1','ZEB2','SNAI1','SNAI2','FOXA1','SOX2','VIM','ERBB2','EGFR','CDKN2A','TBX5','TBX3','TBX6','TCAF1','HOPX','HOXB9','FGF19','CCND1','CCNE1','MYC','MYCN','KRAS','NRAS','NFATC1','CTCF','FOXF1','PCDHB5','CLDN11','TRIM31','HOXB5','WNT16','TERC','TERT','FGFR3','HOXA5','LIF','LIFR','EWSR1','IFITM5','E2F3','ID4','CDKAL1','TNF','VIM-AS1','IRX2','ZEB1','DAB2IP','FOXA1','CCND1','TBX3','MYCL','CDKN2A','TP53','SOX2','ELF3','SOS1','CD3D','NCAM1','CD4','CD8A','IFNG','CD163','MS4A1','CD19','CD79A','SALL1','GZMB','GZMK','CD3D','CD8A','CD4','FOXP3','NCAM1','CD19','CD79A','MS4A1','TNFSF9','CD3D','NCAM1','CD4','CD8A','FOXP3','AIRE','CD79A','SELL','CCR7','PDCD1','CD274','EOMES','NRP1','TOX','TBX21','CD38','AAK1','HAVCR2','LAG3','MAF','MAX','NRP1')),embedding = "UMAP",imputeWeights = getImputeWeights(proj_arrow),threads = 3)
pdf(paste0(workingpath,'/QualityControl/UMAP.MAGIC.pdf'),height=5,width=5)
print(pmagic)
dev.off()

pmagic <- plotEmbedding(ArchRProj = proj_arrow,colorBy = "GeneScoreMatrix",name = unique(c('KCNJ8','EDNRB','HAND2','TNF','DAB2IP','SOX14','APOC2','ZIM3','WNT16','SOX15','FOXA1','MMP9','KRT7','COL6A2','COL6A1','PTPRB','PTPRC','CD34','CD3E','PLVAP','PECAM1','SALL4','AFP','NANOG','POU5F1','MYOG','MYOD1','DES','VIM','LUM','DCN')),embedding = "UMAP",imputeWeights = getImputeWeights(proj_arrow),threads = 3)
pdf(paste0(workingpath,'/QualityControl/UMAP.MAGIC.stage.pdf'),height=5,width=5)
print(pmagic)
dev.off()
