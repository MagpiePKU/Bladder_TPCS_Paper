
idxSample <- BiocGenerics::which(proj_arrow$cancer_type %in% c('DCD','bladder_cancer') & proj_arrow$tissue_type %ni% c('normal_kidney_tissue','normal_prostate_tissue'))
cellsSample <- proj_arrow$cellNames[idxSample]
proj_subtype = proj_arrow[cellsSample, ]

saveArchRProject(ArchRProj = proj_subtype, outputDirectory = "/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/Save/Bladder_epithelial_and_cancer_pure_20201231", load = TRUE) -> proj_subtype
proj_arrow <- proj_subtype
rm(proj_subtype)

workingpath = '/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/Save/Bladder_epithelial_and_cancer_pure_20201231/'
dir.create(paste0(workingpath,'/QualityControl/'))
dir.create(paste0(workingpath,'/Plots/'))

proj_arrow <- addIterativeLSI( ArchRProj = proj_arrow,useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 3,clusterParams = list( resolution = c(0.1, 0.2, 0.4,1), sampleCells = 10000, n.start = 10), varFeatures = 100000, dimsToUse = 1:30,force = T)
# Add cluster based on LSI
proj_arrow <- addClusters(input = proj_arrow,reducedDims = "IterativeLSI",method = "Seurat",name = "seurat_clusters",resolution = 1.2,force=T,dimsToUse = 3:10)
proj_arrow <- addImputeWeights(proj_arrow,reducedDims = "IterativeLSI",sampleCells=10000,nRep=3)
# Add embeding
proj_arrow <- addUMAP(ArchRProj = proj_arrow, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5,metric = "cosine",force = T,dimsToUse = 3:10)
proj_arrow <- addTSNE(ArchRProj = proj_arrow, reducedDims = "IterativeLSI", name = "TSNE",   perplexity = 30,force = T)

p1 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "colData", name = c("plot_T","tissue_type","cancer_type",'seurat_clusters','sex','age','operation','donor',"TSSEnrichment","PromoterRatio","DoubletScore"),   embedding = "UMAP")
pdf(paste0(workingpath,'/QualityControl/UMAP.pdf'),height=8,width=8)
print(p1)
print(p2)
dev.off()
saveArchRProject(ArchRProj = proj_arrow, load = FALSE)


p1 <- plotEmbedding(ArchRProj = proj_arrow, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

ggplot(data=p1$data,mapping=aes(x=x,y=y)) + geom_density2d(linetype='dashed',size=0.5) + geom_point(alpha=0.2) + facet_wrap(~color,ncol=8) + theme_classic() -> p_per_sample_distr
pdf(paste0(workingpath,'/QualityControl/UMAP.density.per.sample.pdf'),height=50,width=56)
print(p_per_sample_distr)
dev.off() 

pmagic <- plotEmbedding(ArchRProj = proj_arrow,colorBy = "GeneScoreMatrix",name = unique(c('KCNJ8','EDNRB','HAND2','TNF','DAB2IP','SOX14','APOC2','ZIM3','WNT16','SOX15','FOXA1','MMP9','KRT7','COL6A2','COL6A1','PTPRB','PTPRC','CD34','CD3E','PLVAP','PECAM1','SALL4','AFP','NANOG','POU5F1','MYOG','MYOD1','DES','VIM','LUM','DCN')),embedding = "UMAP",imputeWeights = getImputeWeights(proj_arrow),threads = 3)
pdf(paste0(workingpath,'/QualityControl/UMAP.MAGIC.stage.pdf'),height=5,width=5)
print(pmagic)
dev.off()

