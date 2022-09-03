rm(list=ls())
setwd('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/')
addArchRThreads(threads = 64)
addArchRGenome("hg38")

inputFiles <- read.delim('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/Processed_Data/ATAC/fragments.list',header=F)
inputFiles <- inputFiles$V1 
inputFiles <- as.character(inputFiles)
inputFiles <- as.list(inputFiles)
names_of_inputfiles <- gsub('-.+','',gsub('.+/','',gsub('/outs.+','',unlist(inputFiles))))
names(inputFiles) <- names_of_inputfiles
print(names_of_inputfiles)

ArrowFiles <- createArrowFiles(
  inputFiles = as.character(as.vector(inputFiles)),
  sampleNames = names(inputFiles),
  minTSS = 3,
  minFrags = 5000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1,threads = 10
)

workingpath='Bladder_20201225'

proj_arrow <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Bladder_20201225",
  copyArrows = TRUE
)

proj_arrow <- filterDoublets(proj_arrow)

paste0("Memory Size = ", round(object.size(proj_arrow) / 10^6, 3), " MB")

# Diagnostics

p1 <- plotGroups( ArchRProj = proj_arrow, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
p2 <- plotFragmentSizes(ArchRProj = proj_arrow)
p3 <- plotTSSEnrichment(ArchRProj = proj_arrow)

dir.create('QualityControl')
pdf('QualityControl/DiagnosticPlot.pdf',height=5,width=15)
(p1|p2|p3)
dev.off()

# Saving temp
dir.create('Save')
saveArchRProject(ArchRProj = proj_arrow, outputDirectory = "Save/Bladder_20201225", load = FALSE)
