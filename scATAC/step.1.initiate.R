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
load('/gpfs/genomedb/singleR/ref_Human_all.RData')

options(bitmapType='cairo') # solve the PNG/X11 issue

# proj_arrow <- loadArchRProject(path="Save/20201018.save.3.ATAC.ArchR.Object") # this is when you load the file instead of re-running

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
  LSIMethod = 1
)

proj_arrow <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "Bladder_20201126",
  copyArrows = TRUE
)

proj_arrow <- filterDoublets(proj_arrow)
