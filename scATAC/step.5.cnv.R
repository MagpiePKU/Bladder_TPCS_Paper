source('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/bin/ArchR_CNV.R')
proj_subtype <- addCNVMatrix(proj_arrow,threads=8,force=T)
CNVMat <- getMatrixFromProject(proj_arrow, "CNVMatrix")
save(CNVMat,file='/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/Save/Bladder_20201225//CNVmat.Rdata')
