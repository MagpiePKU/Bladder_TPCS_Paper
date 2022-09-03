pmagic <- plotEmbedding(ArchRProj = proj_arrow,colorBy = "GeneScoreMatrix",name = unique(c('ZEB1','ZEB2','SNAI1','SNAI2','FOXA1','SOX2','VIM','ERBB2','EGFR','CDKN2A','TBX5','TBX3','TBX6','TCAF1','HOPX','HOXB9','FGF19','CCND1','CCNE1','MYC','MYCN','KRAS','NRAS','NFATC1','CTCF','FOXF1','PCDHB5','CLDN11','TRIM31','HOXB5','WNT16','TERC','TERT','FGFR3','HOXA5','LIF','LIFR','EWSR1','IFITM5','E2F3','ID4','CDKAL1','TNF','VIM-AS1','IRX2'),c('ZEB1','DAB2IP','FOXA1','CCND1','TBX3','MYCL','CDKN2A','TP53','SOX2','ELF3','SOS1')),embedding = "UMAP",imputeWeights = getImputeWeights(proj_arrow))
pdf(paste0(workingpath,'/QualityControl/UMAP.MAGIC.pdf'),height=5,width=5)
print(pmagic)
dev.off()


