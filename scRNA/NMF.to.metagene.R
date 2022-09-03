## Figure 1M
## Liger NMF
## Run on GPU01 c03b01n13

library(liger)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(dplyr)

dir.create('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/Figure1/Figure1M')
setwd('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/Figure1/Figure1M')


# Use Liger to perform iterative NMF on samples, to get NMF-decomposed metagenes

runOriginalLiger <- function (object, k, assay = NA, split.by = "orig.ident", lambda = 5,
                              thresh = 1e-06, max.iters = 30, reduction.name = "iNMF_raw",
                              reduction.key = "riNMF_", nrep = 1, H.init = NULL, W.init = NULL,
                              V.init = NULL, rand.seed = 1, print.obj = FALSE, ...)
{
  if(is.na(assay)){
    assay <- DefaultAssay(object = object)
  }
  if (is.character(x = split.by) && length(x = split.by) ==
      1) {
    split.by <- object[[split.by]]
  }
  split.cells <- split(x = colnames(x = object), f = split.by)
  scale.data <- lapply(X = split.cells, FUN = function(x) {
    return(t(x = GetAssayData(object = object, slot = "scale.data",
                              assay = assay)[, x]))
  })
  out <- liger::optimizeALS(object = scale.data, k = k, lambda = lambda,
                            thresh = thresh, max.iters = max.iters, nrep = nrep,
                            H.init = H.init, W.init = W.init, V.init = V.init, rand.seed = rand.seed,
                            print.obj = print.obj)
  colnames(x = out$W) <- VariableFeatures(object = object)
  return(out)
}



load('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/Figure7/Figure7D/Figure7D.new_combined_obj_sub_tissue.seurat.obj.for.cellchat.Rdata')
options(future.globals.maxSize=2000000000)
merged_dataset <- ScaleData(new_combined_obj_sub_tissue, split.by = "orig.ident", do.center = FALSE)
# merged_dataset <- RunOptimizeALS(merged_dataset, k = 20, lambda = 5, split.by = "orig.ident")
table(merged_dataset$orig.ident) -> tb1
names(tb1[tb1>=100]) -> names_to_retain # use only samples > 100 cells otherwise would hinder the analysis. 
# names_to_retain
merged_dataset_samples_good <- subset(merged_dataset,subset=orig.ident %in% names_to_retain)
# perform liger iNMF
liger_output_for_merged_dataset <- runOriginalLiger(merged_dataset_samples_good, k = 20, lambda = 5, split.by = "orig.ident")

save(liger_output_for_merged_dataset,file='Figure1L.NMF.liger.first.init.k20.l5.output.obj.Rdata')

# perform liger iNMF on early epithelial cluster (cancer)
merged_dataset_ca_only <- subset(merged_dataset_samples_good, subset=name_for_cellchat %in% c('Basal','Progenitor','Ca','TPCS','FPCS'))
table(merged_dataset_ca_only$orig.ident) -> tb1
names(tb1[tb1>=100]) -> names_to_retain
merged_dataset_samples_good_ca_only <- subset(merged_dataset_ca_only,subset=orig.ident %in% names_to_retain)
liger_output_for_merged_dataset_ca_only <- runOriginalLiger(merged_dataset_samples_good_ca_only, k = 20, lambda = 5)
liger_output_for_merged_dataset_ca_only <- runOriginalLiger(merged_dataset_samples_good_ca_only, k = 20, lambda = 5,split='orig.ident')
save(liger_output_for_merged_dataset_ca_only,file='Figure1M.ca.only.NMF.liger.first.init.k20.l5.output.obj.Rdata')

# fetch names for the weight vector.
# liger iNMF output is composed with three vectors:
# H: the decomposed per-cell metagene expression (metagene = NMF component, likewise as PCA PC component)
# W: each gene contribution to each metagene
# V: metagene variance (interesting thing to look at when you want to see how batch effect comes from)

colnames(liger_output_for_merged_dataset$W) <- rownames(GetAssayData(object = merged_dataset_samples_good,assay='RNA'))
colnames(liger_output_for_merged_dataset_ca_only$W) <- rownames(GetAssayData(object = merged_dataset_samples_good_ca_only,assay='RNA'))

# QC
lapply(c(1:20),function(x){
  sort(liger_output_for_merged_dataset$W[x,],decreasing=T) %>% head
})

# Extract cancer specific features of NMF
lapply(c(1:20),function(x){
  sort(liger_output_for_merged_dataset_ca_only$W[x,],decreasing=T) -> test1
  test1 <- test1[!grepl('\\-|\\.|^MT\\d|^RPL|^RPS',names(test1))]
  head(test1,20) %>% names 
}) %>% unlist %>% unique -> ca_only_feature_from_NMF
lapply(c(1:20),function(x){
  sort(liger_output_for_merged_dataset_ca_only$W[x,],decreasing=T) -> test1
  test1 <- test1[!grepl('\\-|\\.|^MT\\d|^RPL|^RPS',names(test1))]
  head(test1,20) %>% names 
}) -> ca_only_feature_list_from_NMF


## Plot NMF loadings

lapply(c(1:10),function(x){
  as.data.frame(liger_output_for_merged_dataset_ca_only$H[[x]])
})  %>%  bind_rows() -> liger_output_for_merged_dataset_ca_only_binded_H_loading_mtx
colnames(liger_output_for_merged_dataset_ca_only_binded_H_loading_mtx) <- paste0('NMF_',c(1:ncol(liger_output_for_merged_dataset_ca_only_binded_H_loading_mtx)))
liger_output_for_merged_dataset_ca_only_binded_H_loading_mtx$cell <- rownames(liger_output_for_merged_dataset_ca_only_binded_H_loading_mtx)
merged_dataset_samples_good_ca_only@meta.data %>% as.data.frame->temp
temp$cell <- rownames(temp)
left_join(temp,liger_output_for_merged_dataset_ca_only_binded_H_loading_mtx,by='cell') -> ca_only_liger
rownames(ca_only_liger) <- c(1:nrow(ca_only_liger))

pdf('Figure1M.Liger.NMF.loading.ca.pdf',height=20,width=20)
pheatmap::pheatmap(ca_only_liger[,grepl("NMF",colnames(ca_only_liger))],annotation_row = ca_only_liger %>% dplyr::select(name_for_cellchat,plot_T,tissue_type,sex,age),cellwidth = 10,cellheight = 0.02,scale = 'row',clustering_method = 'ward.D2',show_rownames = F)
dev.off() 


lapply(c(1:10),function(x){
  as.data.frame(liger_output_for_merged_dataset$H[[x]])
})  %>%  bind_rows() -> liger_output_for_merged_dataset_binded_H_loading_mtx
colnames(liger_output_for_merged_dataset_binded_H_loading_mtx) <- paste0('NMF_',c(1:ncol(liger_output_for_merged_dataset_binded_H_loading_mtx)))
liger_output_for_merged_dataset_binded_H_loading_mtx$cell <- rownames(liger_output_for_merged_dataset_binded_H_loading_mtx)
merged_dataset_samples_good@meta.data %>% as.data.frame->temp
temp$cell <- rownames(temp)
left_join(temp,liger_output_for_merged_dataset_binded_H_loading_mtx,by='cell') -> all_liger
rownames(all_liger) <- c(1:nrow(all_liger))
all_liger$celltype_rough <- gsub('_.+','',all_liger$name_for_cellchat)
all_liger$celltype_rough[all_liger$name_for_cellchat %in% c('Ca','FPCS','TPCS')] <- "Ca"

pdf('Figure1M.Liger.NMF.loading.boxviolin.pdf',height=7,width=40)
for (x in c(1:20)){
  nmfid <- paste0('NMF_',x)
  all_liger$tempNMF <- all_liger[,nmfid]
  ggplot(all_liger,aes(x=name_for_cellchat,y=tempNMF,fill=celltype_rough)) + geom_violin() + geom_boxplot(width=0.3,fill='black',outlier.alpha = 0) + theme_classic() + facet_grid(~celltype_rough,scale='free_x',space = 'free_x',switch = 'x')  + theme(legend.position = 'none',axis.text.x = element_text(angle=90)) + ggtitle(nmfid) -> p
  print(p)
}
dev.off() 

pdf('Figure1M.Liger.ca_only.NMF.loading.boxviolin.pdf',height=2.5,width=4)
for (x in c(1:20)){
  nmfid <- paste0('NMF_',x)
  ca_only_liger$tempNMF <- ca_only_liger[,nmfid]
  ggplot(ca_only_liger,aes(x=name_for_cellchat,y=tempNMF,fill=name_for_cellchat)) + geom_violin(scale='width') + geom_boxplot(width=0.3,fill='black',outlier.alpha = 0) + theme_classic()  + theme(legend.position = 'none',axis.text.x = element_text(angle=90)) + ggtitle(nmfid) -> p
  print(p)
}
dev.off() 

ca_only_liger <- ca_only_liger[,!grepl('\\.|pANN|RNA_snn',colnames(ca_only_liger))]
all_liger <- ca_only_liger[,!grepl('\\.|pANN|RNA_snn',colnames(all_liger))]



# Load ealy_clusters_epi cancer cell seurat obj
load('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021/./Figure4/Figure4F/Figure4F.early_clusters_epi.Rdata') # Stem cell, Basal, Ca, Ca-Stem cell
early_clusters_epi$mark_cluster_by_evolution <- gsub('Basal.+','Basal',early_clusters_epi$mark_cluster_by_evolution)
early_clusters_epi$mark_cluster_by_evolution[grepl('Fast Proliferating',early_clusters_epi$Clusters_detail)] <- 'FPCS'
early_clusters_epi$mark_cluster_by_evolution[grepl('Basal Ca',early_clusters_epi$Clusters_detail)] <- 'Ca'
early_clusters_epi$mark_cluster_by_evolution[grepl('Basal Progenitor',early_clusters_epi$Clusters_detail)] <- 'Progenitor'

# use addmodulescore to compute NMF vector for each cell
AddModuleScore(early_clusters_epi,ca_only_feature_list_from_NMF) -> early_clusters_epi

# Get markers for early_clusters_epi
early_clusters_epi$mark_cluster_by_evolution -> Idents(early_clusters_epi)
FindAllMarkers(early_clusters_epi,logfc.threshold = 0.5,min.pct = 0.3,min.diff.pct = 0.2) -> markers_test_ca_only

# plot the NMF vector loading on each cell (scatter)
pdf('Figure1M.NMF.Cancer.Only.Result.pdf',height=16,width=20)
FeaturePlot(early_clusters_epi,features=paste0('Cluster',c(1:20)),reduction = 'umap',cols=c('lightgrey','red'),pt.size=0.1,order = T) 
dev.off() 

# use RcisTarget to get motif enrichment in each NMF vector
library(RcisTarget)
data(motifAnnotations_hgnc)
motifRankings <- importRankings("/gpfs/genomedb/cisTopic/")
# Motif enrichment analysis:
lapply(c(1:20),function(x){
  sort(liger_output_for_merged_dataset_ca_only$W[x,],decreasing=T) -> test1
  test1 <- test1[test1>0]
  head(test1,100) %>% names 
}) -> ca_only_feature_list_from_NMF_for_cistarget
names(ca_only_feature_list_from_NMF_for_cistarget) <- paste0('NMF_',c(1:20))
lapply(names(ca_only_feature_list_from_NMF_for_cistarget),function(x){
  cat (x,'\n')
  tryCatch({
    t1 <- cisTarget(ca_only_feature_list_from_NMF_for_cistarget[[x]], motifRankings,
              motifAnnot=motifAnnotations_hgnc,verbose = T)
    return(t1)
  },error=function(e){cat ('Cannot work',x,'\n')})
}) -> motifEnrichmentTable_wGenes_list

names(motifEnrichmentTable_wGenes_list) <- paste0('NMF_',c(1:20))
# motifEnrichmentTable_wGenes_list <- motifEnrichmentTable_wGenes_list[motifEnrichmentTable_wGenes_list!=NULL]
lapply(names(motifEnrichmentTable_wGenes_list),function(x){
  motifEnrichmentTable_wGenes_list[[x]] -> y
  if(!is.null(y)){
    y$nmfid <- x   
    return(y)
  }
}) %>% bind_rows -> motifEnrichmentTable_wGenes_df
motifEnrichmentTable_wGenes_df <- motifEnrichmentTable_wGenes_df[motifEnrichmentTable_wGenes_df$NES>3,]
motifEnrichmentTable_wGenes_df <- left_join(motifEnrichmentTable_wGenes_df,motifEnrichmentTable_wGenes_df %>% dplyr::group_by(TF_highConf) %>% dplyr::summarise(TF_observed_n=n()))

motifEnrichmentTable_wGenes_df_specific <- motifEnrichmentTable_wGenes_df[motifEnrichmentTable_wGenes_df$TF_observed_n<=2,]
anotatedTfs  <- lapply(split(motifEnrichmentTable_wGenes_df_specific$TF_highConf,
                             motifEnrichmentTable_wGenes_df_specific$nmfid),
                       function(x) {
                         genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                         genesSplit <- unique(unlist(strsplit(genes, "; ")))
                         return(genesSplit)
                       })
anotatedTfs  

## Draw TF-gene network underlying NMF-features metagene 
for (nmfid in unique(motifEnrichmentTable_wGenes_df_specific$nmfid)){
  signifMotifNames <- motifEnrichmentTable_wGenes_df_specific$motif[motifEnrichmentTable_wGenes_df_specific$nmfid %in% nmfid][1:10]
  signifMotifMatchTF <- gsub(' .+','',motifEnrichmentTable_wGenes_df_specific$TF_highConf[motifEnrichmentTable_wGenes_df_specific$nmfid %in% nmfid][1:10])
  names(signifMotifMatchTF) <- signifMotifNames
  incidenceMatrix <- getSignificantGenes(ca_only_feature_list_from_NMF_for_cistarget[[nmfid]], 
                                         motifRankings,
                                         signifRankingNames=signifMotifNames,
                                         plotCurve=TRUE, maxRank=5000, 
                                         genesFormat="incidMatrix",
                                         method="aprox")$incidMatrix
  library(igraph)
  library(reshape2)
  edges <- melt(incidenceMatrix)
  edges <- edges[which(edges[,3]==1),1:2]
  colnames(edges) <- c("from","to")
  library(visNetwork)
  motifs <- unique(as.character(edges[,1]))
  genes <- unique(as.character(edges[,2]))
  nodes <- data.frame(id=c(motifs, genes),   
                      label=c(motifs, genes),    
                      title=c(motifs, genes), # tooltip 
                      shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                      color=c(rep("red", length(motifs)), rep("cornflowerblue", length(genes))))
  data.frame(label=as.character(signifMotifNames),new_name=as.character(signifMotifMatchTF)) -> temp
  temp$new_name <- as.character(temp$new_name)
  left_join(nodes,temp,by='label') -> nodes
  nodes$label <- as.character(nodes$label)
  nodes$new_name <- as.character(nodes$new_name)
  nodes$label[!is.na(nodes$new_name)] <- nodes$new_name[!is.na(nodes$new_name)]
  visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                          nodesIdSelection = TRUE,width = "300%",height = "300%") -> vn_obj
  visSave(vn_obj, file = paste0("Figure1M.VisNetwork.NMF.Features.after.Rcistopic.TF.gene.network.",nmfid,".html"), background = "white",selfcontained = T)
}


#### Combine NMF and findAllMarkers features, and make a combined 'marker gene' list for each cell type: Progenitor, Basal, Ca, FPCS, TPCS

# Ca:NMF_1, NMF_19, markers_ca
# Progenitor: NMF_7, NMF_14, NMF_20, markers_progenitor
# Basal: markers_basal
# FPCS: NMF_18, NMF_13, NMF_5, markers_FPCS
# TPCS: NMF_16, NMF_11, NMF_6, markers_NMF

geneList <- list(
  'Ca' = unique(c(ca_only_feature_list_from_NMF_for_cistarget[['NMF_1']] %>% head(50),
                ca_only_feature_list_from_NMF_for_cistarget[['NMF_19']] %>% head(50),
                markers_test_ca_only$gene[markers_test_ca_only$cluster %in% 'Ca']%>% head(50))),
  'Progenitor' = unique(c(ca_only_feature_list_from_NMF_for_cistarget[['NMF_7']] %>% head(50),
                  ca_only_feature_list_from_NMF_for_cistarget[['NMF_14']] %>% head(50),
                  ca_only_feature_list_from_NMF_for_cistarget[['NMF_20']] %>% head(50),
                  markers_test_ca_only$gene[markers_test_ca_only$cluster %in% 'Progenitor']%>% head(50)) ),
  'Basal' = unique(c(
                  markers_test_ca_only$gene[markers_test_ca_only$cluster %in% 'Basal'] %>% head(50))),
  'FPCS' = unique(c(ca_only_feature_list_from_NMF_for_cistarget[['NMF_18']] %>% head(50),
                  ca_only_feature_list_from_NMF_for_cistarget[['NMF_13']] %>% head(50),
                  ca_only_feature_list_from_NMF_for_cistarget[['NMF_5']] %>% head(50),
                  markers_test_ca_only$gene[markers_test_ca_only$cluster %in% 'FPCS'] %>% head(50))),
  'TPCS' = unique(c(ca_only_feature_list_from_NMF_for_cistarget[['NMF_16']] %>% head(50),
                  ca_only_feature_list_from_NMF_for_cistarget[['NMF_11']] %>% head(50),
                  ca_only_feature_list_from_NMF_for_cistarget[['NMF_6']] %>% head(50),
                  markers_test_ca_only$gene[markers_test_ca_only$cluster %in% 'TPCS'] %>% head(50))),
  'TPCS_NMF16_TM4' = unique(c(ca_only_feature_list_from_NMF_for_cistarget[['NMF_16']],'TM4SF1')),
  'TPCS_NMF11_TM4' = unique(c(ca_only_feature_list_from_NMF_for_cistarget[['NMF_11']],'TM4SF1')),
  'TPCS_NMF6_TM4' = unique(c(ca_only_feature_list_from_NMF_for_cistarget[['NMF_6']],'TM4SF1'))
                
)


# Motif enrichment analysis:

lapply(names(geneList),function(x){
  cat (x,'\n')
  tryCatch({
    t1 <- cisTarget(geneList[[x]], motifRankings,
                    motifAnnot=motifAnnotations_hgnc,verbose = T)
    return(t1)
  },error=function(e){cat ('Cannot work',x,'\n')})
}) -> motifEnrichmentTable_wGenes_list_combined

names(motifEnrichmentTable_wGenes_list_combined) <- names(geneList)

# Should we iterative with TF? 

# lapply(names(geneList),function(x){
#   cat (x,'\n')
#   tryCatch({
#     t1 <- cisTarget(geneList[[x]], motifRankings,
#                     motifAnnot=motifAnnotations_hgnc,verbose = T)
#     return(t1)
#   },error=function(e){cat ('Cannot work',x,'\n')})
# }) -> motifEnrichmentTable_wGenes_list_combined

names(motifEnrichmentTable_wGenes_list_combined) <- names(geneList)


lapply(names(motifEnrichmentTable_wGenes_list_combined),function(x){
  motifEnrichmentTable_wGenes_list_combined[[x]] -> y
  if(!is.null(y)){
    y$geneset <- x   
    return(y)
  }
}) %>% bind_rows -> motifEnrichmentTable_wGenes_df_combined
motifEnrichmentTable_wGenes_df_combined <- motifEnrichmentTable_wGenes_df_combined[motifEnrichmentTable_wGenes_df_combined]
motifEnrichmentTable_wGenes_df_combined <- left_join(motifEnrichmentTable_wGenes_df_combined,motifEnrichmentTable_wGenes_df_combined %>% dplyr::group_by(TF_highConf,geneset) %>% dplyr::summarise(TF_observed_n=n()))

motifEnrichmentTable_wGenes_df_combined_specific <- motifEnrichmentTable_wGenes_df_combined[motifEnrichmentTable_wGenes_df_combined$TF_observed_n<=2,]
anotatedTfs_combinedMarker  <- lapply(split(motifEnrichmentTable_wGenes_df_combined_specific$TF_highConf,
                             motifEnrichmentTable_wGenes_df_combined_specific$geneset),
                       function(x) {
                         genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                         genesSplit <- unique(unlist(strsplit(genes, "; ")))
                         return(genesSplit)
                       })
anotatedTfs_combinedMarker 


## Draw TF-gene network underlying NMF-features metagene 
for (genesetid in unique(motifEnrichmentTable_wGenes_df_combined_specific$geneset)){
  signifMotifNames <- motifEnrichmentTable_wGenes_df_combined_specific$motif[motifEnrichmentTable_wGenes_df_combined_specific$geneset %in% genesetid][1:10]
  signifMotifMatchTF <- gsub(' .+','',motifEnrichmentTable_wGenes_df_combined_specific$TF_highConf[motifEnrichmentTable_wGenes_df_combined_specific$geneset %in% genesetid][1:10])
  names(signifMotifMatchTF) <- signifMotifNames
  incidenceMatrix <- getSignificantGenes(geneList[[genesetid]], 
                                         motifRankings,
                                         signifRankingNames=signifMotifNames,
                                         plotCurve=TRUE, maxRank=5000, 
                                         genesFormat="incidMatrix",
                                         method="aprox")$incidMatrix
  library(igraph)
  library(reshape2)
  edges <- melt(incidenceMatrix)
  edges <- edges[which(edges[,3]==1),1:2]
  colnames(edges) <- c("from","to")
  library(visNetwork)
  motifs <- unique(as.character(edges[,1]))
  genes <- unique(as.character(edges[,2]))
  nodes <- data.frame(id=c(motifs, genes),   
                      label=c(motifs, genes),    
                      title=c(motifs, genes), # tooltip 
                      shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                      color=c(rep("red", length(motifs)), rep("cornflowerblue", length(genes))))
  data.frame(label=as.character(signifMotifNames),new_name=as.character(signifMotifMatchTF)) -> temp
  temp$new_name <- as.character(temp$new_name)
  left_join(nodes,temp,by='label') -> nodes
  nodes$label <- as.character(nodes$label)
  nodes$new_name <- as.character(nodes$new_name)
  nodes$label[!is.na(nodes$new_name)] <- nodes$new_name[!is.na(nodes$new_name)]
  visNetwork(nodes, edges) %>%visNetwork::visPhysics(enabled = TRUE,stabilization = 2000,solver = 'forceAtlas2Based') %>%  visOptions(highlightNearest = TRUE, 
                                          nodesIdSelection = TRUE,width = "600%",height = "600%") -> vn_obj
  visSave(vn_obj, file = paste0("Figure1M.VisNetwork.combined.Features.after.Rcistopic.TF.gene.network.",genesetid,".html"), background = "white",selfcontained = T)
}

#### 

#### Save


save(list=c("motifEnrichmentTable_wGenes",
            "motifEnrichmentTable_wGenes_df",
            "motifEnrichmentTable_wGenes_df_combined",
            "motifEnrichmentTable_wGenes_df_combined_specific",
            "motifEnrichmentTable_wGenes_df_specific",
            "motifEnrichmentTable_wGenes_list",
            "motifEnrichmentTable_wGenes_list_combined",'ca_only_liger','all_liger'),file='Figure1M.NMF.motifEnrichmentTable.wGenes.df.objs.Rdata')

write.table(motifEnrichmentTable_wGenes_df_combined,file='Figure1/Figure1M/cisTarget.result.for.NMF.TPCS.Ca.20210621.tsv',sep="\t",quote=F,row.names=F)


## TPCS cell cycle pre- and post-chemo



## Study methylation-dependent TF binding

motifEnrichmentTable_wGenes_df %>% dplyr::group_by(nmfid) %>% dplyr::summarise(TFBS_Taipale=sum(grepl('taipale',motif)),MethylatedTFBS_Taipale=sum(grepl('meth$',motif)),MethRepressTFBS_Taipale = sum(grepl('meth_repr$',motif)), Meth_Ratio=sum(grepl('meth$',motif))/sum(grepl('taipale',motif)),MethRepressRatio=sum(grepl('meth_repr$',motif))/sum(grepl('taipale',motif)))


## Study TPCS specific motif, from ATAC

ATAC_motifs <- c(
  'ASCL1','ASCL2','BACH2','BACH1','FOSL2','JDP2','JUN','FOS','FOSB','FOSL1','JUND','JUNB','SMARCC1',
)


motifEnrichmentTable_wGenes_df_combined %>% dplyr::group_by(geneset) %>% dplyr::summarise(TFBS_Taipale=sum(grepl('taipale',motif)),MethylatedTFBS_Taipale=sum(grepl('meth$',motif)),Meth_Ratio=sum(grepl('meth$',motif))/sum(grepl('taipale',motif)))

