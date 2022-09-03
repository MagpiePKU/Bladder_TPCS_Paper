
source('/gpfs/output/SingleCell_Placenta/analyse/bins/environment_setup.R') 
library(slingshot)
library(CytoTRACE)

setwd('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/Paper/TPCS_2021')
load('./New/RNA.full.Data.20210421.ZY.Rdata')
ArchR_matched <- loadArchRProject('/gpfs/output/Bladder_Cancer_Landscape_Project/SingleCell/analysis/ATAC/Save/epithelial_drps_origin_20210218')

ArchR_matched@cellColData%>%as.data.frame -> meta_from_atac
total_cell_for_paper_epi <- subset(total_cell_for_paper,subset=final_annotation_zy %in% c('Urothelial_Im/Um','Urothelial_Basal','Urothelial_Progenitor','Ca_LuminalLike','Ca_FPCS','Ca_Basal','Ca_TPCS'))

total_cell_for_paper_epi@meta.data %>% as.data.frame -> meta_from_RNA
meta_from_RNA$predictedCell_RNA <- meta_from_RNA$cell
total_cell_for_paper_epi$epi_TPCS <- total_cell_for_paper_epi$cell %in% meta_from_atac$predictedCell_RNA[meta_from_atac$simpleType %in% "TPCS"]
total_cell_for_paper_epi$epi_Luminal <- total_cell_for_paper_epi$cell %in% meta_from_atac$predictedCell_RNA[meta_from_atac$simpleType %in% "Luminal"]
total_cell_for_paper_epi$epi_Basal <- total_cell_for_paper_epi$cell %in% meta_from_atac$predictedCell_RNA[meta_from_atac$simpleType %in% "Basal"]
total_cell_for_paper_epi$epi_Normal <- total_cell_for_paper_epi$cell %in% meta_from_atac$predictedCell_RNA[meta_from_atac$simpleType %in% "Normal"]

GetAssayData(total_cell_for_paper_epi,assay='RNA') -> mtx
Var_mtx <- mtx[Seurat::VariableFeatures(total_cell_for_paper_epi),]
Var_mtx <- as.matrix(Var_mtx)
CytoTRACE::CytoTRACE(Var_mtx,ncores=10) -> res

total_cell_for_paper_epi <- Seurat::CellCycleScoring(total_cell_for_paper_epi,s.features = cc.genes.updated.2019$s.genes,g2m.features = cc.genes.updated.2019$g2m.genes)
data.frame(CytoTRACE=res$CytoTRACE,CytoTRACErank=res$CytoTRACErank,cell=names(res$CytoTRACE)) -> df1
left_join(df1,as.data.frame(total_cell_for_paper_epi@meta.data)) -> df1

df1%>% dplyr::group_by(epi_TPCS) %>% dplyr::summarise(trace=mean(CytoTRACE),sd_trace=sd(CytoTRACE))
df1%>% dplyr::group_by(epi_Basal) %>% dplyr::summarise(trace=mean(CytoTRACE),sd_trace=sd(CytoTRACE))
df1%>% dplyr::group_by(epi_Luminal) %>% dplyr::summarise(trace=mean(CytoTRACE),sd_trace=sd(CytoTRACE))
df1%>% dplyr::group_by(epi_Normal) %>% dplyr::summarise(trace=mean(CytoTRACE),sd_trace=sd(CytoTRACE))

save(df1,file='./New/cytotrace.variable.features.only.RNA.full.Data.20210421.ZY.Rdata')

df1$class_simple_for_plot_cytotrace <- ''
df1$class_simple_for_plot_cytotrace[df1$epi_Basal] <- 'Basal'
df1$class_simple_for_plot_cytotrace[df1$epi_Luminal] <- 'Luminal'
df1$class_simple_for_plot_cytotrace[df1$epi_TPCS] <- 'TPCS'

df1$class_simple_for_plot_cytotrace <- factor(df1$class_simple_for_plot_cytotrace,levels=c('TPCS','Luminal','Basal',''))

df1%>% dplyr::group_by(class_simple_for_plot_cytotrace) %>% dplyr::summarise(trace=mean(CytoTRACE),sd_trace=sd(CytoTRACE))

dir.create('result')
pdf('result/CytoTRACE.cells.cancer.pdf',height=3,width=7)
ggplot(df1[df1$class_simple_for_plot_cytotrace %ni% '',],aes(x=class_simple_for_plot_cytotrace,y=CytoTRACE,fill=class_simple_for_plot_cytotrace)) + geom_violin(scale='width') + geom_boxplot(outlier.alpha = 0,width=0.3,fill='black',color='black',size=0.5) + theme_classic() + ylim(c(0,1)) + ggsci::scale_fill_jco()
dev.off() 

df1$class_simple_for_plot_cytotrace <- NULL
df1$class_simple_for_plot_cytotrace <- ''
df1$class_simple_for_plot_cytotrace[df1$epi_Basal | df1$final_annotation_zy %in% 'Ca_Basal'] <- 'Ca_Basal'
df1$class_simple_for_plot_cytotrace[df1$epi_Luminal | df1$final_annotation_zy %in% 'Ca_LuminalLike'] <- 'Ca_LuminalLike'
df1$class_simple_for_plot_cytotrace[df1$final_annotation_zy %in% 'Ca_FPCS'] <- 'Ca_FPCS'
df1$class_simple_for_plot_cytotrace[grepl('Uro',df1$final_annotation_zy)] <- df1$final_annotation_zy[grepl('Uro',df1$final_annotation_zy)]
df1$class_simple_for_plot_cytotrace[df1$epi_TPCS] <- 'Ca_TPCS'
df1%>% dplyr::group_by(class_simple_for_plot_cytotrace) %>% dplyr::summarise(trace=mean(CytoTRACE),sd_trace=sd(CytoTRACE)) %>% arrange(trace) -> df2
df1$class_simple_for_plot_cytotrace <- factor(df1$class_simple_for_plot_cytotrace,levels=rev(df2$class_simple_for_plot_cytotrace))

library(ggpubr)
compare_means(data=df1,formula=CytoTRACE~class_simple_for_plot_cytotrace)

df1$statuse <- as.character(df1$class_simple_for_plot_cytotrace)

pdf('result/CytoTRACE.cells.cancer.pdf',height=6,width=6)
ggplot(df1[df1$class_simple_for_plot_cytotrace %ni% '',],aes(x=class_simple_for_plot_cytotrace,y=CytoTRACE,fill=class_simple_for_plot_cytotrace)) + geom_violin(scale='width') + geom_boxplot(outlier.alpha = 0,width=0.2,fill='black',color='white',size=0.5) + theme_classic() + ylim(c(0,1)) + ggsci::scale_fill_jco() + ggpubr::stat_compare_means(ref.group='Ca_TPCS',hide.ns = F,label = 'p.signif',method = 'wilcox.test') + theme(legend.position = 'none') + theme(legend.position = 'none',axis.text.x = element_text(angle=90,size=20,hjust=1,vjust=0.5),text=element_text(size=20),axis.title.x =element_blank())
dev.off() 