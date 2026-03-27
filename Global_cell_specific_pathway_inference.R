library(dplyr)
library(Seurat)
library(SeuratObject)
library(sctransform)
library(ggplot2)
library(ggpubr)
library(NMF)
library(remotes)
library(SeuratDisk)
library(Matrix)
library(openxlsx)
library(Nebulosa)
library(HGNChelper)
library(scCustomize)
library(AnnotationDbi)
library(CellChat)
library(scRNAtoolVis)
library(cowplot)
library(pheatmap)
library(dittoSeq)
library(genekitr)
library(enrichplot)
library(SCP)
library(decoupleR)
library(readr)


#############################################global and cell specific pathway inference######################################################
#global
Idents(yearstatus) <- "status"
sample <- yearstatus
assay <- "SCT"

network <- decoupleR::get_progeny(organism = "human")
activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)

out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                      activities = activities)
out

#cell-specific
Idents(yearstatus) <- "majority_voting"
yearstatusTNMK <- subset(yearstatus, majority_voting %in% c("CD8 T cells", "CD4 T cells", "Neutrophils", "CD16+ Monocytes", "CD14+ Monocytes", "CD16+ NK cells", "CD16- NK cells"))
sample <- yearstatusTNMK
assay <- "SCT"
network <- decoupleR::get_progeny(organism = "human")
activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)

out <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                      activities = activities,
                                      split.by = "status")
out
