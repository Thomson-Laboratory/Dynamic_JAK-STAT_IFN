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


setwd("~/Desktop/Research/DCreg/Hepatology_P7/GSE285476_RAW")

folders <- c(
  "GSM8703037_liver_S1D1",
  "GSM8703038_PBMC_S1D1",
  "GSM8703039_liver_S2D1",
  "GSM8703040_PBMC_S2D1",
  "GSM8703041_liver_S1D7",
  "GSM8703042_PBMC_S1D7",
  "GSM8703043_liver_S2D7",
  "GSM8703044_PBMC_S2D7",
  "GSM8703045_liver_S3D7",
  "GSM8703046_PBMC_S3D7",
  "GSM8703047_PBMC_HC",
  "GSM8703048_liver_R1D1",
  "GSM8703049_PBMC_R1D1",
  "GSM8703050_liver_R2D1",
  "GSM8703051_PBMC_R2D1",
  "GSM8703052_liver_R1D7",
  "GSM8703053_PBMC_R1D7",
  "GSM8703054_liver_R2D7",
  "GSM8703055_PBMC_R2D7",
  "GSM8703056_liver_R3D7",
  "GSM8703057_PBMC_R3D7",
  "GSM8703058_liver_R1D14",
  "GSM8703059_PBMC_R1D14",
  "GSM8703060_liver_HC"
)

objs <- list()

for (folder in folders) {
  path <- file.path(getwd(), folder)
  
  mat <- Read10X(data.dir = path)
  
  so <- CreateSeuratObject(
    counts       = mat,
    project      = folder,
    min.cells    = 3,
    min.features = 200
  )
  
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
  so[["percent.hb"]] <- PercentageFeatureSet(so, pattern = "^Hba|^Hbb")
  
  so <- subset(so, subset = nFeature_RNA >= 500 & percent.mt <= 25)
  
  # metadata
  so$sample_id  <- folder
  so$orig.ident <- folder
  
  so <- SCTransform(
    so,
    vars.to.regress = c("percent.mt"),
    verbose = FALSE
  )
  
  objs[[folder]] <- so
}

saveRDS(objs, file = "murine_objs.rds")
objs <- read_rds("murine_objs.rds")


objs <- lapply(objs, function(x) {
  if ("RNA" %in% names(x@assays)) {
    x[["RNA"]] <- as(object = x[["RNA"]], Class = "Assay")
  }
  x
})

objs <- lapply(objs, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x
})


for (i in seq_along(objs)) {
  objs[[i]] <- RenameCells(
    objs[[i]],
    add.cell.id = paste0("sample", i)
  )
}

features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)


combined <- RunFastMNN(
  object.list = objs,
  features = features,
  # you can tune the dimensionality if you want, e.g. d = 50
  verbose = FALSE
)

combined <- FindNeighbors(
  combined,
  reduction = "mnn",
  dims = 1:30,
  verbose = FALSE
)

combined <- FindClusters(
  combined,
  resolution = 0.5,
  algorithm = 1,
  verbose = FALSE
)

combined <- RunUMAP(
  combined,
  reduction = "mnn",
  dims = 1:30,
  verbose = FALSE
)

saveRDS(combined, "combined_murine_complete.rds")
combined <- readRDS("combined_murine_complete.rds")


##labeling##
Idents(object = combined) <- "self_anno"
combined$status <- paste0(combined$orig.ident, "x", combined$orig.ident)

# Get the unique titles from the 'self_annotation' column
unique_titles <- unique(combined$status)

# Print unique titles
print(unique_titles)


# Individual scattered cell relabeling
target_cell_types <-  c("GSM8703037_liver_S1D1xGSM8703037_liver_S1D1",
                        "GSM8703038_PBMC_S1D1xGSM8703038_PBMC_S1D1",  
                        "GSM8703039_liver_S2D1xGSM8703039_liver_S2D1",
                        "GSM8703040_PBMC_S2D1xGSM8703040_PBMC_S2D1",  
                         "GSM8703041_liver_S1D7xGSM8703041_liver_S1D7",
                        "GSM8703042_PBMC_S1D7xGSM8703042_PBMC_S1D7",  
                        "GSM8703043_liver_S2D7xGSM8703043_liver_S2D7",
                       "GSM8703044_PBMC_S2D7xGSM8703044_PBMC_S2D7",  
                       "GSM8703045_liver_S3D7xGSM8703045_liver_S3D7",
                       "GSM8703046_PBMC_S3D7xGSM8703046_PBMC_S3D7"  
                     
                  )
new_label <- "Syngeneic" 
combined$status[combined$status %in% target_cell_types] <- new_label

#################################Mouse ISG list##################################
ISG_list <- c( "Procr","Elf1",
               "Ifi27","Trim26","Csf1","Samd9l","Usp18","Psmb9","Psmb8","Uba7","Cxcl10",
               "Eif2ak2","C1s1","Isg15","Irf7","Cxcl11","Nub1","Txnip","Adar","Ripk2","Ifitm3",
               "Gmpr","Lap3","Herc6","Lpar6","Ube2l6","Rtp4","Epsti1","Tmem140","Ifitm1","Bst2",
               "Ifi35","Ifih1","Pnpt1","Ogfr","Parp14","Ccrl2","Mvb12a","Batf2","Trim14",
               "Sp110","Trafd1","Gbp3","Nmi","Parp9",
               "Ifi30","Tdrd7","Parp12","Oas1a","Ddx60","Lamp3",
               "Tent5a","Trim12c","Cnp","Cd47","Mov10","Sell",
               "Cd86","Trim25","Il15","Stat3","Stat2","Stat4","Stat1","Fgl2","Cdkn1a","Xcl1",
               "Wars1","Il15ra","Pml","Nfkbia","Psma3","Psma2","Tnfaip2","Il4ra","Cfb","St8sia4",
               "Ly6e","Trim21","Hif1a","Tnfsf10","Fpr1","Cd38","Irf9","Casp4","Casp3","Gpr18",
               "Myd88","Ripk1","Ciita","Gzma","Casp7","Cmklr1","Psme2","Psme1","Psmb10","Irf4",
               "Upp1","Bpgm","Ifnar2","Ifit3","Pla2g4a","Tnfaip6","Tnfaip3","Tapbp","Socs3",
               "Casp8","Ncoa3","Rnf213","Lcp2","Il18bp","Vamp8","Mthfd2",
               "St3gal5","Nod1","Rbck1","Psmb2","Irf5","Cxcl9",
               "Sspn","Tor1b","Lats2","Socs1","Vamp5",
               "Pfkp","Isoc1","Sppl2a","Eif4e3",
               "Peli1", "Lysmd2","Tmt1b",
               "Nup93","Apol6","Auts2","Marchf1","Cmtr1",
               "Slamf7","Mvp","Cd274","Zbp1","Samhd1",
               "Isg20","Rsad2","Nampt","Dhx58","Ifitm2","Rnf31","Ifi30","Znfx1","Tdrd7",
               "Parp12","P2ry14","Arid5b","Slc25a28","Oasl1","Oas3","Oas2","Ddx60","Rapgef6",
               "Sectm1a","Helz2","Bank1","Rigi","Ifi44","Nlrc5","Xaf1","B2m","Btg1","Cd40",
               "Cd69","Cfh","Plscr1","Serping1","Fas","Fcgr1","Gch1","H2-Aa","H2-DMa","Ifi44l",
               "Ptpn6","Icam1","Irf8","Ido1","Cd74","Il10ra","Casp1","Il2rb","Il6","Il7","Irf1",
               "Irf2","Itgb7","Jak2","Mx2","Nfkb1","Pnp","Pim1","Ptgs2","Ptpn1","Ptpn2","Ccl5",
               "Selp","Sod2","Sri","Tap1","Vcam1","Arl4a","Ifit2","Ccl7","Lgals3bp","Pde4b",
               "Cmpk2")


filter_ISG_list <- function(combined, ISG_list) {
  assay_data_genes <- rownames(GetAssayData(combined, assay = "SCT", slot = "data"))
  return(ISG_list[ISG_list %in% assay_data_genes])
}

ISG_list <- filter_ISG_list(combined, ISG_list)
ISG_scores <- Matrix::colSums(GetAssayData(combined, assay = "SCT", slot = "data")[ISG_list, ])
combined <- AddMetaData(combined, metadata = ISG_scores, col.name = "ISG_scores")


p1 <- FeatureStatPlot(combined, stat.by = "ISG_scores", group.by = "status", multiplegroup_comparisons = TRUE)

dittoDimPlot(combined, "ISG_scores",
             max.color = "red", min.color = "gray90", split.by = "status")

DimPlot(merged, group.by = "status")

###################################pathway activity inference analysis#############
Idents(combined) <- "status"
sample <- combined
assay <- "SCT"

network <- decoupleR::get_progeny(organism = "mouse")
activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 20,
                                   minsize = 5)

out <- SCpubr::do_PathwayActivityHeatmap(sample = sample,
                                         activities = activities)
out

pathway_heatmap


p1 <- CellDimPlot(combined, group.by = "annotation", split.by = "status")
p2 <- FeatureDimPlot(combined, features = "ISG_scores", split.by = "status")

###########################Marker annotation##################
Idents(object = combined) <- "annotation"

# Find all markers after self annotation
markers_self_annotation <- FindAllMarkers(object = combined, min.pct = 0.25, logfc.threshold = 0.25)

#top 10 
top3_genes_per_cluster <- markers_self_annotation %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_max(n = 3, order_by = avg_log2FC)

markers_self_annotation %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

# final heatmap 
get_unique_genes <- function(data, gene_column) {
  genes <- data[[gene_column]]
  unique_genes <- unique(genes)
  return(unique_genes)
}

unique_genes_list <- get_unique_genes(top3, "gene")

DotPlot(combined, features = unique_genes_list, col.min = 0, col.max = 5, cols = c("white", "#8b0000"), dot.scale = 3) + 
  rotate() + 
  RotatedAxis() + 
  theme(text = element_text(size = 9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))
