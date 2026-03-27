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

setwd("~/Desktop/Research/DCreg/Hepatology_P7/GSE256141_RAW")

folders <- c(
  "Patient_1_Pre_TXP",
  "Patient_2_Pre_TXP",
  "Patient_3_Pre_TXP",
  "Patient_4_Pre_TXP",
  "Patient_4_BX2",
  "Patient_5_Pre_TXP",
  "Patient_5_BX1",
  "Patient_6_Pre_TXP",
  "Patient_6_BX2",
  "Patient_6_BX3",
  "Patient_6_BX4",
  "Patient_6_BX5",
  "Patient_7_BX2",
  "Patient_7_BX3",
  "Patient_7_BX4",
  "Patient_8_BX2",
  "Patient_9_BX1",
  "Patient_10_BX1",
  "Patient_11_BX1",
  "Patient_11_BX2",
  "Patient_11_BX3",
  "Patient_12_BX1",
  "Patient_12_BX2",
  "Patient_12_BX3",
  "Patient_12_BX4",
  "Patient_12_BX5",
  "Patient_13_BX1",
  "Patient_13_BX2",
  "Patient_14_BX1",
  "Patient_14_BX4"
)

objs <- list()

for (folder in folders) {
  
  path <- file.path(getwd(), folder)
  
  # Load data
  M.data <- Read10X(data.dir = path)
  
  # Create Seurat object
  M <- CreateSeuratObject(counts = M.data,
                          project = folder,
                          min.cells = 3,
                          min.features = 200)
  
  # % mitochondrial genes
  M[["percent.mt"]] <- PercentageFeatureSet(M, pattern = "^MT-")
  M[["percent.hb"]] <- PercentageFeatureSet(M, pattern = "^HB")
  
  # Subset cells based on thresholds
  M <- subset(M, nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 50)
  
  # SCTransform normalization 
  M <- SCTransform(M, vst.flavor = "v2", verbose = FALSE)
  
  # Store in list for integration
  objs[[folder]] <- M
}


#SCT
resolution <- 1.0
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
objs <- PrepSCTIntegration(object.list = objs, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = objs,
                                  normalization.method = "SCT",
                                  anchor.features = features)


combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:30, verbose = TRUE, n.neighbors = 20)
combined <- RunTSNE(combined, reduction = "pca", dims = 1:30, check_duplicates = FALSE)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30, verbose = TRUE, k.param = 20)
combined <- FindClusters(combined, resolution = resolution)
#####

# Renaming identities 
Idents(object = combined) <- "majority_voting"
combined$status <- paste0(combined$orig.ident, "x", combined$orig.ident)

# Get the unique titles from the 'self_annotation' column
unique_titles <- unique(combined$status)

# Print unique titles
print(unique_titles)

# Individual scattered cell relabeling
target_cell_types <- c("Patient_12_BX3xPatient_12_BX3", "Patient_12_BX4xPatient_12_BX4",
                       "Patient_12_BX5xPatient_12_BX5", "Patient_13_BX2xPatient_13_BX2", "Patient_14_BX4xPatient_14_BX4")  

new_label <- "Resolved_TCMR" 
combined$status[combined$status %in% target_cell_types] <- new_label

combined$status <- sub("x.*", "", combined$status)


combined[["RNA"]] <- JoinLayers(combined[["RNA"]])

#Labeling with celltypes; Access assay for celltypist
raw_counts_combined <- GetAssayData(object = combined, assay = "RNA", layer = "counts")
dense_counts <- as.data.frame(as.matrix(raw_counts_combined))
write.csv(dense_counts, file = "raw_counts_combined.csv", row.names = TRUE)


#Added labels
predicted_labels <- read.csv("predicted_labels.csv", row.names = 1)
combined <- AddMetaData(combined, metadata = predicted_labels$predicted_label, col.name = "predicted_cell_type")
combined <- AddMetaData(combined, metadata = predicted_labels$majority_voting, col.name = "majority_voting")
Idents(object = combined) <- "majority_voting"

CellDimPlot(srt = combined, group.by = "seurat_clusters",
            reduction = "UMAP", label = TRUE)

CellDimPlot(srt = combined, group.by = c("majority_voting", "status"),
            reduction = "UMAP", theme_use = "theme_blank", label = TRUE) + guides(color = guide_legend(ncol = 1, override.aes = list(size=4))) + coord_fixed() 

ISG_list <- c("ADAR", "B2M", "BATF2", "BST2", "C1S", "CASP1", "CASP8", "CCRL2", "CD47", "CD74", 
              "CMPK2", "CNP", "CSF1", "CXCL10", "CXCL11", "DDX60", "DHX58", "EIF2AK2", "ELF1", "EPSTI1", 
              "MVB12A", "TENT5A", "CMTR1", "GBP2", "GBP4", "GMPR", "HERC6", "HLA-C", "IFI27", "IFI30", 
              "IFI35", "IFI44", "IFI44L", "IFIH1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "IL15", 
              "IL4R", "IL7", "IRF1", "IRF2", "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP", 
              "LPAR6", "LY6E", "MOV10", "MX1", "NCOA7", "NMI", "NUB1", "OAS1", "OASL", "OGFR", "PARP12", "PARP14", 
              "PARP9", "PLSCR1", "PNPT1", "HELZ2", "PROCR", "PSMA3", "PSMB8", "PSMB9", "PSME1", "PSME2", 
              "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L", "SELL", "SLC25A28", "SP110", "STAT2", 
              "TAP1", "TDRD7", "TMEM140", "TRAFD1", "TRIM14", "TRIM21", "TRIM25", "TRIM26", "TRIM5", "TXNIP", 
              "UBA7", "UBE2L6", "USP18", "WARS1", "APOL6", "ARID5B", "ARL4A", "AUTS2", "BANK1", "BPGM", "BTG1", 
              "C1R", "CASP3", "CASP4", "CASP7", "CCL2", "CCL5", "CCL7", "CD274", "CD38", "CD40", "CD69", "CD86", 
              "CDKN1A", "CFB", "CFH", "CIITA", "CMKLR1", "CSF2RB", "CXCL9", "DDX58", "EIF4E3", "FAS", "FCGR1A", 
              "FGL2", "FPR1", "GBP6", "GCH1", "GPR18", "GZMA", "HIF1A", "HLA-A", "HLA-B", "HLA-DMA", "HLA-DQA1", 
              "HLA-DRB1", "HLA-G", "ICAM1", "IDO1", "IFIT1", "IFNAR2", "IL10RA", "IL15RA", "IL18BP", "IL2RB", "IL6", 
              "IRF4", "IRF5", "IRF8", "ISOC1", "ITGB7", "JAK2", "KLRK1", "LATS2", "LCP2", "LYSMD2", "MARCH1", 
              "METTL7B", "MT2A", "MTHFD2", "MVP", "MX2", "MYD88", "NAMPT", "NCOA3", "NFKB1", "NFKBIA", 
              "NLRC5", "NOD1", "NUP93", "OAS2", "OAS3", "P2RY14", "PDE4B", "PELI1", "PFKP", "PIM1", 
              "PLA2G4A", "PML", "PNP", "PSMA2", "PSMB10", "PSMB2", "PTGS2", "PTPN1", "PTPN2", "PTPN6", 
              "RAPGEF6", "RBCK1", "RIPK1", "RNF213", "SAMHD1", "SECTM1", "SELP", "SERPING1", 
              "SLAMF7", "SOCS1", "SOCS3", "SOD2", "SPPL2A", "SRI", "SSPN", "ST3GAL5", "ST8SIA4", "STAT1", "STAT3", 
              "STAT4", "TAPBP", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TOR1B", "UPP1", "VAMP5", "VAMP8", 
              "VCAM1", "XAF1", "XCL1", "ZBP1", "ZNFX1", "WARS")


filter_ISG_list <- function(combined, ISG_list) {
  assay_data_genes <- rownames(GetAssayData(combined, assay = "SCT", slot = "data"))
  return(ISG_list[ISG_list %in% assay_data_genes])
}

ISG_list <- filter_ISG_list(combined, ISG_list)
ISG_scores <- Matrix::colSums(GetAssayData(combined, assay = "SCT", slot = "data")[ISG_list, ])
combined <- AddMetaData(combined, metadata = ISG_scores, col.name = "ISG_scores")
p1 <- FeatureStatPlot(combined, stat.by = "ISG_scores", split.by = "status", multiplegroup_comparisons = TRUE)
FeatureDimPlot(combined, features = "ISG_scores", split.by = "status", assay = "SCT")
CellDimPlot(combined, group.by = "majority_voting", split.by = "status")

#####
Idents(combined) <- "status"
sample <- combined
assay <- "SCT"

network <- decoupleR::get_progeny(organism = "human")
activities <- decoupleR::run_wmean(mat = as.matrix(sample@assays[[assay]]@data),
                                   network = network,
                                   .source = "source",
                                   .targe = "target",
                                   .mor = "weight",
                                   times = 100,
                                   minsize = 5)

out5 <- SCpubr::do_PathwayActivityPlot(sample = sample,
                                       activities = activities)

out5


