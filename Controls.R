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
library(readr)

# Function to process each seurat object
process_sample <- function(sample_id) {
  path <- sprintf("%s.raw_feature_bc_matrix.h5", sample_id)
  obj <- CreateSeuratObject(Read10X_h5(path), min.features = 200, project = sample_id)
  obj <- RenameCells(obj, add.cell.id = sample_id)
  
  # Calculate percentages
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  
  # Gene annotations
  mito_genes <- grep("^MT-", rownames(obj), value = TRUE)
  ribo_genes <- grep("^RP[SL]", rownames(obj), value = TRUE)
  hemo_genes <- grep("^HB", rownames(obj), value = TRUE)
  
  # Apply QC filters
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 35 & percent.hb < 0.5)
  
  # SCTransform
  obj <- SCTransform(obj, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE, variable.features.n = NULL)
  new_var_genes <- setdiff(VariableFeatures(obj), c(mito_genes, ribo_genes, hemo_genes, "XIST"))
  VariableFeatures(obj) <- new_var_genes
  
  return(obj)
}

# Sample IDs
sample_ids <- c("001_12M", "002_12M", "004_12M", "007_12M", "001_PRB", "002_PRB", "004_PRB", "007_PRB" )
samples_list <- lapply(sample_ids, process_sample)

# Assign original sample ID and status
statuses <- c("Non-permissive", "Non-permissive", "Permissive", "Permissive", "Non-permissive", "Non-permissive", "Permissive", "Permissive")
names(samples_list) <- sample_ids
for (i in seq_along(samples_list)) {
  samples_list[[i]]$status <- statuses[i]
}


################################Analysis#########################################
# Integration features selection
features <- SelectIntegrationFeatures(object.list = samples_list, nfeatures = 3000)

# Preparing for SCT integration
yearstatus.list <- PrepSCTIntegration(object.list = samples_list, anchor.features = features)

# Finding integration anchors
immune.anchors <- FindIntegrationAnchors(object.list = yearstatus.list, normalization.method = "SCT", anchor.features = features)

# Integrating data
controls <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

# Principal Component Analysis
controls <- RunPCA(controls, features = VariableFeatures(object = controls), npcs = 50, verbose = TRUE)

# Elbow Plot to determine the number of PCs to use
ElbowPlot(object = controls, ndims = 50)

# UMAP for visualization
controls <- RunUMAP(object = controls, dims = 1:30, verbose = TRUE, n.neighbors = 20)

# Finding neighbors
controls <- FindNeighbors(object = controls, dims = 1:30, verbose = TRUE, k.param = 20)

# Clustering
controls <- FindClusters(object = controls, verbose = TRUE, resolution = 1.0)


######################################################cell typist###################################################### 
controls[["RNA"]] <- JoinLayers(controls[["RNA"]])

#Labeling with celltypes; Access assay for celltypist
raw_counts_yearstatus <- GetAssayData(object = controls, assay = "RNA", layer = "counts")
dense_counts <- as.data.frame(as.matrix(raw_counts_yearstatus))
write.csv(dense_counts, file = "raw_counts_controls.csv", row.names = TRUE)

#Added labels
predicted_labels <- read.csv("predicted_labels.csv", row.names = 1)
controls <- AddMetaData(controls, metadata = predicted_labels$predicted_label, col.name = "predicted_cell_type")
controls <- AddMetaData(controls, metadata = predicted_labels$majority_voting, col.name = "majority_voting")
Idents(object = controls) <- "majority_voting"

#################################

controls$id <- paste0(controls$status, "x", controls$orig.ident)
unique_titles <- unique(controls$id)
print(unique_titles)

# Renaming identities 
Idents(object = controls) <- "seurat_clusters"
controls$majority_voting <- paste0(controls$majority_voting, "x", controls$seurat_clusters)

# Get the unique titles from the 'self_annotation' column
unique_titles <- unique(controls$majority_voting)

# Print unique titles
print(unique_titles)

# Individual scattered cell relabeling
target_cell_types <- c(      
  
  
  
  "Neutrophilsx11"            
  
  
  
)


new_label <- "CD4 T cells"                                                                    

controls$majority_voting[controls$majority_voting %in% target_cell_types] <- new_label
controls$majority_voting <- sub("x.*", "", controls$majority_voting)


CellDimPlot(controls, group.by = "majority_voting", theme_use = "theme_blank")
CellDimPlot(controls, group.by = "seurat_clusters", theme_use = "theme_blank")
CellDimPlot(controls, group.by = "predicted.celltype.l2")
CellDimPlot(controls, group.by = "predicted_cell_type", theme_use = "theme_blank")


Idents(object = controls) <- "seurat_clusters"
Subsetting <- subset(x = controls, idents = c("8", "27", "17", "25"))
Subsetting <- subset(x = Subsetting, subset = CD8A > 0.1)

Subsetting$majority_voting <- paste0(Subsetting$majority_voting, "x", Subsetting$predicted.celltype.l2)
unique_titles <- unique(Subsetting$majority_voting)

# Print unique titles
print(unique_titles)
target_cell_types <- c(  "CD4 T cellsxCD4 T" ,               "CD8 T cellsxCD8 T"  ,             
                         "NeutrophilsxNK" ,                  "NeutrophilsxMonocyte"  ,           "NeutrophilsxHepatocyte"       
)


new_label <- "CD8 T cells"
Subsetting$majority_voting[Subsetting$majority_voting %in% target_cell_types] <- new_label
Subsetting$majority_voting <- sub("x.*", "", Subsetting$majority_voting)
controls@meta.data[Cells(Subsetting), "majority_voting"] <- Subsetting$majority_voting


FeatureDimPlot(Subsetting, features = "CD8A", assay = "SCT", split.by = "status")
CellDimPlot(Subsetting, group.by = "majority_voting", label = TRUE)


#Save
saveRDS(controls, file = "controls.rds")
controls <- read_rds("controls.rds")
unique_titles <- unique(controls$id)

# Print unique titles
print(unique_titles)


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


MHC_list <- c("HLA-A", "HLA-B", "HLA-C", "HLA-DMA", 
              "HLA-DMB", "HLA-DPA1", "HLA-DPB1", 
              "HLA-DQA1", "HLA-DQB1", "HLA-DQB3", "HLA-DRA", 
              "HLA-DRB1", "HLA-DRB4", 
              "HLA-DRB5")


filter_ISG_list <- function(controls, MHC_list) {
  assay_data_genes <- rownames(GetAssayData(controls, assay = "SCT", slot = "data"))
  return(MHC_list[MHC_list %in% assay_data_genes])
}

MHC_list <- filter_ISG_list(controls, MHC_list)
MHC_scores <- Matrix::colSums(GetAssayData(controls, assay = "SCT", slot = "data")[MHC_list, ])
controls <- AddMetaData(controls, metadata = MHC_scores, col.name = "MHC_scores")
FeatureDimPlot(controls, features = "MHC_scores", split.by = "status")
FeatureStatPlot(controls, stat.by = "ISG_scores", split.by = "status", comparisons = TRUE, add_box = TRUE)
FeatureDimPlot(controls, features = "STAT1", split.by = "status")

########################

Idents(controls) <- "status"
sample <- controls
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


saveRDS(controls, file = "controls.rds")

##############heatmaps##########


MMSD <- c("IFNG", "IL1B", "IL2", "IL4", "CXCL8", "IL10",
          "IL6", "TNF", "CCL3", "CCL17", "CXCL10", "CCL2", "CSF2", 
          "IL1A", "IL7", "IL15", "LTA", "VEGFA", "FGF2", "TEK", 
          "FLT1", "CRP", "ICAM1", "SAA1", "VCAM1", "IFNB1", "IFNA1", "CXCL9", "FOS", 
          "STAT1", "ISG15", "S100A12")

controls@meta.data$status <- factor(controls@meta.data$status, levels = c('Permissive0', 'Non-permissive0','Permissive', "Non-permissive"))


GroupHeatmap(controls,
             features = MMSD, group.by = c("status","orig.ident"), 
             assay = "SCT", nlabel = 46, cluster_rows = TRUE, 
             cluster_columns = TRUE, 
             height = 5, cell_annotation_params = list(width = unit(5, "mm")))


p1 <- FeatureStatPlot(controls, stat.by = "ISG_scores", split.by = "status", comparisons = TRUE, add_box = TRUE)


write_csv(p1[["data"]], "MHC.csv")
