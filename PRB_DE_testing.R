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


# Split into permissive, non-permissive groups
yearstatus <- NormalizeData(yearstatus, assay = "RNA")
permissive <- subset(yearstatus, subCSDCDset = status == "Permissive")
nonpermissive <- subset(yearstatus, subset = status == "Non-permissive")
permissive0 <- subset(yearstatus, subset = status == "Permissive_0")
nonpermissive0 <- subset(yearstatus, subset = status == "Non-permissive_0")
PversusNP <- subset(yearstatus, status %in% c("Permissive", "Non-permissive"))
PversusNP0 <- subset(yearstatus, status %in% c("Permissive_0", "Non-permissive_0"))

#########################################Differential gene expression and sub-setting at 0 year #####################################
#Subseting CD4 T cells at 0
CD40  <- subset(PversusNP0, subset = majority_voting == "CD4 T cells")
FeatureDimPlot(CD40, features = "CD4", reduction = "UMAP")
CD40 <- subset(x = CD40, subset = CD4 > 0.5)
CD40 <- SCTransform(CD40, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD40 <- RunPCA(CD40, verbose = FALSE)
CD40 <- RunUMAP(CD40, dims = 1:50, verbose = FALSE)
CD40 <- FindNeighbors(CD40, dims = 1:50, verbose = FALSE)
CD40 <- FindClusters(CD40, verbose = FALSE)
CD40markers <- FindAllMarkers(object = CD40, min.pct = 0.25, logfc.threshold = 0.25)

CD40markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val < 0.0001) %>%
  slice_head(n = 7) %>%
  ungroup() -> top7

# final heatmap 
get_unique_genes <- function(data, gene_column) {
  genes <- data[[gene_column]]
  unique_genes <- unique(genes)
  return(unique_genes)
}

unique_genes_list <- get_unique_genes(top7, "gene")

DotPlot(CD40, features = unique_genes_list, col.min = 0, col.max = 5, cols = c("white", "#8b0000"), dot.scale = 3) + 
  rotate() + 
  RotatedAxis() + 
  theme(text = element_text(size = 9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))

#celltypist labeling
raw_counts_CD40 <- GetAssayData(object = CD40, assay = "RNA", layer = "counts")
dense_counts <- as.data.frame(as.matrix(raw_counts_CD40))
write.csv(dense_counts, file = "raw_counts_CD40.csv", row.names = TRUE)

predicted_labels <- read.csv("predicted_labels.csv", row.names = 1)
CD40 <- AddMetaData(CD40, metadata = predicted_labels$predicted_label, col.name = "predicted_cell_type_sub")
CD40 <- AddMetaData(CD40, metadata = predicted_labels$majority_voting, col.name = "majority_voting_sub")

#view results
CellDimPlot(CD40, group.by = c("majority_voting_sub","seurat_clusters"), reduction = "UMAP")
CellStatPlot(CD40, stat.by = "status", group.by = "majority_voting_sub", stat_type = "count", position = "dodge", label = TRUE, palcolor = c("#D2042D", "#0000FF"))
CellDimPlot(CD40, group.by = "majority_voting_sub", split.by = "status", reduction = "UMAP", cells.highlight = TRUE) 

#relabeling CD4 T cells
Idents(object = CD40) <- "majority_voting_sub"
CD40$majority_voting_sub <- paste0(CD40$majority_voting_sub, "_", CD40$seurat_clusters)

unique_titles <- unique(CD40$majority_voting_sub)
print(unique_titles)

target_cell_types <- c("Tcm/Naive helper T cells_1", "Trm helper T cells_1") 
new_label <- "Tem/Temra helper T cells_4"
CD40$majority_voting_sub [CD40$majority_voting_sub %in% target_cell_types] <- new_label
CD40$majority_voting_sub <- sub("_.*", "", CD40$majority_voting_sub)

#####DE expression
CD40$celltype.stim <- paste(CD40$majority_voting, CD40$status, sep = "_")
DefaultAssay(object = CD40) <- "SCT"
Idents(CD40) <- "celltype.stim"
unique_titles <- unique(CD40$celltype.stim)
print(unique_titles)

CD40_T_cells_DE <- FindMarkers(CD40, ident.1 = "CD4 T cells_Permissive_0", ident.2 = "CD4 T cells_Non-permissive_0", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD40_T_cells_DE$Significant <- with(CD40_T_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD40_T_cells_filtered <- subset(CD40_T_cells_DE, Significant == "Yes")
Idents(CD40) <- "majority_voting"
genes1 <- rownames(CD40_T_cells_filtered)

GroupHeatmap(CD40, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 35, 
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD40 <- RunDEtest(srt = CD40, group_by = "status", fc.threshold = 1)
CD40 <- RunEnrichment(
  srt = CD40, group_by = "status", db = c("GO_BP","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD40, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"))

EnrichmentPlot(srt = CD40, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("M5950", "GO:0006119", "GO:1902600", "GO:0009060", "GO:0045333", "GO:0019646", "GO:00427731", "GO:0002181", "GO:1900118", 
                          "GO:1900117", "GO:0097194", "GO:0042274", "GO:0000028"))

#Subsetting CD8 T cells
CD80 <- subset(PversusNP0, subset = majority_voting == "CD8 T cells")
FeatureDimPlot(CD80, features = "CD8A", reduction = "UMAP")
CD80 <- subset(x = CD80, subset = CD8A > 0.5)
CD80 <- SCTransform(CD80, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD80 <- RunPCA(CD80, verbose = FALSE)
CD80 <- RunUMAP(CD80, dims = 1:50, verbose = FALSE)
CD80 <- FindNeighbors(CD80, dims = 1:50, verbose = FALSE)
CD80 <- FindClusters(CD80, verbose = FALSE)
CD80markers <- FindAllMarkers(object = CD80, min.pct = 0.25, logfc.threshold = 0.25)

CD80markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val < 0.00001) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

# final heatmap 
get_unique_genes <- function(data, gene_column) {
  genes <- data[[gene_column]]
  unique_genes <- unique(genes)
  return(unique_genes)
}

unique_genes_list <- get_unique_genes(top10, "gene")

DotPlot(CD80, features = unique_genes_list, col.min = 0, col.max = 5, cols = c("white", "#8b0000"), dot.scale = 3) + 
  rotate() + 
  RotatedAxis() + 
  theme(text = element_text(size = 9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))

raw_counts_CD80 <- GetAssayData(object = CD80, assay = "RNA", layer = "counts")
dense_counts <- as.data.frame(as.matrix(raw_counts_CD80))
write.csv(dense_counts, file = "raw_counts_CD80.csv", row.names = TRUE)

predicted_labels <- read.csv("predicted_labels.csv", row.names = 1)
CD80 <- AddMetaData(CD80, metadata = predicted_labels$predicted_label, col.name = "predicted_cell_type_sub")
CD80 <- AddMetaData(CD80, metadata = predicted_labels$majority_voting, col.name = "majority_voting_sub")

#view results from celltypist
CellDimPlot(CD80, group.by = c("majority_voting_sub","seurat_clusters"), reduction = "UMAP")
CellStatPlot(CD80, stat.by = "status", group.by = "majority_voting_sub", stat_type = "count", position = "dodge", label = TRUE, palcolor = c("#0000FF", "#D2042D"))
CellDimPlot(CD80, group.by = "majority_voting_sub", split.by = "status", reduction = "UMAP", cells.highlight = TRUE) 

# Individual scattered cell relabeling
target_cell_types <- c("Tem/Temra cytotoxic T cells-2", "Tem/Trm cytotoxic T cells-2")   
new_label <- "Cytotoxic CD8 T cells"
CD80$majority_voting_sub [CD80$majority_voting_sub %in% target_cell_types] <- new_label
CD80$majority_voting_sub <- sub("-.*", "", CD80$majority_voting_sub)
Idents(CD80) <- "majority_voting_sub"

#####DE expression
CD80$celltype.stim <- paste(CD80$majority_voting, CD80$status, sep = "_")
DefaultAssay(object = CD80) <- "SCT"
Idents(CD80) <- "celltype.stim"
unique_titles <- unique(CD80$celltype.stim)
print(unique_titles)
CD80_cells_DE <- FindMarkers(CD80, ident.1 = "CD8 T cells_Permissive_0", ident.2 = "CD8 T cells_Non-permissive_0", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD80_cells_DE$Significant <- with(CD80_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD80_cells_filtered <- subset(CD80_cells_DE, Significant == "Yes")
Idents(CD80) <- "majority_voting"
genes1 <- rownames(CD80_cells_filtered)
CD80@meta.data$status <- factor(CD80@meta.data$status, levels = c('Permissive_0', 'Non-permissive_0'))

genes1 <- c("S100A8", "GZMH", "CSK", "GZMB", "ISG15", "LTB", "ADGRE5", "KLF13", 
            "IL2RG", "CTSW", "IKZF1", "HLA-DPB1", "CD3G", "GZMA", "IL32", "CCL5", 
            "CXCR4", "CCL3L1", "TMEM2", "TMEM66", "CD69", "LY6E")

GroupHeatmap(CD80, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 35, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD80 <- RunDEtest(srt = CD80, group_by = "status", fc.threshold = 1)

CD80 <- RunEnrichment(
  srt = CD80, group_by = "status", db = c("GO_BP", "MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD80, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("M5913", "M5891", "M5911", "M5950", "M5921", "M5926", "M5936", "M5911", "M5950", 
                          "M5890", "M5950", "GO:0006754", "GO:0015986", "GO:0009206", "GO:0009145", "GO:0009201", "GO:0006119", "GO:0019646",
                          "GO:0042773", "GO:0042775", "GO:0009060"))

#subsetting CD16 positive NK cells
CD16positiveNK0 <- subset(PversusNP0, subset = majority_voting == "CD16+ NK cells")
FeatureDimPlot(CD16positiveNK0, features = "FCGR3A", reduction = "UMAP")
CD16positiveNK0 <- subset(x = CD16positiveNK0, subset = FCGR3A > 0.5)
CD16positiveNK0 <- SCTransform(CD16positiveNK0, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD16positiveNK0 <- RunPCA(CD16positiveNK0, verbose = FALSE)
CD16positiveNK0 <- RunUMAP(CD16positiveNK0, dims = 1:50, verbose = FALSE)
CD16positiveNK0 <- FindNeighbors(CD16positiveNK0, dims = 1:50, verbose = FALSE)
CD16positiveNK0 <- FindClusters(CD16positiveNK0, verbose = FALSE)
CD16positiveNK0markers <- FindAllMarkers(object = CD16positiveNK0, min.pct = 0.25, logfc.threshold = 0.25)

#####DE expression
CD16positiveNK0$celltype.stim <- paste(CD16positiveNK0$majority_voting, CD16positiveNK0$status, sep = "_")
DefaultAssay(object = CD16positiveNK0) <- "SCT"
Idents(CD16positiveNK0) <- "celltype.stim"
unique_titles <- unique(CD16positiveNK0$celltype.stim)
print(unique_titles)
CD16positiveNK0_cells_DE <- FindMarkers(CD16positiveNK0, ident.1 = "CD16+ NK cells_Permissive_0", ident.2 = "CD16+ NK cells_Non-permissive_0", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD16positiveNK0_cells_DE$Significant <- with(CD16positiveNK0_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD16positiveNK0_cells_filtered <- subset(CD16positiveNK0_cells_DE, Significant == "Yes")
Idents(CD16positiveNK0) <- "majority_voting"
genes1 <- row.names(CD16positiveNK0_cells_filtered)
CD16positiveNK0@meta.data$status <- factor(CD16positiveNK0@meta.data$status, levels = c('Permissive_0', 'Non-permissive_0'))

genes1 <- c("S100B", "TAP2.1", "LSMD1", "IRAK3", "S100A8", "IL10RB", "IL18", "KIF13B", 
            "KARS", "GSK3A", "THEMIS2", "IL2RG", "KLRD1", "FCER1G", "GNLY", 
            "ADGRG1", "CCL4", "CCL3L1", "CCL4L2", "CCL3", "KIR2DL2", "SH2D2A", "XCL2", "SOCS1", "AREG")

GroupHeatmap(CD16positiveNK0, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 30, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD16positiveNK0 <- RunDEtest(srt = CD16positiveNK0, group_by = "status", fc.threshold = 1)
CD16positiveNK0 <- RunEnrichment(
  srt = CD16positiveNK0, group_by = "status", db = c("GO_BP", "GO_CC", "GO_MF","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD16positiveNK0, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("M5906", "M5923", "M5950", "M5890", "M5936", "M5902", "M5913", "M5950", 
                          "M5926", "M5939", "GO:0060213", "GO:0042776", "GO:0006754", "GO:0042744", "GO:0002181", "GO:0046034", "GO:0006754",
                          "GO:0009060", "GO:0006119", "GO:0015986"))

###############subsetting CD16 negative NK cells
CD16negNK0 <- subset(PversusNP0, subset = majority_voting == "CD16- NK cells")
FeatureDimPlot(CD16negNK0, features = "FCGR3A", reduction = "UMAP")
CD16negNK0 <- subset(x = CD16negNK0, subset = FCGR3A < 4.0)
CD16negNK0 <- SCTransform(CD16negNK0, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD16negNK0 <- RunPCA(CD16negNK0, verbose = FALSE)
CD16negNK0 <- RunUMAP(CD16negNK0, dims = 1:50, verbose = FALSE)
CD16negNK0 <- FindNeighbors(CD16negNK0, dims = 1:50, verbose = FALSE)
CD16negNK0 <- FindClusters(CD16negNK0, verbose = FALSE)
CD16negNK0markers <- FindAllMarkers(object = CD16negNK0, min.pct = 0.25, logfc.threshold = 0.25)

#####DE expression
FeatureDimPlot(yearstatus, features = "FCGR3A", reduction = "UMAP")
CD16negNK0$celltype.stim <- paste(CD16negNK0$majority_voting, CD16negNK0$status, sep = "_")
DefaultAssay(object = CD16negNK0) <- "SCT"
Idents(CD16negNK0) <- "celltype.stim"
unique_titles <- unique(CD16negNK0$celltype.stim)
print(unique_titles)
CD16negNK0_cells_DE <- FindMarkers(CD16negNK0, ident.1 = "CD16- NK cells_Permissive_0", ident.2 = "CD16- NK cells_Non-permissive_0", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD16negNK0_cells_DE$Significant <- with(CD16negNK0_cells_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD16negNK0_cells_filtered <- subset(CD16negNK0_cells_DE, Significant == "Yes")
Idents(CD16negNK0) <- "majority_voting"
genes1 <- row.names(CD16negNK0_cells_filtered)
CD16negNK0@meta.data$status <- factor(CD16negNK0@meta.data$status, levels = c('Permissive_0', 'Non-permissive_0'))

genes1 <- c("LAG3", "KDM5D", "IL7R", "FRG1B", "CSF3R", "TNFSF8", "S100A8", "IL18R1", 
            "CD27", "IST1", "FEZ1", "GZMK", "IL2RB", "TIGIT", "CXCR6", "CD160", 
            "PDCD4", "CD63", "IFITM2", "TSC22D3", "KLRF1", "ID2", "S100A4", "IFITM1", "CCL3", "CD226", 
            "IFNG","CCL3L1", "GZMB", "CCL4L2","CX3CR1", "SOCS1", "KLRD1")

GroupHeatmap(CD16negNK0, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 40, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD16negNK0 <- RunDEtest(srt = CD16negNK0, group_by = "status", fc.threshold = 1)
CD16negNK0 <- RunEnrichment(
  srt = CD16negNK0, group_by = "status", db = c("GO_BP", "GO_CC", "GO_MF","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD16negNK0, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("M5932", "M5921", "M5930", "M5950", "M5890", "M5916", "M5913", "M5936", 
                          "M5950", "M5926", "GO:0022618", "GO:0042274", "GO:0042254", "GO:0022613", "GO:0002181", "GO:1990868", "GO:0072677",
                          "GO:0007159", "GO:0002696", "GO:0050867"))

###############subsetting CD14+ MONOCYTES
CD14monocytes0 <- subset(PversusNP0, subset = majority_voting == "CD14+ Monocytes")
FeatureDimPlot(CD14monocytes0, features = "CD14", reduction = "UMAP")
CD14monocytes0 <- SCTransform(CD14monocytes0, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
CD14monocytes0 <- RunPCA(CD14monocytes0, verbose = FALSE)
CD14monocytes0 <- RunUMAP(CD14monocytes0, dims = 1:50, verbose = FALSE)
CD14monocytes0 <- FindNeighbors(CD14monocytes0, dims = 1:50, verbose = FALSE)
CD14monocytes0 <- FindClusters(CD14monocytes0, verbose = FALSE)
CD14monocytes0markers <- FindAllMarkers(object = CD14monocytes0, min.pct = 0.25, logfc.threshold = 0.25)

#####DE expression
CD14monocytes0$celltype.stim <- paste(CD14monocytes0$majority_voting, CD14monocytes0$status, sep = "_")
DefaultAssay(object = CD14monocytes0) <- "SCT"
Idents(CD14monocytes0) <- "celltype.stim"
unique_titles <- unique(CD14monocytes0$celltype.stim)
print(unique_titles)
CD14monocytes0_DE <- FindMarkers(CD14monocytes0, ident.1 = "CD14+ Monocytes_Permissive_0", ident.2 = "CD14+ Monocytes_Non-permissive_0", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
CD14monocytes0_DE$Significant <- with(CD14monocytes0_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
CD14monocytes0_filtered <- subset(CD14monocytes0_DE, Significant == "Yes")
Idents(CD14monocytes0) <- "majority_voting"
genes1 <- row.names(CD14monocytes0_filtered)
genes1 <- c("CXCR2", "CXCR1", "MMP9", "CD47", "CD177", "IL18R1", "IL18RAP", "IL1R2", 
            "S100P", "IFITM2", "CD55", "SLC25A37", "CLEC4D", "TLR4", "SLC2A3", "SLC16A3", 
            "S100A8", "S100A12", "S100A9", "LGALS1", "S100A10", "ID2", "TMEM230", "CD74", "CD36", "HLA-DRB1", 
            "HLA-DRA","AHR", "HLA-DRB5", "SOD1","CCL3", "HLA-DMA", "CLEC4A", "LYZ", "FOS")

CD14monocytes0@meta.data$status <- factor(CD14monocytes0@meta.data$status, levels = c('Permissive_0', 'Non-permissive_0'))

GroupHeatmap(CD14monocytes0, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 40, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF", "#D2042D")), cluster_rows = TRUE)

CD14monocytes0 <- RunDEtest(srt = CD14monocytes0, group_by = "status", fc.threshold = 1)
VolcanoPlot(srt = CD14monocytes0, group_by = "status", nlabel = 20)
CD14monocytes0 <- RunEnrichment(
  srt = CD14monocytes0, group_by = "status", db = c("GO_BP", "MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = CD14monocytes0, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("GO:0002181", "GO:0022613", "GO:0042254", "GO:0022618", "GO:0071826", "GO:0072593", "GO:0042776", "GO:0000302", 
                          "GO:0015986", "GO:2000377", "M5890", "M5932", "M5921", "M5947", "M5913", "M5926", "M5936", "M5950", "M5939", "M5913"))

###############subsetting neutrophils
Neutrophils0 <- subset(PversusNP0, subset = majority_voting == "Neutrophils")
Neutrophils0 <- SCTransform(Neutrophils0, vars.to.regress = c("percent.mt", "percent.ribo", "percent.hb"), verbose = TRUE)
Neutrophils0 <- RunPCA(Neutrophils0, verbose = FALSE)
Neutrophils0 <- RunUMAP(Neutrophils0, dims = 1:50, verbose = FALSE)
Neutrophils0 <- FindNeighbors(Neutrophils0, dims = 1:50, verbose = FALSE)
Neutrophils0 <- FindClusters(Neutrophils0, verbose = FALSE)
Neutrophils0markers <- FindAllMarkers(object = Neutrophils0, min.pct = 0.25, logfc.threshold = 0.25)

#####DE expression
Neutrophils0$celltype.stim <- paste(Neutrophils0$majority_voting, Neutrophils0$status, sep = "_")
DefaultAssay(object = Neutrophils0) <- "SCT"
Idents(Neutrophils0) <- "celltype.stim"
unique_titles <- unique(Neutrophils0$celltype.stim)
print(unique_titles)
Neutrophils0_DE <- FindMarkers(Neutrophils0, ident.1 = "Neutrophils_Permissive_0", ident.2 = "Neutrophils_Non-permissive_0", verbose = TRUE, logfc.threshold = 0.25, assay = "SCT", recorrect_umi = FALSE)
Neutrophils0_DE$Significant <- with(Neutrophils0_DE, ifelse(p_val_adj < 0.05, "Yes", "No"))
Neutrophils0_filtered <- subset(Neutrophils0_DE, Significant == "Yes")
Idents(Neutrophils0) <- "majority_voting"
genes1 <- row.names(Neutrophils0_filtered)
genes1 <- c("IL8", "C5AR1", "SLC12A6", "IL2RG", "SLC11A1", "TREM1", "TNFRSF1B", "S100P", 
            "FOSL2", "CXCR4", "MMP25", "HLA-A", "HLA-B", "CD48", "CST3", "LGALS1", 
            "LYZ", "SLC25A6", "HLA-DRA", "CD74", "HLA-DRB1", "S100A10", "CD300E", "CD36", "LILRB1", "VCAN", "FOS", "CSF3R", "IL18R1", "THBS1", "STAT3", "CD55")

Neutrophils0@meta.data$status <- factor(Neutrophils0@meta.data$status, levels = c('Permissive_0', 'Non-permissive_0'))

GroupHeatmap(Neutrophils0, features = genes1,  group.by = "status", exp_method = "log2fc", 
             slot = "counts", nlabel = 40, features_fontsize = 4,
             heatmap_palette = "material-blue", 
             group_palcolor = list(c("#0000FF","#D2042D")), cluster_rows = TRUE)

Neutrophils0 <- RunDEtest(srt = Neutrophils0, group_by = "status", fc.threshold = 1)
VolcanoPlot(srt = Neutrophils0, group_by = "status", nlabel = 20)
Neutrophils0 <- RunEnrichment(
  srt = Neutrophils0, group_by = "status", db = c("GO_BP", "GO_CC", "GO_MF","MSigDB"), species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05")

EnrichmentPlot(srt = Neutrophils0, group_by = "status", plot_type = "comparison", db = c("GO_BP", "MSigDB"), 
               id_use = c("M5930", "M5946", "M5950", "M5939", "M5926", "M5902", "M5921", "M5913", "M5932", 
                          "M5890", "GO:0030099", "GO:0008154", "GO:0007265", "GO:0007264", "GO:0007015", "GO:0042255", "GO:0042274", "GO:0022613", "GO:0042254", "GO:0002181"))

