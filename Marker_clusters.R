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



######################################markers for clusters#################################
Idents(object = yearstatus) <- "seurat_clusters"

# Find all markers after self annotation
markers_self_annotation <- FindAllMarkers(object = yearstatus, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_self_annotation, file = "markers_self_annotation.csv")

#top 10 
top3_genes_per_cluster <- markers_self_annotation %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_max(n = 3, order_by = avg_log2FC)

write.csv(top3_genes_per_cluster, file = "top3_genes_per_cluster.csv")

# Basic heatmap of clusters without annotation and by cluster
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

DotPlot(yearstatus, features = unique_genes_list, col.min = 0, col.max = 5, cols = c("white", "#8b0000"), dot.scale = 3) + 
  rotate() + 
  RotatedAxis() + 
  theme(text = element_text(size = 9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))


# Find all markers after self annotation
markers_self_annotation <- FindAllMarkers(object = yearstatus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_self_annotation, file = "markers_self_annotation.csv")

#top 10 
top3_genes_per_cluster <- markers_self_annotation %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_max(n = 3, order_by = avg_log2FC)

write.csv(top3_genes_per_cluster, file = "top3_genes_per_self_annotation.csv")

# Basic heatmap of clusters without annotation and by cluster
markers_self_annotation %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  dplyr::filter(p_val < 0.00001) %>%
  slice_head(n = 3) %>%
  ungroup() -> top3

# final heatmap 
get_unique_genes <- function(data, gene_column) {
  genes <- data[[gene_column]]
  unique_genes <- unique(genes)
  return(unique_genes)
}

unique_genes_list <- get_unique_genes(top3, "gene")

DotPlot(yearstatus, features = unique_genes_list, col.min = 0, col.max = 5, cols = c("white", "#8b0000"), dot.scale = 3) + 
  rotate() + 
  RotatedAxis() + 
  theme(text = element_text(size = 9), axis.text.x = element_text(size = 9), axis.text.y = element_text(size = 9))


# Clustered_DotPplot of labeled cells
markersp7 <- Add_Pct_Diff(markersp7,
                          pct.1_name = "pct.1",
                          pct.2_name = "pct.2",
                          overwrite = FALSE)
