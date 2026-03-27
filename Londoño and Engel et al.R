library(dplyr)
library(RColorBrewer)
library("org.Hs.eg.db", character.only = TRUE)
library(AnnotationDbi)
library(pheatmap)
library(decoupleR)
library(edgeR)
library(DESeq2)
library(stringr)
library(tidyr)
library(tibble)
library(limma)

setwd("~/Desktop/Research/DCreg/Hepatology_P7")

counts_path <- "Barcelona_data.csv"
labels_path <- "Barcelona_labels.csv"


counts <- read.csv(counts_path, header = TRUE, check.names = FALSE)
if (!is.numeric(counts[[1]])) {
  rownames(counts) <- counts[[1]]
  counts <- counts[,-1, drop = FALSE]
}
counts <- as.data.frame(counts)

labels <- read.csv(labels_path, header = TRUE, check.names = FALSE)
names(labels)[1:2] <- c("Sample", "Condition")
labels <- labels %>%
  mutate(
    Sample = as.character(Sample),
    Condition = as.character(Condition),
    Condition = str_squish(Condition)       
  )

labels$Condition <- case_when(
  str_to_lower(labels$Condition) %in% c("cabmr", "cabmr") ~ "cABMR",
  str_to_lower(labels$Condition) %in% c("clintcmr", "clintcmr") ~ "clinTCMR",
  str_to_lower(labels$Condition) %in% c("nhr") ~ "NHR",
  TRUE ~ labels$Condition
)


keep <- intersect(labels$Sample, colnames(counts))
if (length(keep) == 0) stop("No overlapping IDs between expression and labels!")
counts <- counts[, keep, drop = FALSE]
labels <- labels[match(keep, labels$Sample), ]
stopifnot(all(labels$Sample == colnames(counts)))

mx <- suppressWarnings(max(as.matrix(counts), na.rm = TRUE))
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = labels,
                              design = ~ Condition)


dds <- estimateSizeFactors(dds)
expr_norm <- counts(dds, normalized = TRUE)
expr_log2 <- log2(expr_norm + 1)

net <- get_progeny(organism = "human", top = 500)
acts <- run_mlm(mat = expr_log2, net = net, .source = 'source', .target = 'target', .mor = 'weight', minsize = 5)
acts <- acts %>% 
  left_join(labels, by = c("condition" = "Sample"))


prot <- acts %>% group_by(Condition, source) %>% summarise(group_score = mean(score, na.rm = TRUE), .groups = 'drop')
mat <- prot %>% pivot_wider(names_from = source, values_from = group_score) %>% tibble::column_to_rownames("Condition") %>% as.matrix()
mat_scaled <- scale(mat)

pheatmap(mat_scaled,
         main = "Pathway activity: Barcelona",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA)

##############################MMH cohort#######################
counts_path <- "count_table_Hannover_validation.csv"
labels_path <- "MHH_validation.csv"


counts <- read.csv(counts_path, header = TRUE, check.names = FALSE)
if (!is.numeric(counts[[1]])) {
  rownames(counts) <- counts[[1]]
  counts <- counts[,-1, drop = FALSE]
}
counts <- as.data.frame(counts)

labels <- read.csv(labels_path, header = TRUE, check.names = FALSE)
names(labels)[1:2] <- c("Sample", "Condition")
labels <- labels %>%
  mutate(
    Sample = as.character(Sample),
    Condition = as.character(Condition),
    Condition = str_squish(Condition)        # trim internal/trailing spaces
  )

keep <- intersect(labels$Sample, colnames(counts))
if (length(keep) == 0) stop("No overlapping sample IDs between expression and labels!")
counts <- counts[, keep, drop = FALSE]
labels <- labels[match(keep, labels$Sample), ]
stopifnot(all(labels$Sample == colnames(counts)))

dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = labels,
                              design = ~ Condition)


dds <- estimateSizeFactors(dds)
expr_norm <- counts(dds, normalized = TRUE)
expr_log2 <- log2(expr_norm + 1)

net <- get_progeny(organism = "human", top = 500)
acts <- run_mlm(mat = expr_log2, net = net, .source = 'source', .target = 'target', .mor = 'weight', minsize = 5)
acts <- acts %>% 
  left_join(labels, by = c("condition" = "Sample"))

prot <- acts %>% group_by(Condition, source) %>% summarise(group_score = mean(score, na.rm = TRUE), .groups = 'drop')
mat <- prot %>% pivot_wider(names_from = source, values_from = group_score) %>% tibble::column_to_rownames("Condition") %>% as.matrix()
mat_scaled <- scale(mat)

pheatmap(mat_scaled,
         main = "Pathway activity: MMH",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA)

#############################Sanchez feuyo - Londoño et al - data##################################

rawcounts <- read.csv(file.choose(), row.names = 1, check.names = FALSE) |> as.matrix()
sampleinfo <- read.csv(file.choose(), row.names = 1, header = TRUE)
sampleinfo$Condition <- factor(sampleinfo$Condition, levels = c("Non-SCR", "SCR", "Other", "Rejection"))
sampleinfo$Sample <- rownames(sampleinfo)

stopifnot(all(sampleinfo$Sample %in% colnames(rawcounts)))
rawcounts <- rawcounts[, sampleinfo$Sample, drop = FALSE]

design <- model.matrix(~ 0 + sampleinfo$Condition)
colnames(design) <- levels(sampleinfo$Condition)

dge <- DGEList(counts = rawcounts)
keep <- filterByExpr(dge, design = design)   # uses group structure & lib sizes
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
contrast.matrix <- makeContrasts(
  SCR_vs_NonSCR       = SCR - `Non-SCR`,
  Other_vs_NonSCR     = Other - `Non-SCR`,
  Rejection_vs_NonSCR = Rejection - `Non-SCR`,
  levels = design
)


kt <- "REFSEQ"
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = rownames(v$E),
  columns = "SYMBOL",
  keytype = kt # "REFSEQ" in your case
)

id_map <- id_map[!is.na(id_map$SYMBOL) & nzchar(id_map$SYMBOL), ]
id_map <- id_map[!duplicated(id_map[[kt]]), ]  # one row per REFSEQ

expr_logcpm <- v$E[rownames(v$E) %in% id_map[[kt]], , drop = FALSE]

mapped_symbols <- id_map$SYMBOL[match(rownames(expr_logcpm), id_map[[kt]])]
keep <- !is.na(mapped_symbols) & nzchar(mapped_symbols)
expr_logcpm <- expr_logcpm[keep, , drop = FALSE]
mapped_symbols <- mapped_symbols[keep]

rownames(expr_logcpm) <- mapped_symbols

expr_sum <- rowsum(expr_logcpm, group = rownames(expr_logcpm), na.rm = TRUE)
u <- rownames(expr_sum)
grp_sizes <- as.numeric(table(rownames(expr_logcpm))[u])  # align sizes to row order
expr_logcpm <- expr_sum / grp_sizes  # average per symbol


net <- get_progeny(organism = "human", top = 500)  # pathways x targets with weights

acts <- run_mlm(
  mat    = expr_logcpm,
  net    = net,
  .source = "source",
  .target = "target",
  .mor    = "weight",
  minsize = 5
)

sample_acts <- acts |>
  rename(Sample = condition) |>
  left_join(sampleinfo, by = "Sample")

group_acts <- sample_acts |>
  group_by(Condition, source) |>
  summarize(group_score = mean(score, na.rm = TRUE), .groups = "drop")

group_acts_mat <- group_acts |>
  pivot_wider(names_from = source, values_from = group_score) |>
  column_to_rownames("Condition") |>
  as.matrix()

group_acts_scaled <- scale(group_acts_mat)

pheatmap(group_acts_scaled,
         main = "Global Pathway Activity per Group (PROGENy)",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         border_color = NA)


