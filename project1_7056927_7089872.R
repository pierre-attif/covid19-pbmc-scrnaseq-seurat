################################################################################
# Single Cell Bioinformatics 2026-27 - Project 1: scRNA-seq analysis
# Weeks 1-4 (Tasks 1-11)
#
# Environment: course conda environment (R 4.0.5, Seurat 4.0.1,
# DoubletFinder 2.0.3, igraph, SingleR, celldex, enrichR; for Week 4 also
# monocle3 1.0.0, SeuratWrappers and CellChat 1.6.1). Run this script from the folder that contains
# the four *.counts.matrices.rds files. Figures are written to ./figures and
# tables to ./tables.
################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(42)

# when run with Rscript, do not write the printed plots to Rplots.pdf
# (all figures are saved as PNG files below)
if (!interactive()) pdf(NULL, width = 20, height = 20)

dir.create("figures", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

sample_names <- c("COVID556", "COVID557", "COVID558", "HIP043")

# helper: save a ggplot/patchwork figure into ./figures
save_fig <- function(plot, name, width = 8, height = 5) {
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = width, height = height, dpi = 150)
}

################################################################################
############################### WEEK 1 #########################################
################################################################################

#### Task 1: load the data #####################################################
# Every .rds file holds a list with three gene x cell UMI count matrices:
# exon (spliced mRNA), intron (unspliced pre-mRNA) and spanning (reads over an
# exon-intron boundary). We use the exon matrix for the downstream analysis.
raw <- list(
  COVID556 = readRDS("GSM4557329_GSM4557329_556_cell.counts.matrices.rds"),
  COVID557 = readRDS("GSM4557330_GSM4557330_557_cell.counts.matrices.rds"),
  COVID558 = readRDS("GSM4557331_GSM4557331_558_cell.counts.matrices.rds"),
  HIP043   = readRDS("GSM4557337_GSM4557337_HIP043_cell.counts.matrices.rds")
)

# the three matrices of a sample contain the same cells (same barcodes)
sapply(raw, function(s) identical(colnames(s$exon), colnames(s$intron)) &&
                        identical(colnames(s$exon), colnames(s$spanning)))

# one Seurat object per sample from the exon counts
# min.cells = 3:     keep genes detected in at least 3 cells
# min.features = 200: keep barcodes with at least 200 detected genes
seu_list <- lapply(sample_names, function(n) {
  CreateSeuratObject(counts = raw[[n]]$exon, project = n,
                     min.cells = 3, min.features = 200)
})
names(seu_list) <- sample_names
seu_list

#### Task 2: add metadata ######################################################
meta <- data.frame(
  sample    = sample_names,
  donor     = c("C2", "C3", "C4", "H4"),
  sex       = c("M", "M", "M", "M"),
  condition = c("COVID", "COVID", "COVID", "Healthy")
)

for (n in sample_names) {
  row <- meta[meta$sample == n, ]
  seu_list[[n]]$sample    <- row$sample
  seu_list[[n]]$donor     <- row$donor
  seu_list[[n]]$sex       <- row$sex
  seu_list[[n]]$condition <- row$condition
}
head(seu_list$COVID556@meta.data)

#### Task 3: raw-data information ##############################################
# 3.1 UMIs per matrix, cells and genes per sample (raw matrices, before any
#     filtering by CreateSeuratObject)
raw_table <- data.frame(
  sample         = sample_names,
  UMIs_exon      = sapply(raw, function(s) sum(s$exon)),
  UMIs_intron    = sapply(raw, function(s) sum(s$intron)),
  UMIs_spanning  = sapply(raw, function(s) sum(s$spanning)),
  cells          = sapply(raw, function(s) ncol(s$exon)),
  genes_exon     = sapply(raw, function(s) nrow(s$exon)),
  genes_intron   = sapply(raw, function(s) nrow(s$intron)),
  genes_spanning = sapply(raw, function(s) nrow(s$spanning)),
  row.names = NULL
)
raw_table
write.csv(raw_table, "tables/3_raw_table.csv", row.names = FALSE)

# cells and genes that remain in the Seurat objects (after min.cells/min.features)
data.frame(sample = sample_names,
           cells  = sapply(seu_list, ncol),
           genes  = sapply(seu_list, nrow))

# metadata columns of the objects
colnames(seu_list$COVID556@meta.data)

# 3.2 the two automatic metadata fields are column-wise summaries of the counts:
#     nCount_RNA = colSums(counts), nFeature_RNA = colSums(counts > 0)
cnt <- GetAssayData(seu_list$COVID556, slot = "counts")
all(Matrix::colSums(cnt)     == seu_list$COVID556$nCount_RNA)
all(Matrix::colSums(cnt > 0) == seu_list$COVID556$nFeature_RNA)
rm(cnt)

# 3.3 histogram of UMIs per cell for each sample (log10 x axis)
meta_all <- do.call(rbind, lapply(seu_list, function(s) s@meta.data))
meta_all$sample <- factor(meta_all$sample, levels = sample_names)

p_umi_w1 <- ggplot(meta_all, aes(x = nCount_RNA)) +
  geom_histogram(bins = 60) +
  scale_x_log10() +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = "UMIs per cell (log10 scale)", y = "Number of cells") +
  theme_bw()
p_umi_w1
save_fig(p_umi_w1, "3_umis_per_cell")

#### Task 4: save the objects ##################################################
saveRDS(seu_list, "week1_seurat_objects.rds")


################################################################################
############################### WEEK 2 #########################################
################################################################################

seu_list <- readRDS("week1_seurat_objects.rds")

#### Task 5.1: quality control and filtering ###################################
# 5.1.1 percentage of mitochondrial UMIs per cell (genes starting with "MT-")
seu_list <- lapply(seu_list, function(s) {
  s$percent.mt <- PercentageFeatureSet(s, pattern = "^MT-")
  s
})

meta_all <- do.call(rbind, lapply(seu_list, function(s) s@meta.data))
meta_all$sample <- factor(meta_all$sample, levels = sample_names)

# 5.1.2 distributions of genes, UMIs and mitochondrial % per cell
p_genes <- ggplot(meta_all, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 60) + scale_x_log10() +
  facet_wrap(~ sample, scales = "free_y") +
  labs(title = "Detected genes per cell", x = "Genes per cell (log10 scale)",
       y = "Number of cells") + theme_bw()
p_umis <- ggplot(meta_all, aes(x = nCount_RNA)) +
  geom_histogram(bins = 60) + scale_x_log10() +
  facet_wrap(~ sample, scales = "free_y") +
  labs(title = "UMIs per cell", x = "UMIs per cell (log10 scale)",
       y = "Number of cells") + theme_bw()
p_mt <- ggplot(meta_all, aes(x = percent.mt)) +
  geom_histogram(bins = 60) +
  facet_wrap(~ sample, scales = "free_y") +
  labs(title = "Mitochondrial UMIs per cell", x = "Mitochondrial UMIs (%)",
       y = "Number of cells") + theme_bw()
p_genes; p_umis; p_mt

# 5.1.3 cut-offs per sample, chosen from the histograms above and from the
#       Week 1 UMI histograms (see report)
cutoffs <- data.frame(
  sample    = sample_names,
  min_genes = c(200,  350,  200,  350),
  max_genes = c(1000, 1800, 1000, 1600),
  max_mt    = c(14,   14,   12,   14)
)
cutoffs

# how many cells each mitochondrial threshold would remove (used to compare a
# 10 % threshold with the chosen ones)
sapply(sample_names, function(n) {
  mt <- seu_list[[n]]$percent.mt
  c(removed_at_10 = round(100 * mean(mt >= 10), 1),
    removed_at_chosen = round(100 * mean(mt >= cutoffs$max_mt[cutoffs$sample == n]), 1))
})

# QC plots with the chosen cut-offs drawn in
cut_long <- data.frame(sample = factor(rep(sample_names, 2), levels = sample_names),
                       value  = c(cutoffs$min_genes, cutoffs$max_genes))
cut_mt   <- data.frame(sample = factor(sample_names, levels = sample_names),
                       value  = cutoffs$max_mt)
p_genes_cut <- p_genes + geom_vline(data = cut_long, aes(xintercept = value),
                                    colour = "red")
p_mt_cut    <- p_mt + geom_vline(data = cut_mt, aes(xintercept = value),
                                 colour = "red")
save_fig(p_genes_cut, "5_1_genes_per_cell")
save_fig(p_umis,      "5_1_umis_per_cell")
save_fig(p_mt_cut,    "5_1_mito_per_cell")

# genes vs UMIs coloured by mitochondrial %, to see which cells each cut-off hits
p_scatter <- ggplot(meta_all, aes(x = nCount_RNA, y = nFeature_RNA, colour = percent.mt)) +
  geom_point(size = 0.2) + scale_x_log10() + scale_y_log10() +
  scale_colour_viridis_c(limits = c(0, 30), oob = scales::squish) +
  facet_wrap(~ sample) +
  geom_hline(data = cut_long, aes(yintercept = value), colour = "red") +
  labs(x = "UMIs per cell", y = "Genes per cell", colour = "% mito") + theme_bw()
save_fig(p_scatter, "5_1_genes_vs_umis", width = 8, height = 6)

# keep a flag of which cells pass QC (needed for the bonus task)
for (n in sample_names) {
  cf <- cutoffs[cutoffs$sample == n, ]
  s  <- seu_list[[n]]
  seu_list[[n]]$qc_pass <- s$nFeature_RNA > cf$min_genes &
                           s$nFeature_RNA < cf$max_genes &
                           s$percent.mt   < cf$max_mt
}

seu_filtered <- lapply(seu_list, function(s) subset(s, subset = qc_pass))

# 5.1.4 cells before and after filtering
qc_table <- data.frame(
  sample = sample_names,
  before = sapply(seu_list, ncol),
  after  = sapply(seu_filtered, ncol),
  row.names = NULL
)
qc_table$removed         <- qc_table$before - qc_table$after
qc_table$percent_removed <- round(100 * qc_table$removed / qc_table$before, 1)
qc_table

#### Bonus: fraction of intronic UMIs per cell #################################
# intron UMIs / (exon + intron + spanning UMIs), compared between the cells
# kept and removed in Task 5.1
intron_frac_df <- do.call(rbind, lapply(sample_names, function(n) {
  cells <- colnames(seu_list[[n]])
  ex <- Matrix::colSums(raw[[n]]$exon[, cells])
  it <- Matrix::colSums(raw[[n]]$intron[, cells])
  sp <- Matrix::colSums(raw[[n]]$spanning[, cells])
  data.frame(sample = n, cell = cells,
             intron_frac = it / (ex + it + sp),
             percent.mt  = seu_list[[n]]$percent.mt,
             nFeature_RNA = seu_list[[n]]$nFeature_RNA,
             status = ifelse(seu_list[[n]]$qc_pass, "kept", "removed"))
}))
intron_frac <- setNames(intron_frac_df$intron_frac,
                        paste(intron_frac_df$sample, intron_frac_df$cell, sep = "_"))

# median intronic fraction of kept / removed cells per sample
bonus_table <- aggregate(intron_frac ~ sample + status, data = intron_frac_df,
                         FUN = median)
bonus_table$intron_frac <- round(bonus_table$intron_frac, 3)
bonus_table

# split the removed cells by the reason they were removed
intron_frac_df$reason <- "kept"
for (n in sample_names) {
  cf <- cutoffs[cutoffs$sample == n, ]
  i  <- intron_frac_df$sample == n & intron_frac_df$status == "removed"
  d  <- intron_frac_df[i, ]
  intron_frac_df$reason[i] <- ifelse(d$percent.mt >= cf$max_mt, "high mito",
                              ifelse(d$nFeature_RNA <= cf$min_genes, "few genes",
                                     "many genes"))
}
aggregate(intron_frac ~ sample + reason, data = intron_frac_df,
          FUN = function(x) round(median(x), 3))

p_bonus <- ggplot(intron_frac_df, aes(x = reason, y = intron_frac, fill = status)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ sample, nrow = 1) +
  labs(x = NULL, y = "Fraction of intronic UMIs") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_bonus, "5_bonus_intron_fraction", width = 10, height = 4)

p_bonus_scatter <- ggplot(intron_frac_df, aes(x = percent.mt, y = intron_frac, colour = status)) +
  geom_point(size = 0.2, alpha = 0.4) + facet_wrap(~ sample, nrow = 1) +
  labs(x = "Mitochondrial UMIs (%)", y = "Fraction of intronic UMIs") + theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
save_fig(p_bonus_scatter, "5_bonus_intron_vs_mito", width = 10, height = 3.5)

rm(raw); invisible(gc())

#### Task 5.2: normalization and feature selection #############################
# LogNormalize: counts / total counts of the cell * 10,000, then log1p
seu_filtered <- lapply(seu_filtered, function(s) {
  NormalizeData(s, normalization.method = "LogNormalize", scale.factor = 10000)
})

# 2000 most variable genes per sample (vst)
seu_filtered <- lapply(seu_filtered, function(s) {
  FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
})

# variable genes of COVID556 with the top 10 labelled
top10 <- head(VariableFeatures(seu_filtered$COVID556), 10)
top10
p_var <- LabelPoints(plot = VariableFeaturePlot(seu_filtered$COVID556),
                     points = top10, repel = TRUE) +
  ggtitle("COVID556: 2000 most variable genes")
p_var
save_fig(p_var, "5_2_variable_features_COVID556", width = 8, height = 5)

# top 10 of the other samples, for comparison
lapply(seu_filtered, function(s) head(VariableFeatures(s), 10))

saveRDS(seu_filtered, "week2_normalized.rds")

#### Task 5.3: doublet removal (DoubletFinder) #################################
library(DoubletFinder)

# DoubletFinder needs scaled data and a PCA
seu_filtered <- lapply(seu_filtered, function(s) {
  s <- ScaleData(s)
  RunPCA(s, verbose = FALSE)
})

# elbow plots to choose the PCs used for DoubletFinder
p_elbow_samples <- wrap_plots(lapply(sample_names, function(n)
  ElbowPlot(seu_filtered[[n]], ndims = 30) + ggtitle(n)), ncol = 2)
p_elbow_samples
save_fig(p_elbow_samples, "5_3_elbow_per_sample", width = 9, height = 6)

# PCs per sample from the elbow plots (the curves flatten after ~7-8 PCs for
# the COVID samples and ~5 PCs for HIP043; we take a few more to be safe)
n_pcs_df <- c(COVID556 = 10, COVID557 = 10, COVID558 = 10, HIP043 = 8)

# expected doublet rate for Seq-Well (Poisson loading, see report):
# lambda = 15,000 cells / 86,000 wells, P(>=2 cells | >=1 cell)
lambda       <- 15000 / 86000
doublet_rate <- (1 - exp(-lambda) * (1 + lambda)) / (1 - exp(-lambda))
round(doublet_rate, 3)

sweep_plots   <- list()
pk_best       <- c()
doublet_info  <- list()
singlets_list <- list()

for (n in sample_names) {
  s   <- seu_filtered[[n]]
  pcs <- 1:n_pcs_df[[n]]

  # clusters are needed to model the fraction of homotypic doublets
  s <- FindNeighbors(s, dims = pcs, verbose = FALSE)
  s <- FindClusters(s, resolution = 0.5, verbose = FALSE)

  # pK sweep: artificial doublets are simulated for every pN/pK combination;
  # the best pK has the highest mean-variance normalised bimodality (BCmetric)
  sweep <- paramSweep_v3(s, PCs = pcs, sct = FALSE)
  stats <- summarizeSweep(sweep, GT = FALSE)
  bcmvn <- find.pK(stats)
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  pk_best[n] <- bcmvn$pK[which.max(bcmvn$BCmetric)]

  sweep_plots[[n]] <- ggplot(bcmvn, aes(x = pK, y = BCmetric)) +
    geom_line() + geom_point() +
    geom_vline(xintercept = pk_best[n], colour = "red", linetype = "dashed") +
    labs(title = paste0(n, " (best pK = ", pk_best[n], ")"),
         x = "pK", y = "BCmetric") + theme_bw()

  # expected number of doublets, reduced by the homotypic fraction
  # (doublets of two cells of the same cluster cannot be detected)
  homotypic <- modelHomotypic(s$seurat_clusters)
  n_exp     <- round(doublet_rate * ncol(s))
  n_exp_adj <- round(n_exp * (1 - homotypic))

  s <- doubletFinder_v3(s, PCs = pcs, pN = 0.25, pK = pk_best[n],
                        nExp = n_exp_adj, reuse.pANN = FALSE, sct = FALSE)

  # the classification column name depends on the parameters
  df_col  <- grep("^DF.classifications", colnames(s@meta.data), value = TRUE)
  s$doublet_class <- s@meta.data[[df_col]]
  print(n)
  print(table(s$doublet_class))

  doublet_info[[n]] <- data.frame(sample = n, pK = pk_best[n],
                                  homotypic = round(homotypic, 3),
                                  n_exp = n_exp, n_exp_adj = n_exp_adj)

  singlets_list[[n]] <- subset(s, subset = doublet_class == "Singlet")
}

p_sweep <- wrap_plots(sweep_plots, ncol = 2)
p_sweep
save_fig(p_sweep, "5_3_pk_sweep", width = 9, height = 6)

doublet_info <- do.call(rbind, doublet_info)
doublet_info

# 5.3.2 extend the QC table by the removed doublets and the remaining cells
qc_table$doublets  <- sapply(sample_names, function(n) ncol(seu_filtered[[n]]) - ncol(singlets_list[[n]]))
qc_table$remaining <- sapply(singlets_list, ncol)
qc_table
write.csv(qc_table, "tables/5_qc_table.csv", row.names = FALSE)

saveRDS(singlets_list, "week2_singlets.rds")
rm(seu_filtered); invisible(gc())

#### Task 5.4: merging and batch correction ####################################
# (A) merge without batch correction; cell names get the sample as prefix
#     because the same barcode can occur in several samples
merged <- merge(singlets_list$COVID556,
                y = list(singlets_list$COVID557, singlets_list$COVID558,
                         singlets_list$HIP043),
                add.cell.ids = sample_names)
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:20)

# (B) Seurat integration (CCA + anchors) on the per-sample normalized objects
singlets_list <- lapply(sample_names, function(n)
  RenameCells(singlets_list[[n]], add.cell.id = n))
names(singlets_list) <- sample_names

features <- SelectIntegrationFeatures(object.list = singlets_list, nfeatures = 2000)
anchors  <- FindIntegrationAnchors(object.list = singlets_list,
                                   anchor.features = features, dims = 1:20)
integrated <- IntegrateData(anchorset = anchors, dims = 1:20)
rm(anchors, singlets_list); invisible(gc())

# PCA/UMAP on the batch-corrected "integrated" assay
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, verbose = FALSE)
integrated <- RunUMAP(integrated, dims = 1:20)

p_54 <- (DimPlot(merged, group.by = "sample") + ggtitle("A: merged, by sample") |
         DimPlot(merged, group.by = "condition") + ggtitle("A: merged, by condition")) /
        (DimPlot(integrated, group.by = "sample") + ggtitle("B: integrated, by sample") |
         DimPlot(integrated, group.by = "condition") + ggtitle("B: integrated, by condition"))
p_54
save_fig(p_54, "5_4_umap_A_vs_B", width = 11, height = 9)

saveRDS(merged, "week2_A_merged.rds")
saveRDS(integrated, "week2_B_integrated.rds")

#### Task 6.1: dimensionality reduction ########################################
p_elbow <- ElbowPlot(merged, ndims = 40) + ggtitle("A: merged") |
           ElbowPlot(integrated, ndims = 40) + ggtitle("B: integrated")
p_elbow
save_fig(p_elbow, "6_1_elbow_A_B", width = 10, height = 4)

# variance explained by each PC (standard deviation squared, relative to the
# variance of the first 50 PCs), to support the choice from the elbow plot
var_expl <- function(obj) {
  sd <- obj[["pca"]]@stdev
  round(100 * sd^2 / sum(sd^2), 2)
}
data.frame(PC = 1:25, A = var_expl(merged)[1:25], B = var_expl(integrated)[1:25])

n_pcs <- 15  # both elbow plots flatten between PC 10 and PC 15

merged     <- RunUMAP(merged, dims = 1:n_pcs)
integrated <- RunUMAP(integrated, dims = 1:n_pcs)

p_61 <- DimPlot(merged, group.by = "sample") + ggtitle(paste0("A: merged (", n_pcs, " PCs)")) |
        DimPlot(integrated, group.by = "sample") + ggtitle(paste0("B: integrated (", n_pcs, " PCs)"))
p_61
save_fig(p_61, "6_1_umap_A_B", width = 11, height = 4.5)

#### Task 6.2: clustering ######################################################
# shared nearest neighbour graph on the selected PCs (not on the UMAP)
merged     <- FindNeighbors(merged, dims = 1:n_pcs)
integrated <- FindNeighbors(integrated, dims = 1:n_pcs)

# 6.2.1 Louvain (Seurat default) for several resolutions; number of clusters
res_l <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
merged     <- FindClusters(merged, resolution = res_l, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = res_l, verbose = FALSE)

n_clusters <- data.frame(
  resolution = res_l,
  A = sapply(res_l, function(r) length(unique(merged@meta.data[[paste0("RNA_snn_res.", r)]]))),
  B = sapply(res_l, function(r) length(unique(integrated@meta.data[[paste0("integrated_snn_res.", r)]])))
)
n_clusters

res_A <- 0.4
res_B <- 0.4
merged$louvain     <- merged@meta.data[[paste0("RNA_snn_res.", res_A)]]
integrated$louvain <- integrated@meta.data[[paste0("integrated_snn_res.", res_B)]]
Idents(merged)     <- "louvain"
Idents(integrated) <- "louvain"

p_621 <- DimPlot(merged, group.by = "louvain", label = TRUE) +
           ggtitle(paste0("A: Louvain, resolution ", res_A)) |
         DimPlot(integrated, group.by = "louvain", label = TRUE) +
           ggtitle(paste0("B: Louvain, resolution ", res_B))
p_621
save_fig(p_621, "6_2_louvain_A_B", width = 11, height = 4.5)

# 6.2.2 Leiden on the same SNN graph of data set A (igraph, no Python needed)
snn <- as(merged@graphs$RNA_snn, "dgCMatrix")
g   <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
leiden <- igraph::cluster_leiden(g, objective_function = "modularity",
                                 weights = igraph::E(g)$weight,
                                 resolution_parameter = res_A, n_iterations = 10)
merged$leiden <- factor(setNames(igraph::membership(leiden),
                                 igraph::V(g)$name)[colnames(merged)])
length(unique(merged$leiden))

louvain_leiden <- table(Louvain = merged$louvain, Leiden = merged$leiden)
louvain_leiden
write.csv(as.data.frame.matrix(louvain_leiden), "tables/6_2_louvain_vs_leiden.csv")
# share of each Louvain cluster that goes into its best-matching Leiden cluster
round(apply(prop.table(louvain_leiden, 1), 1, max), 2)
# adjusted Rand index between the two partitions
igraph::compare(as.integer(merged$louvain), as.integer(merged$leiden), method = "adjusted.rand")

p_leiden <- DimPlot(merged, group.by = "louvain", label = TRUE) + ggtitle("A: Louvain") |
            DimPlot(merged, group.by = "leiden", label = TRUE) + ggtitle("A: Leiden")
save_fig(p_leiden, "6_2_louvain_vs_leiden", width = 11, height = 4.5)

# 6.2.3 cells per Louvain cluster and sample
tab_A <- table(Cluster = merged$louvain, Sample = merged$sample)
tab_B <- table(Cluster = integrated$louvain, Sample = integrated$sample)
tab_A; tab_B
round(100 * prop.table(tab_A, 1), 1)
round(100 * prop.table(tab_B, 1), 1)
# overall share of each sample, as a reference for "dominated"
round(100 * prop.table(table(merged$sample)), 1)
write.csv(as.data.frame.matrix(tab_A), "tables/6_2_cluster_sample_A.csv")
write.csv(as.data.frame.matrix(tab_B), "tables/6_2_cluster_sample_B.csv")

# 6.2.4 a cluster of A dominated by one sample that is also separate in B
dom_A <- apply(prop.table(tab_A, 1), 1, max)
cluster_info <- data.frame(cluster = rownames(tab_A), cells = rowSums(tab_A),
                           dominant_sample = colnames(tab_A)[apply(tab_A, 1, which.max)],
                           share = round(dom_A, 2))
cluster_info

# where do the cells of each A cluster end up in B?
integrated$louvain_A <- merged$louvain[colnames(integrated)]
tab_AB <- table(A = integrated$louvain_A, B = integrated$louvain)

# candidates: large A clusters (>= 500 cells) with >= 75 % cells of one sample;
# for each, the B cluster that receives most of its cells and the share of
# that B cluster coming from the same sample
cand <- cluster_info[cluster_info$cells >= 500 & cluster_info$share >= 0.75, ]
cand$best_B <- sapply(cand$cluster, function(cl) names(which.max(tab_AB[cl, ])))
cand$share_to_B <- sapply(cand$cluster, function(cl) round(max(prop.table(tab_AB, 1)[cl, ]), 2))
cand$B_share_same_sample <- sapply(seq_len(nrow(cand)), function(i) {
  cells_B <- integrated$sample[integrated$louvain == cand$best_B[i]]
  round(mean(cells_B == cand$dominant_sample[i]), 2)
})
cand

# the A cluster that stays a single-sample cluster after integration
cluster_A <- cand$cluster[which.max(cand$B_share_same_sample)]
cluster_B <- cand$best_B[which.max(cand$B_share_same_sample)]
cluster_A; cluster_B
round(prop.table(tab_AB, 1)[cluster_A, ], 2)
# share of B-cluster cells that come from the A cluster, and samples in it
round(prop.table(tab_AB, 2)[cluster_A, cluster_B], 2)
table(integrated$sample[integrated$louvain == cluster_B])

# intronic fraction (bonus task) as an additional quality metric
merged$intron_frac     <- intron_frac[colnames(merged)]
integrated$intron_frac <- intron_frac[colnames(integrated)]

# quality metrics of all clusters
p_qc_A <- VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "intron_frac"),
                  group.by = "louvain", pt.size = 0, ncol = 4, log = FALSE)
p_qc_A
save_fig(p_qc_A, "6_2_qc_per_cluster_A", width = 14, height = 4)

qc_median <- function(obj, cl) {
  m <- obj@meta.data
  m$group <- ifelse(m$louvain == cl, paste("cluster", cl), "all other cells")
  aggregate(cbind(nFeature_RNA, nCount_RNA, percent.mt, intron_frac) ~ group,
            data = m, FUN = median)
}
# medians of all clusters of A
aggregate(cbind(nFeature_RNA, nCount_RNA, percent.mt, intron_frac) ~ louvain,
          data = merged@meta.data, FUN = function(x) round(median(x), 2))
qc_median(merged, cluster_A)
qc_median(integrated, cluster_B)
# the same comparison inside COVID558 only (cluster vs rest of the sample)
qc_median(subset(merged, subset = sample == "COVID558"), cluster_A)

# marker genes of the cluster vs all other cells, sorted by fold change
markers_A <- FindMarkers(merged, ident.1 = cluster_A,
                         min.pct = 0.25, logfc.threshold = 0.25)
markers_A <- markers_A[order(markers_A$avg_log2FC, decreasing = TRUE), ]
head(markers_A[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 25)
tail(markers_A[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 15)
write.csv(head(markers_A, 50), "tables/6_2_markers_cluster_A.csv")

# for B the marker genes are computed on the RNA assay (gene expression)
DefaultAssay(integrated) <- "RNA"
markers_B <- FindMarkers(integrated, ident.1 = cluster_B,
                         min.pct = 0.25, logfc.threshold = 0.25)
markers_B <- markers_B[order(markers_B$avg_log2FC, decreasing = TRUE), ]
head(markers_B[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 25)
DefaultAssay(integrated) <- "integrated"

# typical cell-type genes in the cluster vs the rest
check_genes <- c("CD3D", "IL7R", "NKG7", "GNLY", "MS4A1", "CD79A", "LYZ", "CD14",
                 "S100A8", "S100A9", "FCGR3A", "FCGR3B", "CSF3R", "DEFA3", "MPO",
                 "IGJ", "IGHG1", "HBB", "PPBP", "IFI27", "MALAT1", "MT-CO1")
check_genes <- check_genes[check_genes %in% rownames(merged)]
p_dot <- DotPlot(merged, features = check_genes, group.by = "louvain") +
  RotatedAxis() + ggtitle("A: cell-type and QC genes per cluster")
save_fig(p_dot, "6_2_dotplot_A", width = 10, height = 5)

# remove the low-quality cells (cluster_A of A) from both data sets and repeat
# the clustering
removed_cells <- colnames(merged)[merged$louvain == cluster_A]
length(removed_cells)
table(merged$sample[merged$louvain == cluster_A])

merged2     <- subset(merged, cells = setdiff(colnames(merged), removed_cells))
integrated2 <- subset(integrated, cells = setdiff(colnames(integrated), removed_cells))
merged2$louvain_before     <- merged2$louvain
integrated2$louvain_before <- integrated2$louvain

merged2 <- FindVariableFeatures(merged2, selection.method = "vst", nfeatures = 2000)
merged2 <- ScaleData(merged2)
merged2 <- RunPCA(merged2, verbose = FALSE)
merged2 <- RunUMAP(merged2, dims = 1:n_pcs)
merged2 <- FindNeighbors(merged2, dims = 1:n_pcs)
merged2 <- FindClusters(merged2, resolution = res_A, verbose = FALSE)
merged2$louvain <- merged2$seurat_clusters

DefaultAssay(integrated2) <- "integrated"
integrated2 <- ScaleData(integrated2)
integrated2 <- RunPCA(integrated2, verbose = FALSE)
integrated2 <- RunUMAP(integrated2, dims = 1:n_pcs)
integrated2 <- FindNeighbors(integrated2, dims = 1:n_pcs)
integrated2 <- FindClusters(integrated2, resolution = res_B, verbose = FALSE)
integrated2$louvain <- integrated2$seurat_clusters

p_after <- (DimPlot(merged2, group.by = "louvain", label = TRUE) + ggtitle("A: after removal") |
            DimPlot(integrated2, group.by = "louvain", label = TRUE) + ggtitle("B: after removal")) /
           (DimPlot(merged2, group.by = "sample") + ggtitle("A: by sample") |
            DimPlot(integrated2, group.by = "sample") + ggtitle("B: by sample"))
p_after
save_fig(p_after, "6_2_after_removal", width = 11, height = 9)

tab_A2 <- table(Cluster = merged2$louvain, Sample = merged2$sample)
tab_B2 <- table(Cluster = integrated2$louvain, Sample = integrated2$sample)
tab_A2; tab_B2
write.csv(as.data.frame.matrix(tab_A2), "tables/6_2_cluster_sample_A_after.csv")
write.csv(as.data.frame.matrix(tab_B2), "tables/6_2_cluster_sample_B_after.csv")

# save both data sets with their clusters for Week 3
saveRDS(merged2, "A_merged_clustered.rds")
saveRDS(integrated2, "B_integrated_clustered.rds")

rm(list = setdiff(ls(), c("save_fig", "sample_names")))
invisible(gc())


################################################################################
############################### WEEK 3 #########################################
################################################################################

library(SingleR)
library(celldex)
library(ggrepel)
library(enrichR)

A <- readRDS("A_merged_clustered.rds")
B <- readRDS("B_integrated_clustered.rds")

#### Task 7.1: automatic annotation with SingleR ###############################
# reference: Human Primary Cell Atlas (main labels); SingleR labels every cell
# using the log-normalized expression of the RNA assay
hpca <- HumanPrimaryCellAtlasData()

singler_A <- SingleR(test = GetAssayData(A, assay = "RNA", slot = "data"),
                     ref = hpca, labels = hpca$label.main)
singler_B <- SingleR(test = GetAssayData(B, assay = "RNA", slot = "data"),
                     ref = hpca, labels = hpca$label.main)

# pruned labels: NA for cells whose assignment is ambiguous
A$singler <- singler_A$pruned.labels
B$singler <- singler_B$pruned.labels
sort(table(A$singler, useNA = "ifany"), decreasing = TRUE)
sort(table(B$singler, useNA = "ifany"), decreasing = TRUE)

p_71 <- DimPlot(A, group.by = "singler", label = TRUE, repel = TRUE, label.size = 3) +
          ggtitle("A: merged, SingleR (HPCA)") |
        DimPlot(B, group.by = "singler", label = TRUE, repel = TRUE, label.size = 3) +
          ggtitle("B: integrated, SingleR (HPCA)")
p_71
save_fig(p_71, "7_1_singler_A_B", width = 15, height = 6)

#### Task 7.2: choice of the data set ##########################################
# 7.2.1 do cells of the same cell type from different samples lie together in A?
main_types <- c("T_cells", "NK_cell", "B_cell", "Monocyte")
A_main <- subset(A, cells = colnames(A)[!is.na(A$singler) & A$singler %in% main_types])
p_72 <- DimPlot(A_main, group.by = "sample", split.by = "singler", ncol = 2, pt.size = 0.3) +
  ggtitle("A: cells of the four main SingleR types, coloured by sample")
p_72
save_fig(p_72, "7_2_A_celltype_by_sample", width = 10, height = 9)
rm(A_main)

# SingleR label per sample and per Louvain cluster
table(SingleR = A$singler, Sample = A$sample)
cluster_sample <- table(Cluster = A$louvain, Sample = A$sample)
cluster_sample
round(100 * prop.table(cluster_sample, 1), 1)

# 7.2.2 cell types in the clusters dominated by a single sample (>= 60 %)
dom <- apply(prop.table(cluster_sample, 1), 1, max)
dominated <- names(dom)[dom >= 0.6]
dominated
round(100 * prop.table(table(Cluster = A$louvain, SingleR = A$singler), 1)[dominated, ], 1)
write.csv(as.data.frame.matrix(table(Cluster = A$louvain, SingleR = A$singler)),
          "tables/7_2_cluster_singler_A.csv")

# decision (see report): continue with the merged data set A
obj <- A
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "louvain"
rm(A, B, singler_A, singler_B); invisible(gc())

#### Task 7.3: manual annotation ###############################################
# 7.3.1 marker genes of each cluster vs all other cells
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- do.call(rbind, lapply(split(markers, markers$cluster), function(x)
  head(x[order(x$avg_log2FC, decreasing = TRUE), ], 5)))
top5_table <- aggregate(gene ~ cluster, data = top5, FUN = paste, collapse = ", ")
top5_table
write.csv(markers, "tables/7_3_markers_all_clusters.csv", row.names = FALSE)
write.csv(top5_table, "tables/7_3_top5_markers.csv", row.names = FALSE)

# 7.3.2 marker genes of Table 2 in one dot plot
table2_genes <- c("CD3D", "CD3E", "IL7R", "CD4", "CD8A", "CD8B", "TRDC", "TRGC2",
                  "GNLY", "NKG7", "KLRD1", "KLRF1", "MS4A1", "CD79A", "CD19",
                  "IGJ", "MZB1", "XBP1", "CD14", "LYZ", "S100A8", "S100A9",
                  "FCGR3A", "MS4A7", "LST1", "CD1C", "FCER1A", "CLEC10A",
                  "LILRA4", "CLEC4C", "IL3RA", "FCGR3B", "CXCR2", "CSF3R",
                  "PPBP", "PF4", "HBB", "HBA1", "HBA2", "MKI67", "TOP2A", "CD34")
setdiff(table2_genes, rownames(obj))   # all genes are present

p_dot_t2 <- DotPlot(obj, features = table2_genes, group.by = "louvain") +
  RotatedAxis() + ggtitle("Marker genes of Table 2 per Louvain cluster")
p_dot_t2
save_fig(p_dot_t2, "7_3_dotplot_table2", width = 15, height = 5.5)

# cell type per cluster (from the dot plot and the cluster markers, see report)
cluster_annotation <- c(
  "0"  = "NK",
  "1"  = "T",       # naive T cells (CD4 and CD8, TCF7/LEF1/CCR7)
  "2"  = "CD8",
  "3"  = "CD4",
  "4"  = "B",
  "5"  = "CD14",
  "6"  = "PB",
  "7"  = "Prolif",
  "8"  = "Plt",
  "9"  = "Neu",
  "10" = "RBC",
  "11" = "devNeu"   # developing neutrophils (not in Table 2)
)
obj$cluster_label <- paste0(obj$louvain, " (", cluster_annotation[as.character(obj$louvain)], ")")
obj$cluster_label <- factor(obj$cluster_label,
                            levels = paste0(names(cluster_annotation), " (", cluster_annotation, ")"))
obj$celltype <- factor(cluster_annotation[as.character(obj$louvain)],
                       levels = unique(cluster_annotation))
table(obj$cluster_label)

# 7.3.3 manual annotation on the UMAP, next to the automatic one
p_733 <- DimPlot(obj, group.by = "cluster_label", label = TRUE, repel = TRUE) +
           ggtitle("A: manual annotation") + NoLegend() |
         DimPlot(obj, group.by = "singler", label = TRUE, repel = TRUE, label.size = 3) +
           ggtitle("A: SingleR (HPCA)")
p_733
save_fig(p_733, "7_3_manual_vs_singler", width = 15, height = 6)

manual_vs_auto <- table(Manual = obj$celltype, SingleR = obj$singler)
manual_vs_auto
write.csv(as.data.frame.matrix(manual_vs_auto), "tables/7_3_manual_vs_singler.csv")

# 7.3.4 three marker genes that separate the main cell types
three <- c("CD3D", "MS4A1", "LYZ")
p_vln <- VlnPlot(obj, features = three, group.by = "cluster_label", pt.size = 0, ncol = 1)
p_fp  <- FeaturePlot(obj, features = three, ncol = 3)
p_vln; p_fp
save_fig(p_vln, "7_3_violin_markers", width = 10, height = 9)
save_fig(p_fp,  "7_3_featureplot_markers", width = 15, height = 4.5)

# from now on all clusters of one cell type form one group
Idents(obj) <- "celltype"

#### Task 7.4 (bonus): cell-type proportions ###################################
props <- prop.table(table(Sample = obj$sample, CellType = obj$celltype), 1)
round(100 * props, 1)
write.csv(round(100 * as.data.frame.matrix(props), 2), "tables/7_4_celltype_proportions.csv")

props_df <- as.data.frame(props)
p_props <- ggplot(props_df, aes(x = Sample, y = 100 * Freq, fill = CellType)) +
  geom_col() +
  labs(y = "Share of cells (%)", x = NULL, fill = "Cell type") + theme_bw()
p_props_dodge <- ggplot(props_df, aes(x = CellType, y = 100 * Freq, fill = Sample)) +
  geom_col(position = "dodge") +
  labs(y = "Share of cells (%)", x = NULL) + theme_bw()
p_props | p_props_dodge
save_fig(p_props | p_props_dodge, "7_4_celltype_proportions", width = 14, height = 5)

saveRDS(obj, "week3_annotated.rds")

#### Task 8: differential expression COVID-19 vs healthy #######################
# 8.1 CD14 monocytes: COVID (three patients) vs healthy (HIP043)
de_cd14 <- FindMarkers(obj, ident.1 = "COVID", ident.2 = "Healthy",
                       group.by = "condition", subset.ident = "CD14",
                       logfc.threshold = 0.25, min.pct = 0.1)
de_cd14$gene <- rownames(de_cd14)
de_cd14$direction <- ifelse(de_cd14$p_val_adj < 0.05 & de_cd14$avg_log2FC > 0, "up in COVID",
                     ifelse(de_cd14$p_val_adj < 0.05 & de_cd14$avg_log2FC < 0, "down in COVID",
                            "not significant"))
table(de_cd14$direction)

sig_cd14 <- de_cd14[de_cd14$p_val_adj < 0.05, ]
top_up   <- head(sig_cd14[order(sig_cd14$avg_log2FC, decreasing = TRUE), ], 10)
top_down <- head(sig_cd14[order(sig_cd14$avg_log2FC), ], 10)
top_up[,   c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
top_down[, c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
write.csv(de_cd14, "tables/8_1_de_cd14.csv", row.names = FALSE)

# volcano plot; p-values of 0 are set to the smallest positive value
min_p <- min(de_cd14$p_val_adj[de_cd14$p_val_adj > 0])
de_cd14$log10p <- -log10(pmax(de_cd14$p_val_adj, min_p))
p_volcano <- ggplot(de_cd14, aes(x = avg_log2FC, y = log10p, colour = direction)) +
  geom_point(size = 0.8) +
  geom_text_repel(data = de_cd14[de_cd14$gene %in% c(top_up$gene, top_down$gene), ],
                  aes(label = gene), size = 3, max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = c("up in COVID" = "firebrick", "down in COVID" = "steelblue",
                                 "not significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "CD14 monocytes: COVID-19 vs healthy",
       x = "avg log2 fold change (COVID / healthy)", y = "-log10 adjusted p-value",
       colour = NULL) + theme_bw()
p_volcano
save_fig(p_volcano, "8_1_volcano_cd14", width = 9, height = 6.5)

# 8.2 the same comparison for CD4 T, CD8 T, NK and B cells
de_types <- c("CD14", "CD4", "CD8", "NK", "B")
de_list <- list(CD14 = de_cd14)
for (ct in setdiff(de_types, "CD14")) {
  de <- FindMarkers(obj, ident.1 = "COVID", ident.2 = "Healthy",
                    group.by = "condition", subset.ident = ct,
                    logfc.threshold = 0.25, min.pct = 0.1)
  de$gene <- rownames(de)
  de_list[[ct]] <- de
  write.csv(de, paste0("tables/8_2_de_", ct, ".csv"), row.names = FALSE)
}

de_counts <- data.frame(
  cell_type = de_types,
  cells_COVID   = sapply(de_types, function(ct) sum(obj$celltype == ct & obj$condition == "COVID")),
  cells_healthy = sapply(de_types, function(ct) sum(obj$celltype == ct & obj$condition == "Healthy")),
  up   = sapply(de_list, function(d) sum(d$p_val_adj < 0.05 & d$avg_log2FC > 0)),
  down = sapply(de_list, function(d) sum(d$p_val_adj < 0.05 & d$avg_log2FC < 0)),
  row.names = NULL
)
de_counts

# genes that are significantly up in more than one cell type
up_genes <- lapply(de_list, function(d) d$gene[d$p_val_adj < 0.05 & d$avg_log2FC > 0])
up_count <- sort(table(unlist(up_genes)), decreasing = TRUE)
shared_up <- names(up_count)[up_count >= 3]
length(shared_up)
# log2 fold changes of the most frequently shared genes (NA = not significant)
shared_fc <- sapply(de_list, function(d) {
  s <- d[d$p_val_adj < 0.05 & d$avg_log2FC > 0, ]
  round(s[match(head(shared_up, 30), s$gene), "avg_log2FC"], 2)
})
rownames(shared_fc) <- head(shared_up, 30)
shared_fc

# are the shared genes expressed by these cells themselves? -> where are they
# expressed most strongly (cell types of the COVID samples)?
check_shared <- intersect(c("IGJ", "IGHG1", "IGHM", "IGLC2", "HBB", "HBA2", "IFI27",
                            "IFI6", "ISG15", "MT2A", "S100A8", "S100A9", "LYZ"), shared_up)
check_shared
p_shared <- DotPlot(subset(obj, subset = condition == "COVID"), features = check_shared,
                    group.by = "celltype") + RotatedAxis() +
  ggtitle("Shared up-regulated genes in the cell types of the COVID samples")
p_shared
save_fig(p_shared, "8_2_shared_genes", width = 10, height = 5)

# 8.3 ISGs, HLA class II and cytokines in the CD14 monocytes of each sample
dot_genes <- c("ISG15", "IFI6", "IFI27", "IFITM3", "MX1", "IFI44L",
               "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "CD74",
               "TNF", "IL6", "IL1B", "CCL3", "CCL4", "CXCL2", "CXCL8")
# CXCL8 has its older name IL8 in this annotation
dot_genes[!dot_genes %in% rownames(obj)]
dot_genes[dot_genes == "CXCL8" & !("CXCL8" %in% rownames(obj))] <- "IL8"
stopifnot(all(dot_genes %in% rownames(obj)))

cd14 <- subset(obj, idents = "CD14")
p_dot_cd14 <- DotPlot(cd14, features = dot_genes, group.by = "sample") +
  RotatedAxis() + ggtitle("CD14 monocytes per sample: ISGs | HLA class II | cytokines")
p_dot_cd14
save_fig(p_dot_cd14, "8_3_dotplot_cd14_samples", width = 12, height = 4)
# average expression per sample, for the report
round(AverageExpression(cd14, features = dot_genes, group.by = "sample",
                        assays = "RNA")$RNA, 2)

#### Task 9: pathway analysis (GO, CD14 monocytes) #############################
# up- and down-regulated genes of 8.1 are tested separately with Enrichr
setEnrichrSite("Enrichr")
go_db <- "GO_Biological_Process_2021"
up_cd14   <- sig_cd14$gene[sig_cd14$avg_log2FC > 0]
down_cd14 <- sig_cd14$gene[sig_cd14$avg_log2FC < 0]
length(up_cd14); length(down_cd14)

go_enrich <- list(up   = enrichr(up_cd14,   go_db)[[go_db]],
                  down = enrichr(down_cd14, go_db)[[go_db]])
lapply(go_enrich, function(x) head(x[, c("Term", "Overlap", "P.value", "Adjusted.P.value", "Genes")], 10))
write.csv(go_enrich$up,   "tables/9_go_up_cd14.csv",   row.names = FALSE)
write.csv(go_enrich$down, "tables/9_go_down_cd14.csv", row.names = FALSE)

go_plot <- function(x, title, colour) {
  x <- head(x[order(x$P.value), ], 10)
  x$Term <- factor(x$Term, levels = rev(x$Term))
  ggplot(x, aes(x = -log10(P.value), y = Term)) +
    geom_col(fill = colour) +
    labs(title = title, x = "-log10 p-value", y = NULL) + theme_bw() +
    theme(axis.text.y = element_text(size = 8))
}
p_go <- go_plot(go_enrich$up, "Up in COVID-19 CD14 monocytes", "firebrick") /
        go_plot(go_enrich$down, "Down in COVID-19 CD14 monocytes", "steelblue")
p_go
save_fig(p_go, "9_go_cd14", width = 10, height = 8)

# the same analysis with Seurat's DEenrichRPlot (runs FindMarkers internally)
Idents(cd14) <- "condition"
p_deenrichr <- DEenrichRPlot(cd14, ident.1 = "COVID", ident.2 = "Healthy",
                             enrich.database = go_db, max.genes = 200,
                             num.pathway = 10, balanced = TRUE)
save_fig(p_deenrichr, "9_deenrichrplot_cd14", width = 12, height = 7)


################################################################################
############################### WEEK 4 #########################################
################################################################################
# Requires monocle3 (1.0.0) and CellChat (1.6.1) in addition to the packages
# above; SeuratWrappers provides the Seurat -> monocle3 conversion.

library(monocle3)
library(SeuratWrappers)
library(CellChat)

obj <- readRDS("week3_annotated.rds")
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "celltype"

# convert a Seurat object into a monocle3 cell_data_set that uses the UMAP of
# the Seurat object
to_cds <- function(seu) {
  cds <- as.cell_data_set(seu)
  rowData(cds)$gene_short_name <- rownames(cds)
  cds <- estimate_size_factors(cds)
  cds
}

# graph node closest to a group of cells (function from the monocle3 tutorial)
get_root_node <- function(cds, cells) {
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  igraph::V(principal_graph(cds)[["UMAP"]])$name[
    as.numeric(names(which.max(table(closest_vertex[cells, ]))))]
}

# re-run the dimensionality reduction and clustering on a subset of cells
recluster <- function(seu, resolution = 0.3, n_pcs = 15) {
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:n_pcs)
  seu <- FindNeighbors(seu, dims = 1:n_pcs)
  FindClusters(seu, resolution = resolution)
}

top_markers <- function(seu, n = 10) {
  mk <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  t(sapply(split(mk, mk$cluster), function(x)
    head(x[order(x$avg_log2FC, decreasing = TRUE), "gene"], n)))
}

#### Task 10.1: one trajectory graph for all cells #############################
set.seed(1)
cds_all <- to_cds(obj)
cds_all <- cluster_cells(cds_all, reduction_method = "UMAP")
cds_all <- learn_graph(cds_all)

length(unique(partitions(cds_all)))
table(Partition = partitions(cds_all), CellType = obj$celltype)

p_101 <- plot_cells(cds_all, color_cells_by = "partition", label_groups_by_cluster = FALSE,
                    group_label_size = 4) + ggtitle("All cells: monocle3 partitions") |
         plot_cells(cds_all, color_cells_by = "celltype", label_groups_by_cluster = FALSE,
                    label_cell_groups = TRUE, group_label_size = 3.5) + ggtitle("All cells: cell types")
p_101
save_fig(p_101, "10_1_all_cells_partitions", width = 14, height = 6)
rm(cds_all); invisible(gc())

#### Task 10.2: from B cells to plasmablasts ###################################
# 10.2.1 B cells and plasmablasts: new variable genes, PCA, UMAP, sub-clusters
bcells <- recluster(subset(obj, idents = c("B", "PB")), resolution = 0.3)
table(Subcluster = bcells$seurat_clusters, CellType = bcells$celltype)
table(Subcluster = bcells$seurat_clusters, Sample = bcells$sample)
top_markers(bcells, 12)

b_check <- c("MS4A1", "TCL1A", "IGHD", "FCER2", "IL4R", "FCRL5", "IGHG1", "IGHA1",
             "MZB1", "XBP1", "PRDM1", "MKI67", "CD3D", "GNLY", "NKG7")
p_b_dot <- DotPlot(bcells, features = b_check) + RotatedAxis() +
  ggtitle("Sub-clusters of B cells and plasmablasts")
save_fig(p_b_dot, "10_2_dotplot_subclusters", width = 9, height = 3.5)

b_names <- c("0" = "naive B", "1" = "memory B", "2" = "plasmablasts", "3" = "T/NK cells")
bcells$subcluster <- factor(b_names[as.character(bcells$seurat_clusters)], levels = b_names)
table(bcells$subcluster)

# 10.2.2 trajectory on all these cells, root in the naive B cells
set.seed(1)
cds_b <- to_cds(bcells)
cds_b <- cluster_cells(cds_b, reduction_method = "UMAP")
cds_b <- learn_graph(cds_b, use_partition = FALSE)
cds_b <- order_cells(cds_b, root_pr_nodes =
                       get_root_node(cds_b, colnames(bcells)[bcells$subcluster == "naive B"]))

p_102 <- plot_cells(cds_b, color_cells_by = "subcluster", label_cell_groups = FALSE,
                    label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8) +
           ggtitle("B lineage (all sub-clusters)") |
         plot_cells(cds_b, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                    label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8) +
           ggtitle("Pseudotime (root: naive B)")
p_102
save_fig(p_102, "10_2_trajectory_all_subclusters", width = 13, height = 5.5)
pt_b <- pseudotime(cds_b)
round(tapply(pt_b, bcells$subcluster, median), 2)

# 10.2.3 remove the sub-cluster that is not part of the B lineage and repeat
bcells2 <- subset(bcells, subset = subcluster != "T/NK cells")
bcells2$subcluster <- droplevels(bcells2$subcluster)
bcells2 <- RunPCA(ScaleData(FindVariableFeatures(bcells2, nfeatures = 2000)), verbose = FALSE)
bcells2 <- RunUMAP(bcells2, dims = 1:15)

set.seed(1)
cds_b2 <- to_cds(bcells2)
cds_b2 <- cluster_cells(cds_b2, reduction_method = "UMAP")
cds_b2 <- learn_graph(cds_b2, use_partition = FALSE)
cds_b2 <- order_cells(cds_b2, root_pr_nodes =
                        get_root_node(cds_b2, colnames(bcells2)[bcells2$subcluster == "naive B"]))

p_1023 <- plot_cells(cds_b2, color_cells_by = "subcluster", label_cell_groups = FALSE,
                     label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8) +
            ggtitle("B lineage without T/NK cells") |
          plot_cells(cds_b2, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                     label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.8) +
            ggtitle("Pseudotime (root: naive B)")
p_1023
save_fig(p_1023, "10_2_trajectory_lineage", width = 13, height = 5.5)

# order of the sub-clusters along pseudotime
pt_b2 <- pseudotime(cds_b2)
round(sapply(split(pt_b2, bcells2$subcluster[names(pt_b2)]), median), 2)

b_genes <- c("TCL1A", "MS4A1", "MKI67", "PRDM1", "XBP1", "MZB1")
p_b_genes <- plot_genes_in_pseudotime(cds_b2[b_genes, ], color_cells_by = "subcluster",
                                      min_expr = 0.5, ncol = 2)
p_b_genes
save_fig(p_b_genes, "10_2_genes_in_pseudotime", width = 10, height = 8)

#### Task 10.3: from CD14 to CD16 monocytes ####################################
# 10.3.1 sub-clusters of the monocytes (the CD16 monocytes are part of the
#        CD14 cluster of Week 3)
mono <- recluster(subset(obj, idents = "CD14"), resolution = 0.3)
table(Subcluster = mono$seurat_clusters, Sample = mono$sample)
top_markers(mono, 12)

m_check <- c("CD14", "LYZ", "S100A8", "S100A12", "VCAN", "IFI27", "MX1", "FCGR3A",
             "CDKN1C", "MS4A7", "HLA-DRA", "CD3D", "IL7R", "IRF8", "TCF4", "LILRA4",
             "GZMB", "CD1C", "FCER1A")
p_m_dot <- DotPlot(mono, features = m_check) + RotatedAxis() +
  ggtitle("Sub-clusters of the monocytes")
save_fig(p_m_dot, "10_3_dotplot_subclusters", width = 10, height = 3.5)

m_names <- c("0" = "classical CD14", "1" = "ISG-high CD14", "2" = "T cell-monocyte doublets",
             "3" = "pDC", "4" = "non-classical CD16")
mono$subcluster <- factor(m_names[as.character(mono$seurat_clusters)], levels = m_names)
table(mono$subcluster)

# keep only the monocyte sub-clusters
mono2 <- subset(mono, subset = subcluster %in% c("classical CD14", "ISG-high CD14",
                                                 "non-classical CD16"))
mono2$subcluster <- droplevels(mono2$subcluster)
mono2 <- RunPCA(ScaleData(FindVariableFeatures(mono2, nfeatures = 2000)), verbose = FALSE)
mono2 <- RunUMAP(mono2, dims = 1:15)

# 10.3.2 trajectory with the root in the classical CD14 monocytes
set.seed(1)
cds_m <- to_cds(mono2)
cds_m <- cluster_cells(cds_m, reduction_method = "UMAP")
cds_m <- learn_graph(cds_m, use_partition = FALSE)
cds_m <- order_cells(cds_m, root_pr_nodes =
                       get_root_node(cds_m, colnames(mono2)[mono2$subcluster == "classical CD14"]))

plot_opts <- list(label_cell_groups = FALSE, label_leaves = FALSE,
                  label_branch_points = FALSE, cell_size = 0.8)
p_1032 <- (do.call(plot_cells, c(list(cds_m, color_cells_by = "subcluster"), plot_opts)) +
             ggtitle("Monocytes: sub-clusters") |
           do.call(plot_cells, c(list(cds_m, color_cells_by = "pseudotime"), plot_opts)) +
             ggtitle("Pseudotime (root: classical CD14)") |
           do.call(plot_cells, c(list(cds_m, color_cells_by = "sample"), plot_opts)) +
             ggtitle("Monocytes: samples"))
p_1032
save_fig(p_1032, "10_3_trajectory_root_cd14", width = 18, height = 5.5)

table(Subcluster = mono2$subcluster, Sample = mono2$sample)
pt_m <- pseudotime(cds_m)
round(sapply(split(pt_m, mono2$subcluster[names(pt_m)]), median), 2)

m_genes <- c("S100A8", "VCAN", "HLA-DRA", "FCGR3A", "CDKN1C", "IFI27")
p_m_genes <- plot_genes_in_pseudotime(cds_m[m_genes, ], color_cells_by = "subcluster",
                                      min_expr = 0.5, ncol = 2)
p_m_genes
save_fig(p_m_genes, "10_3_genes_in_pseudotime", width = 10, height = 8)

# 10.3.3 root in the non-classical monocytes instead
cds_m_rev <- order_cells(cds_m, root_pr_nodes =
                           get_root_node(cds_m, colnames(mono2)[mono2$subcluster == "non-classical CD16"]))
pt_m_rev <- pseudotime(cds_m_rev)
round(sapply(split(pt_m_rev, mono2$subcluster[names(pt_m_rev)]), median), 2)
p_1033 <- do.call(plot_cells, c(list(cds_m, color_cells_by = "pseudotime"), plot_opts)) +
            ggtitle("Root: classical CD14") |
          do.call(plot_cells, c(list(cds_m_rev, color_cells_by = "pseudotime"), plot_opts)) +
            ggtitle("Root: non-classical CD16")
p_1033
save_fig(p_1033, "10_3_trajectory_root_comparison", width = 13, height = 5.5)
p_m_genes_rev <- plot_genes_in_pseudotime(cds_m_rev[m_genes, ], color_cells_by = "subcluster",
                                          min_expr = 0.5, ncol = 2)
save_fig(p_m_genes_rev, "10_3_genes_in_pseudotime_root_cd16", width = 10, height = 8)

rm(cds_b, cds_b2, cds_m, cds_m_rev); invisible(gc())

#### Task 11: cell-cell communication COVID-19 vs healthy ######################
# 11.1 CellChat separately for the COVID-19 and the healthy cells, without RBC
obj_cc <- subset(obj, subset = celltype != "RBC")

run_cellchat <- function(seu) {
  seu$celltype <- droplevels(seu$celltype)
  cc <- createCellChat(object = GetAssayData(seu, assay = "RNA", slot = "data"),
                       meta = seu@meta.data, group.by = "celltype")
  cc@DB <- CellChatDB.human            # all categories (secreted, ECM, contact)
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  netAnalysis_computeCentrality(cc, slot.name = "netP")
}

set.seed(1)
cc_covid   <- run_cellchat(subset(obj_cc, subset = condition == "COVID"))
cc_healthy <- run_cellchat(subset(obj_cc, subset = condition == "Healthy"))

# cell groups and their sizes in both conditions
table(obj_cc$celltype, obj_cc$condition)[levels(droplevels(obj_cc$celltype)), ]

# the healthy sample lacks some cell groups -> lift it to the groups of COVID
cc_healthy <- liftCellChat(cc_healthy, group.new = levels(cc_covid@idents))
object.list <- list(Healthy = cc_healthy, COVID = cc_covid)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# total number and strength of the interactions
sapply(object.list, function(x) c(interactions = sum(x@net$count),
                                  strength = round(sum(x@net$weight), 2)))
p_111 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2)) |
         compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
p_111
save_fig(p_111, "11_1_total_interactions", width = 7, height = 4)

# 11.2 differential number and strength between the cell types (COVID - healthy)
# (red: higher in COVID-19, blue: higher in the healthy sample)
png("figures/11_2_diff_count.png", width = 7, height = 7, units = "in", res = 150)
netVisual_diffInteraction(cellchat, weight.scale = TRUE,
                          title.name = "Differential number of interactions")
dev.off()
png("figures/11_2_diff_weight.png", width = 7, height = 7, units = "in", res = 150)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight",
                          title.name = "Differential interaction strength")
dev.off()

# outgoing (sender) and incoming (receiver) change per cell type
diff_count  <- cc_covid@net$count  - cc_healthy@net$count
diff_weight <- cc_covid@net$weight - cc_healthy@net$weight
data.frame(out_count = rowSums(diff_count), in_count = colSums(diff_count),
           out_weight = round(rowSums(diff_weight), 3), in_weight = round(colSums(diff_weight), 3))

# cell groups that send or receive signals in only one condition
active <- function(x) {
  m <- x@net$count
  names(which(rowSums(m) + colSums(m) > 0))
}
setdiff(active(cc_covid), active(cc_healthy))
setdiff(active(cc_healthy), active(cc_covid))

# 11.3 signaling pathways of both conditions
p_rank <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE) |
          rankNet(cellchat, mode = "comparison", stacked = FALSE, do.stat = TRUE)
p_rank
save_fig(p_rank, "11_3_ranknet", width = 12, height = 10)

paths_covid   <- cc_covid@netP$pathways
paths_healthy <- cc_healthy@netP$pathways
setdiff(paths_covid, paths_healthy)   # only in COVID-19
setdiff(paths_healthy, paths_covid)   # only in healthy

# MHC-II pathway: ligand-receptor pairs, senders and receivers
mhc2 <- lapply(object.list, function(x) {
  d <- subsetCommunication(x, signaling = "MHC-II")
  d[order(d$prob, decreasing = TRUE), c("source", "target", "ligand", "receptor", "prob")]
})
lapply(mhc2, function(d) unique(paste(d$ligand, d$receptor, sep = " - ")))
lapply(mhc2, function(d) table(d$source))
lapply(mhc2, function(d) table(d$target))
lapply(mhc2, function(d) round(tapply(d$prob, d$source, sum), 3))

for (n in names(object.list)) {
  png(paste0("figures/11_3_mhc2_circle_", n, ".png"), width = 7, height = 7, units = "in", res = 150)
  netVisual_aggregate(object.list[[n]], signaling = "MHC-II", layout = "circle",
                      signaling.name = paste("MHC-II,", n))
  dev.off()
}

# expression of the MHC-II ligands and receptors in both conditions
p_mhc_dot <- DotPlot(obj_cc, features = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
                                          "HLA-DQA1", "CD4"),
                     group.by = "celltype", split.by = "condition", cols = c("firebrick", "steelblue")) +
  RotatedAxis() + ggtitle("MHC-II ligands and the receptor CD4")
save_fig(p_mhc_dot, "11_3_mhc2_genes", width = 9, height = 7)

saveRDS(cellchat, "week4_cellchat_merged.rds")

sessionInfo()
