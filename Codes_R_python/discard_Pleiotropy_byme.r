library(Seurat)
library(sctransform)
library(anndata)
library(reticulate)
#use_condaenv("tangram-without-squidpy", required = TRUE)
library(SeuratDisk)

library(CellChat)
library(patchwork)
library(BiocNeighbors)
options(stringsAsFactors = FALSE)

library(glue)
library(future)

library(ggplot2)
library(dplyr)
library(forcats)

library(tidyr)

library(tidyverse)
library(rstatix)
library(ggpubr)

library(ggrepel)

library(enrichplot)

library(showtext)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()
# set ggplot2's default font
theme_set(theme_gray(base_family = "Arial"))



glasbey_dark_2 <- c(
  "#a32818", "#0087b1", "#ca0044", "#ffa056", "#eb4d00", "#6b9700", "#528549", "#c8c33f",
  "#91d370", "#4b9793", "#4d230c", "#8300cf", "#9e6e31", "#ac8399", "#c63189", "#015438",
  "#086b83", "#87a8eb", "#6466ef", "#c35dba", "#019e70", "#805059", "#826e8c", "#b3bfda",
  "#ff97b1", "#a793e1", "#698cbd", "#4801cc", "#60006e", "#446966", "#9c5642", "#7bacb5",
  "#cd83bc", "#0054c1", "#7b2f4f", "#fb7c00", "#34bf00", "#ff9c87", "#e1b669", "#526077",
  "#5b3a7c", "#eda5da", "#ef52a3", "#5d7e69", "#c3774f", "#d14867", "#6e00eb", "#1f3400",
  "#c14103", "#6dd4c1", "#46709e", "#a101c3", "#0a8289", "#afa501", "#a55b6b", "#fd77ff",
  "#8a85ae", "#c67ee8", "#9aaa85", "#876bd8", "#01baf6", "#af5dd1", "#59502a", "#b5005e",
  "#7cb569", "#4985ff", "#00c182", "#d195aa", "#a34ba8", "#e205e2", "#16a300", "#382d00",
  "#832f33", "#5d95aa", "#590f00", "#7b4600", "#6e6e31", "#335726", "#4d60b5", "#a19564",
  "#623f28", "#44d457", "#70aacf", "#2d6b4d", "#72af9e", "#fd1500", "#d8b391", "#79893b",
  "#7cc6d8", "#db9036", "#eb605d", "#eb5ed4", "#e47ba7", "#a56b97", "#009744", "#ba5e21",
  "#bcac52", "#87d82f", "#873472", "#aea8d1", "#e28c62", "#d1b1eb", "#36429e", "#3abdc1",
  "#669c4d", "#9e0399", "#4d4d79", "#7b4b85", "#c33431", "#8c6677", "#aa002d", "#7e0175",
  "#01824d", "#724967", "#727790", "#6e0099", "#a0ba52", "#e16e31", "#c46970", "#6d5b95",
  "#a33b74", "#316200", "#87004f", "#335769", "#ba8c7c", "#1859ff", "#909101", "#2b8ad4",
  "#1626ff", "#21d3ff", "#a390af", "#8a6d4f", "#5d213d", "#db03b3", "#6e56ca", "#642821",
  "#ac7700", "#a3bff6", "#b58346", "#9738db", "#b15093", "#7242a3", "#878ed1", "#8970b1",
  "#6baf36", "#5979c8", "#c69eff", "#56831a", "#00d6a7", "#824638", "#11421c", "#59aa75",
  "#905b01", "#f64470", "#ff9703", "#e14231", "#ba91cf", "#34574d", "#f7807c", "#903400",
  "#b3cd00", "#2d9ed3", "#798a9e", "#50807c", "#c136d6", "#eb0552", "#b8ac7e", "#487031",
  "#839564", "#d89c89", "#0064a3", "#4b9077", "#8e6097", "#ff5238", "#a7423b", "#006e70",
  "#97833d", "#dbafc8"
)

# Name glasbey_dark_2
cluster_levels <- as.character(1:31)
colors.use_pre <- glasbey_dark_2[seq_along(cluster_levels)]
names(colors.use_pre) <- cluster_levels

# for E9.5
E9.5_colors.use <- colors.use_pre[!names(colors.use_pre) %in% c("27", "29", "30", "31")]
# for E13.5
E13.5_colors.use <- colors.use_pre

# Parameters
seurat_res <- "0.7"
trim_para <- 0.05



##########################################
#                                        #
#  Load SeuratObject data of each stage  #
#                                        #
##########################################

# Load outputs of SeuratCCA integration of E9.5 and E13.5 datasets
# Leiden clustering resolution = 0.7
# SeuratObjects obtained by splitting "E9.5_bin50_E13.5_bin50_after-CCA-umap_npc=50_res=0.7_Leiden.rds"

E9.5_seurat <- readRDS(file = "/home/project_interconnectivity/result/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA/E9.5_E1S1_bin50_withSeuratCCA_Leiden_res=0.7_clusters_SeuratObject.RDS")
# remove Oi-cluster 27, 29, 30, 31 from E9.5 dataset
E9.5_seurat <- subset(E9.5_seurat, subset = !seurat_clusters %in% c("27","29","30","31"))

E13.5_seurat <- readRDS(file = "/home/project_interconnectivity/result/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA/E13.5_E1S1_bin50_withSeuratCCA_Leiden_res=0.7_clusters_SeuratObject.RDS")

# Set active.ident to seurat_clusters for WhichCells() function
E9.5_seurat <- SetIdent(E9.5_seurat, value = "seurat_clusters")
E13.5_seurat <- SetIdent(E13.5_seurat, value = "seurat_clusters")

# Create output directory for saving results
E9.5_output.dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/Pleiotropy_tau/E9.5/"
E13.5_output.dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/Pleiotropy_tau/E13.5/"
combined_output.dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/Pleiotropy_tau/Combined/"

if(!dir.exists(E9.5_output.dir)) {
  dir.create(E9.5_output.dir, recursive = TRUE)
}

if(!dir.exists(E13.5_output.dir)) {
  dir.create(E13.5_output.dir, recursive = TRUE)
}

if(!dir.exists(combined_output.dir)) {
  dir.create(combined_output.dir, recursive = TRUE)
}



############################################################
#                                                          #
#  Generate pseudo-bulk expression profile per Oi-cluster  #
#                                                          #
############################################################

message(">>>> Apply expression threshold for generating pseudo-bulk expression profile <<<<")
message(">>>> exclude genes expressed in less than ", trim_para*100, "% of spots in each cluster from pseudo-bulk <<<<")
message(">>>> consider these genes as not robustly expressed <<<<")

# keep a copy of original Seurat objects
E9.5_seurat_original <- E9.5_seurat
E13.5_seurat_original <- E13.5_seurat

# preprocess
E9.5_seurat$seurat_clusters  <- droplevels(E9.5_seurat$seurat_clusters)
E13.5_seurat$seurat_clusters <- droplevels(E13.5_seurat$seurat_clusters)


##########
#  E9.5  #
##########

E9.5_cluster_ids <- levels(E9.5_seurat$seurat_clusters)

# get counts matrix for E9.5 dataset
E9.5_counts_all <- LayerData(
  E9.5_seurat,
  assay = "RNA",
  layer = "counts.1"
)

E9.5_pseudo.bulk <- list()

# Generate pseudo-bulk expression profile for each Oi-cluster in E9.5 dataset

for (cl in E9.5_cluster_ids) {
  cells_cl <- WhichCells(E9.5_seurat, idents = cl)
  
  # set expression of genes expressed in less than trim_para of spots to 0
  cl_counts <- E9.5_counts_all[, cells_cl, drop = FALSE]
  
  n_spots <- length(cells_cl)
  gene_expression_proportion <- Matrix::rowSums(cl_counts > 0) / n_spots
  genes_keep <- names(gene_expression_proportion[gene_expression_proportion >= trim_para])
  
  pseudo_bulk <- Matrix::rowSums(cl_counts[genes_keep, , drop = FALSE])
  
  E9.5_pseudo.bulk[[cl]] <- pseudo_bulk
}

E9.5_all_genes <- unique(E9.5_counts_all@Dimnames[[1]])

E9.5_pseudo.bulk_mat <- sapply(
  E9.5_pseudo.bulk,
  function(x) {
    out <- numeric(length(E9.5_all_genes))
    names(out) <- E9.5_all_genes
    out[names(x)] <- x
    out
  }
)

E9.5_pseudo.bulk_df <- as.data.frame(E9.5_pseudo.bulk_mat)
E9.5_pseudo.bulk_df$gene <- rownames(E9.5_pseudo.bulk_df)
E9.5_pseudo.bulk_df <- E9.5_pseudo.bulk_df[
  , c("gene", E9.5_cluster_ids)
]

# save pseudo-bulk expression profile
write.csv(E9.5_pseudo.bulk_df, file = paste0(E9.5_output.dir, "E9.5_pseudo.bulk_Oi-clusters_seurat_res=", seurat_res, "_expression_df_after_trim=", trim_para, ".csv"))

# Obtain remaining genes with expression in each Oi-cluster

E9.5_gene_keep_list <- lapply(E9.5_cluster_ids, function(cl) {
  cells_cl <- WhichCells(E9.5_seurat, idents = cl)
  cl_counts <- E9.5_counts_all[, cells_cl, drop = FALSE]
  rowSums(cl_counts > 0) / ncol(cl_counts) >= trim_para
})

names(E9.5_gene_keep_list) <- E9.5_cluster_ids

E9.5_gene_keep_mask <- do.call(cbind, E9.5_gene_keep_list)
# save E9.5_gene_keep_mask
write.csv(E9.5_gene_keep_mask, file = paste0(E9.5_output.dir, "E9.5_pseudo.bulk_Oi-clusters_seurat_res=", seurat_res, "_gene.keep.mask_after_trim=", trim_para, ".CSV"))

# Check number of remaining genes with expression in each Oi-cluster

E9.5_pseudo_raw <- E9.5_pseudo.bulk_df[, -1]
E9.5_pseudo_raw <- as.matrix(E9.5_pseudo_raw)

E9.5_cluster_expressed_gene_numbers <- colSums(E9.5_pseudo_raw > 0)

E9.5_cluster_expressed_gene_numbers_df <- data.frame(
  cluster = names(E9.5_cluster_expressed_gene_numbers),
  n_expressed_genes = E9.5_cluster_expressed_gene_numbers,
  row.names = NULL
)

# save number of remaining genes with expression as .csv file
write.csv(E9.5_cluster_expressed_gene_numbers_df, file = paste0(E9.5_output.dir, "E9.5_pseudo.bulk_Oi-clusters_seurat_res=", seurat_res, "_Oi-cluster_expressed_gene_numbers_df.CSV"))


###########
#  E13.5  #
###########

E13.5_cluster_ids <- levels(E13.5_seurat$seurat_clusters)

E13.5_counts_all <- LayerData(
  E13.5_seurat,
  assay = "RNA",
  layer = "counts.2"
)

E13.5_pseudo.bulk <- list()

# Generate pseudo-bulk expression profile for each Oi-cluster in E9.5 dataset

for (cl in E13.5_cluster_ids) {
  cells_cl <- WhichCells(E13.5_seurat, idents = cl)
  
  # set expression of genes expressed in less than trim_para of spots to 0
  cl_counts <- E13.5_counts_all[, cells_cl, drop = FALSE]
  
  n_spots <- length(cells_cl)
  gene_expression_proportion <- Matrix::rowSums(cl_counts > 0) / n_spots
  genes_keep <- names(gene_expression_proportion[gene_expression_proportion >= trim_para])
  
  pseudo_bulk <- Matrix::rowSums(cl_counts[genes_keep, , drop = FALSE])
  
  E13.5_pseudo.bulk[[cl]] <- pseudo_bulk
}

E13.5_all_genes <- unique(E13.5_counts_all@Dimnames[[1]])

E13.5_pseudo.bulk_mat <- sapply(
  E13.5_pseudo.bulk,
  function(x) {
    out <- numeric(length(E13.5_all_genes))
    names(out) <- E13.5_all_genes
    out[names(x)] <- x
    out
  }
)

E13.5_pseudo.bulk_df <- as.data.frame(E13.5_pseudo.bulk_mat)
E13.5_pseudo.bulk_df$gene <- rownames(E13.5_pseudo.bulk_df)
E13.5_pseudo.bulk_df <- E13.5_pseudo.bulk_df[
  , c("gene", E13.5_cluster_ids)
]

# save pseudo-bulk expression profile
write.csv(E13.5_pseudo.bulk_df, file = paste0(E13.5_output.dir, "E13.5_pseudo.bulk_Oi-clusters_seurat_res=", seurat_res, "_expression_df_after_trim=", trim_para, ".csv"))

# Obtain remaining genes with expression in each Oi-cluster

E13.5_gene_keep_list <- lapply(E13.5_cluster_ids, function(cl) {
  cells_cl <- WhichCells(E13.5_seurat, idents = cl)
  cl_counts <- E13.5_counts_all[, cells_cl, drop = FALSE]
  rowSums(cl_counts > 0) / ncol(cl_counts) >= trim_para
})

names(E13.5_gene_keep_list) <- E13.5_cluster_ids

E13.5_gene_keep_mask <- do.call(cbind, E13.5_gene_keep_list)
# save E13.5_gene_keep_mask
write.csv(E13.5_gene_keep_mask, file = paste0(E13.5_output.dir, "E13.5_pseudo.bulk_Oi-clusters_seurat_res=", seurat_res, "_gene.keep.mask_after_trim=", trim_para, ".CSV"))

# Check number of remaining genes with expression in each Oi-cluster

E13.5_pseudo_raw <- E13.5_pseudo.bulk_df[, -1]
E13.5_pseudo_raw <- as.matrix(E13.5_pseudo_raw)

E13.5_cluster_expressed_gene_numbers <- colSums(E13.5_pseudo_raw > 0)

E13.5_cluster_expressed_gene_numbers_df <- data.frame(
  cluster = names(E13.5_cluster_expressed_gene_numbers),
  n_expressed_genes = E13.5_cluster_expressed_gene_numbers,
  row.names = NULL
)

# save number of remaining genes with expression as .csv file
write.csv(E13.5_cluster_expressed_gene_numbers_df, file = paste0(E13.5_output.dir, "E13.5_pseudo.bulk_Oi-clusters_seurat_res=", seurat_res, "_Oi-cluster_expressed_gene_numbers_df.CSV"))



############################################################
#                                                          #
#  Normalize each pseudo-bulk to CPM (Counts Per Million)  #
#                                                          #
############################################################

##########
#  E9.5  #
##########

# get read count depth of each Oi-cluster
E9.5_pseudo_bulk_read_counts <- colSums(E9.5_pseudo_raw)
E9.5_pseudo_bulk_read_counts_df <- data.frame(
  cluster = names(E9.5_pseudo_bulk_read_counts),
  n_reads = E9.5_pseudo_bulk_read_counts,
  row.names = NULL
)

E9.5_pseudo_bulk_read_counts_df$cluster <- factor(E9.5_pseudo_bulk_read_counts_df$cluster, levels = as.character(c(1:31)))

# check minimum read count depth across all pseudo-bulk clusters
E9.5_min <- min(E9.5_pseudo_bulk_read_counts)

# CPM-normalization
stopifnot(all(E9.5_pseudo_bulk_read_counts > 0))
E9.5_pseudo_CPM <- sweep(E9.5_pseudo_raw, 2, E9.5_pseudo_bulk_read_counts, "/") * 1e6

E9.5_pseudo_CPM_melt <- reshape2::melt(E9.5_pseudo_CPM)
colnames(E9.5_pseudo_CPM_melt) <- c("gene", "cluster", "expression")
E9.5_pseudo_CPM_melt$cluster <- factor(E9.5_pseudo_CPM_melt$cluster, levels = as.character(c(1:31)))


###########
#  E13.5  #
###########

# get read count depth of each Oi-cluster
E13.5_pseudo_bulk_read_counts <- colSums(E13.5_pseudo_raw)
E13.5_pseudo_bulk_read_counts_df <- data.frame(
  cluster = names(E13.5_pseudo_bulk_read_counts),
  n_reads = E13.5_pseudo_bulk_read_counts,
  row.names = NULL
)

E13.5_pseudo_bulk_read_counts_df$cluster <- factor(E13.5_pseudo_bulk_read_counts_df$cluster, levels = as.character(c(1:31)))

# check minimum read count depth across all pseudo-bulk clusters
E13.5_min <- min(E13.5_pseudo_bulk_read_counts)

# CPM-normalization
stopifnot(all(E13.5_pseudo_bulk_read_counts > 0))
E13.5_pseudo_CPM <- sweep(E13.5_pseudo_raw, 2, E13.5_pseudo_bulk_read_counts, "/") * 1e6

E13.5_pseudo_CPM_melt <- reshape2::melt(E13.5_pseudo_CPM)
colnames(E13.5_pseudo_CPM_melt) <- c("gene", "cluster", "expression")
E13.5_pseudo_CPM_melt$cluster <- factor(E13.5_pseudo_CPM_melt$cluster, levels = as.character(c(1:31)))



#########################
#                       #
#  Calculate Tau index  #
#                       #
#########################

# Calculate the Tau index of tissue specificity

# The following function computes the Tau (tau) index, which measures how specifically a feature is expressed across tissues or conditions.
# Values range from 0 (uniformly expressed) to 1 (specific to a single tissue/condition).

# References:
# Yanai, I., Benjamin, H., Shmoish, M., Chalifa-Caspi, V., Shklar, M., Ophir,
# R., Bar-Even, A., Horn-Saban, S., Safran, M., Domany, E., Lancet, D., &
# Shmueli, O. (2005). Genome-wide midrange transcription profiles reveal
# expression level relationships in human tissue specification.
# Bioinformatics, 21(5), 650-659. \doi{10.1093/bioinformatics/bti042}

tau_index_calc <- function(exp, byRow = TRUE) {
  # transpose so features are always in rows
  if (!byRow) {
    exp <- t(exp)
  }

  # normalize each feature by its maximum value
  max_exp <- apply(exp, 1, max)
  exp_norm <- exp / max_exp

  # Tau = average of (1 - normalized value) across conditions
  tau_index <- Matrix::rowSums(1 - exp_norm) / (ncol(exp) - 1)

  return(tau_index)
}


######################################################################
#  All Oi-clusters (including E13.5-specific Oi-clusters, unpaired)  #
######################################################################

# apply Tau function to E9.5 and E13.5 pseudo-bulk data
E9.5_tau_index <- tau_index_calc(E9.5_pseudo_CPM)
E13.5_tau_index <- tau_index_calc(E13.5_pseudo_CPM)

# combine E9.5 and E13.5
tau_combined_df <- dplyr::bind_rows(
  data.frame(gene = names(E9.5_tau_index), 
             tau = as.numeric(E9.5_tau_index), 
             stage = "E9.5"),
  data.frame(gene = names(E13.5_tau_index), 
             tau = as.numeric(E13.5_tau_index), 
             stage = "E13.5")
) %>%
  dplyr::filter(is.finite(tau))  # remove NA/Inf if any
tau_combined_df$stage <- factor(tau_combined_df$stage, levels = c("E9.5", "E13.5"))


# Clean the tau index data.frame to only common genes shared across E9.5 and E13.5 for comparability

# common genes: intersection of genes with finite Tau in E9.5 (all 27 Oi-clusters) and E13.5 (all 31 Oi-clusters, including 4 E13.5-specific)
common_genes <- intersect(names(E9.5_tau_index), names(E13.5_tau_index))

tau_wide_clean <- tibble::tibble(
  gene  = common_genes,
  `E9.5`  = as.numeric(E9.5_tau_index[common_genes]),
  `E13.5` = as.numeric(E13.5_tau_index[common_genes])
) %>%
  dplyr::filter(is.finite(`E9.5`) & is.finite(`E13.5`))

tau_combined_df_clean <- tau_wide_clean %>%
  tidyr::pivot_longer(cols = c(`E9.5`, `E13.5`),
                      names_to = "stage", values_to = "tau") %>%
  dplyr::mutate(stage = factor(stage, levels = c("E9.5", "E13.5")))

tau_combined_df_clean$stage <- factor(tau_combined_df_clean$stage, levels = c("E9.5", "E13.5"))

### Summary stats of Tau values for quick check

# all genes
tau_summary <- tau_combined_df %>%
  group_by(stage) %>%
  summarise(
    n_genes = n(),
    median_tau = median(tau),
    mean_tau = mean(tau),
    q25 = quantile(tau, 0.25),
    q75 = quantile(tau, 0.75),
    frac_tau_le_0.2 = mean(tau <= 0.2) * 100,
    frac_tau_ge_0.8 = mean(tau >= 0.8) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf(
      "n = %d\nmedian = %.3f (IQR %.3f–%.3f)\nmean = %.3f\nTau ≤ 0.2: %.1f%%\nTau ≥ 0.8: %.1f%%",
      n_genes, median_tau, q25, q75, mean_tau, frac_tau_le_0.2, frac_tau_ge_0.8
    )
  )

# common genes
tau_summary_clean <- tau_combined_df_clean %>%
  group_by(stage) %>%
  summarise(
    n_genes = n(),
    median_tau = median(tau),
    mean_tau = mean(tau),
    q25 = quantile(tau, 0.25),
    q75 = quantile(tau, 0.75),
    frac_tau_le_0.2 = mean(tau <= 0.2) * 100,
    frac_tau_ge_0.8 = mean(tau >= 0.8) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf(
      "n = %d\nmedian = %.3f (IQR %.3f–%.3f)\nmean = %.3f\nTau ≤ 0.2: %.1f%%\nTau ≥ 0.8: %.1f%%",
      n_genes, median_tau, q25, q75, mean_tau, frac_tau_le_0.2, frac_tau_ge_0.8
    )
  )


#######################################################################
#  Shared Oi-clusters (excluding E13.5-specific Oi-clusters, paired)  #
#######################################################################

# apply Tau function to E9.5 and E13.5 pseudo-bulk data
E9.5_tau_index <- tau_index_calc(E9.5_pseudo_CPM)

# remove E13.5-specific clusters for subsequent comparison
rm_cl <- c("27","29","30","31")
E13.5_pseudo_CPM_shared <- E13.5_pseudo_CPM[, !colnames(E13.5_pseudo_CPM) %in% rm_cl, drop = FALSE]

shared_clusters <- intersect(colnames(E9.5_pseudo_CPM), colnames(E13.5_pseudo_CPM_shared))

E9.5_pseudo_CPM_shared  <- E9.5_pseudo_CPM[, shared_clusters, drop = FALSE]
E13.5_pseudo_CPM_shared <- E13.5_pseudo_CPM_shared[, shared_clusters, drop = FALSE]

stopifnot(identical(colnames(E9.5_pseudo_CPM_shared), colnames(E13.5_pseudo_CPM_shared)))

E9.5_tau_index_shared <- tau_index_calc(E9.5_pseudo_CPM_shared)
E13.5_tau_index_shared <- tau_index_calc(E13.5_pseudo_CPM_shared)

# combine E9.5 and E13.5
tau_combined_df_shared <- dplyr::bind_rows(
  data.frame(gene = names(E9.5_tau_index_shared), 
             tau = as.numeric(E9.5_tau_index_shared), 
             stage = "E9.5"),
  data.frame(gene = names(E13.5_tau_index_shared), 
             tau = as.numeric(E13.5_tau_index_shared), 
             stage = "E13.5")
) %>%
  dplyr::filter(is.finite(tau))  # remove NA/Inf if any
tau_combined_df_shared$stage <- factor(tau_combined_df_shared$stage, levels = c("E9.5", "E13.5"))


# Clean the tau index data.frame to only common genes shared across E9.5 and E13.5 for comparability

# common genes: intersection of genes with finite Tau in E9.5 (all 27 Oi-clusters) and E13.5 (shared Oi-clusters, excluding 4 E13.5-specific)
common_genes_shared <- intersect(names(E9.5_tau_index_shared), names(E13.5_tau_index_shared))

tau_wide_clean_shared <- tibble::tibble(
  gene  = common_genes_shared,
  `E9.5`  = as.numeric(E9.5_tau_index_shared[common_genes_shared]),
  `E13.5` = as.numeric(E13.5_tau_index_shared[common_genes_shared])
) %>%
  dplyr::filter(is.finite(`E9.5`) & is.finite(`E13.5`))

cat("common_genes_shared:", length(common_genes_shared), "\n")
cat("kept (finite in both stages):", nrow(tau_wide_clean_shared), "\n")
cat("filtered out:", length(common_genes_shared) - nrow(tau_wide_clean_shared), "\n")

tau_combined_df_clean_shared <- tau_wide_clean_shared %>%
  tidyr::pivot_longer(cols = c(`E9.5`, `E13.5`),
                      names_to = "stage", values_to = "tau") %>%
  dplyr::mutate(stage = factor(stage, levels = c("E9.5", "E13.5")))

tau_combined_df_clean_shared$stage <- factor(tau_combined_df_clean_shared$stage, levels = c("E9.5", "E13.5"))

# Clear any NA in the data
tau_combined_df_clean_shared_filt <- 
  tau_combined_df_clean_shared %>%
  dplyr::filter(is.finite(tau))

tau_combined_df_clean_shared_filt$stage <- factor(tau_combined_df_clean_shared_filt$stage, levels = c("E9.5", "E13.5"))

### Summary stats of Tau values for quick check

# all genes
tau_summary_shared <- tau_combined_df_shared %>%
  group_by(stage) %>%
  summarise(
    n_genes = n(),
    median_tau = median(tau),
    mean_tau = mean(tau),
    q25 = quantile(tau, 0.25),
    q75 = quantile(tau, 0.75),
    frac_tau_le_0.2 = mean(tau <= 0.2) * 100,
    frac_tau_ge_0.8 = mean(tau >= 0.8) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf(
      "n = %d\nmedian = %.3f (IQR %.3f–%.3f)\nmean = %.3f\nTau ≤ 0.2: %.1f%%\nTau ≥ 0.8: %.1f%%",
      n_genes, median_tau, q25, q75, mean_tau, frac_tau_le_0.2, frac_tau_ge_0.8
    )
  )

# common genes
tau_summary_clean_shared <- tau_combined_df_clean_shared %>%
  group_by(stage) %>%
  summarise(
    n_genes = n(),
    median_tau = median(tau),
    mean_tau = mean(tau),
    q25 = quantile(tau, 0.25),
    q75 = quantile(tau, 0.75),
    frac_tau_le_0.2 = mean(tau <= 0.2) * 100,
    frac_tau_ge_0.8 = mean(tau >= 0.8) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf(
      "n = %d\nmedian = %.3f (IQR %.3f–%.3f)\nmean = %.3f\nTau ≤ 0.2: %.1f%%\nTau ≥ 0.8: %.1f%%",
      n_genes, median_tau, q25, q75, mean_tau, frac_tau_le_0.2, frac_tau_ge_0.8
    )
  )



##################################################################################################
#                                                                                                #
#  Detection and cross-stage comparison of pleiotropic genes on common genes/shared Oi-clusters  #
#                                       + Statistical test                                       #
#                                                                                                #
##################################################################################################

# Classify genes by Tau index: 
# those with Tau <= Tau_threshold are considered pleiotropic (broadly expressed). 
# For each cluster, count the expressed genes meeting this criterion.
# Compare the fraction of pleiotropic genes of shared Oi-clusters between stages.

pleiotropy_thresholds <- round(seq(0.1, 0.9, by = 0.1), 1)

for (i in pleiotropy_thresholds) {
  
  pleiotropy_threshold <- i
  
  # create new output.dir with pleiotropy_threshold
  this.E9.5_output.dir <- paste0(E9.5_output.dir, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_shared.clusters/")
  this.E13.5_output.dir <- paste0(E13.5_output.dir, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_shared.clusters/")
  this.combined_output.dir <- paste0(combined_output.dir, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_shared.clusters/")
  
  if (!dir.exists(this.E9.5_output.dir)) {
    dir.create(this.E9.5_output.dir, recursive = TRUE)
  }
  if (!dir.exists(this.E13.5_output.dir)) {
    dir.create(this.E13.5_output.dir, recursive = TRUE)
  }
  if (!dir.exists(this.combined_output.dir)) {
    dir.create(this.combined_output.dir, recursive = TRUE)
  }
  
  # subset the E13.5 data to shared Oi-clusters
  rm_cl <- c("27","29","30","31")
  E13.5_pseudo_raw_shared_cl <- E13.5_pseudo_raw[, !colnames(E13.5_pseudo_raw) %in% rm_cl, drop = FALSE]
  
  # used gene set: common genes between E9.5 (all 27 cl) and E13.5 (shared 27 cl) with finite Tau in both
  used_gene_set <- unique(tau_combined_df_clean_shared_filt$gene)
  used_gene_set <- intersect(used_gene_set, rownames(E9.5_pseudo_raw))
  used_gene_set <- intersect(used_gene_set, rownames(E13.5_pseudo_raw_shared_cl))
  
  # pleiotropic gene in each stage
  pleiotropic_gene_set_E9.5 <- tau_combined_df_clean_shared_filt %>%
    dplyr::filter(stage == "E9.5", tau <= pleiotropy_threshold) %>%
    dplyr::pull(gene) %>% unique() %>%
    intersect(used_gene_set)
  
  pleiotropic_gene_set_E13.5 <- tau_combined_df_clean_shared_filt %>%
    dplyr::filter(stage == "E13.5", tau <= pleiotropy_threshold) %>%
    dplyr::pull(gene) %>% unique() %>%
    intersect(used_gene_set)
  
  
  ##########
  #  E9.5  #
  ##########
  
  # subset the data to only used gene set (common genes)
  E9.5_pseudo_raw_used_subset <- E9.5_pseudo_raw[used_gene_set, , drop = FALSE]
  E9.5_cluster_n_expressed_genes = colSums(E9.5_pseudo_raw_used_subset > 0)
  
  # subset the data to only pleiotropic gene set
  E9.5_pseudo_raw_pleiotropic_subset <- E9.5_pseudo_raw[pleiotropic_gene_set_E9.5, , drop = FALSE]
  E9.5_cluster_n_pleiotropic_genes = colSums(E9.5_pseudo_raw_pleiotropic_subset > 0)
  
  # percentage per Oi-cluster
  E9.5_pleiotropic_perc <- E9.5_cluster_n_pleiotropic_genes / E9.5_cluster_n_expressed_genes
  E9.5_pleiotropic_perc[E9.5_cluster_n_expressed_genes == 0] <- NA_real_
  
  E9.5_pleiotropic_perc_df <- data.frame(
    stage = "E9.5",
    cluster = names(E9.5_pleiotropic_perc),
    pleiotropic_perc = as.numeric(E9.5_pleiotropic_perc),
    row.names = NULL
  )
  
  ###########
  #  E13.5  #
  ###########
  
  # subset the data to only used gene set (common genes)
  E13.5_pseudo_raw_used_subset <- E13.5_pseudo_raw_shared_cl[used_gene_set, , drop = FALSE]
  E13.5_cluster_n_expressed_genes = colSums(E13.5_pseudo_raw_used_subset > 0)
  
  # subset the data to only pleiotropic gene set
  E13.5_pseudo_raw_pleiotropic_subset <- E13.5_pseudo_raw_shared_cl[pleiotropic_gene_set_E13.5, , drop = FALSE]
  E13.5_cluster_n_pleiotropic_genes = colSums(E13.5_pseudo_raw_pleiotropic_subset > 0)
  
  # percentage per Oi-luster
  E13.5_pleiotropic_perc <- E13.5_cluster_n_pleiotropic_genes / E13.5_cluster_n_expressed_genes
  E13.5_pleiotropic_perc[E13.5_cluster_n_expressed_genes == 0] <- NA_real_
  
  E13.5_pleiotropic_perc_df <- data.frame(
    stage = "E13.5",
    cluster = names(E13.5_pleiotropic_perc),
    pleiotropic_perc = as.numeric(E13.5_pleiotropic_perc),
    row.names = NULL
  )
  
  
  
  ###############################################
  #                                             #
  #  Cross-stage comparison + Statistical test  #
  #                                             #
  ###############################################

  # Combine E9.5 and E13.5 results for comparison

  combined_pleiotropic_perc_df <- rbind(
    E9.5_pleiotropic_perc_df,
    E13.5_pleiotropic_perc_df
  )
  
  combined_pleiotropic_perc_df$cluster <- factor(combined_pleiotropic_perc_df$cluster, levels = as.character(c(1:31)))
  combined_pleiotropic_perc_df$stage <- factor(combined_pleiotropic_perc_df$stage, levels = c("E9.5", "E13.5"))

  # save temporary results as CSV
  write.csv(
    combined_pleiotropic_perc_df,
    file = paste0(this.combined_output.dir, "E9.5+E13.5_pseudo.bulk_cluster_percentage_of_pleiotropic_genes_Tau.threshold=", pleiotropy_threshold, "_combined_pleiotropic_perc_df.csv"),
    row.names = FALSE
  )
  

  ### Wilcoxon Signed-Ranks Test to compare the Tau index of shared Oi-clusters between E9.5 and E13.5 (paired)

  # arrange the data for paired comparison
  pleiotropic_perc_long <- combined_pleiotropic_perc_df %>%
    dplyr::select(stage, cluster, pleiotropic_perc) %>%
    dplyr::distinct() %>%
    mutate(stage = factor(stage, levels = c("E9.5", "E13.5")))
  
  pleiotropic_perc_wide <- pleiotropic_perc_long %>%
    tidyr::pivot_wider(names_from = stage, values_from = pleiotropic_perc)
  
  pleiotropic_perc_wide_complete <- pleiotropic_perc_wide %>%
    tidyr::drop_na(`E9.5`, `E13.5`)
  
  combined_pleiotropic_perc_df_paired <- pleiotropic_perc_wide_complete %>%
    tidyr::pivot_longer(cols = c(`E9.5`, `E13.5`), names_to = "stage", values_to = "pleiotropic_perc") %>%
    mutate(stage = factor(stage, levels = c("E9.5", "E13.5"))) %>%
    arrange(cluster, stage)
  
  # perform Wilcoxon Signed-Rank Test (Paired)
  paired_stat.test <- combined_pleiotropic_perc_df_paired %>%
    wilcox_test(pleiotropic_perc ~ stage, paired = TRUE) %>%
    add_significance()
  
  paired_stat.test
  
  # compute effect size
  paired_effect.size <- combined_pleiotropic_perc_df_paired %>%
    wilcox_effsize(pleiotropic_perc ~ stage, paired = TRUE) %>%
    pull(effsize)
  
  
  ### Visualization -- paired dotplot with label
  
  # construct subtitle for plotting
  paired_r.val <- round(paired_effect.size, 3)
  paired_p.label <- signif(paired_stat.test$p, 5)
  paired_subtitle_expr <- bquote(italic(p) == .(paired_p.label) ~ ", " ~ italic(r) == .(paired_r.val))
  
  paired_stat.test <- paired_stat.test %>%
    add_xy_position(x = "stage")

  # paired dotplot with label
  combined_pleiotropic_perc_dotplot_paired_test_with.label <- ggplot(
    combined_pleiotropic_perc_df_paired,
    aes(
      x = stage,
      y = pleiotropic_perc,
      group = cluster
    )
  ) +
    geom_line(
      color = "grey30",
      linewidth = 0.6,
      alpha = 0.8
    ) +
    geom_point(
      aes(color = stage),
      size = 2.5,
      show.legend = FALSE
    ) +
    # E9.5 labels (left)
    geom_text_repel(
      data = subset(
        combined_pleiotropic_perc_df_paired,
        subset = stage == "E9.5"
      ),
      aes(label = cluster),
      nudge_x = -0.25,
      direction = "y",
      hjust = 0.5,
      size = 3,
      color = "black",
      segment.color = "grey80",
      # max.overlaps = Inf
    ) +
    # E13.5 labels (right)
    geom_text_repel(
      data = subset(
        combined_pleiotropic_perc_df_paired,
        subset = stage == "E13.5"
      ),
      aes(label = cluster),
      nudge_x = 0.25,
      direction = "y",
      hjust = 0.5,
      size = 3,
      color = "black",
      segment.color = "grey80",
      # max.overlaps = Inf
    ) +
    scale_color_manual(
      values = c("#eb4d00", "#87a8eb")
    ) +
    stat_pvalue_manual(
      paired_stat.test, tip.length = 0
    ) +
    labs(
      x = "Stage",
      y = "Fraction of pleiotropic genes per Oi-clusters",
      title = paste0("Distribution of percentage of pleiotropic genes (Tau threshold=", pleiotropy_threshold,") (paired)"),
      subtitle = paired_subtitle_expr
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 14, family = "Arial"),
      plot.title = element_text(hjust = 0.5, size = 8),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black")
    )
  
  combined_pleiotropic_perc_dotplot_paired_test_with.label
  
  ggsave(
    filename = paste0(this.combined_output.dir, "E9.5+E13.5_pseudo.bulk_distribution_of_cl_percentage_of_pleio_genes_Tau.threshold=", pleiotropy_threshold, "_dotplot_paired_wilcox_label.pdf"),
    combined_pleiotropic_perc_dotplot_paired_test_with.label,
    width = 5,
    height = 7
  )
}



###############################################################################################
#                                                                                             #
#  Detection and cross-stage comparison of pleiotropic genes on common genes/all Oi-clusters  #
#                                      + Statistical test                                     #
#                                                                                             #
###############################################################################################

# Classify genes by Tau index: 
# those with Tau <= Tau_threshold are considered pleiotropic (broadly expressed). 
# For each cluster, count the expressed genes meeting this criterion.
# Compare the fraction of pleiotropic genes of shared Oi-clusters between stages.

pleiotropy_thresholds <- round(seq(0.1, 0.9, by = 0.1), 1)


for (i in pleiotropy_thresholds) {

  pleiotropy_threshold <- i

  # create new output.dir with pleiotropy_threshold
  this.E9.5_output.dir <- paste0(E9.5_output.dir, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_all.clusters/")
  this.E13.5_output.dir <- paste0(E13.5_output.dir, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_all.clusters/")
  this.combined_output.dir <- paste0(combined_output.dir, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_all.clusters/")

  if (!dir.exists(this.E9.5_output.dir)) {
    dir.create(this.E9.5_output.dir, recursive = TRUE)
  }
  if (!dir.exists(this.E13.5_output.dir)) {
    dir.create(this.E13.5_output.dir, recursive = TRUE)
  }
  if (!dir.exists(this.combined_output.dir)) {
    dir.create(this.combined_output.dir, recursive = TRUE)
  }

  # used gene set: common genes between E9.5 (all 27 cl) and E13.5 (all 31 cl) with finite Tau in both
  tau_combined_df_clean_filt <-
    tau_combined_df_clean %>%
    dplyr::filter(is.finite(tau))

  tau_combined_df_clean_filt$stage <- factor(tau_combined_df_clean_filt$stage, levels = c("E9.5", "E13.5"))

  used_gene_set <- unique(tau_combined_df_clean_filt$gene)
  used_gene_set <- intersect(used_gene_set, rownames(E9.5_pseudo_raw))
  used_gene_set <- intersect(used_gene_set, rownames(E13.5_pseudo_raw))

  # pleiotropic gene in each stage
  pleiotropic_gene_set_E9.5 <- tau_combined_df_clean_filt %>%
    dplyr::filter(stage == "E9.5", tau <= pleiotropy_threshold) %>%
    dplyr::pull(gene) %>% unique() %>%
    intersect(used_gene_set)

  pleiotropic_gene_set_E13.5 <- tau_combined_df_clean_filt %>%
    dplyr::filter(stage == "E13.5", tau <= pleiotropy_threshold) %>%
    dplyr::pull(gene) %>% unique() %>%
    intersect(used_gene_set)


  ############
  #   E9.5   #
  ############

  # subset the data to only used gene set (common genes)
  E9.5_pseudo_raw_used_subset <- E9.5_pseudo_raw[used_gene_set, , drop = FALSE]
  E9.5_cluster_n_expressed_genes = colSums(E9.5_pseudo_raw_used_subset > 0)

  # subset the data to only pleiotropic gene set
  E9.5_pseudo_raw_pleiotropic_subset <- E9.5_pseudo_raw[pleiotropic_gene_set_E9.5, , drop = FALSE]
  E9.5_cluster_n_pleiotropic_genes = colSums(E9.5_pseudo_raw_pleiotropic_subset > 0)

  # percentage per Oi-cluster
  E9.5_pleiotropic_perc <- E9.5_cluster_n_pleiotropic_genes / E9.5_cluster_n_expressed_genes
  E9.5_pleiotropic_perc[E9.5_cluster_n_expressed_genes == 0] <- NA_real_

  E9.5_pleiotropic_perc_df <- data.frame(
    stage = "E9.5",
    cluster = names(E9.5_pleiotropic_perc),
    pleiotropic_perc = as.numeric(E9.5_pleiotropic_perc),
    row.names = NULL
  )

  #############
  #   E13.5   #
  #############

  # subset the data to only used gene set (common genes)
  E13.5_pseudo_raw_used_subset <- E13.5_pseudo_raw[used_gene_set, , drop = FALSE]
  E13.5_cluster_n_expressed_genes = colSums(E13.5_pseudo_raw_used_subset > 0)

  # subset the data to only pleiotropic gene set
  E13.5_pseudo_raw_pleiotropic_subset <- E13.5_pseudo_raw[pleiotropic_gene_set_E13.5, , drop = FALSE]
  E13.5_cluster_n_pleiotropic_genes = colSums(E13.5_pseudo_raw_pleiotropic_subset > 0)

  # percentage per Oi-cluster
  E13.5_pleiotropic_perc <- E13.5_cluster_n_pleiotropic_genes / E13.5_cluster_n_expressed_genes
  E13.5_pleiotropic_perc[E13.5_cluster_n_expressed_genes == 0] <- NA_real_

  E13.5_pleiotropic_perc_df <- data.frame(
    stage = "E13.5",
    cluster = names(E13.5_pleiotropic_perc),
    pleiotropic_perc = as.numeric(E13.5_pleiotropic_perc),
    row.names = NULL
  )



  ###############################################
  #                                             #
  #  Cross-stage comparison + Statistical test  #
  #                                             #
  ###############################################
  
  # Combine E9.5 and E13.5 results for comparison

  combined_pleiotropic_perc_df <- rbind(
    E9.5_pleiotropic_perc_df,
    E13.5_pleiotropic_perc_df
  )

  combined_pleiotropic_perc_df$cluster <- factor(combined_pleiotropic_perc_df$cluster, levels = as.character(c(1:31)))
  combined_pleiotropic_perc_df$stage <- factor(combined_pleiotropic_perc_df$stage, levels = c("E9.5", "E13.5"))

  # save temporary results as CSV
  write.csv(
    combined_pleiotropic_perc_df,
    file = paste0(this.combined_output.dir, "E9.5+E13.5_pseudo.bulk_cluster_percentage_of_pleiotropic_genes_Tau.threshold=", pleiotropy_threshold, "_combined_pleiotropic_perc_df.csv"),
    row.names = FALSE
  )


  ### Mann-Whitney U Test (Wilcoxon Rank Sum Test) to compare the Tau index of all Oi-clusters between E9.5 (27 Oi-clusters) and E13.5 (31 Oi-clusters) (unpaired)

  # perform Mann-Whitney U Test
  unpaired_stat.test <- combined_pleiotropic_perc_df %>%
    wilcox_test(pleiotropic_perc ~ stage) %>%
    add_significance()

  unpaired_stat.test

  # compute effect size
  unpaired_effect.size <- combined_pleiotropic_perc_df %>%
    wilcox_effsize(pleiotropic_perc ~ stage) %>%
    pull(effsize)


  ### Visualization -- unpaired boxplot with label

  # construct subtitle for plotting
  unpaired_r.val <- round(unpaired_effect.size, 3)
  unpaired_p.label <- signif(unpaired_stat.test$p, 5)
  unpaired_subtitle_expr <- bquote(italic(p) == .(unpaired_p.label) ~ ", " ~ italic(r) == .(unpaired_r.val))

  unpaired_stat.test <- unpaired_stat.test %>%
    add_xy_position(x = "stage")

  # unpaired boxplot with label
  combined_pleiotropic_perc_boxplot_with.label <- ggplot(
    combined_pleiotropic_perc_df,
    aes(
      x = stage,
      y = pleiotropic_perc,
      fill = stage
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.37,
    ) +
    geom_jitter(
      color = "black",
      width = 0,
      size = 1.8,
      alpha = 0.8,
      show.legend = FALSE
    ) +
    # E9.5 labels (left)
    geom_text_repel(
      data = subset(
        combined_pleiotropic_perc_df,
        subset = stage == "E9.5"
      ),
      aes(label = cluster),
      nudge_x = -0.25,
      direction = "y",
      hjust = 0.5,
      size = 4,
      color = "black",
      segment.color = "grey80",
      max.overlaps = Inf,
      force = 2,
      force_pull = 0.5,
      box.padding = 0.6,
      point.padding = 0.4
    ) +
    # E13.5 labels (right)
    geom_text_repel(
      data = subset(
        combined_pleiotropic_perc_df,
        subset = stage == "E13.5"
      ),
      aes(label = cluster),
      nudge_x = 0.25,
      direction = "y",
      hjust = 0.5,
      size = 4,
      color = "black",
      segment.color = "grey80",
      max.overlaps = Inf,
      force = 2,
      force_pull = 0.5,
      box.padding = 0.6,
      point.padding = 0.4
    ) +
    scale_fill_manual(
      values = c("#eb4d00", "#87a8eb")
    ) +
    stat_pvalue_manual(
      unpaired_stat.test, tip.length = 0
    ) +
    labs(
      x = "Stage",
      y = "Fraction of pleiotropic genes per Oi-clusters",
      title = paste0("Distribution of percentage of pleiotropic genes (Tau threshold=", pleiotropy_threshold, ") (Unpaired)"),
      subtitle = unpaired_subtitle_expr
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 14, family = "Arial"),
      plot.title = element_text(hjust = 0.5, size = 8),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      axis.text = element_text(color = "black"),
      axis.line = element_line(color = "black")
    )

  combined_pleiotropic_perc_boxplot_with.label

  ggsave(
    filename = paste0(this.combined_output.dir, "E9.5+E13.5_pseudo.bulk_distribution_of_cluster_percentage_of_pleiotropic_genes_Tau.threshold=", pleiotropy_threshold, "_boxplot_unpaired_wilcox_with.label.pdf"),
    combined_pleiotropic_perc_boxplot_with.label,
    width = 5,
    height = 7
  )
}



###################################################################
#                                                                 #
#   Re-plot ALL τ = 0.1 figures (4 analyses × all plots) because the original statistics title could't show correctly          #
#   - reads from existing CSVs                                    #
#   - fixes **** position above max data point                    #
#   - saves with _FIXED suffix (does NOT overwrite originals)     #
#                                                                 #
###################################################################

pleiotropy_threshold <- 0.1

base_combined <- "/home/project_interconnectivity/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/Pseudo-bulk_gene_expression_analysis/Pseudo-bulk_GO_saveTemp_260301/Combined/"

# Pre-compute cluster sizes (needed for Spearman scatter plots)
E9.5_cluster_size <- sapply(E9.5_cluster_ids, function(cl) {
  length(WhichCells(E9.5_seurat, idents = cl))
})
E9.5_cluster_size_df <- data.frame(
  cluster = names(E9.5_cluster_size),
  n_spots = as.integer(E9.5_cluster_size),
  stringsAsFactors = FALSE
)

E13.5_cluster_size <- sapply(E13.5_cluster_ids, function(cl) {
  length(WhichCells(E13.5_seurat, idents = cl))
})
E13.5_cluster_size_df <- data.frame(
  cluster = names(E13.5_cluster_size),
  n_spots = as.integer(E13.5_cluster_size),
  stringsAsFactors = FALSE
)


# =========================================================
# Helper function: regenerate all plots for one analysis
# =========================================================
plot_pleio_analysis_FIXED <- function(combined_output.dir,
                                      threshold,
                                      has_paired = TRUE) {
  
  csv_path <- paste0(combined_output.dir,
                     "E9.5+E13.5_pseudo.bulk_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                     threshold, "_combined_pleiotropic_perc_df.csv")
  
  if (!file.exists(csv_path)) {
    warning("CSV not found, skipping: ", csv_path)
    return(invisible(NULL))
  }
  
  df <- read.csv(csv_path, stringsAsFactors = FALSE) %>%
    mutate(
      cluster = factor(as.character(cluster), levels = as.character(1:31)),
      stage   = factor(stage, levels = c("E9.5", "E13.5"))
    )
  
  ymax_global <- max(df$pleiotropic_perc, na.rm = TRUE)
  
  ##########
  # Barplot
  ##########
  p_barplot <- ggplot(df, aes(x = cluster, y = pleiotropic_perc, fill = stage)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = c("#eb4d00", "#87a8eb")) +
    labs(x = "Oi-clusters (res=0.7)",
         y = "Fraction of pleiotropic genes",
         title = paste0("Percentage of pleiotropic genes (Tau threshold=", threshold, ")")) +
    theme_minimal() +
    theme(text = element_text(size = 14, family = "Arial"),
          strip.text = element_text(size = 16),
          panel.grid = element_blank(),
          axis.line.x = element_line("black"),
          axis.line.y = element_line("black"),
          axis.text = element_text(colour = "black"))
  
  ggsave(paste0(combined_output.dir,
                "E9.5+E13.5_pseudo.bulk_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                threshold, "_barplot_FIXED.pdf"),
         p_barplot, width = 8, height = 4)
  
  ###################
  # Unpaired test
  ###################
  unpaired_stat.test <- df %>%
    wilcox_test(pleiotropic_perc ~ stage) %>%
    add_significance()
  unpaired_effect.size <- df %>%
    wilcox_effsize(pleiotropic_perc ~ stage) %>%
    pull(effsize)
  
  unpaired_r.val   <- round(unpaired_effect.size, 3)
  unpaired_p.label <- signif(unpaired_stat.test$p, 5)
  unpaired_subtitle_expr <- bquote(italic(P) == .(unpaired_p.label) ~ ", " ~
                                     italic(r) == .(unpaired_r.val))
  
  unpaired_stat.test <- unpaired_stat.test %>% add_xy_position(x = "stage")
  unpaired_stat.test$y.position <- ymax_global * 1.08   # ★ FIX
  
  # Unpaired boxplot (no labels)
  p_unpaired_box <- ggplot(df, aes(x = stage, y = pleiotropic_perc, fill = stage)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(color = "black", width = 0.15, size = 1.8, alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = c("#eb4d00", "#87a8eb")) +
    stat_pvalue_manual(unpaired_stat.test, tip.length = 0) +
    labs(x = "Stage",
         y = "Fraction of pleiotropic genes per Oi-cluster",
         title = paste0("Unpaired-distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ")"),
         subtitle = unpaired_subtitle_expr) +
    theme_minimal() +
    theme(text = element_text(size = 14, family = "Arial"),
          plot.title = element_text(hjust = 0.5, size = 8),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          panel.grid = element_blank(), panel.background = element_blank(),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black"))
  
  ggsave(paste0(combined_output.dir,
                "E9.5+E13.5_pseudo.bulk_distribution_of_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                threshold, "_boxplot_unpaired_wilcox_FIXED.pdf"),
         p_unpaired_box, width = 5, height = 7)
  
  # Unpaired boxplot with labels
  p_unpaired_box_label <- ggplot(df, aes(x = stage, y = pleiotropic_perc, fill = stage)) +
    geom_boxplot(outlier.shape = NA, width = 0.37) +
    geom_jitter(color = "black", width = 0, size = 1.8, alpha = 0.8, show.legend = FALSE) +
    geom_text_repel(data = subset(df, stage == "E9.5"),
                    aes(label = cluster), nudge_x = -0.25, direction = "y",
                    hjust = 0.5, size = 4, color = "black", segment.color = "grey80",
                    max.overlaps = Inf, force = 2, force_pull = 0.5,
                    box.padding = 0.6, point.padding = 0.4) +
    geom_text_repel(data = subset(df, stage == "E13.5"),
                    aes(label = cluster), nudge_x = 0.25, direction = "y",
                    hjust = 0.5, size = 4, color = "black", segment.color = "grey80",
                    max.overlaps = Inf, force = 2, force_pull = 0.5,
                    box.padding = 0.6, point.padding = 0.4) +
    scale_fill_manual(values = c("#eb4d00", "#87a8eb")) +
    stat_pvalue_manual(unpaired_stat.test, tip.length = 0) +
    labs(x = "Stage",
         y = "Fraction of pleiotropic genes per Oi-cluster",
         title = paste0("Unpaired-distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ")"),
         subtitle = unpaired_subtitle_expr) +
    theme_minimal() +
    theme(text = element_text(size = 14, family = "Arial"),
          plot.title = element_text(hjust = 0.5, size = 8),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          panel.grid = element_blank(), panel.background = element_blank(),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black"))
  
  ggsave(paste0(combined_output.dir,
                "E9.5+E13.5_pseudo.bulk_distribution_of_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                threshold, "_boxplot_unpaired_wilcox_with.label_FIXED.pdf"),
         p_unpaired_box_label, width = 5, height = 7)
  
  # Unpaired dotplot with labels
  p_unpaired_dot_label <- ggplot(df, aes(x = stage, y = pleiotropic_perc)) +
    geom_point(aes(color = stage), size = 2.5, show.legend = FALSE) +
    geom_text_repel(data = subset(df, stage == "E9.5"),
                    aes(label = cluster), nudge_x = -0.25, direction = "y",
                    hjust = 0.5, size = 4, color = "black", segment.color = "grey80",
                    max.overlaps = Inf, force = 2, force_pull = 0.5,
                    box.padding = 0.6, point.padding = 0.4) +
    geom_text_repel(data = subset(df, stage == "E13.5"),
                    aes(label = cluster), nudge_x = 0.25, direction = "y",
                    hjust = 0.5, size = 4, color = "black", segment.color = "grey80",
                    max.overlaps = Inf, force = 2, force_pull = 0.5,
                    box.padding = 0.6, point.padding = 0.4) +
    scale_fill_manual(values = c("#eb4d00", "#87a8eb")) +
    scale_color_manual(values = c("#eb4d00", "#87a8eb")) +
    stat_pvalue_manual(unpaired_stat.test, tip.length = 0) +
    labs(x = "Stage",
         y = "Fraction of pleiotropic genes per Oi-cluster",
         title = paste0("Unpaired-distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ")"),
         subtitle = unpaired_subtitle_expr) +
    theme_minimal() +
    theme(text = element_text(size = 14, family = "Arial"),
          plot.title = element_text(hjust = 0.5, size = 8),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          panel.grid = element_blank(), panel.background = element_blank(),
          legend.position = "none",
          axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black"))
  
  ggsave(paste0(combined_output.dir,
                "E9.5+E13.5_pseudo.bulk_distribution_of_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                threshold, "_dotplot_unpaired_wilcox_with.label_FIXED.pdf"),
         p_unpaired_dot_label, width = 5, height = 7)
  
  #################
  # Paired test
  #################
  if (has_paired) {
    df_long <- df %>%
      dplyr::select(stage, cluster, pleiotropic_perc) %>%
      dplyr::distinct()
    df_wide <- df_long %>%
      tidyr::pivot_wider(names_from = stage, values_from = pleiotropic_perc) %>%
      tidyr::drop_na(`E9.5`, `E13.5`)
    df_paired <- df_wide %>%
      tidyr::pivot_longer(cols = c(`E9.5`, `E13.5`),
                          names_to = "stage", values_to = "pleiotropic_perc") %>%
      mutate(stage = factor(stage, levels = c("E9.5", "E13.5"))) %>%
      arrange(cluster, stage)
    
    paired_stat.test <- df_paired %>%
      wilcox_test(pleiotropic_perc ~ stage, paired = TRUE) %>%
      add_significance()
    paired_effect.size <- df_paired %>%
      wilcox_effsize(pleiotropic_perc ~ stage, paired = TRUE) %>%
      pull(effsize)
    
    paired_r.val   <- round(paired_effect.size, 3)
    paired_p.label <- signif(paired_stat.test$p, 5)
    paired_subtitle_expr <- bquote(italic(P) == .(paired_p.label) ~ ", " ~
                                     italic(r) == .(paired_r.val))
    
    paired_stat.test <- paired_stat.test %>% add_xy_position(x = "stage")
    paired_stat.test$y.position <- ymax_global * 1.08   # ★ FIX
    
    # Paired boxplot
    p_paired_box <- ggplot(df_paired, aes(x = stage, y = pleiotropic_perc, fill = stage)) +
      geom_boxplot(outlier.shape = NA, width = 0.5) +
      geom_jitter(color = "black", width = 0.15, size = 1.8, alpha = 0.8, show.legend = FALSE) +
      scale_fill_manual(values = c("#eb4d00", "#87a8eb")) +
      stat_pvalue_manual(paired_stat.test, tip.length = 0) +
      labs(x = "Stage",
           y = "Fraction of pleiotropic genes per Oi-cluster",
           title = paste0("Paired-distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ")"),
           subtitle = paired_subtitle_expr) +
      theme_minimal() +
      theme(text = element_text(size = 14, family = "Arial"),
            plot.title = element_text(hjust = 0.5, size = 8),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            panel.grid = element_blank(), panel.background = element_blank(),
            legend.position = "none",
            axis.text = element_text(color = "black"),
            axis.line = element_line(color = "black"))
    
    ggsave(paste0(combined_output.dir,
                  "E9.5+E13.5_pseudo.bulk_distribution_of_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                  threshold, "_boxplot_paired_wilcox_FIXED.pdf"),
           p_paired_box, width = 5, height = 7)
    
    # Paired dotplot (no labels)
    p_paired_dot <- ggplot(df_paired, aes(x = stage, y = pleiotropic_perc, group = cluster)) +
      geom_line(color = "grey30", linewidth = 0.6, alpha = 0.8) +
      geom_point(aes(color = stage), size = 2.5, show.legend = FALSE) +
      scale_fill_manual(values = c("#eb4d00", "#87a8eb")) +
      scale_color_manual(values = c("#eb4d00", "#87a8eb")) +
      stat_pvalue_manual(paired_stat.test, tip.length = 0) +
      labs(x = "Stage",
           y = "Fraction of pleiotropic genes per Oi-cluster",
           title = paste0("Paired-distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ")"),
           subtitle = paired_subtitle_expr) +
      theme_minimal() +
      theme(text = element_text(size = 14, family = "Arial"),
            plot.title = element_text(hjust = 0.5, size = 8),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            panel.grid = element_blank(), panel.background = element_blank(),
            legend.position = "none",
            axis.text = element_text(color = "black"),
            axis.line = element_line(color = "black"))
    
    ggsave(paste0(combined_output.dir,
                  "E9.5+E13.5_pseudo.bulk_distribution_of_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                  threshold, "_dotplot_paired_wilcox_FIXED.pdf"),
           p_paired_dot, width = 5, height = 7)
    
    # Paired dotplot with labels
    p_paired_dot_label <- ggplot(df_paired, aes(x = stage, y = pleiotropic_perc, group = cluster)) +
      geom_line(color = "grey30", linewidth = 0.6, alpha = 0.8) +
      geom_point(aes(color = stage), size = 2.5, show.legend = FALSE) +
      geom_text_repel(data = subset(df_paired, stage == "E9.5"),
                      aes(label = cluster), nudge_x = -0.25, direction = "y",
                      hjust = 0.5, size = 3, color = "black", segment.color = "grey80") +
      geom_text_repel(data = subset(df_paired, stage == "E13.5"),
                      aes(label = cluster), nudge_x = 0.25, direction = "y",
                      hjust = 0.5, size = 3, color = "black", segment.color = "grey80") +
      scale_fill_manual(values = c("#eb4d00", "#87a8eb")) +
      scale_color_manual(values = c("#eb4d00", "#87a8eb")) +
      stat_pvalue_manual(paired_stat.test, tip.length = 0) +
      labs(x = "Stage",
           y = "Fraction of pleiotropic genes per Oi-cluster",
           title = paste0("Paired-distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ")"),
           subtitle = paired_subtitle_expr) +
      theme_minimal() +
      theme(text = element_text(size = 14, family = "Arial"),
            plot.title = element_text(hjust = 0.5, size = 8),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            panel.grid = element_blank(), panel.background = element_blank(),
            legend.position = "none",
            axis.text = element_text(color = "black"),
            axis.line = element_line(color = "black"))
    
    ggsave(paste0(combined_output.dir,
                  "E9.5+E13.5_pseudo.bulk_distribution_of_cl_percentage_of_pleio_genes_Tau.threshold=",
                  threshold, "_dotplot_paired_wilcox_label_FIXED.pdf"),
           p_paired_dot_label, width = 5, height = 7)
  }
  
  ###################################
  # Spearman scatter plots
  ###################################
  E9.5_scatter_df <- df %>%
    filter(stage == "E9.5") %>%
    mutate(cluster = as.character(cluster)) %>%
    left_join(E9.5_cluster_size_df, by = "cluster")
  
  E13.5_scatter_df <- df %>%
    filter(stage == "E13.5") %>%
    mutate(cluster = as.character(cluster)) %>%
    inner_join(E13.5_cluster_size_df, by = "cluster")
  
  for (stage_name in c("E9.5", "E13.5")) {
    scatter_df <- if (stage_name == "E9.5") E9.5_scatter_df else E13.5_scatter_df
    
    scatter_p <- ggplot(scatter_df, aes(x = pleiotropic_perc, y = n_spots)) +
      geom_point() +
      geom_smooth(method = "lm", color = "darkgreen", formula = y ~ x, se = TRUE) +
      ggpubr::stat_cor(method = "spearman",
                       label.x.npc = "left", label.y.npc = "top",
                       size = 4, color = "black") +
      labs(x = "Fraction of pleiotropic genes",
           y = "Numbers of spatial spots",
           title = paste0(stage_name, "_Spearman correlation between numbers of spatial spots and fraction of pleiotropic genes")) +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 8),
            axis.line = element_line(color = "black"),
            axis.text = element_text(color = "black"),
            axis.title = element_text(face = "bold"),
            legend.position = "none")
    
    ggsave(file.path(combined_output.dir,
                     glue("{stage_name}_numbers of spatial spots_fraction of pleiotropic genes_scatter plot_FIXED.pdf")),
           plot = scatter_p, width = 7, height = 5)
  }
  
  message("Done: ", combined_output.dir)
  invisible(NULL)
}


# =========================================================
# Run all 4 analyses at τ = 0.1
# =========================================================

# 1. common.genes / shared.clusters (has paired)
plot_pleio_analysis_FIXED(
  combined_output.dir = paste0(base_combined, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_shared.clusters/"),
  threshold = pleiotropy_threshold,
  has_paired = TRUE
)

# 2. common.genes / all.clusters (no paired)
plot_pleio_analysis_FIXED(
  combined_output.dir = paste0(base_combined, "tau_pleiotropy.threshold=", pleiotropy_threshold, "/common.genes_all.clusters/"),
  threshold = pleiotropy_threshold,
  has_paired = FALSE
)
