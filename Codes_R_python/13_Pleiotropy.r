library(Seurat)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rstatix)
library(ggpubr)
library(ggrepel)
options(stringsAsFactors = FALSE)

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
cluster_levels  <- as.character(1:31)
colors.use_pre  <- glasbey_dark_2[seq_along(cluster_levels)]
names(colors.use_pre) <- cluster_levels

# for E9.5 (excludes E13.5-specific clusters 27, 29, 30, 31)
E9.5_colors.use  <- colors.use_pre[!names(colors.use_pre) %in% c("27", "29", "30", "31")]
# for E13.5 (all clusters)
E13.5_colors.use <- colors.use_pre

# Stage colours reused across all comparison plots
STAGE_COLORS    <- c("E9.5" = "#eb4d00", "E13.5" = "#87a8eb")
E13.5_SPECIFIC  <- c("27", "29", "30", "31")
CLUSTER_LEVELS  <- as.character(1:31)


# Parameters

seurat_res <- "0.7"
trim_para  <- 0.05

# thresholds at which a gene is called "pleiotropic" (Tau <= threshold)
pleiotropy_thresholds <- round(seq(0.1, 0.9, by = 0.1), 1)


#  Output paths

base_dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/Pleiotropy_tau/"

E9.5_output.dir     <- paste0(base_dir, "E9.5/")
E13.5_output.dir    <- paste0(base_dir, "E13.5/")
combined_output.dir <- paste0(base_dir, "Combined/")

for (d in c(E9.5_output.dir, E13.5_output.dir, combined_output.dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


##########################################
#                                        #
#  Load SeuratObject data of each stage  #
#                                        #
##########################################

# Load outputs of SeuratCCA integration of E9.5 and E13.5 datasets
# Leiden clustering resolution = 0.7
# SeuratObjects obtained by splitting
# "E9.5_bin50_E13.5_bin50_after-CCA-umap_npc=50_res=0.7_Leiden.rds"

E9.5_seurat <- readRDS(file = "/home/project_interconnectivity/result/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA/E9.5_E1S1_bin50_withSeuratCCA_Leiden_res=0.7_clusters_SeuratObject.RDS")
# remove Oi-cluster 27, 29, 30, 31 from E9.5 dataset (E13.5-specific)
E9.5_seurat <- subset(E9.5_seurat, subset = !seurat_clusters %in% E13.5_SPECIFIC)

E13.5_seurat <- readRDS(file = "/home/project_interconnectivity/result/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA/E13.5_E1S1_bin50_withSeuratCCA_Leiden_res=0.7_clusters_SeuratObject.RDS")

# Set active.ident to seurat_clusters for WhichCells()
E9.5_seurat  <- SetIdent(E9.5_seurat,  value = "seurat_clusters")
E13.5_seurat <- SetIdent(E13.5_seurat, value = "seurat_clusters")


############################################################
#                                                          #
#  Generate pseudo-bulk expression profile per Oi-cluster  #
#                                                          #
############################################################

message(">>>> Apply expression threshold for generating pseudo-bulk expression profile <<<<")
message(">>>> exclude genes expressed in less than ", trim_para * 100, "% of spots in each cluster from pseudo-bulk <<<<")
message(">>>> consider these genes as not robustly expressed <<<<")

# keep a copy of original Seurat objects
E9.5_seurat_original  <- E9.5_seurat
E13.5_seurat_original <- E13.5_seurat

# preprocess
E9.5_seurat$seurat_clusters  <- droplevels(E9.5_seurat$seurat_clusters)
E13.5_seurat$seurat_clusters <- droplevels(E13.5_seurat$seurat_clusters)


##################################################################################
#  Helper: build pseudo-bulk profile + per-cluster gene-keep mask for one stage  #
##################################################################################

# Returns a list with:
#   pseudo_df   : data.frame (gene + one column per cluster) of summed counts
#                 over genes passing the trim threshold within each cluster
#   pseudo_raw  : numeric matrix (genes x clusters) version of pseudo_df
#   keep_mask   : logical matrix (genes x clusters); TRUE if gene passes trim

build_pseudobulk <- function(seurat_obj, counts_layer, trim) {

  cluster_ids <- levels(seurat_obj$seurat_clusters)
  counts_all  <- LayerData(seurat_obj, assay = "RNA", layer = counts_layer)
  all_genes   <- unique(counts_all@Dimnames[[1]])

  # per-cluster: trim low-prevalence genes, then sum counts -> pseudo-bulk
  pseudo_list <- list()
  keep_list   <- list()
  for (cl in cluster_ids) {
    cells_cl  <- WhichCells(seurat_obj, idents = cl)
    cl_counts <- counts_all[, cells_cl, drop = FALSE]
    n_spots   <- length(cells_cl)

    gene_prop  <- Matrix::rowSums(cl_counts > 0) / n_spots
    genes_keep <- names(gene_prop[gene_prop >= trim])

    pseudo_list[[cl]] <- Matrix::rowSums(cl_counts[genes_keep, , drop = FALSE])
    keep_list[[cl]]   <- gene_prop >= trim
  }

  # assemble dense gene x cluster matrix (missing genes -> 0)
  pseudo_mat <- sapply(pseudo_list, function(x) {
    out <- numeric(length(all_genes))
    names(out) <- all_genes
    out[names(x)] <- x
    out
  })

  pseudo_df <- as.data.frame(pseudo_mat)
  pseudo_df$gene <- rownames(pseudo_df)
  pseudo_df <- pseudo_df[, c("gene", cluster_ids)]

  pseudo_raw <- as.matrix(pseudo_df[, -1])

  keep_mask <- do.call(cbind, keep_list)
  colnames(keep_mask) <- cluster_ids

  list(
    cluster_ids = cluster_ids,
    pseudo_df   = pseudo_df,
    pseudo_raw  = pseudo_raw,
    keep_mask   = keep_mask
  )
}

################################################################################
#  Helper: write per-stage pseudo-bulk QC files (profile, keep-mask, gene #s)  #
################################################################################

save_pseudobulk_qc <- function(pb, stage, out_dir) {
  prefix <- paste0(out_dir, stage, "_pseudo.bulk_Oi-clusters_seurat_res=", seurat_res)

  # pseudo-bulk expression profile
  write.csv(pb$pseudo_df,
            file = paste0(prefix, "_expression_df_after_trim=", trim_para, ".csv"))

  # per-cluster gene-keep mask
  write.csv(pb$keep_mask,
            file = paste0(prefix, "_gene.keep.mask_after_trim=", trim_para, ".CSV"))

  # number of expressed genes per cluster
  n_expressed_df <- data.frame(
    cluster           = colnames(pb$pseudo_raw),
    n_expressed_genes = colSums(pb$pseudo_raw > 0),
    row.names         = NULL
  )
  write.csv(n_expressed_df,
            file = paste0(prefix, "_Oi-cluster_expressed_gene_numbers_df.CSV"))
}


### E9.5

E9.5_pb          <- build_pseudobulk(E9.5_seurat,  counts_layer = "counts.1", trim = trim_para)
E9.5_cluster_ids <- E9.5_pb$cluster_ids
E9.5_pseudo_raw  <- E9.5_pb$pseudo_raw
save_pseudobulk_qc(E9.5_pb, "E9.5", E9.5_output.dir)

### E13.5
E13.5_pb          <- build_pseudobulk(E13.5_seurat, counts_layer = "counts.2", trim = trim_para)
E13.5_cluster_ids <- E13.5_pb$cluster_ids
E13.5_pseudo_raw  <- E13.5_pb$pseudo_raw
save_pseudobulk_qc(E13.5_pb, "E13.5", E13.5_output.dir)


###############################################################################
#                                                                             #
#  Normalize each pseudo-bulk to CPM (Counts Per Million)                     #
#                                                                             #
###############################################################################

#################################################################
#  Helper: CPM-normalize a (genes x clusters) raw count matrix  #
#################################################################

cpm_normalize <- function(pseudo_raw) {
  read_counts <- colSums(pseudo_raw)
  stopifnot(all(read_counts > 0))
  sweep(pseudo_raw, 2, read_counts, "/") * 1e6
}

E9.5_pseudo_CPM  <- cpm_normalize(E9.5_pseudo_raw)
E13.5_pseudo_CPM <- cpm_normalize(E13.5_pseudo_raw)


###############################################################################
#                                                                             #
#  Calculate Tau index of tissue specificity                                  #
#                                                                             #
###############################################################################

# The Tau (tau) index measures how specifically a feature is expressed across
# tissues/conditions. Values range from 0 (uniformly expressed) to 1 (specific
# to a single tissue/condition).
#
# Reference:
#   Yanai I. et al. (2005) Genome-wide midrange transcription profiles reveal
#   expression level relationships in human tissue specification.
#   Bioinformatics 21(5):650-659. doi:10.1093/bioinformatics/bti042

tau_index_calc <- function(exp, byRow = TRUE) {
  # transpose so features are always in rows
  if (!byRow) exp <- t(exp)

  # normalize each feature by its maximum value
  max_exp  <- apply(exp, 1, max)
  exp_norm <- exp / max_exp

  # Tau = average of (1 - normalized value) across conditions
  Matrix::rowSums(1 - exp_norm) / (ncol(exp) - 1)
}

################################################################
#  Helper: per-stage Tau summary stats with a formatted label  #
################################################################

tau_summary_stats <- function(df) {
  df %>%
    group_by(stage) %>%
    summarise(
      n_genes         = n(),
      median_tau      = median(tau),
      mean_tau        = mean(tau),
      q25             = quantile(tau, 0.25),
      q75             = quantile(tau, 0.75),
      frac_tau_le_0.2 = mean(tau <= 0.2) * 100,
      frac_tau_ge_0.8 = mean(tau >= 0.8) * 100,
      .groups = "drop"
    ) %>%
    mutate(
      label = sprintf(
        "n = %d\nmedian = %.3f (IQR %.3f-%.3f)\nmean = %.3f\nTau <= 0.2: %.1f%%\nTau >= 0.8: %.1f%%",
        n_genes, median_tau, q25, q75, mean_tau, frac_tau_le_0.2, frac_tau_ge_0.8
      )
    )
}

###############################################################################
#  Helper: build long + clean(common-gene) Tau data.frames for a stage pair,  #
#  plus matching summary-stat tables.                                         #
###############################################################################

# Returns a list with:
#   long          : gene/tau/stage for all genes with finite Tau in each stage
#   clean         : gene/tau/stage restricted to genes finite in BOTH stages
#   summary       : per-stage summary stats on `long`
#   summary_clean : per-stage summary stats on `clean`
build_tau_tables <- function(tau_E9.5, tau_E13.5) {

  long <- dplyr::bind_rows(
    data.frame(gene = names(tau_E9.5),  tau = as.numeric(tau_E9.5),  stage = "E9.5"),
    data.frame(gene = names(tau_E13.5), tau = as.numeric(tau_E13.5), stage = "E13.5")
  ) %>%
    dplyr::filter(is.finite(tau))
  long$stage <- factor(long$stage, levels = c("E9.5", "E13.5"))

  # common genes: finite Tau in both stages
  common_genes <- intersect(names(tau_E9.5), names(tau_E13.5))
  wide_clean <- tibble::tibble(
    gene    = common_genes,
    `E9.5`  = as.numeric(tau_E9.5[common_genes]),
    `E13.5` = as.numeric(tau_E13.5[common_genes])
  ) %>%
    dplyr::filter(is.finite(`E9.5`) & is.finite(`E13.5`))

  clean <- wide_clean %>%
    tidyr::pivot_longer(cols = c(`E9.5`, `E13.5`),
                        names_to = "stage", values_to = "tau") %>%
    dplyr::mutate(stage = factor(stage, levels = c("E9.5", "E13.5"))) %>%
    dplyr::filter(is.finite(tau))

  list(
    long          = long,
    clean         = clean,
    summary       = tau_summary_stats(long),
    summary_clean = tau_summary_stats(clean)
  )
}


###  (A) All Oi-clusters (including E13.5-specific clusters; unpaired scope)

E9.5_tau_index  <- tau_index_calc(E9.5_pseudo_CPM)
E13.5_tau_index <- tau_index_calc(E13.5_pseudo_CPM)

tau_all <- build_tau_tables(E9.5_tau_index, E13.5_tau_index)

# expose with the original variable names for downstream readability
tau_combined_df        <- tau_all$long
tau_combined_df_clean  <- tau_all$clean
tau_summary            <- tau_all$summary
tau_summary_clean      <- tau_all$summary_clean


###  (B) Shared Oi-clusters (excluding E13.5-specific clusters; paired scope)

E13.5_pseudo_CPM_shared <- E13.5_pseudo_CPM[, !colnames(E13.5_pseudo_CPM) %in% E13.5_SPECIFIC, drop = FALSE]
shared_clusters         <- intersect(colnames(E9.5_pseudo_CPM), colnames(E13.5_pseudo_CPM_shared))

E9.5_pseudo_CPM_shared  <- E9.5_pseudo_CPM[, shared_clusters, drop = FALSE]
E13.5_pseudo_CPM_shared <- E13.5_pseudo_CPM_shared[, shared_clusters, drop = FALSE]
stopifnot(identical(colnames(E9.5_pseudo_CPM_shared), colnames(E13.5_pseudo_CPM_shared)))

E9.5_tau_index_shared  <- tau_index_calc(E9.5_pseudo_CPM_shared)
E13.5_tau_index_shared <- tau_index_calc(E13.5_pseudo_CPM_shared)

tau_shared <- build_tau_tables(E9.5_tau_index_shared, E13.5_tau_index_shared)

tau_combined_df_shared            <- tau_shared$long
tau_combined_df_clean_shared      <- tau_shared$clean   # already finite-filtered
tau_combined_df_clean_shared_filt <- tau_shared$clean
tau_summary_shared                <- tau_shared$summary
tau_summary_clean_shared          <- tau_shared$summary_clean

cat("common_genes_shared (finite in both stages):",
    nrow(dplyr::distinct(tau_combined_df_clean_shared, gene)), "\n")


#########################################################
#                                                       #
#  Pleiotropic-gene detection + cross-stage comparison  #
#                                                       #
#########################################################

#  For a chosen pleiotropy threshold a gene is "pleiotropic" if Tau <=
#  threshold. Within each cluster we compute the fraction of expressed
#  (common) genes that are pleiotropic, then compare that fraction between
#  stages.

#####################################################################
#  Helper: per-cluster fraction of pleiotropic genes for one stage  #
#####################################################################

#   pseudo_raw       : genes x clusters raw pseudo-bulk for the stage
#   used_gene_set    : common-gene universe (rows to consider as "expressed")
#   pleiotropic_genes: genes called pleiotropic in this stage

pleiotropic_fraction <- function(pseudo_raw, used_gene_set, pleiotropic_genes, stage) {
  used_subset        <- pseudo_raw[used_gene_set, , drop = FALSE]
  n_expressed        <- colSums(used_subset > 0)

  pleio_subset       <- pseudo_raw[pleiotropic_genes, , drop = FALSE]
  n_pleiotropic      <- colSums(pleio_subset > 0)

  perc               <- n_pleiotropic / n_expressed
  perc[n_expressed == 0] <- NA_real_

  data.frame(
    stage            = stage,
    cluster          = names(perc),
    pleiotropic_perc = as.numeric(perc),
    row.names        = NULL
  )
}

################################################################################
#  Helper: compute combined per-cluster pleiotropic fractions for both stages  #
################################################################################

#   tau_clean       : clean (common-gene) long Tau table for this scope
#   E9.5_raw / E13.5_raw : pseudo-bulk matrices restricted to the scope's clusters
#   threshold       : pleiotropy threshold

compute_pleiotropic_perc <- function(tau_clean, E9.5_raw, E13.5_raw, threshold) {

  used_gene_set <- unique(tau_clean$gene) %>%
    intersect(rownames(E9.5_raw)) %>%
    intersect(rownames(E13.5_raw))

  pleio_E9.5 <- tau_clean %>%
    dplyr::filter(stage == "E9.5", tau <= threshold) %>%
    dplyr::pull(gene) %>% unique() %>% intersect(used_gene_set)

  pleio_E13.5 <- tau_clean %>%
    dplyr::filter(stage == "E13.5", tau <= threshold) %>%
    dplyr::pull(gene) %>% unique() %>% intersect(used_gene_set)

  df <- rbind(
    pleiotropic_fraction(E9.5_raw,  used_gene_set, pleio_E9.5,  "E9.5"),
    pleiotropic_fraction(E13.5_raw, used_gene_set, pleio_E13.5, "E13.5")
  )
  df$cluster <- factor(df$cluster, levels = CLUSTER_LEVELS)
  df$stage   <- factor(df$stage,   levels = c("E9.5", "E13.5"))
  df
}

##################################################
#  Shared ggplot theme for the comparison plots  #
##################################################

comparison_theme <- function() {
  theme_minimal() +
    theme(
      text            = element_text(size = 14, family = "Arial"),
      plot.title      = element_text(hjust = 0.5, size = 8),
      plot.subtitle   = element_text(hjust = 0.5, size = 8),
      panel.grid      = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      axis.text       = element_text(color = "black"),
      axis.line       = element_line(color = "black")
    )
}


#  Plot: (A) shared / paired -> paired dotplot with labels

plot_paired_dotplot_label <- function(df, threshold, out_dir) {

  # arrange data for a paired comparison (clusters present in both stages)
  df_paired <- df %>%
    dplyr::select(stage, cluster, pleiotropic_perc) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(names_from = stage, values_from = pleiotropic_perc) %>%
    tidyr::drop_na(`E9.5`, `E13.5`) %>%
    tidyr::pivot_longer(cols = c(`E9.5`, `E13.5`),
                        names_to = "stage", values_to = "pleiotropic_perc") %>%
    mutate(stage = factor(stage, levels = c("E9.5", "E13.5"))) %>%
    arrange(cluster, stage)

  # Wilcoxon signed-rank (paired) + effect size
  paired_stat.test <- df_paired %>%
    wilcox_test(pleiotropic_perc ~ stage, paired = TRUE) %>%
    add_significance()
  paired_eff <- df_paired %>%
    wilcox_effsize(pleiotropic_perc ~ stage, paired = TRUE) %>%
    pull(effsize)

  paired_subtitle <- bquote(italic(p) == .(signif(paired_stat.test$p, 5)) ~ ", " ~
                              italic(r) == .(round(paired_eff, 3)))

  paired_stat.test <- paired_stat.test %>% add_xy_position(x = "stage")
  paired_stat.test$y.position <- max(df_paired$pleiotropic_perc, na.rm = TRUE) * 1.08  # manual fix

  p <- ggplot(df_paired, aes(x = stage, y = pleiotropic_perc, group = cluster)) +
    geom_line(color = "grey30", linewidth = 0.6, alpha = 0.8) +
    geom_point(aes(color = stage), size = 2.5, show.legend = FALSE) +
    # E9.5 labels (left)
    geom_text_repel(data = subset(df_paired, stage == "E9.5"),
                    aes(label = cluster), nudge_x = -0.25, direction = "y",
                    hjust = 0.5, size = 3, color = "black", segment.color = "grey80") +
    # E13.5 labels (right)
    geom_text_repel(data = subset(df_paired, stage == "E13.5"),
                    aes(label = cluster), nudge_x = 0.25, direction = "y",
                    hjust = 0.5, size = 3, color = "black", segment.color = "grey80") +
    scale_color_manual(values = STAGE_COLORS) +
    stat_pvalue_manual(paired_stat.test, tip.length = 0) +
    labs(x = "Stage",
         y = "Fraction of pleiotropic genes per Oi-clusters",
         title = paste0("Distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ") (paired)"),
         subtitle = paired_subtitle) +
    comparison_theme()

  ggsave(
    filename = paste0(out_dir,
                      "E9.5+E13.5_pseudo.bulk_distribution_of_cl_percentage_of_pleio_genes_Tau.threshold=",
                      threshold, "_dotplot_paired_wilcox_label.pdf"),
    p, width = 5, height = 7
  )
  invisible(NULL)
}


#  Plot: (B) all / unpaired -> unpaired boxplot with labels

plot_unpaired_boxplot_label <- function(df, threshold, out_dir) {

  # Mann-Whitney U (unpaired) + effect size
  unpaired_stat.test <- df %>%
    wilcox_test(pleiotropic_perc ~ stage) %>%
    add_significance()
  unpaired_eff <- df %>%
    wilcox_effsize(pleiotropic_perc ~ stage) %>%
    pull(effsize)

  unpaired_subtitle <- bquote(italic(p) == .(signif(unpaired_stat.test$p, 5)) ~ ", " ~
                                italic(r) == .(round(unpaired_eff, 3)))

  unpaired_stat.test <- unpaired_stat.test %>% add_xy_position(x = "stage")
  unpaired_stat.test$y.position <- max(df$pleiotropic_perc, na.rm = TRUE) * 1.08  # manual fix

  p <- ggplot(df, aes(x = stage, y = pleiotropic_perc, fill = stage)) +
    geom_boxplot(outlier.shape = NA, width = 0.37) +
    geom_jitter(color = "black", width = 0, size = 1.8, alpha = 0.8, show.legend = FALSE) +
    # E9.5 labels (left)
    geom_text_repel(data = subset(df, stage == "E9.5"),
                    aes(label = cluster), nudge_x = -0.25, direction = "y",
                    hjust = 0.5, size = 4, color = "black", segment.color = "grey80",
                    max.overlaps = Inf, force = 2, force_pull = 0.5,
                    box.padding = 0.6, point.padding = 0.4) +
    # E13.5 labels (right)
    geom_text_repel(data = subset(df, stage == "E13.5"),
                    aes(label = cluster), nudge_x = 0.25, direction = "y",
                    hjust = 0.5, size = 4, color = "black", segment.color = "grey80",
                    max.overlaps = Inf, force = 2, force_pull = 0.5,
                    box.padding = 0.6, point.padding = 0.4) +
    scale_fill_manual(values = STAGE_COLORS) +
    stat_pvalue_manual(unpaired_stat.test, tip.length = 0) +
    labs(x = "Stage",
         y = "Fraction of pleiotropic genes per Oi-clusters",
         title = paste0("Distribution of percentage of pleiotropic genes (Tau threshold=", threshold, ") (Unpaired)"),
         subtitle = unpaired_subtitle) +
    comparison_theme()

  ggsave(
    filename = paste0(out_dir,
                      "E9.5+E13.5_pseudo.bulk_distribution_of_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                      threshold, "_boxplot_unpaired_wilcox_with.label.pdf"),
    p, width = 5, height = 7
  )
  invisible(NULL)
}


##########################################################
#  Driver: run one analysis scope across all thresholds  #
##########################################################

#   scope_subdir : sub-dir name ("common.genes_shared.clusters" / "...all.clusters")
#   tau_clean    : clean long Tau table for the scope
#   E9.5_raw / E13.5_raw : pseudo-bulk matrices restricted to the scope's clusters
#   plot_fun     : plotting function for the scope (paired or unpaired variant)

run_analysis_scope <- function(scope_subdir, tau_clean, E9.5_raw, E13.5_raw, plot_fun) {

  for (threshold in pleiotropy_thresholds) {

    this.E9.5_output.dir     <- paste0(E9.5_output.dir,     "tau_pleiotropy.threshold=", threshold, "/", scope_subdir, "/")
    this.E13.5_output.dir    <- paste0(E13.5_output.dir,    "tau_pleiotropy.threshold=", threshold, "/", scope_subdir, "/")
    this.combined_output.dir <- paste0(combined_output.dir, "tau_pleiotropy.threshold=", threshold, "/", scope_subdir, "/")
    for (d in c(this.E9.5_output.dir, this.E13.5_output.dir, this.combined_output.dir)) {
      if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    }

    # per-cluster pleiotropic fractions, both stages
    combined_pleiotropic_perc_df <- compute_pleiotropic_perc(
      tau_clean = tau_clean,
      E9.5_raw  = E9.5_raw,
      E13.5_raw = E13.5_raw,
      threshold = threshold
    )

    # save the per-cluster fractions
    write.csv(
      combined_pleiotropic_perc_df,
      file = paste0(this.combined_output.dir,
                    "E9.5+E13.5_pseudo.bulk_cluster_percentage_of_pleiotropic_genes_Tau.threshold=",
                    threshold, "_combined_pleiotropic_perc_df.csv"),
      row.names = FALSE
    )

    # the single figure for this scope (with the y.position fix applied)
    plot_fun(combined_pleiotropic_perc_df, threshold, this.combined_output.dir)

    message("Done: ", this.combined_output.dir)
  }
  invisible(NULL)
}


###############################################################################
#                                                                             #
#  Run both analysis scopes across all pleiotropy thresholds                  #
#                                                                             #
###############################################################################

# (A) common genes / shared Oi-clusters -> paired Wilcoxon
# E13.5 restricted to shared (non-E13.5-specific) clusters
E13.5_pseudo_raw_shared <- E13.5_pseudo_raw[, !colnames(E13.5_pseudo_raw) %in% E13.5_SPECIFIC, drop = FALSE]

run_analysis_scope(
  scope_subdir = "common.genes_shared.clusters",
  tau_clean    = tau_combined_df_clean_shared_filt,
  E9.5_raw     = E9.5_pseudo_raw,
  E13.5_raw    = E13.5_pseudo_raw_shared,
  plot_fun     = plot_paired_dotplot_label
)

# (B) common genes / all Oi-clusters -> unpaired Wilcoxon
run_analysis_scope(
  scope_subdir = "common.genes_all.clusters",
  tau_clean    = tau_combined_df_clean,
  E9.5_raw     = E9.5_pseudo_raw,
  E13.5_raw    = E13.5_pseudo_raw,
  plot_fun     = plot_unpaired_boxplot_label
)
