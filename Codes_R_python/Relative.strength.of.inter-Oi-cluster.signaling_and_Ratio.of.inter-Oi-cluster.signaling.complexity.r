library(Seurat)
library(sctransform)
library(anndata)
library(reticulate)
library(SeuratDisk)

library(CellChat)
library(patchwork)
library(BiocNeighbors)
options(stringsAsFactors = FALSE)

library(future)
library(glue)
library(dplyr)
library(forcats)
library(tidyr)
library(tidyverse)

library(rstatix)
library(ggpubr)

library(ggrepel)
library(ggplot2)
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
E9.5_colors.use <- colors.use_pre[!names(colors.use_pre) %in% c("27", "31")]
# for E13.5
E13.5_colors.use <- colors.use_pre



#############################################
#                                           #
#  Read-in CellChat results of both stages  #
#                                           #
#############################################

# Load CellChat and CellChat_netP data

# paths to CellChat results
seurat_res <- "0.7"
trim_para <- 0.05
# E9.5
E9.5_output.dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/E9.5_E1S1/population.size=TRUE_trim=0.05/"
# E13.5
E13.5_output.dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/E13.5_E1S1/population.size=TRUE_trim=0.05/"

# load E9.5_E1S1_bin50 (SeuratCCA_integrated, Leiden-resolution=0.7, trim=0.05)
E9.5_bin50.cellchat <- readRDS(file.path(E9.5_output.dir, "E9.5_E1S1_bin50_seuratCCA_clusters_res=0.7_after_aggregateNet_population.size=TRUE_trim=0.05.RDS"))
E9.5_df_netP <- read.csv2(file.path(E9.5_output.dir, "E9.5_E1S1_bin50_seuratCCA_clusters_res=0.7_after_aggregateNet_population.size=TRUE_trim=0.05_df.netP.CSV"))

E9.5_seurat_clusters <- as.character(unique(E9.5_bin50.cellchat@idents))


# load E13.5_E1S1_bin50 (SeuratCCA_integrated, Leiden-resolution=0.7, trim=0.05)
E13.5_bin50.cellchat <- readRDS(file.path(E13.5_output.dir, "E13.5_E1S1_bin50_seuratCCA_clusters_res=0.7_after_aggregateNet_population.size=TRUE_trim=0.05.RDS"))
E13.5_df_netP <- read.csv2(file.path(E13.5_output.dir, "E13.5_E1S1_bin50_seuratCCA_clusters_res=0.7_after_aggregateNet_population.size=TRUE_trim=0.05_df.netP.CSV"))

E13.5_seurat_clusters <- as.character(unique(E13.5_bin50.cellchat@idents))

# Set new output.dir for inter-Oi-cluster signaling analysis results
# E9.5
E9.5_output.dir <- "/home/project_interconnectivity/result/clusterfilter_inter_signaling_results_onesided/E9.5_E1S1/"
# E13.5
E13.5_output.dir <- "/home/project_interconnectivity/result/clusterfilter_inter_signaling_results_onesided/E13.5_E1S1/"

if (!dir.exists(E9.5_output.dir)) {
  dir.create(E9.5_output.dir, recursive = TRUE)
}

if (!dir.exists(E13.5_output.dir)) {
  dir.create(E13.5_output.dir, recursive = TRUE)
}



#################################
#                               #
#  Obtain signaling properties  #
#                               #
#################################

### E9.5

E9.5_signaling_list <- lapply(E9.5_seurat_clusters, function(cl) {
  
  # incoming
  incoming_all   <- subset(E9.5_df_netP,  target == cl)
  incoming_woSelf <- subset(incoming_all, source != cl)
  
  # outgoing
  outgoing_all   <- subset(E9.5_df_netP,  source == cl)
  outgoing_woSelf <- subset(outgoing_all, target != cl)
  
  # all_related
  all_all <- subset(E9.5_df_netP, (source == cl) | (target == cl))
  
  # inter-cluster signaling
  all_woSelf <- subset(all_all, !(source == cl & target == cl))
  
  # self_signaling
  Self <- subset(all_all, (source == cl & target == cl))
  
  # calculate ratio value
  ratio_value <- ifelse(nrow(all_all) == 0,
                        NA,
                        nrow(all_woSelf) / nrow(all_all))
  
  data.frame(
    cluster                    = cl,
    incoming_signaling_wSelf   = nrow(incoming_all),
    outgoing_signaling_wSelf   = nrow(outgoing_all),
    incoming_signaling_woSelf  = nrow(incoming_woSelf),
    outgoing_signaling_woSelf  = nrow(outgoing_woSelf),
    all_signaling_wSelf        = nrow(all_all),
    all_signaling_woSelf       = nrow(all_woSelf),
    self_signaling             = nrow(Self),
    ratio                      = ratio_value,
    stringsAsFactors = FALSE
  )
})

E9.5_signaling_complexity <- do.call(rbind, E9.5_signaling_list)


### E13.5

E13.5_signaling_list <- lapply(E13.5_seurat_clusters, function(cl) {
  
  # incoming
  incoming_all   <- subset(E13.5_df_netP,  target == cl)
  incoming_woSelf <- subset(incoming_all, source != cl)
  
  # outgoing
  outgoing_all   <- subset(E13.5_df_netP,  source == cl)
  outgoing_woSelf <- subset(outgoing_all, target != cl)
  
  # all_related
  all_all <- subset(E13.5_df_netP, (source == cl) | (target == cl))
  
  # inter-cluster signaling
  all_woSelf <- subset(all_all, !(source == cl & target == cl))
  
  # self_signaling
  Self <- subset(all_all, (source == cl & target == cl))
  
  # calculate ratio value
  ratio_value <- ifelse(nrow(all_all) == 0,
                        NA,
                        nrow(all_woSelf) / nrow(all_all))
  
  data.frame(
    cluster                    = cl,
    incoming_signaling_wSelf   = nrow(incoming_all),
    outgoing_signaling_wSelf   = nrow(outgoing_all),
    incoming_signaling_woSelf  = nrow(incoming_woSelf),
    outgoing_signaling_woSelf  = nrow(outgoing_woSelf),
    all_signaling_wSelf        = nrow(all_all),
    all_signaling_woSelf       = nrow(all_woSelf),
    self_signaling             = nrow(Self),
    ratio                      = ratio_value,
    stringsAsFactors = FALSE
  )
})

E13.5_signaling_complexity <- do.call(rbind, E13.5_signaling_list)


check_df <- E13.5_signaling_complexity %>%
  mutate(
    sum_components =
      incoming_signaling_woSelf +
      outgoing_signaling_woSelf +
      self_signaling,
    logic_ok = (sum_components == all_signaling_wSelf),
    diff = sum_components - all_signaling_wSelf
  )

check_df
all(check_df$logic_ok)

# create a new directory for saving the comparison results
save.path <- "/home/project_interconnectivity/result/clusterfilter_inter_signaling_results_onesided/E9.5+E13.5_comparison/"
if (!dir.exists(save.path)) {
  dir.create(save.path)
}



####################################################################################################
#                                                                                                  #
#  Relative strength of inter-Oi-cluster signaling: Cross-stage comparison + statistical analysis  #
#                                                                                                  #
####################################################################################################

### E9.5

# remove Oi-clusters 27, 29, 30, 31 from E9.5_seurat_clusters in advance
E9.5_seurat_clusters <- E9.5_seurat_clusters[!E9.5_seurat_clusters %in% c("27", "29", "30", "31")]

E9.5_ratio_value_list <- lapply(E9.5_seurat_clusters, function(cl) {
  # cluster specific signaling
  cluster_subset <- E9.5_df_netP %>%
    filter(source == cl | target == cl)
  
  # intra-cluster signaling
  intra_subset <- E9.5_df_netP %>%
    filter(source == cl & target == cl)
  
  # inter-cluster signaling
  inter_subset <- E9.5_df_netP %>%
    filter((source == cl | target == cl) & source != target)
  
  # calculate ratio value
  ratio_value <- ifelse(sum(cluster_subset$prob) == 0,
                        NA,
                        sum(inter_subset$prob) / sum(cluster_subset$prob))
  
  data.frame(
    cluster = cl,
    E9.5_intra.cluster_signaling = sum(intra_subset$prob),
    E9.5_inter.cluster_signaling = sum(inter_subset$prob),
    E9.5_all_signaling = sum(cluster_subset$prob),
    E9.5_ratio = ratio_value
  )
})

E9.5_signaling_ratio_value <- do.call(rbind, E9.5_ratio_value_list)


### E13.5

E13.5_ratio_value_list <- lapply(E13.5_seurat_clusters, function(cl) {
  # cluster specific signaling
  cluster_subset <- E13.5_df_netP %>%
    filter(source == cl | target == cl)
  
  # intra-cluster signaling
  intra_subset <- E13.5_df_netP %>%
    filter(source == cl & target == cl)
  
  # inter-cluster signaling
  inter_subset <- E13.5_df_netP %>%
    filter((source == cl | target == cl) & source != target)
  
  # calculate ratio value
  ratio_value <- ifelse(sum(cluster_subset$prob) == 0,
                        NA,
                        sum(inter_subset$prob) / sum(cluster_subset$prob))
  
  data.frame(
    cluster = cl,
    E13.5_intra.cluster_signaling = sum(intra_subset$prob),
    E13.5_inter.cluster_signaling = sum(inter_subset$prob),
    E13.5_all_signaling = sum(cluster_subset$prob),
    E13.5_ratio = ratio_value
  )
})

E13.5_signaling_ratio_value <- do.call(rbind, E13.5_ratio_value_list)


# combine E9.5 and E13.5 data frames for inter-cluster signaling strength and ratio values
combined_cluster_signaling_strength <- left_join(E13.5_signaling_ratio_value, E9.5_signaling_ratio_value, by = "cluster")
combined_cluster_signaling_strength$cluster <- factor(combined_cluster_signaling_strength$cluster, levels = as.character(c(1:31)))

# save temporary results as CSV
save.relative.strength <- paste0(save.path, "Relative_strength_of_inter-Oi-cluster_signaling_cross-stage_comparison/")
if (!dir.exists(save.relative.strength)) {
  dir.create(save.relative.strength)
}

write.csv(
  combined_cluster_signaling_strength,
  file = paste0(save.relative.strength, "Combined_E9.5+E13.5_Oi-cluster_relative_strength_of_inter-Oi-cluster_signaling.CSV"),
  row.names = FALSE
)


###################################################
#  Cross-stage comparison + statistical analysis  #
###################################################

# create a new directory for saving the cross-stage comparison results
combined_output.dir <- paste0(save.relative.strength, "Cross-stage_comparison_statistical.analysis/")
if (!dir.exists(combined_output.dir)) {
  dir.create(combined_output.dir, recursive = TRUE)
}

any(E9.5_df_netP$source %in% c("27", "29", "30", "31")) |
  any(E9.5_df_netP$target %in% c("27", "29", "30", "31"))

# reshape the combined data frame to long format
combined_cluster_relative_strength_long <- combined_cluster_signaling_strength %>%
  pivot_longer(
    cols = c("E9.5_ratio", "E13.5_ratio"),
    names_to = "stage",
    values_to = "cluster_inter_strength_ratio"
  )

combined_cluster_relative_strength_long$stage <- sub(
  "_ratio$",
  "",
  combined_cluster_relative_strength_long$stage
)

combined_cluster_relative_strength_long$stage <- factor(combined_cluster_relative_strength_long$stage, levels = c("E9.5", "E13.5"))


### Wilcoxon Signed-Ranks Test to compare the relative strength values between E9.5 and E13.5 (paired)

# clear data from Oi-clusters 27, 29, 30, 31 in both E9.5 and E13.5
combined_cluster_relative_strength_long_paired_filt <- 
  combined_cluster_relative_strength_long %>%
  subset(subset = !cluster %in% c("27", "29", "30", "31"))

combined_cluster_relative_strength_long_paired_filt <-
  combined_cluster_relative_strength_long_paired_filt %>%
  arrange(cluster, stage)

combined_cluster_relative_strength_long_paired_filt$stage <-
  factor(combined_cluster_relative_strength_long_paired_filt$stage,
         levels = c("E9.5", "E13.5"))

combined_cluster_relative_strength_long_paired_filt

# perform Wilcoxon Signed-Ranks Test (Paired, one-sided: H0: E9.5 <= E13.5, H1: E9.5 > E13.5)
paired_stat.test <- combined_cluster_relative_strength_long_paired_filt %>%
  wilcox_test(cluster_inter_strength_ratio ~ stage, paired = TRUE, alternative = "greater") %>%
  add_significance()

paired_stat.test

# compute effect size
wilcox_effsize(combined_cluster_relative_strength_long_paired_filt, cluster_inter_strength_ratio ~ stage, paired = TRUE)

paired_effect.size <- combined_cluster_relative_strength_long_paired_filt %>%
  wilcox_effsize(cluster_inter_strength_ratio ~ stage, paired = TRUE) %>%
  pull(effsize)

# visualization -- dotplot with labels

# construct subtitle for plotting
paired_r.val <- round(paired_effect.size, 3)
paired_p.label <- get_test_label(paired_stat.test, detailed = TRUE)
paired_subtitle_expr <- bquote(.(paired_p.label) * ", " * italic(r) == .(paired_r.val))

paired_stat.test <- paired_stat.test %>%
  add_xy_position(x = "stage")

# dotplot with labels
combined_paired_test_dotplot_with.label <- ggplot(
  combined_cluster_relative_strength_long_paired_filt,
  aes(
    x = stage,
    y = cluster_inter_strength_ratio,
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
      combined_cluster_relative_strength_long_paired_filt,
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
      combined_cluster_relative_strength_long_paired_filt,
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
  scale_color_manual(
    values = c("#eb4d00", "#87a8eb")
  ) +
  stat_pvalue_manual(
    paired_stat.test, tip.length = 0
  ) +
  labs(
    x = "Stage",
    y = "Ratio",
    title = "Relative strength of inter-Oi-cluster signaling (E9.5 vs E13.5)_paired",
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

combined_paired_test_dotplot_with.label

ggsave(
  filename = paste0(combined_output.dir, "E9.5+E13.5_Oi-cluster_relative_strength_of_inter-Oi-cluster_signaling_paired_wilcox_dotplot.with.label.pdf"),
  combined_paired_test_dotplot_with.label,
  width = 5,
  height = 7
)


### Mann-Whitney U Test (Wilcoxon Rank Sum Test) to compare the relative strength values between E9.5 and E13.5 (unpaired)

# clear the NA in the original data
combined_cluster_relative_strength_long_filt <- 
  combined_cluster_relative_strength_long %>%
  dplyr::filter(!is.na(cluster_inter_strength_ratio))

combined_cluster_relative_strength_long_filt$stage <- factor(combined_cluster_relative_strength_long_filt$stage, levels = c("E9.5", "E13.5"))

# perform Mann-Whitney U Test (one-sided: H0: E9.5 <= E13.5, H1: E9.5 > E13.5)
unpaired_stat.test <- combined_cluster_relative_strength_long_filt %>%
  wilcox_test(cluster_inter_strength_ratio ~ stage, alternative = "greater") %>%
  add_significance()

unpaired_stat.test

# compute effect size
wilcox_effsize(combined_cluster_relative_strength_long_filt, cluster_inter_strength_ratio ~ stage)

unpaired_effect.size <- combined_cluster_relative_strength_long_filt %>%
  wilcox_effsize(cluster_inter_strength_ratio ~ stage) %>%
  pull(effsize)

# visualization -- boxplot with labels

# construct subtitle for plotting
unpaired_r.val <- round(unpaired_effect.size, 3)
unpaired_p.label <- get_test_label(unpaired_stat.test, detailed = TRUE)
unpaired_subtitle_expr <- bquote(.(unpaired_p.label) * ", " * italic(r) == .(unpaired_r.val))

unpaired_stat.test <- unpaired_stat.test %>%
  add_xy_position(x = "stage")

# boxplot with labels
combined_unpaired_test_boxplot_with.label <- ggplot(
  combined_cluster_relative_strength_long_filt,
  aes(
    x = stage,
    y = cluster_inter_strength_ratio,
    fill = stage
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.37
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
      combined_cluster_relative_strength_long_filt,
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
      combined_cluster_relative_strength_long_filt,
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
    y = "Ratio",
    title = "Relative strength of inter-Oi-cluster signaling (E9.5 vs E13.5)_unpaired",
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

combined_unpaired_test_boxplot_with.label

ggsave(
  filename = paste0(combined_output.dir, "E9.5+E13.5_Oi-cluster_relative_strength_of_inter-Oi-cluster_signaling_unpaired_wilcox_boxplot.with.label.pdf"),
  combined_unpaired_test_boxplot_with.label,
  width = 5,
  height = 7
)



###################################################################################################
#                                                                                                 #
#  Ratio of inter–Oi-cluster signaling complexity: Cross-stage comparison + statistical analysis  #
#                                                                                                 #
###################################################################################################

### E9.5

# remove Oi-clusters 27, 29, 30, 31 from E9.5_seurat_clusters in advance
E9.5_seurat_clusters <- E9.5_seurat_clusters[!E9.5_seurat_clusters %in% c("27", "29", "30", "31")]

E9.5_pw_ratio_list <- lapply(E9.5_seurat_clusters, function(cl) {
  
  # incoming
  incoming_all   <- subset(E9.5_df_netP,  target == cl)
  incoming_woSelf <- subset(incoming_all, source != cl)
  
  # outgoing
  outgoing_all   <- subset(E9.5_df_netP,  source == cl)
  outgoing_woSelf <- subset(outgoing_all, target != cl)
  
  # all_related
  all_all <- subset(E9.5_df_netP, (source == cl) | (target == cl))
  
  # inter-cluster signaling
  all_woSelf <- subset(all_all, !(source == cl & target == cl))
  
  # self_signaling
  Self <- subset(all_all, (source == cl & target == cl))
  
  # calculate ratio value
  ratio_value <- ifelse(nrow(all_all) == 0,
                        NA,
                        nrow(all_woSelf) / nrow(all_all))
  
  data.frame(
    cluster                    = cl,
    # E9.5_incoming_signaling_wSelf   = nrow(incoming_all),
    # E9.5_outgoing_signaling_wSelf   = nrow(outgoing_all),
    # E9.5_incoming_signaling_woSelf  = nrow(incoming_woSelf),
    # E9.5_outgoing_signaling_woSelf  = nrow(outgoing_woSelf),
    # E9.5_all_signaling_wSelf        = nrow(all_all),
    # E9.5_all_signaling_woSelf       = nrow(all_woSelf),
    # E9.5_self_signaling             = nrow(Self),
    E9.5_ratio                      = ratio_value,
    stringsAsFactors = FALSE
  )
})

E9.5_pw_ratio <- do.call(rbind, E9.5_pw_ratio_list)

E9.5_pw_all_list <- lapply(E9.5_seurat_clusters, function(cl) {
  
  # incoming
  incoming_all   <- subset(E9.5_df_netP,  target == cl)
  incoming_woSelf <- subset(incoming_all, source != cl)
  
  # outgoing
  outgoing_all   <- subset(E9.5_df_netP,  source == cl)
  outgoing_woSelf <- subset(outgoing_all, target != cl)
  
  # all_related
  all_all <- subset(E9.5_df_netP, (source == cl) | (target == cl))
  
  # inter-cluster signaling
  all_woSelf <- subset(all_all, !(source == cl & target == cl))
  
  # self_signaling
  Self <- subset(all_all, (source == cl & target == cl))
  
  # calculate ratio value
  ratio_value <- ifelse(nrow(all_all) == 0,
                        NA,
                        nrow(all_woSelf) / nrow(all_all))
  
  data.frame(
    cluster                    = cl,
    # E9.5_incoming_signaling_wSelf   = nrow(incoming_all),
    # E9.5_outgoing_signaling_wSelf   = nrow(outgoing_all),
    # E9.5_incoming_signaling_woSelf  = nrow(incoming_woSelf),
    # E9.5_outgoing_signaling_woSelf  = nrow(outgoing_woSelf),
    E9.5_all_signaling_wSelf        = nrow(all_all),
    E9.5_all_signaling_woSelf       = nrow(all_woSelf),
    E9.5_self_signaling             = nrow(Self),
    E9.5_ratio                      = ratio_value,
    stringsAsFactors = FALSE
  )
})

E9.5_pw_all <- do.call(rbind, E9.5_pw_all_list)


### E13.5

E13.5_pw_ratio_list <- lapply(E13.5_seurat_clusters, function(cl) {
  
  # incoming
  incoming_all   <- subset(E13.5_df_netP,  target == cl)
  incoming_woSelf <- subset(incoming_all, source != cl)
  
  # outgoing
  outgoing_all   <- subset(E13.5_df_netP,  source == cl)
  outgoing_woSelf <- subset(outgoing_all, target != cl)
  
  # all_related
  all_all <- subset(E13.5_df_netP, (source == cl) | (target == cl))
  
  # inter-cluster signaling
  all_woSelf <- subset(all_all, !(source == cl & target == cl))
  
  # self_signaling
  Self <- subset(all_all, (source == cl & target == cl))
  
  # calculate ratio value
  ratio_value <- ifelse(nrow(all_all) == 0,
                        NA,
                        nrow(all_woSelf) / nrow(all_all))
  
  data.frame(
    cluster                    = cl,
    # E13.5_incoming_signaling_wSelf   = nrow(incoming_all),
    # E13.5_outgoing_signaling_wSelf   = nrow(outgoing_all),
    # E13.5_incoming_signaling_woSelf  = nrow(incoming_woSelf),
    # E13.5_outgoing_signaling_woSelf  = nrow(outgoing_woSelf),
    # E13.5_all_signaling_wSelf        = nrow(all_all),
    # E13.5_all_signaling_woSelf       = nrow(all_woSelf),
    # E13.5_self_signaling             = nrow(Self),
    E13.5_ratio                      = ratio_value,
    stringsAsFactors = FALSE
  )
})

E13.5_pw_ratio <- do.call(rbind, E13.5_pw_ratio_list)


E13.5_pw_all_list <- lapply(E13.5_seurat_clusters, function(cl) {
  
  # incoming
  incoming_all   <- subset(E13.5_df_netP,  target == cl)
  incoming_woSelf <- subset(incoming_all, source != cl)
  
  # outgoing
  outgoing_all   <- subset(E13.5_df_netP,  source == cl)
  outgoing_woSelf <- subset(outgoing_all, target != cl)
  
  # all_related
  all_all <- subset(E13.5_df_netP, (source == cl) | (target == cl))
  
  # inter-cluster signaling
  all_woSelf <- subset(all_all, !(source == cl & target == cl))
  
  # self_signaling
  Self <- subset(all_all, (source == cl & target == cl))
  
  # calculate ratio value
  ratio_value <- ifelse(nrow(all_all) == 0,
                        NA,
                        nrow(all_woSelf) / nrow(all_all))
  
  data.frame(
    cluster                    = cl,
    # E13.5_incoming_signaling_wSelf   = nrow(incoming_all),
    # E13.5_outgoing_signaling_wSelf   = nrow(outgoing_all),
    # E13.5_incoming_signaling_woSelf  = nrow(incoming_woSelf),
    # E13.5_outgoing_signaling_woSelf  = nrow(outgoing_woSelf),
    E13.5_all_signaling_wSelf        = nrow(all_all),
    E13.5_all_signaling_woSelf       = nrow(all_woSelf),
    E13.5_self_signaling             = nrow(Self),
    E13.5_ratio                      = ratio_value,
    stringsAsFactors = FALSE
  )
})

E13.5_pw_all <- do.call(rbind, E13.5_pw_all_list)


# combine E9.5 and E13.5 data frames for ratio and all signaling properties
combined_cluster_pw_ratio <- left_join(E13.5_pw_ratio, E9.5_pw_ratio, by = "cluster")
combined_cluster_pw_ratio$cluster <- factor(combined_cluster_pw_ratio$cluster, levels = as.character(c(1:31)))

combined_cluster_pw_all <- left_join(E13.5_pw_all, E9.5_pw_all, by = "cluster")
combined_cluster_pw_all$cluster <- factor(combined_cluster_pw_all$cluster, levels = as.character(c(1:31)))

# save temporary results as CSV
save.path.pathway.number <- paste0(save.path, "Ratio_of_inter-Oi-cluster_signaling_complexity_cross-stage_comparison/")
if (!dir.exists(save.path.pathway.number)) {
  dir.create(save.path.pathway.number)
}

write.csv(
  combined_cluster_pw_ratio,
  file = paste0(save.path.pathway.number, "Combined_E9.5+E13.5_Oi-cluster_inter-Oi-cluster_signaling_complexity_ratio.CSV"),
  row.names = FALSE
)

write.csv(
  combined_cluster_pw_all,
  file = paste0(save.path.pathway.number, "Combined_E9.5+E13.5_Oi-cluster_inter-Oi-cluster_signaling_complexity_all.CSV"),
  row.names = FALSE
)


###################################################
#  Cross-stage comparison + statistical analysis  #
###################################################

# create a new directory for saving the cross-stage comparison results
combined_output.dir <- paste0(save.path.pathway.number, "Cross-stage_comparison_statistical.analysis/")
if (!dir.exists(combined_output.dir)) {
  dir.create(combined_output.dir)
}

any(E9.5_df_netP$source %in% c("27", "29", "30", "31")) |
  any(E9.5_df_netP$target %in% c("27", "29", "30", "31"))

# reshape the combined data frame to long format
combined_cluster_pw_ratio_long <- combined_cluster_pw_ratio %>%
  pivot_longer(
    cols = c("E9.5_ratio", "E13.5_ratio"),
    names_to = "stage",
    values_to = "cluster_pathway_ratio"
  )

combined_cluster_pw_ratio_long$stage <- sub(
  "_ratio$",
  "",
  combined_cluster_pw_ratio_long$stage
)

combined_cluster_pw_ratio_long$stage <- factor(combined_cluster_pw_ratio_long$stage, levels = c("E9.5", "E13.5"))


### Wilcoxon Signed-Ranks Test to compare the ratio values between E9.5 and E13.5 (paired)

# clear data from Oi-clusters 27, 29, 30, 31 in both E9.5 and E13.5
combined_cluster_pw_ratio_long_paired_filt <- 
  combined_cluster_pw_ratio_long %>%
  subset(subset = !cluster %in% c("27", "29", "30", "31"))

combined_cluster_pw_ratio_long_paired_filt <- 
  combined_cluster_pw_ratio_long_paired_filt %>%
  arrange(cluster, stage)

combined_cluster_pw_ratio_long_paired_filt$stage <-
  factor(combined_cluster_pw_ratio_long_paired_filt$stage,
         levels = c("E9.5", "E13.5"))

combined_cluster_pw_ratio_long_paired_filt

# perform Wilcoxon Signed-Ranks Test (Paired, one-sided: H0: E9.5 <= E13.5, H1: E9.5 > E13.5)
paired_stat.test <- combined_cluster_pw_ratio_long_paired_filt %>%
  wilcox_test(cluster_pathway_ratio ~ stage, paired = TRUE, alternative = "greater") %>%
  add_significance()

paired_stat.test

# compute effect size
wilcox_effsize(combined_cluster_pw_ratio_long_paired_filt, cluster_pathway_ratio ~ stage, paired = TRUE)

paired_effect.size <- combined_cluster_pw_ratio_long_paired_filt %>%
  wilcox_effsize(cluster_pathway_ratio ~ stage, paired = TRUE) %>%
  pull(effsize)

# visualization -- dotplot with labels

# construct subtitle for plotting
paired_r.val <- round(paired_effect.size, 3)
paired_p.label <- get_test_label(paired_stat.test, detailed = TRUE)
paired_subtitle_expr <- bquote(.(paired_p.label) * ", " * italic(r) == .(paired_r.val))

paired_stat.test <- paired_stat.test %>%
  add_xy_position(x = "stage")

# dotplot with labels
combined_paired_test_dotplot_with.label <- ggplot(
  combined_cluster_pw_ratio_long_paired_filt,
  aes(
    x = stage,
    y = cluster_pathway_ratio,
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
      combined_cluster_pw_ratio_long_paired_filt,
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
      combined_cluster_pw_ratio_long_paired_filt,
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
  scale_color_manual(
    values = c("#eb4d00", "#87a8eb")
  ) +
  stat_pvalue_manual(
    paired_stat.test, tip.length = 0
  ) +
  labs(
    x = "Stage",
    y = "Ratio",
    title = "Ratio of inter-Oi-cluster signaling complexity (E9.5 vs E13.5)_paired",
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

combined_paired_test_dotplot_with.label

ggsave(
  filename = paste0(combined_output.dir, "E9.5+E13.5_Oi-cluster_inter-Oi-cluster_signaling_complexity_ratio_paired_wilcox_dotplot.with.label.pdf"),
  combined_paired_test_dotplot_with.label,
  width = 5,
  height = 7
)


### Mann-Whitney U Test (Wilcoxon Rank Sum Test) to compare the ratio values between E9.5 and E13.5 (unpaired)

# clear the NA in the original data
combined_cluster_pw_ratio_long_filt <- 
  combined_cluster_pw_ratio_long %>%
  dplyr::filter(!is.na(cluster_pathway_ratio))

combined_cluster_pw_ratio_long_filt$stage <- factor(combined_cluster_pw_ratio_long_filt$stage, levels = c("E9.5", "E13.5"))

# perform Mann-Whitney U Test (one-sided: H0: E9.5 <= E13.5, H1: E9.5 > E13.5)
unpaired_stat.test <- combined_cluster_pw_ratio_long_filt %>%
  wilcox_test(cluster_pathway_ratio ~ stage, alternative = "greater") %>%
  add_significance()

unpaired_stat.test

# compute effect size
wilcox_effsize(combined_cluster_pw_ratio_long_filt, cluster_pathway_ratio ~ stage)

unpaired_effect.size <- combined_cluster_pw_ratio_long_filt %>%
  wilcox_effsize(cluster_pathway_ratio ~ stage) %>%
  pull(effsize)

# visualization -- boxplot with labels

# construct subtitle for plotting
unpaired_r.val <- round(unpaired_effect.size, 3)
unpaired_p.label <- get_test_label(unpaired_stat.test, detailed = TRUE)
unpaired_subtitle_expr <- bquote(.(unpaired_p.label) * ", " * italic(r) == .(unpaired_r.val))

unpaired_stat.test <- unpaired_stat.test %>%
  add_xy_position(x = "stage")

# boxplot with label
combined_unpaired_test_boxplot_with.label <- ggplot(
  combined_cluster_pw_ratio_long_filt,
  aes(
    x = stage,
    y = cluster_pathway_ratio,
    fill = stage
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.37
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
      combined_cluster_pw_ratio_long_filt,
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
      combined_cluster_pw_ratio_long_filt,
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
    y = "Ratio",
    title = "Ratio of inter-Oi-cluster signaling complexity (E9.5 vs E13.5)_unpaired",
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

combined_unpaired_test_boxplot_with.label

ggsave(
  filename = paste0(combined_output.dir, "E9.5+E13.5_Oi-cluster_inter-Oi-cluster_signaling_complexity_ratio_unpaired_wilcox_boxplot.with.label.pdf"),
  combined_unpaired_test_boxplot_with.label,
  width = 5,
  height = 7
)

###################################################################
#  Post-hoc one-sided test (H0: E13.5 <= E9.5, H1: E13.5 > E9.5)  #
###################################################################

### Wilcoxon Signed-Ranks Test (paired, reversed one-sided: H0: E13.5 <= E9.5, H1: E13.5 > E9.5)

# perform Wilcoxon Signed-Ranks Test
paired_stat.test_rev <- combined_cluster_pw_ratio_long_paired_filt %>%
  wilcox_test(cluster_pathway_ratio ~ stage, paired = TRUE, alternative = "less") %>%
  add_significance()

paired_stat.test_rev

# compute effect size
paired_effect.size_rev <- combined_cluster_pw_ratio_long_paired_filt %>%
  wilcox_effsize(cluster_pathway_ratio ~ stage, paired = TRUE) %>%
  pull(effsize)

# visualization -- dotplot with labels

# construct subtitle for plotting
paired_r.val_rev <- round(paired_effect.size_rev, 3)
paired_p.label_rev <- get_test_label(paired_stat.test_rev, detailed = TRUE)
paired_subtitle_expr_rev <- bquote(.(paired_p.label_rev) * ", " * italic(r) == .(paired_r.val_rev))

paired_stat.test_rev <- paired_stat.test_rev %>%
  add_xy_position(x = "stage")

# dotplot with labels
combined_paired_test_dotplot_with.label_rev <- ggplot(
  combined_cluster_pw_ratio_long_paired_filt,
  aes(
    x = stage,
    y = cluster_pathway_ratio,
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
      combined_cluster_pw_ratio_long_paired_filt,
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
      combined_cluster_pw_ratio_long_paired_filt,
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
  scale_color_manual(
    values = c("#eb4d00", "#87a8eb")
  ) +
  stat_pvalue_manual(
    paired_stat.test_rev, tip.length = 0
  ) +
  labs(
    x = "Stage",
    y = "Ratio",
    title = "Ratio of inter-Oi-cluster signaling complexity (E9.5 vs E13.5)_paired_reversed",
    subtitle = paired_subtitle_expr_rev
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

combined_paired_test_dotplot_with.label_rev

ggsave(
  filename = paste0(combined_output.dir, "E9.5+E13.5_Oi-cluster_inter-to-total_pathway_number_ratio_clusterfilter_dotplot_paired_wilcox_with.label_reversed.pdf"),
  combined_paired_test_dotplot_with.label_rev,
  width = 5,
  height = 7
)


### Mann-Whitney U Test (Wilcoxon Rank Sum Test) (unpaired, reversed one-sided: H0: E13.5 <= E9.5, H1: E13.5 > E9.5)

# perform Mann-Whitney U Test
unpaired_stat.test_rev <- combined_cluster_pw_ratio_long_filt %>%
  wilcox_test(cluster_pathway_ratio ~ stage, alternative = "less") %>%
  add_significance()

unpaired_stat.test_rev

# Compute effect size
unpaired_effect.size_rev <- combined_cluster_pw_ratio_long_filt %>%
  wilcox_effsize(cluster_pathway_ratio ~ stage) %>%
  pull(effsize)

# visualization -- boxplot with labels

# construct subtitle for plotting
unpaired_r.val_rev <- round(unpaired_effect.size_rev, 3)
unpaired_p.label_rev <- get_test_label(unpaired_stat.test_rev, detailed = TRUE)
unpaired_subtitle_expr_rev <- bquote(.(unpaired_p.label_rev) * ", " * italic(r) == .(unpaired_r.val_rev))

unpaired_stat.test_rev <- unpaired_stat.test_rev %>%
  add_xy_position(x = "stage")

# boxplot with label
combined_unpaired_test_boxplot_with.label_rev <- ggplot(
  combined_cluster_pw_ratio_long_filt,
  aes(
    x = stage,
    y = cluster_pathway_ratio,
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
      combined_cluster_pw_ratio_long_filt,
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
      combined_cluster_pw_ratio_long_filt,
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
    unpaired_stat.test_rev, tip.length = 0
  ) +
  labs(
    x = "Stage",
    y = "Ratio",
    title = "Ratio of inter-Oi-cluster signaling complexity (E9.5 vs E13.5)_unpaired_reversed",
    subtitle = unpaired_subtitle_expr_rev
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

combined_unpaired_test_boxplot_with.label_rev

ggsave(
  filename = paste0(combined_output.dir, "E13.5+E9.5_Oi-cluster_inter-Oi-cluster_signaling_complexity_ratio_unpaired_wilcox_boxplot.with.label_reversed.pdf"),
  combined_unpaired_test_boxplot_with.label_rev,
  width = 5,
  height = 7
)