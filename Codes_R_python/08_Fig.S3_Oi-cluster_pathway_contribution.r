library(glue)
library(future)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
options(stringsAsFactors = FALSE)

library(showtext)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()
# set ggplot2's default font
theme_set(theme_gray(base_family = "Arial"))



###################
#  Color palette  #
###################

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

colors_E9.5  <- colors.use_pre[!names(colors.use_pre) %in% c("27", "29", "30", "31")]
colors_E13.5 <- colors.use_pre



################
#  Parameters  #
################

seurat_res <- "0.7"
trim_para  <- "0.05"

base_dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter"

# stages to analyze
stages <- c("E9.5", "E13.5")

# paths to the CellChat inferred signaling network (pathway level), per stage
df_netP_path <- setNames(
  lapply(stages, function(st) {
    file.path(
      base_dir,
      glue("{st}_E1S1"),
      glue("population.size=TRUE_trim={trim_para}"),
      glue("{st}_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.netP.CSV")
    )
  }),
  stages
)



####################################################################################################
#                                                                                                  #
#  Analysis per stage: extract intra- & inter-Oi-cluster signaling and keep the results in memory  #
#                                                                                                  #
####################################################################################################

# Containers for per-stage, per-cluster results
# intra_by_cluster[[stage]][[cluster]] -> data.frame (pathway_name, prob, ...)
# inter_by_cluster[[stage]][[cluster]] -> data.frame (pathway_name, total_prob)
intra_by_cluster <- setNames(vector("list", length(stages)), stages)
inter_by_cluster <- setNames(vector("list", length(stages)), stages)

for (stage_used in stages) {

  ###############
  #  Load data  #
  ###############

  df_netP <- read.csv2(
    df_netP_path[[stage_used]]
  )

  output_dir <- file.path(
    base_dir,
    "Check_pathway_contribution",
    glue("{stage_used}_E1S1"),
    glue("population.size=TRUE_trim={trim_para}")
  )

  output_dir_intra <- file.path(output_dir, "intra-signaling")
  output_dir_inter <- file.path(output_dir, "inter-signaling")

  if (!dir.exists(output_dir_intra)) dir.create(output_dir_intra, recursive = TRUE)
  if (!dir.exists(output_dir_inter)) dir.create(output_dir_inter, recursive = TRUE)


  ########################################
  #                                      #
  #  Extract intra-Oi-cluster signaling  #
  #                                      #
  ########################################

  # all intra-Oi-cluster signaling
  df_intra_all <- df_netP %>%
    filter(source == target) %>%
    arrange(desc(prob))

  # save the arranged data
  write.csv(
    df_intra_all,
    file = file.path(
      output_dir_intra,
      "intra-Oi-cluster_pathway_contribution.CSV"
    )
  )

  # per-Oi-cluster
  all_cluster <- sort(unique(c(df_netP$source, df_netP$target)))

  intra_list <- setNames(vector("list", length(all_cluster)), as.character(all_cluster))

  for (cl in all_cluster) {

    df_intra <- df_netP %>%
      filter(source == cl, target == cl) %>%
      arrange(desc(prob))

    intra_list[[as.character(cl)]] <- df_intra

    write.csv(
      df_intra,
      file = file.path(
        output_dir_intra,
        glue("Oi-cluster-{cl}_intra-Oi-cluster_pathway_contribution.CSV")
      )
    )
  }

  intra_by_cluster[[stage_used]] <- intra_list


  ########################################
  #                                      #
  #  Extract inter-Oi-cluster signaling  #
  #                                      #
  ########################################

  # per-Oi-cluster
  inter_list <- setNames(vector("list", length(all_cluster)), as.character(all_cluster))

  for (cl in all_cluster) {

    df_inter_summary <- df_netP %>%
      filter(source != target) %>%
      filter(source == cl | target == cl) %>%
      group_by(pathway_name) %>%
      summarise(total_prob = sum(prob), .groups = "drop") %>%
      arrange(desc(total_prob))

    inter_list[[as.character(cl)]] <- df_inter_summary

    write.csv(
      df_inter_summary,
      file.path(
        output_dir_inter,
        glue("Oi-cluster-{cl}_inter-Oi-cluster_pathway_contribution.CSV")
      )
    )
  }

  inter_by_cluster[[stage_used]] <- inter_list
}



#####################################
#                                   #
#  Visualization -- dodged barplot  #
#                                   #
#####################################

# Shared Oi-clusters for plotting (excluding E13.5-specific Oi-clusters 27, 29, 30, 31)
Oi_cluster_id <- as.character(setdiff(1:31, c(27, 29, 30, 31)))

for (i in Oi_cluster_id) {
  checked_Oi_cluster <- i
  checked_type        <- "intra-signaling" # "inter-signaling" or "intra-signaling"
  checked_type_simple <- "intra" # "inter" or "intra"

  # path for saving output plots
  plot_output_dir <- file.path(
    base_dir,
    "Check_pathway_contribution",
    "Combined_E9.5_E13.5",
    glue("Oi-cluster_{checked_Oi_cluster}"),
    glue("{checked_type}")
  )

  if (!dir.exists(plot_output_dir)) {
    dir.create(plot_output_dir, recursive = TRUE)
  }


  ###############
  #  Load data  #
  ###############

  # read directly from the variables generated above
  if (checked_type_simple == "intra") {
    E95  <- intra_by_cluster[["E9.5"]][[checked_Oi_cluster]]
    E135 <- intra_by_cluster[["E13.5"]][[checked_Oi_cluster]]
  } else {
    E95  <- inter_by_cluster[["E9.5"]][[checked_Oi_cluster]]
    E135 <- inter_by_cluster[["E13.5"]][[checked_Oi_cluster]]
  }

  # skip clusters absent in a stage
  if (is.null(E95))  E95  <- data.frame()
  if (is.null(E135)) E135 <- data.frame()

  if (checked_type_simple == "intra") {
    full <- full_join(
      rename(E95, E9.5 = prob),
      rename(E135, E13.5 = prob),
      by = "pathway_name"
    )
  } else {
    full <- full_join(
      rename(E95, E9.5 = total_prob),
      rename(E135, E13.5 = total_prob),
      by = "pathway_name"
    )
  }

  full <- full %>%
    select(c("pathway_name", "E9.5", "E13.5"))

  # change the naming
  full$pathway_name <- tools::toTitleCase(tolower(full$pathway_name))
  # correct "ncWnt"
  full$pathway_name <-
    gsub("^Ncwnt$", "ncWnt", full$pathway_name)

  full_long <- full %>%
    pivot_longer(
      cols = -pathway_name,
      names_to = "Stage",
      values_to = "Prob"
    )

  pathway_level <- full_long %>%
    filter(Stage == "E9.5") %>%
    arrange(desc(Prob)) %>%
    pull(pathway_name)

  full_long$Stage <- factor(full_long$Stage, levels = c("E9.5", "E13.5"))
  full_long$pathway_name <- factor(full_long$pathway_name, levels = pathway_level)

  ### Dodged barplot
  barplot_p <- ggplot(full_long, aes(x = pathway_name, y = Prob, fill = Stage)) +
    geom_col(
      position = position_dodge(width = 0.75),
      width = 0.65,
      linewidth = 0.4
    ) +
    scale_fill_manual(
      values = c(
        "E9.5" = "#eb4d00",
        "E13.5" = "#87a8eb"
      )
    ) +
    labs(
      x = "Signaling Pathways",
      y = "Signaling strength",
      title = glue("Oi-cluster {checked_Oi_cluster} {checked_type} signaling strength by pathway")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 8),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )

  barplot_p

  width_flex <- 4 * length(colnames(full_long))
  ggsave(
    file.path(plot_output_dir, "Signaling_strength_boxplot_stage_dodged.pdf"),
    plot = barplot_p,
    width = width_flex, height = 4
  )
}