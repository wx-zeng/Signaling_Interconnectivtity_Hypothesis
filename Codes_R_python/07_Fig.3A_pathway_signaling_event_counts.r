library(glue)
library(dplyr)
library(forcats)
library(ggplot2)
library(patchwork)
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

# Parameters
seurat_res <- "0.7"
trim_para  <- "0.05"

base_dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter"

E9.5_output_dir <- file.path(
  base_dir,
  "E9.5_E1S1",
  glue("population.size=TRUE_trim={trim_para}")
)

E13.5_output_dir <- file.path(
  base_dir,
  "E13.5_E1S1",
  glue("population.size=TRUE_trim={trim_para}")
)



##################################################################
#                                                                #
#  Load the CellChat inferred signaling network (pathway level)  #
#                                                                #
##################################################################

### E9.5

E9.5_df.netp <- read.csv2(
  file.path(
    E9.5_output_dir,
    glue("E9.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.netP.CSV")
  )
)

E9.5_df.netp$source <- as.character(E9.5_df.netp$source)
E9.5_df.netp$target <- as.character(E9.5_df.netp$target)


### E13.5

E13.5_df.netp <- read.csv2(
  file.path(
    E13.5_output_dir,
    glue("E13.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.netP.CSV")
  )
)

E13.5_df.netp$source <- as.character(E13.5_df.netp$source)
E13.5_df.netp$target <- as.character(E13.5_df.netp$target)


### Correct pathway naming to initial capital

# save the original data
E9.5_df.netp_original <- E9.5_df.netp
E13.5_df.netp_original <- E13.5_df.netp

# change the naming
E9.5_df.netp$pathway_name <- tools::toTitleCase(tolower(E9.5_df.netp$pathway_name))
E13.5_df.netp$pathway_name <- tools::toTitleCase(tolower(E13.5_df.netp$pathway_name))

# correct "ncWnt"
E9.5_df.netp$pathway_name <- 
  gsub("^Ncwnt$", "ncWnt", E9.5_df.netp$pathway_name)
E13.5_df.netp$pathway_name <- 
  gsub("^Ncwnt$", "ncWnt", E13.5_df.netp$pathway_name)



###################################
#                                 #
#  Pathway signaling event count  #
#                                 #
###################################

### Calculate signaling event count of each pathway

# E9.5
E9.5_pleiotropy_netp_event <- E9.5_df.netp %>%
  distinct(pathway_name, source, target)

E9.5_pleiotropy_netp <- E9.5_df.netp %>%
  distinct(pathway_name, source, target) %>%
  count(pathway_name, name = "n_signaling_events")


# E13.5
E13.5_pleiotropy_netp_event <- E13.5_df.netp %>%
  distinct(pathway_name, source, target)

E13.5_pleiotropy_netp <- E13.5_df.netp %>%
  distinct(pathway_name, source, target) %>%
  count(pathway_name, name = "n_signaling_events")


# Create new path for saving results
output_path <- file.path(
    base_dir,
    glue("Pathway_signaling_event_counts/")
)

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}


### Visualization -- barplots

x_max_global <- max(E9.5_pleiotropy_netp$n_signaling_events, E13.5_pleiotropy_netp$n_signaling_events)

# Reversed horizontal bar plot (for E9.5)

E9.5_plot_df <- E9.5_pleiotropy_netp %>%
  arrange(desc(n_signaling_events)) %>%
  mutate(
    rank_desc = row_number(),
    rank_asc = n() - row_number() + 1
  )

p1 <- E9.5_plot_df %>%
  ggplot(aes(y = fct_reorder(pathway_name, n_signaling_events), 
             x = n_signaling_events)) +
  geom_col(fill = "#eb4d00", width = 0.6) +
  # geom_text(
  #   aes(label = label),
  #   hjust = +0.6,
  #   size = 3.4
  # ) +
  scale_x_reverse(
    limits = c(0, x_max_global),
    expand = expansion(mult = c(0.1, 0))
  ) +
  scale_y_discrete(position = "right") +
  labs(
    y = NULL,
    x = "Signaling event counts",
    title = "E9.5"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line("black"),
    axis.text.y = element_text(colour = "black", hjust = 0)
  )

# Horizontal bar plot (for E13.5)

E13.5_plot_df <- E13.5_pleiotropy_netp %>%
  arrange(desc(n_signaling_events)) %>%
  mutate(
    rank_desc = row_number(),
    rank_asc = n() - row_number() + 1
  )

p2 <- E13.5_plot_df %>%
  ggplot(aes(y = fct_reorder(pathway_name, n_signaling_events), 
             x = n_signaling_events)) +
  geom_col(fill = "#87a8eb", width = 0.6) +
  # geom_text(
  #   aes(label = label),
  #   hjust = -0.3,
  #   size = 3.4
  # ) +
  scale_x_continuous(
    limits = c(0, x_max_global),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    y = NULL,
    x = "Signaling event counts",
    title = "E13.5"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.line = element_line("black"),
    axis.text = element_text(colour = "black")
  )

# Combine two plots
p_combined <- p1 | p2

# Save plots

ggsave(
  file.path(output_path, "Pathway_signaling_event_counts_barplot_combined_thinner_6.pdf"),
  p_combined,
  width = 6,
  height = 10,
  units = "in"
)
# pdf portrait (6,10)

p1
ggsave(
  file.path(output_path, "Pathway_signaling_event_counts_barplot_E9.5.pdf"),
  p1,
  width = 6,
  height = 10,
  units = "in"
)
# pdf portrait (6,10)

p2
ggsave(
  file.path(output_path, "Pathway_signaling_event_counts_barplot_E13.5.pdf"),
  p2,
  width = 6,
  height = 10,
  units = "in"
)
# pdf portrait (6,10)