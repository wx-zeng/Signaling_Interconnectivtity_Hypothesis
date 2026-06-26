options(stringsAsFactors = FALSE)
library(anndata)
library(tidyverse)
library(glue)

library(igraph)
library(assortnet)

library(viridisLite) # heatmap
library(ggplot2)
library(showtext)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()
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

cluster_levels <- as.character(1:31)
colors_all <- glasbey_dark_2[seq_along(cluster_levels)]
names(colors_all) <- cluster_levels

colors_E9.5  <- colors_all[!names(colors_all) %in% c("27","29","30","31")]
colors_E13.5 <- colors_all

# Parameters

stage_used  <- "E9.5"
seurat_res <- "0.7"
trim_para  <- "0.05"

base_dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter"

# path to E9.5/E13.5 bin 50 datasets after Oi-cluster identification
h5ad_path <- list(
  E9.5  = "/home/project_interconnectivity/Mouse_spatialtranscriptome_mosta_bgi_2022/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA_v2/res=0.7/E9.5_E1S1_bin50.h5ad",
  E13.5 = "/home/project_interconnectivity/Mouse_spatialtranscriptome_mosta_bgi_2022/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA_v2/res=0.7/E13.5_E1S1_bin50.h5ad"
)

output_dir <- file.path(
  base_dir,
  glue("{stage_used}_E1S1"),
  glue("population.size=TRUE_trim={trim_para}")
)

colors_use <- if (stage_used == "E9.5") colors_E9.5 else colors_E13.5

# Values of the relative strength of inter-Oi-cluster signaling (used for node heat coloring)
relative_inter_strength_result_path <- "/home/project_interconnectivity/result/clusterfilter_inter_signaling_results_onesided/E9.5+E13.5_comparison/Relative_strength_of_inter-Oi-cluster_signaling_cross-stage_comparison/Combined_E9.5+E13.5_Oi-cluster_relative_strength_of_inter-Oi-cluster_signaling.CSV"



############################
#                          #
#  Functions for analysis  #
#                          #
############################

# function to load E9.5/E13.5 bin 50 datasets after Oi-cluster identification
load_spatial_data <- function(h5ad_file) {
  obj <- read_h5ad(h5ad_file)
  meta <- obj$obs
  meta$labels <- meta$seurat_clusters
  
  center <- meta %>%
    group_by(labels) %>%
    summarise(
      x_mean = mean(x_center),
      y_mean = mean(y_center),
      .groups = "drop"
    )
  
  cluster_size <- meta %>%
    count(seurat_clusters) %>%
    deframe()
  
  storage.mode(cluster_size) <- "double"
  
  list(meta = meta, center = center, cluster_size = cluster_size)
}

# function to plot Oi-clusters in spatial context
plot_spatial_clusters <- function(meta, colors, title) {
  ggplot(
    data.frame(
      x = meta$y_center,
      y = meta$x_center,
      cluster = meta$seurat_clusters
    ),
    aes(x = x, y = y, color = cluster)
  ) +
    geom_point(size = 0.5) +
    scale_color_manual(values = colors) +
    coord_fixed() +
    theme_void() +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
}

# function to map values to thermal colors
map_to_thermal <- function(values, n = 256) { 
  cols <- plasma(n) 
  cols[floor(values * (n - 1)) + 1] 
}

# function to build undirected graph with edge weights based on total strength of interactions between Oi-clusters, and node colors based on Oi-cluster colors
build_graph_strength <- function(df_netP, center_coor, cluster_size, colors, flip_y = FALSE) {
  
  df_undirected <- df_netP %>%
    mutate(
      node1 = pmin(source, target),
      node2 = pmax(source, target)
    ) %>%
    group_by(node1, node2) %>%
    summarise(
      total_strength = sum(prob),
      .groups = "drop"
    )
  
  g <- graph_from_data_frame(
    df_undirected %>% dplyr::select(node1, node2, total_strength),
    directed = FALSE
  )
  
  E(g)$weight <- df_undirected$total_strength
  V(g)$color  <- colors[V(g)$name]
  V(g)$size   <- 3
  E(g)$width  <- E(g)$weight * 120000
  E(g)$edge.color <- "gray20"
  
  center_coor <- center_coor %>%
    filter(labels %in% V(g)$name) %>%
    arrange(match(labels, V(g)$name))
  
  layout0 <- cbind(
    center_coor$y_mean,
    if (flip_y) -center_coor$x_mean else center_coor$x_mean
  )
  
  layout_fr <- layout_with_fr(
    g, 
    coords = layout0, 
    niter = 100,
    grid = "nogrid"
  )
  
  list(graph = g, layout = layout_fr)
}

# function to build undirected graph with edge weights based on total strength of interactions between Oi-clusters, 
# and node colors based on heatmap colors mapped from the values of the relative strength of inter–Oi-cluster signaling
build_graph_strength_heat <- function(df_netP, center_coor, ratio_vec, flip_y = FALSE) {
  
  df_undirected <- df_netP %>%
    mutate(
      node1 = pmin(source, target),
      node2 = pmax(source, target)
    ) %>%
    group_by(node1, node2) %>%
    summarise(
      total_strength = sum(prob),
      .groups = "drop"
    )
  
  g <- graph_from_data_frame(
    df_undirected %>% dplyr::select(node1, node2, total_strength),
    directed = FALSE
  )
  
  E(g)$weight <- df_undirected$total_strength
  V(g)$size   <- 12
  E(g)$width  <- E(g)$weight * 120000
  E(g)$edge.color <- "gray20"
  
  # heatmap color mapping for node coloration
  node_values <- ratio_vec[V(g)$name]
  V(g)$color <- map_to_thermal(node_values)
  
  center_coor <- center_coor %>%
    filter(labels %in% V(g)$name) %>%
    arrange(match(labels, V(g)$name))
  
  layout0 <- cbind(
    center_coor$y_mean,
    if (flip_y) -center_coor$x_mean else center_coor$x_mean
  )
  
  layout_fr <- layout_with_fr(
    g, 
    coords = layout0, 
    niter = 100,
    grid = "nogrid"
  )
  
  list(graph = g, layout = layout_fr)
}



##################
#                #
#  Run analysis  #
#                #
##################

# Load spatial bin 50 data
spatial <- load_spatial_data(h5ad_path[[stage_used]])

# Plot Oi-clusters in spatial context
spatial_p <- plot_spatial_clusters(
  spatial$meta,
  colors_use,
  glue("{stage_used}_E1S1_bin50 | Seurat res={seurat_res}")
)
print(spatial_p)

# Load the CellChat inferred signaling network (pathway level)
df_netP <- read.csv2(
  file.path(
    output_dir,
    glue("{stage_used}_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.netP.CSV")
  )
)

# Load the values of the relative strength of inter–Oi-cluster signaling of Oi-clusters
relative_inter_strength_all <- read.csv(
  file = relative_inter_strength_result_path
)

relative_inter_strength_E9.5 <- relative_inter_strength_all[, "E9.5_ratio"]
names(relative_inter_strength_E9.5) <- rownames(relative_inter_strength_all)

# check
print(relative_inter_strength_E9.5)

# Build graph -- signaling strength as weight
graph_res_strength <- build_graph_strength(
  df_netP,
  spatial$center,
  spatial$cluster_size,
  colors_use,
  flip_y = (stage_used == "E13.5")
)

graph_res_strength_heat <- build_graph_strength_heat(
  df_netP,
  spatial$center,
  relative_inter_strength_E9.5,
  flip_y = (stage_used == "E13.5")
)


###########
#  Plots  #
###########

# create output directory for saving results
plot_dir_strength <- file.path(output_dir, "assortativity_analysis", "signaling.strength_as_weight")

if (!dir.exists(plot_dir_strength)) {
  dir.create(plot_dir_strength, recursive = TRUE, showWarnings = FALSE)
}


### Network plot with signaling strength as weight and normal node

lay <- graph_res_strength$layout
# center
lay[,1] <- lay[,1] - mean(lay[,1])
lay[,2] <- lay[,2] - mean(lay[,2])

# use a single scale factor for both axes
rng <- max(diff(range(lay[,1])), diff(range(lay[,2])))
enlarge_factor <- 2
lay_scaled <- (lay / rng) * enlarge_factor

png(
  file.path(plot_dir_strength, glue("{stage_used}_res={seurat_res}_trim={trim_para}_network_strength.as.weight_node.size=fix_normal_node.png")),
  width = 4500, height = 5000,
  res = 300,
  bg = "transparent"
)

plot(
  graph_res_strength$graph,
  layout = lay_scaled,
  edge.loop.angle = 0,
  vertex.frame.color = NA,
  vertex.label = V(graph_res_strength$graph)$name,
  vertex.label.family = "Arial",
  vertex.label.cex = 2.5,
  vertex.label.color = "black",
  vertex.label.dist = 0.8, # 0 = on top of node; increase to push outward
  vertex.label.degree = 0, # direction (radians); 0 = to the right
  asp = 1,
  rescale = FALSE
)

dev.off()


### Network plot with signaling strength as weight and node color mapped from the relative strength of inter–Oi-cluster signaling

lay <- graph_res_strength_heat$layout
# center
lay[,1] <- lay[,1] - mean(lay[,1])
lay[,2] <- lay[,2] - mean(lay[,2])

# use a single scale factor for both axes
rng <- max(diff(range(lay[,1])), diff(range(lay[,2])))
enlarge_factor <- 2
lay_scaled <- (lay / rng) * enlarge_factor

png(
  file.path(plot_dir_strength, glue("{stage_used}_res={seurat_res}_trim={trim_para}_network_strength.as.weight_node.size=fix_node.color=heat.png")),
  width = 4500, height = 5000,
  res = 300,
  bg = "transparent"
)

plot(
  graph_res_strength_heat$graph,
  layout = lay_scaled,
  edge.loop.angle = 0,
  vertex.frame.color = NA,
  vertex.label = V(graph_res_strength_heat$graph)$name,
  vertex.label.family = "Arial",
  vertex.label.cex = 2.5,
  vertex.label.color = "black",
  vertex.label.dist = 0, # 0 = on top of node; increase to push outward
  vertex.label.degree = 0, # direction (radians); 0 = to the right
  asp = 1,
  rescale = FALSE
)

dev.off()



#############################################################################
#                                                                           #
#  Calculate assotativity coefficient using assortnet::assortment.discrete  #
#                                                                           #
#############################################################################

# signaling strength as weight
adj_mat_strength <- as_adjacency_matrix(graph_res_strength$graph, sparse = FALSE, attr = "weight")
assortment.discrete(adj_mat_strength, types = V(graph_res_strength$graph)$name, weighted = TRUE, SE = TRUE)
assort_coefficient_strength <- assortment.discrete(adj_mat_strength, types = V(graph_res_strength$graph)$name, weighted = TRUE)$r
cat("Assortativity coefficient of", stage_used, "is", assort_coefficient_strength, "\n")