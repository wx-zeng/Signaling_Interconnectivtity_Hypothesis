library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(scCustomize)
options(stringsAsFactors = FALSE)



# load the integrated E9.5 + E13.5 bin 50 Seurat object; Oi-clusters already annotated
aggr_obj <- readRDS("/home/project_interconnectivity/Mouse_spatialtranscriptome_mosta_bgi_2022/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA_v2/res=0.7/E9.5_bin50_E13.5_bin50_after-CCA-umap_npc=50_res=0.7_leiden.rds")



######################################
#                                    #
#   Find DEGs using FindAllMarkers   #
#                                    #
######################################

DefaultAssay(aggr_obj) <- "RNA"

aggr_obj <- JoinLayers(
  aggr_obj,
  assay = "RNA",
  layers = "counts"
)

print(aggr_obj@assays$RNA$counts[1:10, 1:10])

# use Log-Normalized raw counts for detecting DEGs
aggr_obj[["RNA.LogNormalized"]] <- aggr_obj[["RNA"]]
DefaultAssay(aggr_obj) <- "RNA.LogNormalized"

aggr_obj <- NormalizeData(
  aggr_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


save_path <- "/home/project_interconnectivity/result/Oi-Cluster_seurat_DEG_FindAllMarkers/"

if (!dir.exists(save_path)){
  dir.create(save_path, recursive = TRUE)
}



######################################
#                                    #
#   Find DEGs using FindAllMarkers   #
#                                    #
######################################

all_markers <- FindAllMarkers(
  object = aggr_obj,
  assay = "RNA.LogNormalized",
  #slot = "data", # default
  group.by = "seurat_clusters",
  only.pos = TRUE,
  min.pct = 0.01,
  logfc.threshold = 0.1,
  test.use = "wilcox",
  verbose = TRUE,
  random.seed = 1
)

# preview
all_markers %>% head(20)
all_markers %>% tail(20)


# select top 30 DEGs based on avg_log2FC for each Oi-cluster
method_select_top_markers <- "fc"   # avg_log2FC

if (method_select_top_markers == "fc") {
  # make sure pval < 0.05
  pval_cutoff <- 0.05
  
  all_markers_rows_ordered <- all_markers_OnlyPos %>%
    filter(p_val_adj < pval_cutoff) %>%
    group_by(cluster) %>%  
    slice_max(n = 30, order_by = avg_log2FC)
  
  top_markers <- all_markers_rows_ordered %>%
    pull(gene)
  
  names(top_markers) <- all_markers_rows_ordered$cluster
} else if (method_select_top_markers == "pval") {
  all_markers_rows_ordered <- all_markers_OnlyPos %>%
    group_by(cluster) %>%  
    slice_min(n = 20, order_by = p_val_adj) 
  
  top_markers <- all_markers_rows_ordered %>%
    pull(gene)
  
  names(top_markers) <- all_markers_rows_ordered$cluster
} else {
  stop("Invalid method to select top markers")
}

# preview
top_markers %>% head(20)
top_markers %>% tail(20)

{
  message('> Preview of the markers:')
  print(all_markers %>% head())
  message('...')
  print(all_markers %>% tail())
  message('\n\n')
  
  message('> Preview of the top markers:')
  print(top_markers %>% head())
  message('...')
  print(top_markers %>% tail())
  message('\n\n')
}

marker_matrix <- all_markers_rows_ordered %>%
  group_by(cluster) %>%
  mutate(rank = row_number()) %>%
  select(cluster, rank, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene) %>%
  arrange(rank) %>%
  select(-rank)

save_fic_fc <- file.path(save_path, "selected_by_fc")
if(!dir.exists(save_fic_fc)) {
  dir.create(save_fic_fc)
}

write.csv(all_markers, file = file.path(save_fic_fc, "all_markers.csv"))
write.csv(all_markers_OnlyPos, file = file.path(save_fic_fc, "all_markers_OnlyPos.csv"))
write.csv(all_markers_rows_ordered, file = file.path(save_fic_fc, "Top30_markers_by_fc_RowOrdered.csv"))
write.csv(marker_matrix, file = file.path(save_fic_fc, "Top30_markers_by_fc_OnlyGene.csv"))