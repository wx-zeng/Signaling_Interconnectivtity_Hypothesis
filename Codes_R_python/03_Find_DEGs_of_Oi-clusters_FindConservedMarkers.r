library(dplyr)
library(Seurat)
library(patchwork)
library(tidyr)
library(scCustomize)
options(stringsAsFactors = FALSE)



# load the integrated E9.5 + E13.5 bin 50 Seurat object; Oi-clusters already annotated
aggr_obj <- readRDS("/home/project_interconnectivity/Mouse_spatialtranscriptome_mosta_bgi_2022/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA_v2/res=0.7/E9.5_bin50_E13.5_bin50_after-CCA-umap_npc=50_res=0.7_leiden.rds")



############################################
#                                          #
#   Find DEGs using FindConservedMarkers   #
#                                          #
############################################

aggr_obj[['RNA']] <- JoinLayers(aggr_obj[['RNA']])
aggr_obj[['RNA']]$data <- aggr_obj[['RNA']]$counts
Idents(aggr_obj) <- "seurat_clusters"
cluster_ID <- sort(as.numeric(unique(as.character(aggr_obj@meta.data$seurat_clusters))))

origin_marker_list <- list()
filtered_marker_list <- list()
top_marker_list <- list()

save_path <- "/home/project_interconnectivity/result/Oi-Cluster_seurat_DEG_FindConservedMarkers/"

if (!dir.exists(save_path)){
  dir.create(save_path, recursive = TRUE)
}


for (i in cluster_ID) {
  
  cluster_now <- as.character(i)
  
  message("\n----------------------------")
  message("Now Processing Oi-cluster: ", cluster_now)
  message("----------------------------\n")

  save_path_detailed <- file.path(save_path, paste0("All.conserved.markers_Oi-cluster_", cluster_now, ".csv"))
  
  conserved_markers <- FindConservedMarkers(
    joined_obj,
    ident.1 = i,
    ident.2 = NULL, # default 
    grouping.var = "orig.ident", # E9.5 or E13.5
    verbose = TRUE
  )
  
  saveRDS(conserved_markers, file = file.path(save_path, paste0("All.conserved.markers_Oi-cluster_", cluster_now, ".rds")))
  write.csv(conserved_markers, file = save_path_detailed, row.names = TRUE)
  
  origin_markers <- conserved_markers %>%
    tibble::rownames_to_column("gene") %>%
    filter(!is.na(gene)) %>%
    mutate(cluster = cluster_now)
  
  # filter markers that pass the significant test
  colnames_p_val <- grep("_p_val_adj$", colnames(origin_markers), value = TRUE)
  filtered_markers <- origin_markers %>%
    filter(if_all(all_of(colnames_p_val), ~ . < 0.05))
  
  # top 30 based on raw FindConservedMarkers() result
  top_markers <- filtered_markers %>%
    slice_head(n = 30)
  
  origin_marker_list[[cluster_now]] <- origin_markers
  filtered_marker_list[[cluster_now]] <- filtered_markers
  top_marker_list[[cluster_now]] <- top_markers
}


origin_marker_all <- dplyr::bind_rows(origin_marker_list)
origin_marker_all <- origin_marker_all %>% arrange(as.numeric(cluster))

filtered_marker_all <- dplyr::bind_rows(filtered_marker_list)
filtered_marker_all <- filtered_marker_all %>% arrange(as.numeric(cluster))

top_marker_all <- dplyr::bind_rows(top_marker_list)
top_marker_all <- top_marker_all %>% arrange(as.numeric(cluster))

write.csv(origin_marker_all, file = file.path(save_path, "Oi-Cluster_All_original_conserved_markers_all_clusters.csv"), row.names = FALSE)
write.csv(filtered_marker_all, file = file.path(save_path, "Oi-Cluster_Filtered_conserved_markers_all_clusters.csv"), row.names = FALSE)
write.csv(top_marker_all, file = file.path(save_path, "Oi-Cluster_Top30_conserved_markers_all_clusters.csv"), row.names = FALSE)