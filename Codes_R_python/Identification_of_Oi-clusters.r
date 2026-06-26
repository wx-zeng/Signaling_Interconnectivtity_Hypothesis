library(dplyr)
library(tidyr)
library(glue)

library(Seurat)
library(sctransform)
library(patchwork)
# library(glmGamPoi) # for faster SCTransform
library(anndata)
library(SeuratDisk)
library(scCustomize)

library(CellChat)

options(stringsAsFactors = FALSE)



color_palette<- c(
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

# Output directory for Seurat integration
output_dir <- "/home/project_interconnectivity/result/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}



########################################################
#                                                      #
#  Load bin 50 h5ad data and convert to SeurateObject  #
#                                                      #
########################################################

### E9.5_bin50 data

E9.5_bin50_data <- read_h5ad("/home/project_interconnectivity/result/Mm_E9.5/E9.5_E1S1_bin50/E9.5_bin50_for_Seurat.h5ad")

# access count data
E9.5_bin50_counts <- t(as.matrix(E9.5_bin50_data$X))
# preview the count data
E9.5_bin50_counts[1:5, 1:5] # raw counts data

# access meta data
E9.5_bin50.meta <- E9.5_bin50_data$obs
E9.5_bin50.meta$labels <- E9.5_bin50.meta[['mclust']]
E9.5_bin50.meta$samples <- "E9.5_E1S1_bin50"
E9.5_bin50.meta$samples <- factor(E9.5_bin50.meta$samples)

unique(E9.5_bin50.meta$labels)

E9.5_seurat.data <- E9.5_bin50_counts
E9.5_seurat <- CreateSeuratObject(E9.5_seurat.data, meta.data = E9.5_bin50.meta)
E9.5_seurat$orig.ident <- E9.5_seurat$samples # set orig.ident to samples

### E13.5_bin50 data

E13.5_bin50_data <- read_h5ad("/home/project_interconnectivity/result/Mm_E13.5/E13.5_E1S1_bin50/E13.5_bin50_for_Seurat.h5ad")

# access count data
E13.5_bin50_counts <- t(as.matrix(E13.5_bin50_data$X))
# preview the count data
E13.5_bin50_counts[1:5, 1:5] # raw counts data

# access meta data
E13.5_bin50.meta <- E13.5_bin50_data$obs
E13.5_bin50.meta$labels <- E13.5_bin50.meta[['mclust']]
E13.5_bin50.meta$samples <- "E13.5_E1S1_bin50"
E13.5_bin50.meta$samples <- factor(E13.5_bin50.meta$samples)

unique(E13.5_bin50.meta$labels)

E13.5_seurat.data <- E13.5_bin50_counts
E13.5_seurat <- CreateSeuratObject(E13.5_seurat.data, meta.data = E13.5_bin50.meta)
E13.5_seurat$orig.ident <- E13.5_seurat$samples # set orig.ident to samples



# merge E9.5 and E13.5 data
aggr_obj <- merge(
  x = E9.5_seurat,
  y = E13.5_seurat,
  add.cell.ids = c("E9.5", "E13.5")
)

aggr_obj$orig.ident <- factor(aggr_obj$orig.ident)

gene.name <- rownames(aggr_obj)
mito.gene <- grep("mt-", gene.name, value = TRUE)

# SCTransform on merged data
aggr_obj <- PercentageFeatureSet(aggr_obj, pattern = "mt-", col.name = "percent.mt")

options(future.globals.maxSize = 32 * 1024^3)
aggr_obj_SCT <- SCTransform(aggr_obj, vars.to.regress = "percent.mt", vst.flavor = "v2", verbose = TRUE)

# check which genes are still remained in the SeuratObject
genes_all <- rownames(aggr_obj[["RNA"]])
genes_retained_sct <- rownames(aggr_obj_SCT[["SCT"]])
writeLines(genes_all, file.path(output_dir, "E9.5_bin50_E13.5_bin50_all_genes.txt"))
writeLines(genes_retained_sct, file.path(output_dir, "E9.5_bin50_E13.5_bin50_retained_genes_after-SCT.txt"))

aggr_obj <- aggr_obj_SCT
rm(aggr_obj_SCT) # remove the SCT object to save memory

# Run PCA
n_PCs_to_use <- 50
aggr_obj <- RunPCA(aggr_obj, npcs = n_PCs_to_use, verbose = TRUE)

# save PCA results
save.path <- file.path(output_dir, "E9.5_bin50_E13.5_bin50_after-SCT_pca.rds")
saveRDS(aggr_obj, file = save.path)



######################################################
#                                                    #
#  Find Neighbors & clusters (on unintegrated data)  #
#                                                    #
######################################################

aggr_obj_PCA <- aggr_obj
aggr_obj_NonIntegrated <- FindNeighbors(
    aggr_obj_PCA,
    reduction = "pca",
    dims = 1:n_PCs_to_use,
    verbose = TRUE,
)
aggr_obj_NonIntegrated <- FindClusters(
    aggr_obj_NonIntegrated,
    resolution = 0.8,
    algorithm = 4, # 4 for Leiden
    random.seed = 1, # default is 1
    verbose = TRUE,
)
aggr_obj_NonIntegrated <- RunUMAP(
    aggr_obj_NonIntegrated,
    dims = 1:n_PCs_to_use,
    reduction = "pca",
    reduction.name = "umap.unintegrated",
    verbose = TRUE
)

# save UMAP DimPlots
p <- DimPlot(aggr_obj_NonIntegrated,
        reduction = "umap.unintegrated",
        group.by = c("orig.ident", "seurat_clusters"),
        label = FALSE,
        label.size = 3) +
    ggtitle("SCTransform | Non integrated | UMAP of E9.5_bin50+E13.5_bin50 data")
ggsave.path <- file.path(output_dir, "E9.5_bin50_E13.5_bin50_after-SCT_umap.pdf")
ggsave(ggsave.path, plot = p, width = 12, height = 7)

# save the Seurat object after SCTransform
save.path <- file.path(output_dir, "E9.5_bin50_E13.5_bin50_after-SCT_umap.rds")
saveRDS(aggr_obj, file = save.path)



##########################################################
#                                                        #
#  Integrate E9.5 and E13.5 bin50 data using Seurat CCA  #
#                                                        #
##########################################################

# Create a new Assay in the Seurat object for integrated data
aggr_obj[["SCT.integrated"]] <- aggr_obj[["SCT"]]
DefaultAssay(aggr_obj) <- "SCT.integrated"

integrated.dr <- "CCAintegration" # Name for the integrated reduction

aggr_obj <- IntegrateLayers(
    aggr_obj,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = integrated.dr,
    normalization.method = "SCT",
    verbose = TRUE,
)

save.path <- file.path(output_dir, "E9.5_bin50_E13.5_bin50_after-CCA.rds")
saveRDS(aggr_obj, file = save.path)



####################################################
#                                                  #
#  Find Neighbors & clusters (on integrated data)  #
#                                                  #
####################################################

cl_algorithm <- 4 # 4 for Leiden

aggr_obj <- FindNeighbors(
    aggr_obj,
    dims = 1:n_PCs_to_use,
    reduction = integrated.dr,
    assay = "SCT.integrated",
    verbose = TRUE,
)

if (cl_algorithm == 1) {
    cl_algorithm_name <- "louvain"
} else if (cl_algorithm == 4) {
    cl_algorithm_name <- "leiden"
} else {
    stop("Unknown clustering algorithm")
}

# Run UMAP on integrated data (independent of clustering resolution)
aggr_obj <- RunUMAP(
    aggr_obj,
    reduction = integrated.dr,
    dims = 1:n_PCs_to_use,
    reduction.name = "umap.integrated",
    verbose = TRUE,
)


# Export different FindClusters results with different clustering resolutions (from 0.5 to 1.0, step = 0.1)
# find_cl_res: clustering resolution parameter for FindClusters function
for (find_cl_res in seq(0.5, 1.0, by = 0.1)) {

    message(paste0("\n\n--- FindClusters: resolution = ", find_cl_res, " ---\n\n"))

    aggr_obj_cl <- FindClusters(
        aggr_obj,
        resolution = find_cl_res,
        algorithm = cl_algorithm,
        graph.name = "SCT_snn",
        random.seed = 1, # default is 1
        verbose = TRUE,
    )

    # print cluster numbers
    cluster_numbers <- table(aggr_obj_cl$seurat_clusters)
    message("Cluster numbers after FindClusters:")
    print(cluster_numbers)

    save.path <- file.path(output_dir, paste0("E9.5_bin50_E13.5_bin50_after-CCA-umap_npc=", n_PCs_to_use, "_res=", find_cl_res, "_", cl_algorithm_name, ".rds"))
    saveRDS(aggr_obj_cl, file = save.path)

    # save UMAP DimPlots: samples + clusters
    p <- DimPlot(
        aggr_obj_cl,
        reduction = "umap.integrated",
        group.by = c("samples", "seurat_clusters"),
        label = FALSE,
        label.size = 3) +
        ggtitle(paste0("UMAP of E9.5_bin50+E13.5_bin50 after CCA integration (npc=", n_PCs_to_use, ", res=", find_cl_res, ", ", cl_algorithm_name, ")"))

    ggsave.path <- file.path(output_dir, paste0("E9.5_bin50_E13.5_bin50_umap_after-CCAIntegration_npc=", n_PCs_to_use, "_res=", find_cl_res, "_", cl_algorithm_name, ".pdf"))
    ggsave(ggsave.path, plot = p, width = 15, height = 6)

    # save UMAP DimPlots: samples only
    p <- DimPlot(
        aggr_obj_cl,
        reduction = "umap.integrated",
        group.by = c("samples"),
        label = FALSE,
        label.size = 3) +
        ggtitle(paste0("UMAP of E9.5_bin50+E13.5_bin50 after CCA integration (npc=", n_PCs_to_use, ", res=", find_cl_res, ", ", cl_algorithm_name, ")"))

    ggsave.path <- file.path(output_dir, paste0("E9.5_bin50_E13.5_bin50_umap_after-CCAIntegration_npc=", n_PCs_to_use, "_res=", find_cl_res, "_", cl_algorithm_name, "_integration.only.pdf"))
    ggsave(ggsave.path, plot = p, width = 8, height = 6)

    # save UMAP DimPlots: clusters only
    p <- DimPlot(
        aggr_obj_cl,
        reduction = "umap.integrated",
        group.by = c("seurat_clusters"),
        label = FALSE,
        label.size = 3) +
        ggtitle(paste0("UMAP of E9.5_bin50+E13.5_bin50 after CCA integration (npc=", n_PCs_to_use, ", res=", find_cl_res, ", ", cl_algorithm_name, ")"))

    ggsave.path <- file.path(output_dir, paste0("E9.5_bin50_E13.5_bin50_umap_after-CCAIntegration_npc=", n_PCs_to_use, "_res=", find_cl_res, "_", cl_algorithm_name, "_cluster.only.pdf"))
    ggsave(ggsave.path, plot = p, width = 8, height = 6)

    # Save the Seurat cluster as a data frame
    seurat_clusters_df <- data.frame(
        cell = colnames(aggr_obj_cl),
        seurat_clusters = aggr_obj_cl$seurat_clusters,
        samples = aggr_obj_cl$samples
    )
    write.csv(
        seurat_clusters_df,
        file = file.path(output_dir, paste0("E9.5_bin50_E13.5_bin50_after-CCA_npc=", n_PCs_to_use, "_res=", find_cl_res, "_", cl_algorithm_name, "_seurat_clusters.csv")),
        row.names = FALSE
    )
}