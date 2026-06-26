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


stage_used <- "E9.5" # E9.5 or E13.5

if (stage_used == "E9.5") {
  output.dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/E9.5_E1S1/"
} else if (stage_used == "E13.5") {
  output.dir <- "/home/project_interconnectivity/result/MOSTA2022_Mm_E9.5+E13.5_E1S1_bin50_CCACluster_res=0.7_CCC_analysis_clusterfilter/E13.5_E1S1/"
}

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE)
}



####################
#                  #
#  Load h5ad data  #
#                  #
####################

# load E9.5 & E13.5 bin 50 datasets after Oi-cluster identification
if (stage_used == "E9.5") {
  sp_bin50_data_seurat <- read_h5ad("/home/project_interconnectivity/Mouse_spatialtranscriptome_mosta_bgi_2022/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA_v2/res=0.7/E9.5_E1S1_bin50.h5ad")
} else if (stage_used == "E13.5") {
  sp_bin50_data_seurat <- read_h5ad("/home/project_interconnectivity/Mouse_spatialtranscriptome_mosta_bgi_2022/Seurat_integration_E9.5+E13.5/E9.5_bin50+E13.5_bin50/SCTransform+CCA_v2/res=0.7/E13.5_E1S1_bin50.h5ad")
}



##############################################
#                                            #
#  Prepare input data for CellChat analysis  #
#                                            #
##############################################

sp_bin50_raw_counts <- t(as.matrix(sp_bin50_data_seurat$X))
sp_bin50_counts <- sp_bin50_raw_counts
# check the count data
sp_bin50_counts[1:10, 1:10]

# normalize the count data
bin50.library.size <- Matrix::colSums(sp_bin50_counts)
bin50.data.input <- as(log1p(Matrix::t(Matrix::t(sp_bin50_counts)/bin50.library.size) * 10000), "dgCMatrix")

# access meta data
bin50.meta <- sp_bin50_data_seurat$obs

# use Oi-clusters as the labels of spatial bins for CellChat analysis (stored in "seurat_clusters")
bin50.meta$labels <- bin50.meta[["seurat_clusters"]]

unique(bin50.meta$labels)
bin50.meta$samples <- "sample1"
bin50.meta$samples <- factor(bin50.meta$samples)

# load spatial coordinates
bin50.coordinates <- as.data.frame(bin50.meta[, c("x_center", "y_center")])
# reverse the x-axis
bin50.coordinates$x_center <- -bin50.coordinates$x_center

# prepare spatial factors of spatial distance
bin50.ratio <- 0.5 # for stereo-seq data (500nm is the distance between two DNBs)
bin50.tol <- 12.5
bin50.spatial.factors <- data.frame(ratio = bin50.ratio, tol = bin50.tol)



##############################
#                            #
#  Create a CellChat object  #
#                            #
##############################

bin50.cellchat <- createCellChat(object = bin50.data.input, 
                                 meta = bin50.meta, 
                                 group.by = "labels",
                                 datatype = "spatial",
                                 coordinates = bin50.coordinates,
                                 spatial.factors = bin50.spatial.factors)
bin50.cellchat



##################################################
#                                                #
#  Set the ligand-receptor interaction database  #
#                                                #
##################################################

CellChatDB <- CellChatDB.mouse
dplyr::glimpse(CellChatDB$interaction)

# subset to only "Secreted signaling" for the analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")

bin50.cellchat@DB <- CellChatDB.use

# check CellChatDB.mouse coverage in the spatial data before running CellChat
gene_db <- bin50.cellchat@DB[["geneInfo"]][["Symbol"]] # all genes in CellChatDB.mouse
gene_sp_data <- bin50.cellchat@data@Dimnames[[1]] # all genes in MOSTA spatial data

overlap_gene <- intersect(gene_db, gene_sp_data)
length(overlap_gene)

overlap_ratio_sp_data <- length(overlap_gene) / length(gene_sp_data)
overlap_ratio_db <- length(overlap_gene) / length(gene_db)



#############################################################
#                                                           #
#  Preprocessing the expression data for CellChat analysis  #
#                                                           #
#############################################################

# subset the expression data of signaling genes for saving computation cost
bin50.cellchat <- subsetData(bin50.cellchat)

signaling_genes_sp_data <- bin50.cellchat@data.signaling@Dimnames[[1]]

future::plan("multicore", workers = 6)
options(future.globals.maxSize = 20 * 1024^3) 

bin50.cellchat <- identifyOverExpressedGenes(bin50.cellchat)
bin50.cellchat <- identifyOverExpressedInteractions(bin50.cellchat, variable.both = F) # variable.both = F: either ligand or receptor must be in variable features (but both must be expressed)
# Results:
# E9.5: The number of highly variable ligand-receptor pairs used for signaling inference is 1128
# E13.5: The number of highly variable ligand-receptor pairs used for signaling inference is 1128



####################################################################
#                                                                  #
#  Infer signaling networks and compute communication probability  #
#                                                                  #
####################################################################

bin50.cellchat@images$coordinates <- as.matrix(bin50.cellchat@images$coordinates) # to solve the error of 'queryKNN' when running computeCommunProb()

# Set trim threshold: 0.1, 0.05, 0.01
trim_para <- 0.05

message(">>>>>>>>>> Run computeCommunProb() on Oi-clusters, trim=", trim_para)

output.dir <- paste0(output.dir, "population.size=TRUE_trim=", trim_para, "/")
if (!dir.exists(output.dir)) {
  dir.create(output.dir)
}

bin50.cellchat <- computeCommunProb(bin50.cellchat,
                                    type              = "truncatedMean",
                                    trim              = trim_para,
                                    distance.use      = TRUE,
                                    interaction.range = 250,
                                    scale.distance    = 0.1,
                                    contact.dependent = FALSE,
                                    population.size   = TRUE)

# population size filtering threshold: min.cells = 10
bin50.cellchat <- filterCommunication(bin50.cellchat, min.cells = 10)

#  Infer the cell-cell communication at a signaling pathway level
bin50.cellchat <- computeCommunProbPathway(bin50.cellchat)

#  Extract the inferred signaling network as a data frame
seurat_res <- "0.7"

# at the level of ligand-receptor pairs
df.net <- subsetCommunication(bin50.cellchat)

if (stage_used == "E9.5") {
    write.csv2(df.net, file = file.path(output.dir, glue("E9.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.net.CSV")))
} else if (stage_used == "E13.5") {
    write.csv2(df.net, file = file.path(output.dir, glue("E13.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.net.CSV")))
}

# at the level of signaling pathways
df.netp <- subsetCommunication(bin50.cellchat, slot.name = "netP")

if (stage_used == "E9.5") {
    write.csv2(df.netp, file = file.path(output.dir, glue("E9.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.netP.CSV")))
} else if (stage_used == "E13.5") {
    write.csv2(df.netp, file = file.path(output.dir, glue("E13.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_df.netP.CSV")))
}

#  Calculate the aggregated cell-cell communication network
bin50.cellchat <- aggregateNet(bin50.cellchat)

if (stage_used == "E9.5") {
    write.csv2(bin50.cellchat@net$weight, file = file.path(output.dir, glue("E9.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_net.weight.CSV")))
} else if (stage_used == "E13.5") {
    write.csv2(bin50.cellchat@net$weight, file = file.path(output.dir, glue("E13.5_E1S1_bin50_seuratCCA_clusters_res={seurat_res}_after_aggregateNet_population.size=TRUE_trim={trim_para}_net.weight.CSV")))
}



#############################################################################################
#                                                                                           #
#  Extract detected over expressed genes & signaling pathways & LR pairs in bin50.cellchat  #
#                                                                                           #
#############################################################################################

bin50_over_expr_genes <- unique(bin50.cellchat@var.features$features)
write.csv(bin50_over_expr_genes, 
          file = file.path(output.dir, "overExpr_genes_list.csv"), 
          row.names = FALSE)

pathways.show <- unique(bin50.cellchat@netP$pathways)
pathways.show.df <- data.frame(Pathway = unlist(pathways.show))
write.csv(pathways.show.df,
          file = file.path(output.dir, "Detected_signaling.pathways.csv"),
          row.names = FALSE)

pairLR.show <- unique(bin50.cellchat@net$LRs)
pairLR.show.df <- data.frame(LRpair = unlist(pairLR.show))
write.csv(pairLR.show.df,
          file = file.path(output.dir, "Detected_LRpairs.csv"),
          row.names = FALSE)


#####################################################################
#                                                                   #
#  Visualization of the aggregated cell-cell communication network  #
#                                                                   #
#####################################################################

# name glasbey_dark_2
cluster_levels <- as.character(1:31)
colors.use_pre <- glasbey_dark_2[seq_along(cluster_levels)]
names(colors.use_pre) <- cluster_levels

if (stage_used == "E9.5") {
  colors.use <- colors.use_pre[!names(colors.use_pre) %in% c("27", "31")]
} else if (stage_used == "E13.5") {
  colors.use <- colors.use_pre
}

# visualize the aggregated cell-cell communication network

bin50.cellchat@meta$labels <- as.character(bin50.cellchat@meta$labels)
groupSize <- table(bin50.cellchat@idents)

### Circle plot of [all-in-one] number of interactions and total interaction weights
pdf(file = file.path(output.dir, "all-in-one_Circle plot_number of interaction + interaction weights.pdf"), width = 10, height = 5)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(bin50.cellchat@net$count, color.use = colors.use, vertex.weight = groupSize, edge.width.max = 4, weight.scale = T, label.edge= F, title.name = "Number of interactions_Circle Plot")
netVisual_circle(bin50.cellchat@net$weight, color.use = colors.use, vertex.weight = groupSize, edge.width.max = 4, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength_Circle Plot")
dev.off()

### Heatmap of [all-in-one] number of interactions and total interaction weights
pdf(file = file.path(output.dir, "all-in-one_Heatmap_number of interaction + interaction weights.pdf"), width = 15, height = 8)
ht1 <- netVisual_heatmap(bin50.cellchat, color.use = colors.use, measure = "count", color.heatmap = "Blues", font.size = 11, font.size.title = 13, slot.name = "net", title.name = "Number of interactions between cell clusters")
# the same as using "net": netVisual_heatmap(bin20.cellchat, measure = "count", color.heatmap = "Blues", slot.name = "netP", title.name = "Number of interactions (netP)")

ht2 <- netVisual_heatmap(bin50.cellchat, color.use = colors.use, measure = "weight", color.heatmap = "Blues", font.size = 11, font.size.title = 13, title.name = "Interaction strength between cell clusters")
ht1 + ht2
dev.off()



###############################################################################
#                                                                             #
#  Identify the top-contributing signaling pathways for each Oi-cluster pair  #
#                                                                             #
###############################################################################

for (i in unique(bin50.cellchat@idents)){
  inquired.source <- i
  inquired.netp <- subsetCommunication(bin50.cellchat, sources.use = inquired.source, slot.name = "netP")
  inquired.net <- subsetCommunication(bin50.cellchat, sources.use = inquired.source, slot.name = "net")

  cluster_id <- paste0("source=", i)
  op_prefix <- file.path(output.dir, "all_clusters_as_source_signaling/")
  if (!dir.exists(op_prefix)) {
    dir.create(op_prefix)
  }
  net_loc <- file.path(op_prefix, paste0(cluster_id,"_df.net.csv"))
  netp_loc <- file.path(op_prefix, paste0(cluster_id,"_df.netP.csv"))

  write.csv2(inquired.net, file = net_loc)
  write.csv2(inquired.netp, file = netp_loc)
}


for (i in unique(bin50.cellchat@idents)){
  inquired.target <- i
  inquired.netp <- subsetCommunication(bin50.cellchat, targets.use = inquired.target, slot.name = "netP")
  inquired.net <- subsetCommunication(bin50.cellchat, targets.use = inquired.target, slot.name = "net")

  cluster_id <- paste0("target=", i)
  op_prefix <- file.path(output.dir, "all_clusters_as_target_signaling/")
  if (!dir.exists(op_prefix)) {
    dir.create(op_prefix)
  }
  net_loc <- file.path(op_prefix, paste0(cluster_id,"_df.net.csv"))
  netp_loc <- file.path(op_prefix, paste0(cluster_id,"_df.netP.csv"))

  write.csv2(inquired.net, file = net_loc)
  write.csv2(inquired.netp, file = netp_loc)
}



############################################################################
#                                                                          #
#  Visualization of cell-cell communication network by signaling pathways  #
#                                                                          #
############################################################################

pathways.show <- unique(bin50.cellchat@netP$pathways)

### Circle plot
for (i in 1:length(pathways.show)){
  pathway = pathways.show[i]
  op_prefix <- file.path(output.dir, "pathway_plots/")
  if (!dir.exists(op_prefix)) {
    dir.create(op_prefix)
  }
  pdf(file = file.path(op_prefix, paste0(pathway, "_Circle_plot.pdf")), width = 8, height = 8)
  par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
  print(
    netVisual_aggregate(bin50.cellchat, 
                        color.use = colors.use,
                        vertex.size.max = 8,
                        signaling = pathway, 
                        vertex.weight = NULL,
                        layout = "circle")
  )
  dev.off()
  while (dev.cur() > 1) dev.off()
}

### Spatial plot
if (stage_used == "E9.5") {
  point_size_set = 1
} else if (stage_used == "E13.5") {
  point_size_set = 0.2
}

for (i in 1:length(pathways.show)) {
  pathway <- pathways.show[i]
  
  # Check if this pathway has any interaction (non-zero weights)
  if (sum(bin50.cellchat@netP$prob[,,pathway]) == 0) {
    message("Skipping pathway with no interactions: ", pathway)
    next
  }
  
  op_prefix <- file.path(output.dir, "pathway_plots_new/")
  pdf(
    file = file.path(op_prefix, paste0(pathway, "_Spatial_plot.pdf")),
    width = 8, height = 8
  )
  par(
    mfrow=c(1,1),
    xpd = TRUE
  )
  print(
    netVisual_aggregate(bin50.cellchat, 
                        color.use = colors.use, 
                        signaling = pathway,
                        vertex.size.max = 3,
                        edge.width.max = 2,
                        vertex.label.cex = 0,
                        alpha.image = 0.2,
                        point.size = point_size_set,
                        alpha.edge = 10,
                        layout = "spatial")
  )
  dev.off()
}

### Spatial plot with labels (Fig. 3B)
if (stage_used == "E9.5") {
  alpha_image_set = 0.2
} else if (stage_used == "E13.5") {
  alpha_image_set = 0.15
}

for (i in 1:length(pathways.show)) {
  pathway <- pathways.show[i]
  
  # Check if this pathway has any interaction (non-zero weights)
  if (sum(bin50.cellchat@netP$prob[,,pathway]) == 0) {
    message("Skipping pathway with no interactions: ", pathway)
    next
  }
  
  op_prefix <- file.path(output.dir, "pathway_plots_new.with.label")
  pdf(
    file = file.path(op_prefix, paste0(pathway, "_Spatial_plot_larger.label.pdf")),
    width = 8, height = 8
  )
  par(
    mfrow=c(1,1),
    xpd = TRUE
  )
  print(
    netVisual_aggregate(bin50.cellchat, 
                        color.use = colors.use, 
                        signaling = pathway,
                        vertex.size.max = 3,
                        edge.width.max = 2,
                        vertex.label.cex = 5,                        
                        alpha.image = alpha_image_set,
                        point.size = point_size_set,                        
                        alpha.edge = 10,
                        layout = "spatial")
  )
  dev.off()
}