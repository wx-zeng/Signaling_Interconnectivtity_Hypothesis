# Code for analyzing signaling network interconnectivity among organ primordia in mouse embryonic development

Analysis code accompanying the manuscript "Mid-embryonic conservation cannot be explained solely by high inductive signaling network interconnectivity among organ primordia" by Zeng et al. 2026.
Using the Stereo-seq spatial transcriptome atlas (MOSTA, BGI 2022), the pipeline identifies Organ-identity clusters ("Oi-clusters") as proxies for stage-comparable organ primordia, infers their inductive signaling networks with CellChat, and compares the network interconnectivity between two developmental stages (**E9.5** and **E13.5**, sections E1S1, bin50).

## Overview

Throughout, **Oi-clusters** are stored as `seurat_clusters`. Clusters **27, 29, 30, 31** are E13.5-specific and are excluded from most cross-stage comparisons (E9.5 has clusters 1–26, 28).

## Scripts

**01_Identification_of_Oi-clusters.r**
Load E9.5 and E13.5 bin50 `.h5ad` data, build Seurat objects, merge, `SCTransform` (v2), PCA, CCA integration, then Leiden clustering (npc=50) across resolutions 0.5–1.0 to define Oi-clusters. Resolution **0.7** is used downstream.

**02_Find_DEGs_of_Oi-clusters_FindAllMarkers.r**
Marker genes (DEGs) per Oi-cluster via `FindAllMarkers` (Wilcoxon) on log-normalized RNA counts; export top 30 by `avg_log2FC`.

**03_Find_DEGs_of_Oi-clusters_FindConservedMarkers.r**
Conserved markers per Oi-cluster via `FindConservedMarkers`, grouped by stage (`orig.ident`).

**04_CellChat_analysis.r**
Core CellChat run per stage: build spatial CellChat object, use the mouse Secreted Signaling DB, compute communication probability (`truncatedMean`, distance-aware, `population.size=TRUE`, `trim=0.05`), aggregate networks, export `df.net` / `df.netP` / weights, over-expressed genes, pathways, LR pairs, and produce circle / heatmap / spatial pathway plots (Fig. 3B).

**05_Fig_S2_E13_5_number_of_retained_genes_using_different_trim_parameters.ipynb**
Fig. S2 (E13.5): genes retained per Oi-cluster after intersecting over-expressed genes with non-zero detection fractions at trim = 0.1 / 0.05 / 0.01.

**06_Fig_S2_E9_5_number_of_retained_genes_using_different_trim_parameters.ipynb**
Fig. S2 (E9.5): same as 05 for E9.5.

**07_Fig_3A_pathway_signaling_event_counts.r**
Fig. 3A: per-pathway signaling-event counts (distinct source→target pairs) for each stage, mirrored bar plots.

**08_Fig_S3_Oi-cluster_pathway_contribution.r**
Fig. S3: per-Oi-cluster intra- vs inter-cluster pathway contributions, plotted as stage-dodged bar plots.

**09_Relative_strength_of_inter-Oi-cluster_signaling_and_Ratio_of_inter-Oi-cluster_signaling_complexity.r**
Per-cluster inter/intra signaling strength and complexity ratios; cross-stage comparison with paired (Wilcoxon signed-rank) and unpaired (Mann–Whitney U) tests plus effect sizes.

**10_E9_5_assortativity.r**
E9.5 signaling network: build weighted undirected graph (edge = total interaction strength), spatially anchored layout, node coloring, and assortativity coefficient (`assortnet::assortment.discrete`).

**11_E13_5_assortativity_w_and_wo_E13_5-specific_Oi-clusters.R**
E13.5 assortativity, with a `omit_clusters` switch to include/exclude the E13.5-specific clusters (27/29/30/31).

**12_Assortativity_coefficient_cross-stage_comparison_z-test.r**
Two-sided z-test comparing E9.5 vs E13.5 assortativity coefficients (Newman 2003 jackknife SEs).

**13_Pleiotropy.r**
Pseudo-bulk per Oi-cluster, CPM normalization, Tau tissue-specificity index (Yanai et al. 2005), per-cluster fraction of "pleiotropic" genes (Tau ≤ threshold) across thresholds 0.1–0.9, compared between stages (paired for shared clusters, unpaired for all clusters).

## Requirements

**R** (Seurat v5, SeuratDisk, sctransform, CellChat, anndata, reticulate, BiocNeighbors, igraph, assortnet, rstatix, ggpubr, ggrepel, tidyverse, glue, scCustomize, showtext, viridisLite).

**Python** (notebooks 05–06): `anndata`, `numpy`, `pandas`, `seaborn`, `matplotlib`.
