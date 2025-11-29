# ==============================================================================
# DATA LOADING & PREPROCESSING
# ==============================================================================

#' Load Single-Cell RNA-seq Data
#' 
#' Loads 10X Genomics Cell Ranger output and creates a Seurat object with
#' optional sample metadata mapping.
#' 
#' @param in_dir Character. Path to Cell Ranger directory (containing 
#'   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) or HDF5 file (.h5)
#' @param sample_map_df Data frame with columns: suffix, sample_name, group.
#'   Maps Cell Ranger sample indices to meaningful names and experimental groups.
#'   Default: NULL
#' @param min_cells Integer. Include features detected in at least this many cells.
#'   Default: 3
#' @param min_features Integer. Include cells with at least this many features.
#'   Default: 200
#'   
#' @return A Seurat object with metadata columns:
#' \describe{
#'   \item{replicate}{Sample/replicate name}
#'   \item{group}{Experimental group/condition}
#' }
#' 
#' @details
#' Cell barcodes from Cell Ranger aggregation contain numeric suffixes (e.g., -1, -2)
#' that identify the sample of origin. This function extracts these suffixes and maps
#' them to sample names and experimental groups via \code{sample_map_df}.
#' 
#' If \code{sample_map_df} is NULL, samples are named "Sample_1", "Sample_2", etc.,
#' and all assigned to group "Unknown".
#' 
#' @examples
#' \dontrun{
#' # Without sample mapping
#' seu <- load_sc_data("/data/cellranger/filtered_feature_bc_matrix/")
#' 
#' # With sample mapping
#' sample_map <- data.frame(
#'   suffix = c("1", "2", "3", "4"),
#'   sample_name = c("Ctrl_Rep1", "Ctrl_Rep2", "Treat_Rep1", "Treat_Rep2"),
#'   group = c("Control", "Control", "Treatment", "Treatment")
#' )
#' seu <- load_sc_data("/data/cellranger/", sample_map)
#' }
#' 
#' @export
#' @importFrom Seurat Read10X Read10X_h5 CreateSeuratObject
load_sc_data <- function(in_dir, sample_map_df = NULL, min_cells = 3, min_features = 200) {
  message("Loading Cell Ranger data...")
  if (file.exists(in_dir) && grepl("\\.h5$", in_dir)) {
    counts <- Read10X_h5(in_dir)
  } else if (dir.exists(in_dir)) {
    counts <- Read10X(in_dir)
  } else {
    stop("Input path must be a directory or .h5 file")
  }
  
  if (is.list(counts)) counts <- counts[[1]] 
  
  seu <- CreateSeuratObject(counts = counts, project = "scRNA_Analysis",
                            min.cells = min_cells, min.features = min_features)
  
  sample_indices <- sub(".*-", "", colnames(seu))
  if (!is.null(sample_map_df)) {
    suffix_to_name <- setNames(sample_map_df$sample_name, sample_map_df$suffix)
    suffix_to_group <- setNames(sample_map_df$group, sample_map_df$suffix)
    seu$replicate <- unname(suffix_to_name[sample_indices])
    seu$group <- unname(suffix_to_group[sample_indices])
    seu$replicate <- ifelse(is.na(seu$replicate), paste0("Sample_", sample_indices), seu$replicate)
    seu$group <- ifelse(is.na(seu$group), "Unknown", seu$group)
  } else {
    seu$replicate <- paste0("Sample_", sample_indices)
    seu$group <- "Unknown"
  }
  seu$group <- factor(seu$group)
  return(seu)
}

#' Run Quality Control and Filtering
#' 
#' Performs QC metric calculation, generates diagnostic plots, and filters cells
#' based on feature count and mitochondrial percentage thresholds.
#' 
#' @param seu Seurat object from \code{\link{load_sc_data}}
#' @param dirs Named list of output directories from \code{\link{setup_output_dirs}}
#' @param min_features Integer. Minimum features per cell threshold
#' @param max_features Integer. Maximum features per cell threshold
#' @param max_mito_pct Numeric. Maximum mitochondrial percentage threshold
#' @param plot_format Character. "pdf", "png", or "both"
#' @param width Numeric. Plot width in inches
#' @param height Numeric. Plot height in inches
#'   
#' @return Filtered Seurat object with QC metrics in metadata
#' 
#' @details
#' Calculates and plots:
#' \itemize{
#'   \item nFeature_RNA: Number of genes detected per cell
#'   \item nCount_RNA: Total UMI counts per cell
#'   \item percent.mt: Percentage of reads mapping to mitochondrial genes
#' }
#' 
#' Mitochondrial genes are identified by pattern "^mt-" (mouse) or "^MT-" (human).
#' 
#' QC plots saved to \code{dirs$qc}:
#' \itemize{
#'   \item QC_Metrics_Violin.pdf: Violin plots of QC metrics by group
#'   \item QC_Count_vs_Feature.pdf: Scatter plot of sequencing depth vs complexity
#' }
#' 
#' @examples
#' \dontrun{
#' dirs <- setup_output_dirs("/results/")
#' seu <- load_sc_data("/data/cellranger/")
#' seu <- run_qc(seu, dirs, 
#'               min_features = 200, 
#'               max_features = 6000, 
#'               max_mito_pct = 15,
#'               plot_format = "pdf", 
#'               width = 14, 
#'               height = 12)
#' }
#' 
#' @export
#' @importFrom Seurat PercentageFeatureSet VlnPlot FeatureScatter subset
#' @importFrom ggplot2 ggtitle
run_qc <- function(seu, dirs, min_features, max_features, max_mito_pct, plot_format, width, height) {
  message("Running QC...")
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-|^MT-")
  
  p_qc <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  group.by = "group", pt.size = 0, ncol = 3) + .theme_bw_custom()
  .save_plot(p_qc, "QC_Metrics_Violin", dirs$qc, width = 16, height = 6, format = plot_format)
  
  p_scatter <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
                              group.by = "group") + 
    .theme_bw_custom() + ggtitle("Sequencing Depth vs. Complexity")
  .save_plot(p_scatter, "QC_Count_vs_Feature", dirs$qc, width = 10, height = 8, format = plot_format)
  
  seu <- subset(seu, subset = nFeature_RNA >= min_features & 
                  nFeature_RNA <= max_features & percent.mt <= max_mito_pct)
  message(sprintf("After QC: %d cells retained", ncol(seu)))
  return(seu)
}

#' Preprocess scRNA-seq Data
#' 
#' Performs normalization, feature selection, scaling, PCA, clustering,
#' and dimensionality reduction (UMAP and t-SNE).
#' 
#' @param seu Filtered Seurat object from \code{\link{run_qc}}
#' @param n_pcs Integer. Number of principal components to compute and use.
#'   Default: 30
#' @param resolution Numeric. Clustering resolution for \code{FindClusters}.
#'   Higher values yield more clusters. Default: 0.5
#'   
#' @return Preprocessed Seurat object with:
#' \itemize{
#'   \item Normalized and scaled data
#'   \item PCA, UMAP, and t-SNE embeddings
#'   \item Cluster assignments in \code{seurat_clusters}
#' }
#' 
#' @details
#' Processing steps:
#' \enumerate{
#'   \item Log-normalization (scale factor 10,000)
#'   \item Identify 2,000 highly variable features (VST method)
#'   \item Scale data (z-score transformation)
#'   \item PCA on variable features
#'   \item Shared nearest neighbor (SNN) graph construction
#'   \item Louvain clustering
#'   \item UMAP dimensionality reduction
#'   \item t-SNE dimensionality reduction
#' }
#' 
#' @examples
#' \dontrun{
#' seu <- load_sc_data("/data/cellranger/")
#' seu <- run_qc(seu, dirs, 200, 6000, 15, "pdf", 14, 12)
#' seu <- preprocess_data(seu, n_pcs = 30, resolution = 0.8)
#' }
#' 
#' @export
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA
#' @importFrom Seurat FindNeighbors FindClusters RunUMAP RunTSNE
preprocess_data <- function(seu, n_pcs = 30, resolution = 0.5) {
  message("Normalizing and scaling...")
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu)
  message("Running PCA, UMAP, tSNE...")
  seu <- RunPCA(seu, npcs = n_pcs, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:n_pcs)
  seu <- FindClusters(seu, resolution = resolution)
  seu <- RunUMAP(seu, dims = 1:n_pcs, verbose = FALSE)
  seu <- RunTSNE(seu, dims = 1:n_pcs, verbose = FALSE)
  return(seu)
}

#' Annotate Cells by Marker Genes
#' 
#' Assigns cell type labels based on marker gene expression patterns using
#' a simple scoring approach.
#' 
#' @param seu Preprocessed Seurat object from \code{\link{preprocess_data}}
#' @param marker_list Named list where names are cell types and values are
#'   character vectors of marker gene symbols
#'   
#' @return Seurat object with new metadata columns:
#' \describe{
#'   \item{cell_type_marker}{Per-cell annotation based on marker expression}
#'   \item{all_cell_types}{Final cell type labels (consensus per cluster)}
#' }
#' 
#' @details
#' Annotation algorithm:
#' \enumerate{
#'   \item For each cell, count how many markers are detected (count > 0) for each cell type
#'   \item Assign cell to type with highest marker count
#'   \item If no markers detected, label as "Unknown"
#'   \item If multiple types tie, label as "Ambiguous"
#'   \item Refine by cluster: assign each cluster the most common cell type within it
#' }
#' 
#' The final \code{all_cell_types} column is used as the primary cell type identity
#' throughout downstream analysis.
#' 
#' @examples
#' \dontrun{
#' # Using built-in markers
#' seu <- annotate_by_markers(seu, default_mouse_markers)
#' 
#' # Using custom markers
#' my_markers <- list(
#'   "T_cells" = c("Cd3d", "Cd3e"),
#'   "B_cells" = c("Cd79a", "Ms4a1"),
#'   "Macrophages" = c("Cd68", "Itgam")
#' )
#' seu <- annotate_by_markers(seu, my_markers)
#' 
#' # Check results
#' table(seu$all_cell_types)
#' }
#' 
#' @export
#' @importFrom Seurat DefaultAssay LayerData WhichCells Idents
#' @importFrom Matrix colSums
annotate_by_markers <- function(seu, marker_list) {
  message("Annotating cells...")
  DefaultAssay(seu) <- "RNA"
  mat <- LayerData(seu[["RNA"]], layer = "counts")
  
  if (is.null(mat) || !length(marker_list)) {
    seu$all_cell_types <- "Unknown"
    return(seu)
  }
  
  sets <- names(marker_list)
  score <- matrix(0L, nrow = ncol(mat), ncol = length(sets), dimnames = list(colnames(mat), sets))
  
  for (ct in sets) {
    present <- intersect(marker_list[[ct]], rownames(mat))
    if (length(present) > 0) score[, ct] <- as.integer(Matrix::colSums(mat[present, , drop = FALSE] > 0))
  }
  
  best <- apply(score, 1, function(x) {
    if (all(x == 0)) return("Unknown")
    w <- which(x == max(x))
    if (length(w) == 1) colnames(score)[w] else "Ambiguous"
  })
  
  seu$cell_type_marker <- unname(best)
  seu$all_cell_types <- seu$cell_type_marker
  
  if ("seurat_clusters" %in% names(seu@meta.data)) {
    cluster_labels <- sapply(levels(seu$seurat_clusters), function(cl) {
      cells <- WhichCells(seu, idents = cl)
      tab <- table(seu$all_cell_types[cells])
      tab <- tab[!names(tab) %in% c("Unknown", "Ambiguous")]
      if (length(tab) == 0) return("Unknown")
      names(tab)[which.max(tab)]
    })
    final_labels <- cluster_labels[as.character(seu$seurat_clusters)]
    names(final_labels) <- colnames(seu)
    seu$all_cell_types <- factor(final_labels)
  }
  Idents(seu) <- "all_cell_types"
  return(seu)
}

# ==============================================================================
# VISUALIZATION FUNCTIONS
# ==============================================================================

#' Plot Global Dimensionality Reductions
#' 
#' Generates UMAP plots colored by clusters, cell types, and experimental groups
#' for the complete dataset.
#' 
#' @param seu Annotated Seurat object
#' @param dirs Named list of output directories
#' @param format Character. "pdf", "png", or "both"
#' @param width Numeric. Plot width in inches
#' @param height Numeric. Plot height in inches
#' @param use_ggrepel Logical. Use ggrepel for non-overlapping labels
#'   
#' @return NULL (invisible). Saves plots to \code{dirs$global}:
#' \itemize{
#'   \item Global_UMAP_Clusters.pdf: Colored by Seurat clusters
#'   \item Global_UMAP_CellTypes.pdf: Colored by annotated cell types
#'   \item Global_UMAP_Groups.pdf: Colored by experimental groups
#' }
#' 
#' @examples
#' \dontrun{
#' plot_global_reductions(seu, dirs, "pdf", 14, 12, use_ggrepel = TRUE)
#' }
#' 
#' @export
#' @importFrom Seurat DimPlot
plot_global_reductions <- function(seu, dirs, format, width, height, use_ggrepel) {
  p1 <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = use_ggrepel) + .theme_bw_custom()
  .save_plot(p1, "Global_UMAP_Clusters", dirs$global, width, height, format)
  
  p2 <- DimPlot(seu, reduction = "umap", group.by = "all_cell_types", label = TRUE, repel = use_ggrepel) + .theme_bw_custom()
  .save_plot(p2, "Global_UMAP_CellTypes", dirs$global, width, height, format)
  
  p3 <- DimPlot(seu, reduction = "umap", group.by = "group", label = FALSE) + .theme_bw_custom()
  .save_plot(p3, "Global_UMAP_Groups", dirs$global, width, height, format)
}

#' Plot Per-Replicate Visualizations
#' 
#' Creates separate UMAP plots for each biological replicate.
#' 
#' @param seu Annotated Seurat object
#' @param dirs Named list of output directories
#' @param format Character. "pdf", "png", or "both"
#' @param width Numeric. Plot width in inches
#' @param height Numeric. Plot height in inches
#' @param use_ggrepel Logical. Use ggrepel for labels
#'   
#' @return NULL (invisible). Saves plots to subdirectories within \code{dirs$replicates}
#' 
#' @details
#' Skips replicates with fewer than 50 cells. Each replicate gets its own subdirectory
#' with annotated UMAP plot.
#' 
#' @examples
#' \dontrun{
#' plot_per_replicate(seu, dirs, "pdf", 14, 12, use_ggrepel = TRUE)
#' }
#' 
#' @export
#' @importFrom Seurat DimPlot subset
#' @importFrom ggplot2 ggtitle
plot_per_replicate <- function(seu, dirs, format, width, height, use_ggrepel) {
  reps <- unique(seu$replicate)
  for (r in reps) {
    r_dir <- file.path(dirs$replicates, gsub("[^[:alnum:]]", "_", r))
    dir.create(r_dir, showWarnings = FALSE, recursive = TRUE)
    seu_sub <- subset(seu, replicate == r)
    if (ncol(seu_sub) < 50) next
    
    p <- DimPlot(seu_sub, reduction = "umap", group.by = "all_cell_types", label = TRUE, repel = use_ggrepel) +
      .theme_bw_custom() + ggtitle(paste("UMAP -", r))
    .save_plot(p, "UMAP_Annotated", r_dir, width, height, format)
  }
}

#' Plot Per-Sample Visualizations
#' 
#' Creates UMAP plots split by replicate for each experimental group.
#' 
#' @param seu Annotated Seurat object
#' @param dirs Named list of output directories
#' @param format Character. "pdf", "png", or "both"
#' @param width Numeric. Plot width in inches
#' @param height Numeric. Plot height in inches
#' @param use_ggrepel Logical. Use ggrepel for labels
#'   
#' @return NULL (invisible). Saves plots to subdirectories within \code{dirs$samples}
#' 
#' @details
#' Skips groups with fewer than 30 cells. Creates faceted UMAP plots showing
#' all replicates within each experimental group.
#' 
#' @examples
#' \dontrun{
#' plot_per_sample(seu, dirs, "pdf", 14, 12, use_ggrepel = TRUE)
#' }
#' 
#' @export
#' @importFrom Seurat DimPlot subset
#' @importFrom ggplot2 ggtitle
plot_per_sample <- function(seu, dirs, format, width, height, use_ggrepel) {
  grps <- unique(seu$group)
  grps <- grps[grps != "Unknown"]
  for (g in grps) {
    g_dir <- file.path(dirs$samples, gsub("[^[:alnum:]]", "_", g))
    dir.create(g_dir, showWarnings = FALSE, recursive = TRUE)
    seu_sub <- subset(seu, group == g)
    if (ncol(seu_sub) < 30) next
    
    p <- DimPlot(seu_sub, reduction = "umap", group.by = "all_cell_types", split.by = "replicate", label = FALSE) +
      .theme_bw_custom() + ggtitle(paste("UMAP -", g))
    dims <- .auto_dims(length(unique(seu_sub$replicate)), width, height)
    .save_plot(p, "UMAP_Annotated_by_Replicate", g_dir, dims$width * 1.8, dims$height, format)
  }
}

#' Find and Plot Marker Genes
#' 
#' Identifies marker genes for each cell type and generates a heatmap of
#' top markers.
#' 
#' @param seu Annotated Seurat object
#' @param dirs Named list of output directories
#' @param format Character. "pdf", "png", or "both"
#' @param width Numeric. Base plot width
#' @param height Numeric. Base plot height
#'   
#' @return Data frame of marker genes with columns:
#' \describe{
#'   \item{cluster}{Cell type identity}
#'   \item{gene}{Gene symbol}
#'   \item{avg_log2FC}{Average log2 fold-change}
#'   \item{pct.1}{Percent cells expressing in cluster}
#'   \item{pct.2}{Percent cells expressing in other clusters}
#'   \item{p_val_adj}{Adjusted p-value}
#' }
#' 
#' @details
#' Uses Wilcoxon rank-sum test to identify genes with:
#' \itemize{
#'   \item Positive fold-change only
#'   \item min.pct = 0.1 (expressed in ≥10% of cells)
#'   \item logfc.threshold = 0.25
#' }
#' 
#' Outputs saved to \code{dirs$markers}:
#' \itemize{
#'   \item All_Cluster_Markers.csv: Complete marker table
#'   \item Heatmap_Top10_Markers.pdf: Heatmap of top 10 markers per cell type
#' }
#' 
#' @examples
#' \dontrun{
#' markers <- find_and_plot_markers(seu, dirs, "pdf", 14, 12)
#' head(markers)
#' }
#' 
#' @export
#' @importFrom Seurat FindAllMarkers ScaleData DoHeatmap Idents
#' @importFrom data.table fwrite
#' @importFrom dplyr %>% group_by slice_max
#' @importFrom ggplot2 ggtitle theme element_blank
find_and_plot_markers <- function(seu, dirs, format, width, height) {
  message("Finding markers...")
  Idents(seu) <- "all_cell_types"
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  fwrite(markers, file.path(dirs$markers, "All_Cluster_Markers.csv"))
  
  top10 <- markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
  if (nrow(top10) > 0) {
    feats <- intersect(unique(top10$gene), rownames(seu))
    if (length(feats) >= 2) {
      seu <- ScaleData(seu, features = feats, verbose = FALSE)
      # Calculate dynamic dimensions
      n_types <- length(unique(seu$all_cell_types))
      n_feats <- length(feats)
      final_w <- min(max(16, n_types * 0.8 + 8), 40) * 1.8
      final_h <- min(max(12, n_feats * 0.15 + 6), 50) * 2
      
      p <- DoHeatmap(seu, features = feats, group.by = "all_cell_types", size = 4.5) + 
        .theme_bw_custom() + ggtitle("Top 10 Markers") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      .save_plot(p, "Heatmap_Top10_Markers", dirs$markers, final_w, final_h, format)
    }
  }
  return(markers)
}

# ==============================================================================
# ENRICHMENT & DIFFERENTIAL EXPRESSION
# ==============================================================================

#' Internal Enrichment Analysis Helper
#' 
#' Performs over-representation analysis (ORA) and gene set enrichment analysis
#' (GSEA) on differential expression results.
#' 
#' @param de_df Data frame of differential expression results
#' @param all_genes_universe Character vector of all genes in the dataset
#' @param outdir Output directory path
#' @param file_stem Prefix for output files
#' @param use_offline Logical. Use offline GMT files
#' @param gmt_cols List of GMT collections (term2gene data frames)
#' @param format Character. "pdf", "png", or "both"
#'   
#' @return NULL (invisible). Saves enrichment results and plots
#' 
#' @details
#' Performs:
#' \itemize{
#'   \item ORA on upregulated genes (log2FC > 0.25, adj.p < 0.1)
#'   \item ORA on downregulated genes (log2FC < -0.25, adj.p < 0.1)
#'   \item GSEA on ranked gene list (all genes by log2FC)
#'   \item GO enrichment (BP and MF ontologies)
#'   \item Custom GMT enrichment (if provided)
#' }
#' 
#' Converts gene symbols to Entrez IDs for GO analysis using org.Mm.eg.db.
#' 
#' @keywords internal
#' @importFrom clusterProfiler enrichGO gseGO enricher GSEA read.gmt
#' @importFrom enrichplot dotplot barplot
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom AnnotationDbi mapIds
#' @importFrom data.table fwrite
#' @importFrom stats na.omit setNames
.run_enrichment_internal <- function(de_df, all_genes_universe, outdir, file_stem, use_offline, gmt_cols, format) {
  if (is.null(de_df) || nrow(de_df) == 0) return()
  de_df <- de_df[!is.na(de_df$p_val_adj) & !is.na(de_df$avg_log2FC), , drop = FALSE]
  
  # Symbol to Entrez conversion
  .sym2ent <- function(s) {
    s <- unique(as.character(s))
    if (!length(s)) return(character(0))
    suppressMessages(tryCatch(AnnotationDbi::mapIds(org.Mm.eg.db, keys = s, column = "ENTREZID", 
                                                    keytype = "SYMBOL", multiVals = "first"), 
                              error = function(e) setNames(rep(NA, length(s)), s))) %>% stats::na.omit()
  }
  
  universe_entrez <- .sym2ent(all_genes_universe)
  up_syms <- rownames(de_df)[de_df$avg_log2FC > 0.25 & de_df$p_val_adj < 0.1]
  down_syms <- rownames(de_df)[de_df$avg_log2FC < -0.25 & de_df$p_val_adj < 0.1]
  up_entrez <- .sym2ent(up_syms)
  down_entrez <- .sym2ent(down_syms)
  
  # Ranked lists for GSEA
  ranked_sym <- sort(setNames(de_df$avg_log2FC, rownames(de_df)), decreasing = TRUE)
  ent_map <- .sym2ent(names(ranked_sym))
  ranked_entrez <- ranked_sym[names(ent_map)]
  names(ranked_entrez) <- ent_map
  ranked_entrez <- ranked_entrez[!is.na(names(ranked_entrez))]
  ranked_entrez <- sort(tapply(ranked_entrez, names(ranked_entrez), max), decreasing = TRUE)

  # Output directories
  ora_dir <- file.path(outdir, "ORA")
  gsea_dir <- file.path(outdir, "GSEA")
  dir.create(ora_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Helper to save enrichment plots
  .save_enr <- function(res, name, dir, title) {
    if (is.null(res) || nrow(as.data.frame(res)) == 0) return()
    tryCatch({
      p <- dotplot(res, showCategory = 10) + ggtitle(title)
      .save_plot(p, paste0(name, "_dotplot"), dir, 12, 8, format)
      p2 <- barplot(res, showCategory = 10) + ggtitle(title)
      .save_plot(p2, paste0(name, "_barplot"), dir, 12, 8, format)
    }, error = function(e) NULL)
  }
  
  # ORA (GO)
  for (ont in c("BP", "MF")) {
    if (length(up_entrez) > 5) {
      ego <- tryCatch(enrichGO(up_entrez, universe = universe_entrez, OrgDb = org.Mm.eg.db, ont = ont, keyType = "ENTREZID"), error = function(e) NULL)
      if (!is.null(ego)) { fwrite(as.data.frame(ego), file.path(ora_dir, paste0(file_stem, "_UP_GO_", ont, ".csv")))
      .save_enr(ego, paste0(file_stem, "_UP_GO_", ont), ora_dir, paste("ORA UP", ont)) }
    }
    if (length(down_entrez) > 5) {
      ego <- tryCatch(enrichGO(down_entrez, universe = universe_entrez, OrgDb = org.Mm.eg.db, ont = ont, keyType = "ENTREZID"), error = function(e) NULL)
      if (!is.null(ego)) { fwrite(as.data.frame(ego), file.path(ora_dir, paste0(file_stem, "_DOWN_GO_", ont, ".csv")))
      .save_enr(ego, paste0(file_stem, "_DOWN_GO_", ont), ora_dir, paste("ORA DOWN", ont)) }
    }
  }
  
  # GMT ORA
  if (use_offline && length(gmt_cols) > 0) {
    for (nm in names(gmt_cols)) {
      t2g <- gmt_cols[[nm]]$t2g
      if (length(up_syms) > 5) {
        enr <- tryCatch(clusterProfiler::enricher(up_syms, TERM2GENE = t2g), error = function(e) NULL)
        if (!is.null(enr)) { fwrite(as.data.frame(enr), file.path(ora_dir, paste0(file_stem, "_UP_GMT_", nm, ".csv")))
        .save_enr(enr, paste0(file_stem, "_UP_GMT_", nm), ora_dir, paste("ORA UP", nm)) }
      }
      if (length(down_syms) > 5) {
        enr <- tryCatch(clusterProfiler::enricher(down_syms, TERM2GENE = t2g), error = function(e) NULL)
        if (!is.null(enr)) { fwrite(as.data.frame(enr), file.path(ora_dir, paste0(file_stem, "_DOWN_GMT_", nm, ".csv")))
        .save_enr(enr, paste0(file_stem, "_DOWN_GMT_", nm), ora_dir, paste("ORA DOWN", nm)) }
      }
    }
  }
  
  # GSEA (GO)
  if (length(ranked_entrez) > 20) {
    for (ont in c("BP", "MF")) {
      ggo <- tryCatch(gseGO(ranked_entrez, OrgDb = org.Mm.eg.db, ont = ont, keyType = "ENTREZID", seed = TRUE), error = function(e) NULL)
      if (!is.null(ggo)) {
        fwrite(as.data.frame(ggo), file.path(gsea_dir, paste0(file_stem, "_GSEA_GO_", ont, ".csv")))
        .save_enr(ggo, paste0(file_stem, "_GSEA_GO_", ont), gsea_dir, paste("GSEA", ont))
      }
    }
  }
  
  # GSEA (GMT)
  if (use_offline && length(gmt_cols) > 0 && length(ranked_sym) > 20) {
    for (nm in names(gmt_cols)) {
      t2g <- gmt_cols[[nm]]$t2g
      gmt_res <- tryCatch(clusterProfiler::GSEA(ranked_sym, TERM2GENE = t2g, seed = TRUE), error = function(e) NULL)
      if (!is.null(gmt_res)) {
        fwrite(as.data.frame(gmt_res), file.path(gsea_dir, paste0(file_stem, "_GSEA_GMT_", nm, ".csv")))
        .save_enr(gmt_res, paste0(file_stem, "_GSEA_GMT_", nm), gsea_dir, paste("GSEA", nm))
      }
    }
  }
}

#' Run Differential Expression Analysis
#' 
#' Performs pairwise differential expression testing between experimental groups
#' for each cell type, followed by pathway enrichment analysis.
#' 
#' @param seu Annotated Seurat object
#' @param dirs Named list of output directories
#' @param conditions_str Character. Comma-separated comparisons (e.g., "GroupA_vs_GroupB,GroupC_vs_GroupD")
#' @param format Character. "pdf", "png", or "both"
#' @param use_offline_gmt Logical. Use local GMT files for enrichment
#' @param gmt_dir Character. Directory containing GMT files
#' @param gmt_globs Character. Comma-separated GMT filenames
#'   
#' @return NULL (invisible). Saves results to:
#' \itemize{
#'   \item \code{dirs$de}: DE results and volcano plots per cell type
#'   \item \code{dirs$enrichment}: ORA and GSEA results per cell type
#' }
#' 
#' @details
#' For each comparison and cell type:
#' \enumerate{
#'   \item Subset to cells from comparison groups
#'   \item Require ≥10 cells per group for testing
#'   \item Perform Wilcoxon rank-sum test (logfc.threshold = 0)
#'   \item Generate volcano plots
#'   \item Run GO enrichment (BP, MF)
#'   \item Run GMT enrichment (if enabled)
#'   \item Perform GSEA on ranked genes
#' }
#' 
#' Output structure per comparison:
#' \preformatted{
#' DE/GroupA_vs_GroupB/
#'   CellType1/
#'     DE_results.csv
#'     Volcano_Plot.pdf
#'   CellType2/
#'     ...
#' 
#' Enrichment/GroupA_vs_GroupB/
#'   CellType1/
#'     ORA/
#'       CellType1_GroupA_vs_GroupB_UP_GO_BP.csv
#'       CellType1_GroupA_vs_GroupB_DOWN_GO_BP.csv
#'       ...
#'     GSEA/
#'       CellType1_GroupA_vs_GroupB_GSEA_GO_BP.csv
#'       ...
#' }
#' 
#' @section GMT Files:
#' Download GMT files from MSigDB (\url{https://www.gsea-msigdb.org/gsea/msigdb/}):
#' \itemize{
#'   \item m2.cp.wikipathways: WikiPathways
#'   \item m2.cp.kegg: KEGG pathways
#'   \item m5.go.bp: GO Biological Process
#'   \item h.all: Hallmark gene sets
#' }
#' 
#' @examples
#' \dontrun{
#' # Single comparison without GMT enrichment
#' run_de_analysis(seu, dirs, "Treatment_vs_Control", "pdf", FALSE, "", "")
#' 
#' # Multiple comparisons with GMT enrichment
#' run_de_analysis(
#'   seu, dirs,
#'   conditions_str = "Mutant_vs_WT,Rescue_vs_Mutant",
#'   format = "both",
#'   use_offline_gmt = TRUE,
#'   gmt_dir = "/reference/msigdb/mouse/",
#'   gmt_globs = "m2.cp.wikipathways.v2024.1.Mm.symbols.gmt,m2.cp.kegg.v2024.1.Mm.symbols.gmt"
#' )
#' }
#' 
#' @export
#' @importFrom Seurat FindMarkers subset
#' @importFrom data.table fwrite
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom clusterProfiler read.gmt
#' @importFrom ggplot2 ggtitle
run_de_analysis <- function(seu, dirs, conditions_str, format, use_offline_gmt, gmt_dir, gmt_globs) {
  message("Running Differential Expression...")
  comparisons <- strsplit(conditions_str, ",")[[1]]
  
  # Load GMTs if needed
  gmt_cols <- list()
  if (use_offline_gmt && dir.exists(gmt_dir)) {
    gfiles <- trimws(unlist(strsplit(gmt_globs, ",")))
    for (gf in gfiles) {
      gp <- file.path(gmt_dir, gf)
      if (file.exists(gp)) {
        df <- tryCatch(clusterProfiler::read.gmt(gp), error = function(e) NULL)
        if (!is.null(df)) gmt_cols[[gf]] <- list(t2g = df[, c("term", "gene")])
      }
    }
  }
  
  for (comp in comparisons) {
    parts <- strsplit(comp, "_vs_")[[1]]
    if (length(parts) != 2) next
    g1 <- parts[1]; g2 <- parts[2]
    message(sprintf("  %s vs %s", g1, g2))
    
    comp_dir <- file.path(dirs$de, comp)
    dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
    
    cts <- unique(seu$all_cell_types)
    cts <- cts[!cts %in% c("Unknown", "Ambiguous")]
    
    for (ct in cts) {
      seu_sub <- subset(seu, all_cell_types == ct)
      if (sum(seu_sub$group == g1) < 10 || sum(seu_sub$group == g2) < 10) next
      
      tryCatch({
        # NOTE: logfc.threshold = 0 to capture ALL genes for GSEA
        de_res <- FindMarkers(seu_sub, ident.1 = g1, ident.2 = g2, group.by = "group",
                              test.use = "wilcox", logfc.threshold = 0, min.pct = 0.05)
        
        if (nrow(de_res) > 0) {
          de_res$gene <- rownames(de_res)
          ct_safe <- gsub("[^[:alnum:]]", "_", ct)
          ct_dir <- file.path(comp_dir, ct_safe)
          dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
          fwrite(de_res, file.path(ct_dir, "DE_results.csv"))
          
          # Volcano plot
          up_n <- sum(de_res$avg_log2FC > 0.25 & de_res$p_val_adj < 0.1, na.rm = TRUE)
          down_n <- sum(de_res$avg_log2FC < -0.25 & de_res$p_val_adj < 0.1, na.rm = TRUE)
          
          p_vol <- EnhancedVolcano(de_res, lab = rownames(de_res), x = 'avg_log2FC', y = 'p_val_adj',
                                   title = paste(ct, comp), subtitle = sprintf("Up: %d | Down: %d", up_n, down_n),
                                   pCutoff = 0.1, FCcutoff = 0.25, labSize = 3) + .theme_bw_custom()
          .save_plot(p_vol, "Volcano_Plot", ct_dir, 14, 12, format)
          
          # Enrichment analysis
          enrich_dir <- file.path(dirs$enrichment, comp, ct_safe)
          dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)
          .run_enrichment_internal(de_res, rownames(seu_sub), enrich_dir, paste(ct_safe, comp, sep="_"),
                                   use_offline_gmt, gmt_cols, format)
        }
      }, error = function(e) message(sprintf("DE failed for %s: %s", ct, e$message)))
    }
  }
}

# ==============================================================================
# TRAJECTORY ANALYSIS
# ==============================================================================

#' Run Pseudotime Trajectory Analysis
#' 
#' Performs trajectory inference using Slingshot to identify developmental
#' or differentiation lineages within each experimental group.
#' 
#' @param seu Annotated Seurat object
#' @param dirs Named list of output directories
#' @param format Character. "pdf", "png", or "both"
#'   
#' @return NULL (invisible). Saves results to \code{dirs$trajectory}:
#' \itemize{
#'   \item SCE_<GroupName>.rds: SingleCellExperiment with lineages
#'   \item Trajectory_<GroupName>.pdf: UMAP with trajectory curves
#' }
#' 
#' @details
#' For each experimental group:
#' \enumerate{
#'   \item Convert Seurat to SingleCellExperiment
#'   \item Transfer PCA and UMAP embeddings
#'   \item Run Slingshot trajectory inference using cell type labels as clusters
#'   \item Plot trajectories on UMAP
#'   \item Save SCE object for downstream analysis
#' }
#' 
#' Skips groups with fewer than 30 cells.
#' 
#' @section Downstream Analysis:
#' Load saved SCE objects for further exploration:
#' \preformatted{
#' library(slingshot)
#' sce <- readRDS("Trajectory/SCE_Treatment.rds")
#' 
#' # Extract pseudotime
#' pseudotime <- slingPseudotime(sce)
#' 
#' # Get lineages
#' curves <- slingCurves(sce)
#' 
#' # Analyze genes along pseudotime
#' library(tradeSeq)
#' # ... downstream trajectory analysis
#' }
#' 
#' @examples
#' \dontrun{
#' run_trajectory_analysis(seu, dirs, format = "pdf")
#' 
#' # Load results
#' sce <- readRDS("output/Trajectory/SCE_Control.rds")
#' pt <- slingPseudotime(sce)
#' head(pt)
#' }
#' 
#' @export
#' @importFrom SeuratObject as.SingleCellExperiment
#' @importFrom Seurat GetAssayData Embeddings subset
#' @importFrom SingleCellExperiment logcounts reducedDim colLabels
#' @importFrom slingshot slingshot SlingshotDataSet
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette png pdf dev.off
#' @importFrom graphics plot lines
run_trajectory_analysis <- function(seu, dirs, format) {
  message("Running Trajectory...")
  grps <- unique(seu$group)
  grps <- grps[grps != "Unknown"]
  
  for (g in grps) {
    seu_sub <- subset(seu, group == g)
    if (ncol(seu_sub) < 30) next
    tryCatch({
      sce <- as.SingleCellExperiment(seu_sub)
      logcounts(sce) <- as.matrix(GetAssayData(seu_sub, layer = "data"))
      reducedDim(sce, "PCA") <- Embeddings(seu_sub, "pca")
      reducedDim(sce, "UMAP") <- Embeddings(seu_sub, "umap")
      colLabels(sce) <- seu_sub$all_cell_types
      sce <- slingshot(sce, clusterLabels = colLabels(sce), reducedDim = 'PCA')
      
      saveRDS(sce, file.path(dirs$trajectory, paste0("SCE_", gsub("[^[:alnum:]]", "_", g), ".rds")))
      
      # Plotting
      cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(unique(colLabels(sce))))
      plot_file <- file.path(dirs$trajectory, paste0("Trajectory_", gsub("[^[:alnum:]]", "_", g)))
      
      if (format %in% c("png", "both")) {
        png(paste0(plot_file, ".png"), width = 12, height = 10, units = "in", res = 300)
        plot(reducedDim(sce, "UMAP"), col = cols[as.factor(colLabels(sce))], pch=16, cex=0.5, main=g)
        lines(SlingshotDataSet(sce), lwd=2, col='black')
        dev.off()
      }
      if (format %in% c("pdf", "both")) {
        pdf(paste0(plot_file, ".pdf"), width = 12, height = 10)
        plot(reducedDim(sce, "UMAP"), col = cols[as.factor(colLabels(sce))], pch=16, cex=0.5, main=g)
        lines(SlingshotDataSet(sce), lwd=2, col='black')
        dev.off()
      }
    }, error = function(e) message("Trajectory failed: ", e$message))
  }
}
