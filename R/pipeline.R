#' Run Full scRNAbassoon Pipeline
#' 
#' Executes the complete single-cell RNA-seq analysis pipeline including QC,
#' normalization, dimensionality reduction, clustering, cell type annotation,
#' differential expression, pathway enrichment, and trajectory analysis.
#' 
#' @param in_dir Character. Path to Cell Ranger output directory 
#'   (containing barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz) or 
#'   path to HDF5 file (.h5). Required.
#' @param out_dir Character. Output directory path. Default: "./output"
#' @param markers_file Character. Path to custom markers file (CSV or RDS).
#'   If empty, uses \code{\link{default_mouse_markers}}. Default: ""
#' @param sample_map_file Character. Path to sample mapping CSV with columns:
#'   suffix, sample_name, group. Default: ""
#' @param conditions Character. Comma-separated differential expression 
#'   comparisons (e.g., "GroupA_vs_GroupB,GroupC_vs_GroupD"). Default: "Eftud2_vs_Control"
#' @param no_compare Logical. If TRUE, skip differential expression analysis.
#'   Default: FALSE
#' @param n_pcs Integer. Number of principal components to use. Default: 30
#' @param resolution Numeric. Clustering resolution (higher = more clusters).
#'   Default: 0.5
#' @param min_cells Integer. Minimum number of cells expressing a gene.
#'   Default: 3
#' @param min_features Integer. Minimum features (genes) per cell. Default: 200
#' @param max_features Integer. Maximum features per cell. Default: 6000
#' @param max_mito_pct Numeric. Maximum mitochondrial percentage. Default: 15
#' @param plot_format Character. Output format: "pdf", "png", or "both".
#'   Default: "pdf"
#' @param width Numeric. Base plot width in inches. Default: 14
#' @param height Numeric. Base plot height in inches. Default: 12
#' @param use_ggrepel Logical. Use ggrepel for non-overlapping labels.
#'   Default: TRUE
#' @param zip_outputs Logical. Create ZIP bundle of all results. Default: TRUE
#' @param run_trajectory Logical. Run Slingshot trajectory analysis.
#'   Default: TRUE
#' @param use_offline_gmt Logical. Use local GMT files for enrichment.
#'   Default: FALSE
#' @param gmt_dir Character. Directory containing GMT files. Default: ""
#' @param gmt_globs Character. Comma-separated GMT filenames to use. Default: ""
#'   
#' @return NULL (invisible). Analysis results are saved to \code{out_dir}.
#'   
#' @details
#' The pipeline performs the following steps:
#' \enumerate{
#'   \item Create organized output directory structure
#'   \item Load 10X Cell Ranger data
#'   \item Quality control filtering (features, UMI counts, mitochondrial %)
#'   \item Normalization, scaling, and PCA
#'   \item UMAP and t-SNE dimensionality reduction
#'   \item Clustering and marker-based cell type annotation
#'   \item Generate global, per-replicate, and per-sample visualizations
#'   \item Find marker genes for each cluster
#'   \item Differential expression between conditions (optional)
#'   \item Pathway enrichment (GO terms and custom GMT files)
#'   \item Pseudotime trajectory analysis (optional)
#'   \item Save complete Seurat object
#'   \item Create ZIP bundle (optional)
#' }
#' 
#' @section Output Structure:
#' Results are organized into subdirectories:
#' \itemize{
#'   \item \code{QC/}: Quality control plots
#'   \item \code{Global/}: Complete Seurat object and global UMAPs
#'   \item \code{Replicates/}: Per-replicate visualizations
#'   \item \code{Samples/}: Per-condition visualizations
#'   \item \code{Markers/}: Cluster marker genes and heatmaps
#'   \item \code{DE/}: Differential expression results by cell type
#'   \item \code{Enrichment/}: ORA and GSEA pathway analysis
#'   \item \code{Trajectory/}: Pseudotime trajectories (if enabled)
#'   \item \code{Complete_Analysis.zip}: Full results archive (if enabled)
#' }
#' 
#' @section Parallel Processing:
#' Configure parallel processing before running:
#' \preformatted{
#' library(future)
#' plan("multisession", workers = 8)
#' options(future.globals.maxSize = 100 * 1024^3)
#' }
#' 
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' run_bassoon_pipeline(
#'   in_dir = "/data/cellranger/filtered_feature_bc_matrix/",
#'   out_dir = "/results/my_analysis/"
#' )
#' 
#' # Full analysis with custom parameters
#' library(future)
#' plan("multisession", workers = 8)
#' 
#' run_bassoon_pipeline(
#'   in_dir = "/data/cellranger/aggregated/",
#'   out_dir = "/results/experiment1/",
#'   sample_map_file = "/metadata/sample_map.csv",
#'   markers_file = "/markers/custom_markers.rds",
#'   conditions = "Treatment_vs_Control,KO_vs_WT",
#'   n_pcs = 30,
#'   resolution = 0.8,
#'   min_features = 500,
#'   max_features = 7000,
#'   max_mito_pct = 10,
#'   plot_format = "both",
#'   use_offline_gmt = TRUE,
#'   gmt_dir = "/reference/msigdb/mouse/",
#'   gmt_globs = "m2.cp.wikipathways.v2024.1.Mm.symbols.gmt"
#' )
#' }
#' 
#' @seealso 
#' \code{\link{load_sc_data}}, \code{\link{run_qc}}, 
#' \code{\link{preprocess_data}}, \code{\link{annotate_by_markers}},
#' \code{\link{run_de_analysis}}, \code{\link{run_trajectory_analysis}}
#' 
#' @export
#' @importFrom data.table fread
run_bassoon_pipeline <- function(in_dir, out_dir, markers_file = "", sample_map_file = "",
                                 conditions = "Eftud2_vs_Control", no_compare = FALSE,
                                 n_pcs = 30, resolution = 0.5, min_cells = 3,
                                 min_features = 200, max_features = 6000, max_mito_pct = 15,
                                 plot_format = "pdf", width = 14, height = 12,
                                 use_ggrepel = TRUE, zip_outputs = TRUE, run_trajectory = TRUE,
                                 use_offline_gmt = FALSE, gmt_dir = "", gmt_globs = "") {
  
  message("=== Starting scRNAbassoon Pipeline ===")
  
  # Setup
  dirs <- setup_output_dirs(out_dir)
  
  # Markers
  marker_list <- default_mouse_markers
  if (nzchar(markers_file) && file.exists(markers_file)) {
    ext <- tolower(tools::file_ext(markers_file))
    if (ext == "rds") {
      marker_list <- readRDS(markers_file)
    } else {
      dt <- fread(markers_file)
      cn <- tolower(names(dt))
      sc <- intersect(c("set", "group", "celltype"), cn)[1]
      gc <- intersect(c("gene", "symbol", "feature"), cn)[1]
      if (!is.na(sc) && !is.na(gc)) marker_list <- split(dt[[gc]], dt[[sc]])
    }
  }
  
  # Sample Map
  map_df <- NULL
  if (nzchar(sample_map_file) && file.exists(sample_map_file)) {
    map_df <- fread(sample_map_file)
    colnames(map_df) <- tolower(colnames(map_df))
  }
  
  # Execution
  seu <- load_sc_data(in_dir, map_df, min_cells, min_features)
  seu <- run_qc(seu, dirs, min_features, max_features, max_mito_pct, plot_format, width, height)
  seu <- preprocess_data(seu, n_pcs, resolution)
  seu <- annotate_by_markers(seu, marker_list)
  
  plot_global_reductions(seu, dirs, plot_format, width, height, use_ggrepel)
  
  if (length(unique(seu$replicate)) > length(unique(seu$group))) {
    plot_per_replicate(seu, dirs, plot_format, width, height, use_ggrepel)
  }
  
  plot_per_sample(seu, dirs, plot_format, width, height, use_ggrepel)
  find_and_plot_markers(seu, dirs, plot_format, width, height)
  
  if (!no_compare && length(unique(seu$group)) >= 2) {
    run_de_analysis(seu, dirs, conditions, plot_format, use_offline_gmt, gmt_dir, gmt_globs)
  }
  
  if (run_trajectory) run_trajectory_analysis(seu, dirs, plot_format)
  
  saveRDS(seu, file.path(dirs$global, "Complete_Seurat_Object.rds"))
  
  if (zip_outputs) .create_zip_bundle(out_dir)
  
  message("=== Pipeline Complete ===")
}
