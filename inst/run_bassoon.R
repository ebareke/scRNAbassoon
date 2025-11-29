#!/usr/bin/env Rscript

#' @title scRNAbassoon Command-Line Interface
#' @description Command-line wrapper for running the complete scRNAbassoon
#'   single-cell RNA-seq analysis pipeline with customizable parameters.
#' @details
#' This script provides a user-friendly command-line interface to the
#' scRNAbassoon package. It handles argument parsing, parallel processing
#' configuration, and executes the full analysis pipeline.
#'
#' Usage:
#'   Rscript run_bassoon.R --in_dir <path> [options]
#'
#' Required argument:
#'   --in_dir        Input Cell Ranger directory or H5 file
#'
#' For full help:
#'   Rscript run_bassoon.R --help

# ==============================================================================
# LIBRARY LOADING
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(scRNAbassoon)
  library(future)
})

# ==============================================================================
# PARALLEL PROCESSING CONFIGURATION
# ==============================================================================

# Configure future for parallel processing
# Adjust workers based on your system capabilities
plan("multisession", workers = 18)

# Set memory limit for large datasets
# 360 GB allocated - adjust based on available RAM
options(future.globals.maxSize = 360 * 1024^3)

message("Parallel processing configured: 18 workers, 360 GB memory limit")

# ==============================================================================
# COMMAND-LINE ARGUMENT SPECIFICATION
# ==============================================================================

option_list <- list(
  
  # ============================================================================
  # REQUIRED ARGUMENTS
  # ============================================================================
  
  make_option(
    c("--in_dir"),
    type = "character",
    default = NULL,
    help = "Input Cell Ranger directory (filtered_feature_bc_matrix/) or HDF5 file (.h5). REQUIRED.",
    metavar = "PATH"
  ),
  
  # ============================================================================
  # OUTPUT CONFIGURATION
  # ============================================================================
  
  make_option(
    c("--out_dir"),
    type = "character",
    default = "./output",
    help = "Output directory for all results [default: %default]",
    metavar = "PATH"
  ),
  
  make_option(
    c("--plot_format"),
    type = "character",
    default = "pdf",
    help = "Plot output format: 'pdf', 'png', or 'both' [default: %default]",
    metavar = "FORMAT"
  ),
  
  make_option(
    c("--width"),
    type = "integer",
    default = 14,
    help = "Base plot width in inches [default: %default]",
    metavar = "NUMBER"
  ),
  
  make_option(
    c("--height"),
    type = "integer",
    default = 12,
    help = "Base plot height in inches [default: %default]",
    metavar = "NUMBER"
  ),
  
  make_option(
    c("--use_ggrepel"),
    type = "logical",
    default = TRUE,
    help = "Use ggrepel for non-overlapping labels [default: %default]",
    metavar = "TRUE/FALSE"
  ),
  
  make_option(
    c("--zip_outputs"),
    type = "logical",
    default = TRUE,
    help = "Create ZIP bundle of all results [default: %default]",
    metavar = "TRUE/FALSE"
  ),
  
  # ============================================================================
  # INPUT DATA CONFIGURATION
  # ============================================================================
  
  make_option(
    c("--markers"),
    type = "character",
    default = "",
    help = "Path to custom markers file (CSV or RDS). If not provided, uses built-in mouse markers.",
    metavar = "PATH"
  ),
  
  make_option(
    c("--sample_map"),
    type = "character",
    default = "",
    help = "Path to sample mapping CSV with columns: suffix, sample_name, group",
    metavar = "PATH"
  ),
  
  # ============================================================================
  # QUALITY CONTROL PARAMETERS
  # ============================================================================
  
  make_option(
    c("--min_cells"),
    type = "integer",
    default = 3,
    help = "Minimum number of cells expressing a gene [default: %default]",
    metavar = "NUMBER"
  ),
  
  make_option(
    c("--min_features"),
    type = "integer",
    default = 200,
    help = "Minimum features (genes) per cell [default: %default]",
    metavar = "NUMBER"
  ),
  
  make_option(
    c("--max_features"),
    type = "integer",
    default = 6000,
    help = "Maximum features per cell (doublet filter) [default: %default]",
    metavar = "NUMBER"
  ),
  
  make_option(
    c("--max_mito_pct"),
    type = "double",
    default = 15,
    help = "Maximum mitochondrial percentage [default: %default]",
    metavar = "NUMBER"
  ),
  
  # ============================================================================
  # ANALYSIS PARAMETERS
  # ============================================================================
  
  make_option(
    c("--n_pcs"),
    type = "integer",
    default = 30,
    help = "Number of principal components to compute and use [default: %default]",
    metavar = "NUMBER"
  ),
  
  make_option(
    c("--resolution"),
    type = "double",
    default = 0.5,
    help = "Clustering resolution (higher = more clusters) [default: %default]",
    metavar = "NUMBER"
  ),
  
  # ============================================================================
  # DIFFERENTIAL EXPRESSION
  # ============================================================================
  
  make_option(
    c("--conditions"),
    type = "character",
    default = "",
    help = "Comma-separated DE comparisons (e.g., 'GroupA_vs_GroupB,GroupC_vs_GroupD')",
    metavar = "STRING"
  ),
  
  make_option(
    c("--no_compare"),
    type = "logical",
    default = FALSE,
    help = "Skip differential expression analysis [default: %default]",
    metavar = "TRUE/FALSE"
  ),
  
  # ============================================================================
  # ENRICHMENT ANALYSIS
  # ============================================================================
  
  make_option(
    c("--use_offline_gmt"),
    type = "logical",
    default = FALSE,
    help = "Use local GMT files for pathway enrichment [default: %default]",
    metavar = "TRUE/FALSE"
  ),
  
  make_option(
    c("--gmt_dir"),
    type = "character",
    default = "",
    help = "Directory containing GMT files (required if --use_offline_gmt is TRUE)",
    metavar = "PATH"
  ),
  
  make_option(
    c("--gmt_globs"),
    type = "character",
    default = "",
    help = "Comma-separated GMT filenames to use (e.g., 'msigdb.v2024.1.Mm.symbols.gmt')",
    metavar = "STRING"
  ),
  
  # ============================================================================
  # OPTIONAL ANALYSES
  # ============================================================================
  
  make_option(
    c("--run_trajectory"),
    type = "logical",
    default = TRUE,
    help = "Run Slingshot trajectory analysis [default: %default]",
    metavar = "TRUE/FALSE"
  )
)

# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

# Create parser
parser <- OptionParser(
  option_list = option_list,
  description = "\nscRNAbassoon: Optimized Single-Cell RNA-seq Analysis Pipeline\n",
  epilogue = paste0(
    "\nExample usage:\n",
    "  # Basic analysis with default parameters\n",
    "  Rscript run_bassoon.R --in_dir /data/cellranger/filtered_feature_bc_matrix/\n\n",
    "  # Full analysis with custom parameters\n",
    "  Rscript run_bassoon.R \\\n",
    "    --in_dir /data/cellranger/aggregated/ \\\n",
    "    --out_dir /results/experiment1/ \\\n",
    "    --sample_map /metadata/samples.csv \\\n",
    "    --conditions 'Treatment_vs_Control,KO_vs_WT' \\\n",
    "    --n_pcs 30 \\\n",
    "    --resolution 0.8 \\\n",
    "    --max_mito_pct 10 \\\n",
    "    --plot_format both \\\n",
    "    --use_offline_gmt TRUE \\\n",
    "    --gmt_dir /reference/msigdb/mouse/ \\\n",
    "    --gmt_globs 'm2.cp.wikipathways.v2024.1.Mm.symbols.gmt'\n\n",
    "For more information, visit: https://github.com/yourusername/scRNAbassoon\n"
  )
)

# Parse arguments
opt <- parse_args(parser)

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

# Check for required argument
if (is.null(opt$in_dir)) {
  cat("\nError: --in_dir is required.\n\n")
  print_help(parser)
  quit(status = 1)
}

# Validate input directory/file exists
if (!file.exists(opt$in_dir) && !dir.exists(opt$in_dir)) {
  cat(sprintf("\nError: Input path does not exist: %s\n\n", opt$in_dir))
  quit(status = 1)
}

# Validate GMT configuration
if (opt$use_offline_gmt) {
  if (opt$gmt_dir == "") {
    cat("\nError: --gmt_dir must be specified when --use_offline_gmt is TRUE\n\n")
    quit(status = 1)
  }
  if (!dir.exists(opt$gmt_dir)) {
    cat(sprintf("\nError: GMT directory does not exist: %s\n\n", opt$gmt_dir))
    quit(status = 1)
  }
  if (opt$gmt_globs == "") {
    cat("\nWarning: --gmt_globs is empty. No GMT enrichment will be performed.\n")
  }
}

# Validate plot format
if (!opt$plot_format %in% c("pdf", "png", "both")) {
  cat(sprintf("\nError: Invalid plot format '%s'. Must be 'pdf', 'png', or 'both'\n\n", opt$plot_format))
  quit(status = 1)
}

# Validate numeric parameters
if (opt$min_features >= opt$max_features) {
  cat("\nError: --min_features must be less than --max_features\n\n")
  quit(status = 1)
}

if (opt$max_mito_pct < 0 || opt$max_mito_pct > 100) {
  cat("\nError: --max_mito_pct must be between 0 and 100\n\n")
  quit(status = 1)
}

if (opt$resolution <= 0) {
  cat("\nError: --resolution must be greater than 0\n\n")
  quit(status = 1)
}

# ==============================================================================
# DISPLAY CONFIGURATION
# ==============================================================================

cat("\n")
cat("================================================================================\n")
cat("                        scRNAbassoon Pipeline v1.0.0                            \n")
cat("================================================================================\n")
cat("\n")
cat("Configuration:\n")
cat(sprintf("  Input:              %s\n", opt$in_dir))
cat(sprintf("  Output:             %s\n", opt$out_dir))
cat(sprintf("  Sample map:         %s\n", ifelse(opt$sample_map == "", "None (auto-detect)", opt$sample_map)))
cat(sprintf("  Markers:            %s\n", ifelse(opt$markers == "", "Built-in mouse markers", opt$markers)))
cat("\n")
cat("Quality Control:\n")
cat(sprintf("  Min cells/gene:     %d\n", opt$min_cells))
cat(sprintf("  Min features/cell:  %d\n", opt$min_features))
cat(sprintf("  Max features/cell:  %d\n", opt$max_features))
cat(sprintf("  Max mito %%:         %.1f\n", opt$max_mito_pct))
cat("\n")
cat("Analysis Parameters:\n")
cat(sprintf("  PCs:                %d\n", opt$n_pcs))
cat(sprintf("  Resolution:         %.2f\n", opt$resolution))
cat(sprintf("  DE comparisons:     %s\n", ifelse(opt$conditions == "", "None", opt$conditions)))
cat(sprintf("  Trajectory:         %s\n", opt$run_trajectory))
cat("\n")
cat("Output Options:\n")
cat(sprintf("  Plot format:        %s\n", opt$plot_format))
cat(sprintf("  Plot size:          %d x %d inches\n", opt$width, opt$height))
cat(sprintf("  ZIP bundle:         %s\n", opt$zip_outputs))
cat("\n")
if (opt$use_offline_gmt) {
  cat("Enrichment:\n")
  cat(sprintf("  GMT directory:      %s\n", opt$gmt_dir))
  cat(sprintf("  GMT files:          %s\n", opt$gmt_globs))
  cat("\n")
}
cat("================================================================================\n")
cat("\n")

# ==============================================================================
# PIPELINE EXECUTION
# ==============================================================================

# Wrap in tryCatch for error handling
result <- tryCatch(
  {
    run_bassoon_pipeline(
      in_dir = opt$in_dir,
      out_dir = opt$out_dir,
      markers_file = opt$markers,
      sample_map_file = opt$sample_map,
      conditions = opt$conditions,
      no_compare = opt$no_compare,
      n_pcs = opt$n_pcs,
      resolution = opt$resolution,
      min_cells = opt$min_cells,
      min_features = opt$min_features,
      max_features = opt$max_features,
      max_mito_pct = opt$max_mito_pct,
      plot_format = opt$plot_format,
      width = opt$width,
      height = opt$height,
      use_ggrepel = opt$use_ggrepel,
      zip_outputs = opt$zip_outputs,
      run_trajectory = opt$run_trajectory,
      use_offline_gmt = opt$use_offline_gmt,
      gmt_dir = opt$gmt_dir,
      gmt_globs = opt$gmt_globs
    )
    
    # Success
    cat("\n")
    cat("================================================================================\n")
    cat("                           PIPELINE COMPLETED SUCCESSFULLY                      \n")
    cat("================================================================================\n")
    cat("\n")
    cat(sprintf("Results saved to: %s\n", opt$out_dir))
    cat("\n")
    
    return(0)
  },
  error = function(e) {
    cat("\n")
    cat("================================================================================\n")
    cat("                              PIPELINE FAILED                                   \n")
    cat("================================================================================\n")
    cat("\n")
    cat("Error message:\n")
    cat(sprintf("  %s\n", e$message))
    cat("\n")
    cat("Traceback:\n")
    print(sys.calls())
    cat("\n")
    cat("Please check your input parameters and data, then try again.\n")
    cat("For help: Rscript run_bassoon.R --help\n")
    cat("\n")
    
    return(1)
  }
)

# Exit with appropriate status
quit(status = result)
