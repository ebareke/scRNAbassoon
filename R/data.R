#' Default Mouse Markers for Annotation
#' 
#' A curated list of gene markers used for automated cell type annotation
#' in mouse embryonic single-cell RNA-seq data. Covers developmental stages
#' from E7.5 to E9.5 and includes markers for various tissue types including
#' neural crest, mesoderm, endoderm, and ectoderm derivatives.
#' 
#' @format A named list where each element contains a character vector of gene symbols:
#' \describe{
#'   \item{E7.5-Gut}{Early gut markers: Apela, Krt8}
#'   \item{E7.5-Caudal_lateral_epiblast}{Posterior epiblast markers}
#'   \item{E7.5-Rostral_neuroectoderm}{Anterior neural markers}
#'   \item{E8.25-Neural_crest}{Neural crest markers: Sox10, Dlx2, Foxd3, Sox9}
#'   \item{E8.25-Forebrain_midbrain}{Forebrain and midbrain markers}
#'   \item{E9.5-Skeletal_muscle_progenitors}{Muscle progenitor markers}
#'   \item{E9.5-Hepatocytes}{Liver cell markers}
#'   \item{...}{And many more developmental cell types}
#' }
#' 
#' @details
#' This marker set is particularly useful for:
#' \itemize{
#'   \item Mouse embryonic development studies (E7.5-E9.5)
#'   \item Gastrulation and neurulation analysis
#'   \item Organogenesis studies
#'   \item Neural crest and mesoderm differentiation
#' }
#' 
#' For other species or developmental stages, users should provide custom
#' marker lists using the \code{markers_file} parameter in 
#' \code{\link{run_bassoon_pipeline}}.
#' 
#' @source Curated from published single-cell atlases of mouse embryonic development
#' 
#' @examples
#' # View all available cell types
#' names(default_mouse_markers)
#' 
#' # Check markers for a specific cell type
#' default_mouse_markers$`E8.25-Neural_crest`
#' 
#' # Use in annotation
#' \dontrun{
#' seu <- annotate_by_markers(seu, default_mouse_markers)
#' }
#' 
#' @export
default_mouse_markers <- list(
  "E7.5-Gut" = c("Apela","Krt8"),
  "E7.5-Caudal_lateral_epiblast" = c("Nkx1-2","Cdx2","Gbx2"),
  "E7.5-Rostral_neuroectoderm" = c("Six3","Hesx1"),
  "E7.5-Caudal_neuroectoderm" = c("Gbx2","Hes3"),
  "E7.5-Surface_ectoderm" = c("Foxg1","Trp63","Grhl2","Grhl3"),
  "E7.5-Extraembryonic_mesoderm" = c("Bmp4","Cdx2","Hoxa10"),
  "E7.5-Paraxial_mesoderm_A" = c("Tbx1","Pax3"),
  "E7.5-Paraxial_mesoderm_B" = c("Tbx6","Dll1"),
  "E7.5-Splanchnic_mesoderm" = c("Tcf21","Isl1","Gata4"),
  "E7.75-Paraxial_mesoderm_C" = c("Cdx1","Cdx2","Hes7"),
  "E7.75-Intermediate_mesoderm" = c("Osr1","Lhx1","Pax2","Pax8"),
  "E7.75-First_heart_field" = c("Tbx5","Hcn4","Gata4"),
  "E7.75-Allantois" = c("Tbx4","Hoxa11"),
  "E7.75-Endothelium" = c("Kdr","Pecam1","Cdh5"),
  "E7.75-Primitive_erythroid_cells" = c("Hba-a1","Hbb-y","Hba-x","Hbb-bh1"),
  "E7.75-Somatic_mesoderm" = c("Lhx1","Prrx1","Lix1","Msx1"),
  "E8-Spinal_cord" = c("Foxb1","Pax6","Crabp2"),
  "E8.25-Neuromesodermal_progenitors" = c("Cdx4","Epha5","Hes3"),
  "E8.25-Neural_crest" = c("Sox10","Dlx2","Foxd3","Sox9"),
  "E8.25-Forebrain_midbrain" = c("Pax2","Sox2","Igfbp2","Otx2","Pcsk1n"),
  "E8.25-Hindbrain" = c("Pax2","Sox2","Crabp1","Fst","Hoxa2"),
  "E8.25-Amniochorionic_mesoderm_A" = c("Hlx","Postn"),
  "E8.25-Amniochorionic_mesoderm_B" = c("Bmp2","Postn"),
  "E8.25-Second_heart_field" = c("Isl1","Tbx1"),
  "E8.5_ab-Pre-epidermal_keratinocytes" = c("Tfap2b","Trp63","Egfr"),
  "E8.5_ab-Placodal_area" = c("Six1","Eya1"),
  "E8.5_ab-Fusing_epithelium" = c("Itgb1","Rhoa"),
  "E8.5_ab-Anterior_floor_plate" = c("Foxa2","Shh","Ntn1","Bmp7"),
  "E8.5_ab-Posterior_floor_plate" = c("Foxa2","Shh","Ntn1"),
  "E9.5-Chondrocyte_and_osteoblast_progenitors" = c("Pax1","Pax9"),
  "E9.5-Osteoblast_progenitors_A" = c("Runx2","Twist2","Prrx1","Pax3"),
  "E9.5-Osteoblast_progenitors_B" = c("Runx2","Twist2","Prrx1","Meis2"),
  "E9.5-Skeletal_muscle_progenitors" = c("Fap","Sim1"),
  "E9.5-Limb_mesenchyme_progenitors" = c("Cpa2","Msx1","Fgf10","Wnt5a","Lmx1b"),
  "E9.5-Branchial_arch_epithelium" = c("Has2","Pax1","Pitx2"),
  "E9.5-Olfactory_epithelium" = c("Pax3","Six3"),
  "E9.5-Otic_epithelium" = c("Pax2","Lmx1a"),
  "E9.5-Renal_epithelium" = c("Sim1","Pax2","Mecom"),
  "E9.5-Gut_and_lung_epithelium" = c("Hoxa10","Gata4","Prox1","Foxp2"),
  "E9.5-Hepatocytes" = c("A1cf","Afp","Alb","Apoa2","Afp29","Pik3c2g","Hoga1"),
  "E9.5-Mesencephalon_MHB" = c("Pax5","En1","Fgf8","Pax3","Pax7","Dmbx1"),
  "E9.5-Neuron_progenitor_cells" = c("Mybl1","Prmt8","Gadd45g","Cdkn1c","Btg2"),
  "E9.5-Di_telencephalon" = c("Emx2","Pax6","Sox5","Sox6","Wnt8b"),
  "E9.5-Spinal_cord_ventral" = c("Foxb1","Pax6","Crabp2"),
  "E9.5-Spinal_cord_dorsal" = c("Pth2r","Fabp7","Pax3","Fzd10","Hes5"),
  "E9.5-Neural_crest_PNS_neurons" = c("Ppp1r1c","Syt13","Shox2","Ptprr","Pcbp3"),
  "E9.5-Neural_crest_PNS_glia" = c("Sox10","Dlx2","Foxd3","Pax3"),
  "E9.5-Retinal_primordium" = c("Pax2","Nlgn1","Six3"),
  "E9.5-Roof_plate" = c("Lmx1a","Msx1"),
  "E9.5-Apical_ectodermal_ridge" = c("Fgf8","Msx2","Rspo2"),
  "E9.5-Motor_neurons" = c("Mnx1","Isl2","Lhx3","Lhx4"),
  "GFP" = c("GFP"),
  "RFP" = c("RFP")
)
