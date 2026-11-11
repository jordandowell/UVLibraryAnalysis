## ========================================================================== ##
## ORGANIZE RAW-DAD OUTPUTS INTO AN INTERPRETATION-ORIENTED DIRECTORY TREE
## ========================================================================== ##
## This script runs only after the alignment, annotation, final-figure, and
## validation stages have succeeded. It changes file locations, never file
## contents. Exact basenames remain stable so VALIDATE_OUTPUTS.R can locate
## them recursively and older analysis code can still identify each product.
## ========================================================================== ##

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) {
  stop("Usage: Rscript 09_organize_raw_dad_outputs.R <policy-output-root>")
}
OUTPUT_ROOT <- arguments[1]
if (!dir.exists(OUTPUT_ROOT)) stop("Output root not found: ", OUTPUT_ROOT)

alignment_root <- file.path(OUTPUT_ROOT, "01_Alignment_QC")
annotation_root <- file.path(OUTPUT_ROOT, "02_Library_Annotation")
final_root <- file.path(OUTPUT_ROOT, "03_Final_Results")
if (!all(dir.exists(c(alignment_root, annotation_root, final_root)))) {
  stop("The raw-DAD output root is missing one or more numbered stages.")
}

# Stage each replacement beside its destination, verify the copy, and report
# the exact locked file if Windows/Excel/a PDF viewer/OneDrive prevents the old
# target from being removed. A partially completed organizer is safe to rerun:
# successfully moved files stay organized and remaining flat files are retried.
replace_organized_file <- function(source, target) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  staged <- paste0(target, ".organizing-", Sys.getpid())
  if (file.exists(staged)) unlink(staged, force = TRUE)
  copied <- file.copy(source, staged, overwrite = TRUE, copy.mode = TRUE)
  if (!copied || !file.exists(staged) ||
      unname(tools::md5sum(source)) != unname(tools::md5sum(staged))) {
    if (file.exists(staged)) unlink(staged, force = TRUE)
    stop("Could not stage output for organization: ", source)
  }

  if (file.exists(target)) {
    for (attempt in seq_len(4L)) {
      unlink(target, force = TRUE)
      if (!file.exists(target)) break
      Sys.sleep(0.25)
    }
  }
  if (file.exists(target)) {
    unlink(staged, force = TRUE)
    stop(
      "Cannot replace organized output: ", target,
      "\nThe file is probably open in Excel, a PDF/RStudio viewer, or locked ",
      "by cloud sync. Close it (or pause syncing), then rerun only:\n  ",
      "Rscript 09_organize_raw_dad_outputs.R ", OUTPUT_ROOT
    )
  }

  installed <- file.rename(staged, target)
  if (!installed) {
    installed <- file.copy(staged, target, overwrite = FALSE, copy.mode = TRUE)
    if (installed) unlink(staged, force = TRUE)
  }
  if (!installed || !file.exists(target) ||
      unname(tools::md5sum(source)) != unname(tools::md5sum(target))) {
    if (file.exists(staged)) unlink(staged, force = TRUE)
    stop("Could not install organized output: ", target)
  }
  unlink(source, force = TRUE)
  if (file.exists(source)) {
    stop(
      "Organized copy was created, but the flat source is locked: ", source,
      "\nClose the file and rerun the organizer."
    )
  }
  invisible(TRUE)
}

move_flat_files <- function(stage_root, classify_file) {
  flat_files <- list.files(stage_root, full.names = TRUE, recursive = FALSE)
  flat_files <- flat_files[file.info(flat_files)$isdir %in% FALSE]
  # Folder guides belong permanently at the stage root. Excluding them makes
  # this organizer idempotent: rerunning it cannot move yesterday's guide into
  # a model-QC subfolder before writing today's replacement.
  flat_files <- flat_files[basename(flat_files) != "00_FOLDER_GUIDE.csv"]
  if (!length(flat_files)) return(invisible(character()))
  destinations <- vapply(basename(flat_files), classify_file, character(1))
  for (destination in unique(destinations)) {
    dir.create(file.path(stage_root, destination), recursive = TRUE,
               showWarnings = FALSE)
  }
  for (index in seq_along(flat_files)) {
    target <- file.path(stage_root, destinations[index], basename(flat_files[index]))
    replace_organized_file(flat_files[index], target)
  }
  invisible(destinations)
}

alignment_category <- function(filename) {
  if (grepl("policy_manifest|alignment_parameters|input_audit|input_exclusions|wavelength_grid", filename)) {
    return("00_Run_Manifest")
  }
  if (grepl("peak_picking|peak_detection_summary|peak_detections", filename)) {
    return("01_Peak_Picking")
  }
  if (grepl("retention_time|rt_landmark", filename)) return("02_RT_Correction")
  if (grepl("rare_fragment|split_feature|neighboring_feature|competition_family|ambiguous_unforced|provisional_rare", filename)) {
    return("04_Reconciliation_and_Review")
  }
  if (grepl(
    "sample_by_aligned_feature_(integrated|assignment)|sample_by_feature_with_summed",
    filename
  )) {
    return("05_Quantitative_Matrices")
  }
  if (grepl("representative_spectra|aligned_reference_spectra", filename)) {
    return("06_Representative_Spectra")
  }
  if (grepl("^genus_", filename)) return("07_Annotated_Overviews")
  if (grepl("raw_dad_alignment_objects[.]rds$", filename)) return("99_Reproducibility")
  if (grepl("aligned_feature|aligned_peak_map|feature_prevalence", filename)) {
    return("03_Feature_Alignment")
  }
  "99_Reproducibility"
}

annotation_category <- function(filename) {
  if (grepl("dictionary", filename)) return("00_Dictionaries")
  if (grepl("analysis_parameters|library_input_validation|wavelength_grid_alignment", filename)) {
    return("00_Run_Manifest")
  }
  if (grepl("comprehensive_results|compound_first_interpretation|ensemble_.*probabilities|peak_annotation_metadata|parallel_label_interpretations|original_data_with_peak_annotations|practical_.*decision_summary|hybrid_label_type_summary", filename)) {
    return("01_Primary_Annotations")
  }
  if (grepl("pseudochromatogram|signal_composition_mirror", filename)) {
    return("06_Annotation_Overviews")
  }
  if (grepl("ambiguous|coelution|mixture|unsupported_family", filename)) {
    return("05_Ambiguity_and_Coelution")
  }
  if (grepl("consensus_cluster|consensus_matrix|consensus_dendrogram|consensus_coclustering|representative_spectral_overlays", filename)) {
    return("04_Unsupervised_Clustering")
  }
  if (grepl("bootstrap|insertion|ablation|wavelength_omission|wavelength_block|seed_sensitivity", filename)) {
    return("03_Robustness_and_Sensitivity")
  }
  if (grepl("^library_|class_distance|class_silhouette|spectral_quality|distance_method|analogue_set", filename)) {
    return("02_Library_Model_QC")
  }
  "02_Library_Model_QC"
}

final_category <- function(filename) {
  if (grepl("READ_FIRST", filename)) return("00_READ_FIRST")
  if (grepl("upper_bound", filename)) return("03_Sensitivity_and_Review")
  if (grepl("_strict[.]csv$", filename)) return("01_Primary_Tables")
  if (grepl("matrix_strict[.]pdf$|genus_.*pseudochromatogram[.]pdf$", filename)) {
    return("02_Primary_Figures")
  }
  "03_Sensitivity_and_Review"
}

move_flat_files(alignment_root, alignment_category)
move_flat_files(annotation_root, annotation_category)
move_flat_files(final_root, final_category)

write_stage_guide <- function(stage_root, rows) {
  target <- file.path(stage_root, "00_FOLDER_GUIDE.csv")
  staged <- tempfile("folder-guide-", tmpdir = stage_root, fileext = ".csv")
  write.csv(rows, staged, row.names = FALSE)
  replace_organized_file(staged, target)
}

write_stage_guide(alignment_root, data.frame(
  folder = c(
    "00_Run_Manifest", "01_Peak_Picking", "02_RT_Correction",
    "03_Feature_Alignment", "04_Reconciliation_and_Review",
    "05_Quantitative_Matrices", "06_Representative_Spectra",
    "07_Annotated_Overviews", "99_Reproducibility"
  ),
  use = c(
    "Confirm inputs, policy, parameters, and wavelength handling.",
    "Audit candidate detection and sample-level peak-picking performance.",
    "Inspect RT landmarks and the magnitude of sample corrections.",
    "Inspect homologous-feature identities, support, and sample maps.",
    paste(
      "Review reconciled detections, probable split-alignment pairs, and",
      "unresolved neighboring features."
    ),
    paste(
      "Access NA-free confident abundance matrices, uncertainty flags,",
      "feature-wise upper bounds, and family sensitivity tables."
    ),
    "Inspect the spectra used for library comparison, including uncertainty.",
    "Inspect annotated genus-level class, chromophore, and hybrid summaries.",
    "Retain serialized objects needed to reproduce or extend the run."
  ), stringsAsFactors = FALSE
))

write_stage_guide(annotation_root, data.frame(
  folder = c(
    "00_Run_Manifest", "00_Dictionaries", "01_Primary_Annotations",
    "02_Library_Model_QC", "03_Robustness_and_Sensitivity",
    "04_Unsupervised_Clustering", "05_Ambiguity_and_Coelution",
    "06_Annotation_Overviews"
  ),
  use = c(
    "Confirm the library fingerprint and analysis parameters.",
    "Interpret class and chromophore names and their compound membership.",
    "Start here for compound-first analogue and final label evidence.",
    "Assess whether the fixed library can support the requested distinctions.",
    "Assess stability across bootstraps, metrics, linkage, and wavelengths.",
    "Explore sample-feature structure without using class labels to fit it.",
    "Investigate uncertain analogues, tied families, and possible co-elution.",
    "View annotation-specific pseudochromatograms and their scaling audits."
  ), stringsAsFactors = FALSE
))

write_stage_guide(final_root, data.frame(
  folder = c("00_READ_FIRST", "01_Primary_Tables", "02_Primary_Figures",
             "03_Sensitivity_and_Review"),
  use = c(
    "Open READ_FIRST_output_order.csv before interpreting other files.",
    paste(
      "Primary NA-free confident matrix, matching uncertainty flags, and",
      "annotated per-feature table."
    ),
    "Primary sample-by-peak matrix and genus chemical overview figures.",
    paste(
      "Feature-wise upper bounds, complete matrices, and excluded/review",
      "features for sensitivity analysis."
    )
  ), stringsAsFactors = FALSE
))

start_here <- data.frame(
  interpretation_order = 1:7,
  location = c(
    "03_Final_Results/00_READ_FIRST/READ_FIRST_output_order.csv",
    "03_Final_Results/01_Primary_Tables",
    "03_Final_Results/02_Primary_Figures",
    "03_Final_Results/03_Sensitivity_and_Review",
    "01_Alignment_QC/00_FOLDER_GUIDE.csv",
    "02_Library_Annotation/00_FOLDER_GUIDE.csv",
    "00_Library_QC"
  ),
  purpose = c(
    "Read the file-by-file interpretation order.",
    paste(
      "Use the strict confident-abundance table for primary analyses and its",
      "matching flag table to locate unresolved candidate signal."
    ),
    "Inspect detection patterns and class/chromophore overviews.",
    "Determine whether conclusions change when review features are retained.",
    "Trace peak detection, RT correction, alignment, and quantitative evidence.",
    "Trace library analogue, uncertainty, and unsupervised spectral evidence.",
    "Audit raw-standard extraction, compound admission, fallback, and validation."
  ), stringsAsFactors = FALSE
)
start_here_target <- file.path(OUTPUT_ROOT, "00_START_HERE.csv")
start_here_staged <- tempfile(
  "start-here-", tmpdir = OUTPUT_ROOT, fileext = ".csv"
)
write.csv(start_here, start_here_staged, row.names = FALSE)
replace_organized_file(start_here_staged, start_here_target)

message("Organized raw-DAD outputs: ", OUTPUT_ROOT)
