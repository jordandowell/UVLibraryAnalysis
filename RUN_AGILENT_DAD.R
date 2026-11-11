## ========================================================================== ##
## SINGLE ENTRY POINT FOR RAW AGILENT HPLC-DAD DATA
## ========================================================================== ##
## Usage (the alignment policy is deliberately required)
##   Rscript RUN_AGILENT_DAD.R Rhododendron \
##     --alignment-policy=retained_uniqueness
##   Rscript RUN_AGILENT_DAD.R Rhododendron \
##     --alignment-policy=expanded_rescue
##
## Required input
##   data/<Genus>_Agilent/<sample_name>.D/DAD1.uv
##
## Optional maintainer flag
##   --rebuild-library  Re-extract the active library from raw standard `.D`
##                      files before the sample analysis. Ordinary users do
##                      not need the private/raw standard acquisitions.
##
## This route does not read <Genus>_MaxAbsorbances.csv. It:
##   1. detects peaks independently in every raw DAD cube;
##   2. performs presence/absence-aware RT-spectral alignment;
##   3. exports an integrated-DAD-area sample-by-feature matrix;
##   4. annotates robust feature spectra with the fixed library; and
##   5. draws genus-level prevalence, class, chromophore, and hybrid
##      pseudochromatograms using the aligned quantitative results.
##
## Alignment review figures
##   `sample_by_aligned_peak_map.pdf` shows assigned peaks directly: the same
##   aligned feature forms one vertical consensus-RT column across samples.
##   The paginated identity atlas uses a stable feature color, while the older
##   detection matrix uses color for QC/review status rather than identity.
##
## Retention time is used only by steps 1-3. The library model in step 4 uses
## wavelength shape alone.
## ========================================================================== ##

runner <- file.path(R.home("bin"), "Rscript")
arguments <- commandArgs(trailingOnly = TRUE)
rebuild_library_arguments <- tolower(arguments) == "--rebuild-library"
if (sum(rebuild_library_arguments) > 1L) {
  stop("Supply --rebuild-library at most once.")
}
rebuild_library <- any(rebuild_library_arguments)
arguments <- arguments[!rebuild_library_arguments]

# The repository ships a validated, compact active spectral library so ordinary
# users do not need hundreds of megabytes of raw authenticated-standard files.
# Maintainers can explicitly rebuild it after adding/replacing standards. Raw
# spectra remain auxiliary QC evidence and are never appended as duplicate
# training rows.
if (rebuild_library) {
  message("Rebuilding the active library from raw standard acquisitions.")
  raw_library_status <- system2(runner, "00_build_raw_dad_library.R")
  if (!identical(raw_library_status, 0L)) {
    stop("Raw authenticated-standard library build failed.")
  }
} else {
  message(
    "Using the packaged validated spectral library. Supply --rebuild-library ",
    "only when raw authenticated-standard acquisitions are available."
  )
}

source("00_library_inputs.R")
library_validation <- load_library_inputs()
message(
  "Validated library: ",
  length(unique(library_validation$spectra$Compound)), " compounds; ",
  length(unique(library_validation$class_dictionary$class_name)), " classes; ",
  length(unique(
    library_validation$chromophore_mapping$fine_chromophore_code
  )), " fine chromophore groups."
)

# Alignment policy is a scientific choice, not a hidden tuning parameter. The
# runner therefore refuses to proceed until the user names one of the two
# presets. Policy-specific output folders prevent one choice from overwriting
# the other or silently feeding the wrong representative spectra downstream.
policy_arguments <- grep(
  "^--alignment-policy=", arguments, value = TRUE, ignore.case = TRUE
)
if (length(policy_arguments) != 1L) {
  stop(
    "Choose exactly one raw-alignment policy:\n",
    "  --alignment-policy=retained_uniqueness\n",
    "  --alignment-policy=expanded_rescue\n",
    "Example: Rscript RUN_AGILENT_DAD.R Rhododendron ",
    "--alignment-policy=retained_uniqueness"
  )
}
alignment_policy <- tolower(sub(
  "^--alignment-policy=", "", policy_arguments, ignore.case = TRUE
))
valid_alignment_policies <- c("retained_uniqueness", "expanded_rescue")
if (!alignment_policy %in% valid_alignment_policies) {
  stop(
    "Unknown --alignment-policy value: ", alignment_policy, ". Use ",
    "retained_uniqueness or expanded_rescue."
  )
}
policy_folder_label <- c(
  retained_uniqueness = "RetainedUniqueness",
  expanded_rescue = "ExpandedRescue"
)[alignment_policy]
policy_display_name <- c(
  retained_uniqueness = "Retained uniqueness",
  expanded_rescue = "Expanded rescue"
)[alignment_policy]
datasets <- arguments[!grepl(
  "^--alignment-policy=", arguments, ignore.case = TRUE
)]
if (!length(datasets)) datasets <- "Rhododendron"
run_stage <- function(script, arguments = character(),
                      environment = character(), label) {
  if (length(environment)) {
    previous <- Sys.getenv(names(environment), unset = NA_character_)
    on.exit({
      missing_before <- is.na(previous)
      if (any(missing_before)) Sys.unsetenv(names(previous)[missing_before])
      if (any(!missing_before)) {
        do.call(Sys.setenv, as.list(previous[!missing_before]))
      }
    }, add = TRUE)
    do.call(Sys.setenv, as.list(environment))
  }
  message(label)
  status <- system2(runner, c(script, arguments))
  if (!identical(status, 0L)) stop(label, " failed.")
}

# R and Windows may refuse to replace a file that is open in Excel, a PDF
# viewer, the RStudio data viewer, or a cloud-sync preview. Detect that before
# the expensive raw-data analysis rather than failing at the final organizer.
remove_existing_output <- function(paths, description) {
  paths <- paths[file.exists(paths) | dir.exists(paths)]
  if (!length(paths)) return(invisible(TRUE))
  for (path in paths) {
    unlink(path, recursive = dir.exists(path), force = TRUE)
  }
  remaining <- paths[file.exists(paths) | dir.exists(paths)]
  if (length(remaining)) {
    stop(
      "Cannot replace existing ", description, ":\n  ",
      paste(remaining, collapse = "\n  "),
      "\nClose open CSV/PDF files and RStudio viewers, pause cloud syncing if ",
      "needed, and rerun. No sample processing has started."
    )
  }
  invisible(TRUE)
}

for (source_dataset in datasets) {
  raw_dir <- file.path("data", paste0(source_dataset, "_Agilent"))
  if (!dir.exists(raw_dir)) stop("Raw Agilent folder not found: ", raw_dir)

  legacy_folders <- c(
    paste0(source_dataset, "_RawDAD_Alignment"),
    paste0(source_dataset, "_RawDAD")
  )
  if (any(dir.exists(legacy_folders))) {
    message(
      "Legacy unsuffixed raw-DAD outputs were found and will be ignored: ",
      paste(legacy_folders[dir.exists(legacy_folders)], collapse = ", "),
      ". Only the policy-qualified folders created by this run are current."
    )
  }

  sample_directories <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
  sample_directories <- sample_directories[grepl(
    "[.]D$", sample_directories, ignore.case = TRUE
  )]
  if (!length(sample_directories)) {
    stop("No .D sample folders were found in: ", raw_dir)
  }

  # One policy-qualified root keeps the collaborator handoff coherent while
  # preserving the scientific separation between peak alignment, library-only
  # annotation, and final analysis-ready deliverables.
  output_root <- paste0(
    source_dataset, "_RawDAD_", policy_folder_label
  )
  alignment_dir <- file.path(output_root, "01_Alignment_QC")
  annotation_dataset <- paste0(
    source_dataset, "_RawDAD_", policy_folder_label
  )
  annotation_slug <- tolower(gsub("[^A-Za-z0-9]+", "_", annotation_dataset))
  annotation_dir <- file.path(output_root, "02_Library_Annotation")
  final_dir <- file.path(output_root, "03_Final_Results")
  library_qc_dir <- file.path(output_root, "00_Library_QC")

  # Snapshot the exact active library evidence used by this genus run. The
  # central library folder remains authoritative, while this policy-specific
  # copy makes collaborator handoffs self-contained and fingerprint-auditable.
  dir.create(library_qc_dir, recursive = TRUE, showWarnings = FALSE)
  library_qc_files <- c(
    "library/Absorbance_StandardLibrary_RawDAD.csv",
    "library/raw_dad_library_manifest.csv",
    "library/raw_dad_identity_corrections.csv",
    "library/standard_curve_manifest.csv",
    list.files("library/raw_dad_derived", full.names = TRUE,
               recursive = FALSE)
  )
  library_qc_files <- library_qc_files[file.exists(library_qc_files)]
  copied_library_qc <- file.copy(
    library_qc_files,
    file.path(library_qc_dir, basename(library_qc_files)),
    overwrite = TRUE
  )
  if (!all(copied_library_qc)) {
    stop("Could not copy active raw-standard library QC into: ", library_qc_dir)
  }

  # A new production run starts from an empty purpose-based tree. This avoids
  # validating stale conditional figures from an earlier run while the three
  # analysis stages write fresh files at the stage roots. The organizer
  # recreates the hierarchy only after the flat outputs pass validation.
  organized_subfolders <- c(
    file.path(alignment_dir, c(
      "00_Run_Manifest", "01_Peak_Picking", "02_RT_Correction",
      "03_Feature_Alignment", "04_Reconciliation_and_Review",
      "05_Quantitative_Matrices", "06_Representative_Spectra",
      "07_Annotated_Overviews", "99_Reproducibility"
    )),
    file.path(annotation_dir, c(
      "00_Run_Manifest", "00_Dictionaries", "01_Primary_Annotations",
      "02_Library_Model_QC", "03_Robustness_and_Sensitivity",
      "04_Unsupervised_Clustering", "05_Ambiguity_and_Coelution",
      "06_Annotation_Overviews"
    )),
    file.path(final_dir, c(
      "00_READ_FIRST", "01_Primary_Tables", "02_Primary_Figures",
      "03_Sensitivity_and_Review"
    ))
  )
  remove_existing_output(
    organized_subfolders,
    "organized output folders (a file may be open or cloud-synced)"
  )
  remove_existing_output(c(
    file.path(output_root, "00_START_HERE.csv"),
    file.path(alignment_dir, "00_FOLDER_GUIDE.csv"),
    file.path(annotation_dir, "00_FOLDER_GUIDE.csv"),
    file.path(final_dir, "00_FOLDER_GUIDE.csv")
  ), "run guides")

  # Earlier development versions wrote three unannotated pseudochromatograms
  # and a duplicate matrix alias into the alignment-QC folder. They are not
  # part of the current handoff: the annotated versions live in Final Results,
  # and the strict-primary matrix has one authoritative filename. Removing
  # only these explicitly named files prevents stale outputs from confusing a
  # collaborator when an existing policy folder is rerun.
  obsolete_alignment_outputs <- file.path(alignment_dir, c(
    "sample_by_aligned_feature_integrated_dad_area_conservative.csv",
    "genus_prevalence_pseudochromatogram.pdf",
    "genus_abundance_pseudochromatogram.pdf",
    "genus_abundance_pseudochromatogram_signal_composition_mirror.pdf"
  ))
  remove_existing_output(obsolete_alignment_outputs, "obsolete outputs")

  source_slug <- tolower(gsub("[^A-Za-z0-9]+", "_", source_dataset))
  reference_file <- file.path(
    alignment_dir,
    paste0(source_slug, "_rawdad_aligned_reference_spectra.csv")
  )

  alignment_environment <- c(
    UV_DATASET_NAME = source_dataset,
    UV_AGILENT_DIR = normalizePath(raw_dir),
    UV_RAW_ALIGNMENT_OUTPUT_DIR = alignment_dir,
    UV_RAW_ALIGNMENT_POLICY = alignment_policy,
    UV_RAW_ALIGNMENT_POLICY_DISPLAY_NAME = policy_display_name,
    UV_RAW_ANNOTATION_DATASET = annotation_dataset,
    UV_RAW_ANNOTATION_OUTPUT_DIR = annotation_dir,
    UV_RAW_FINAL_OUTPUT_DIR = final_dir
  )
  run_stage(
    "01_raw_agilent_dad_alignment.R", source_dataset,
    alignment_environment,
    paste0(
      "Aligning raw Agilent DAD peaks: ", source_dataset,
      " [", policy_display_name, "]"
    )
  )
  if (!file.exists(reference_file)) {
    stop("Raw alignment did not create the reference spectra: ", reference_file)
  }
  policy_manifest_file <- file.path(
    alignment_dir, "raw_alignment_policy_manifest.csv"
  )
  if (!file.exists(policy_manifest_file)) {
    stop("Raw alignment did not create its policy manifest: ",
         policy_manifest_file)
  }
  policy_manifest <- read.csv(
    policy_manifest_file, stringsAsFactors = FALSE, check.names = FALSE
  )
  if (nrow(policy_manifest) != 1L ||
      !identical(policy_manifest$alignment_policy, alignment_policy) ||
      !identical(
        policy_manifest$downstream_annotation_dataset, annotation_dataset
      )) {
    stop(
      "Alignment-policy manifest does not match the requested downstream ",
      "annotation route: ", policy_manifest_file
    )
  }

  annotation_environment <- c(
    alignment_environment,
    UV_DATASET_NAME = annotation_dataset,
    UV_DATASET_FILE = normalizePath(reference_file),
    UV_OUTPUT_DIR = annotation_dir,
    UV_SPECTRAL_LIBRARY_FILE = normalizePath(
      "library/Absorbance_StandardLibrary_RawDAD.csv"
    ),
    UV_ANALYSIS_INPUT_MODE = "raw_agilent_dad_joint_alignment"
  )
  run_stage(
    "03_comprehensive_unsupervised_analysis.R", character(),
    annotation_environment,
    paste0("Annotating aligned features from library standards only: ",
           source_dataset)
  )
  expected_annotation <- file.path(
    annotation_dir, paste0(annotation_slug, "_comprehensive_results.csv")
  )
  if (!file.exists(expected_annotation)) {
    stop("Library annotation did not create: ", expected_annotation)
  }
  annotation_parameters <- read.csv(
    file.path(annotation_dir, "analysis_parameters.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!"source_alignment_policy" %in% names(annotation_parameters) ||
      !identical(
        annotation_parameters$source_alignment_policy, alignment_policy
      )) {
    stop(
      "Downstream annotation did not preserve alignment policy: ",
      alignment_policy
    )
  }

  run_stage(
    "08_raw_dad_genus_pseudochromatograms.R", source_dataset,
    annotation_environment,
    paste0("Drawing annotated genus pseudochromatograms: ", source_dataset)
  )
  run_stage(
    "VALIDATE_OUTPUTS.R", annotation_dataset,
    annotation_environment,
    label = paste0(
      "Validating policy-matched alignment and annotation outputs: ",
      source_dataset, " [", policy_display_name, "]"
    )
  )
  run_stage(
    "09_organize_raw_dad_outputs.R", output_root,
    alignment_environment,
    paste0("Organizing collaborator outputs: ", source_dataset)
  )
  # Validate a second time after organization. The validator resolves exact
  # basenames recursively, so this confirms that nesting did not break any
  # scientific handoff or move a file into the wrong policy folder.
  run_stage(
    "VALIDATE_OUTPUTS.R", annotation_dataset,
    alignment_environment,
    label = paste0(
      "Validating organized outputs: ", source_dataset,
      " [", policy_display_name, "]"
    )
  )

  message(
    "Raw Agilent DAD pipeline complete. Alignment policy: ",
    policy_display_name
  )
  message(
    "Complete quantitative matrix: ",
    file.path(alignment_dir, "05_Quantitative_Matrices",
              "sample_by_aligned_feature_integrated_dad_area.csv")
  )
  message(
    "Strict primary quantitative matrix: ",
    file.path(
      alignment_dir, "05_Quantitative_Matrices",
      "sample_by_aligned_feature_integrated_dad_area_strict_primary.csv"
    )
  )
  message(
    "Family-summed sensitivity matrix: ",
    file.path(
      alignment_dir, "05_Quantitative_Matrices",
      "sample_by_feature_with_summed_unresolved_families_sensitivity.csv"
    )
  )
  message(
    "Alignment policy manifest: ",
    file.path(alignment_dir, "00_Run_Manifest",
              "raw_alignment_policy_manifest.csv")
  )
  message(
    "Library interpretation: ", file.path(
      annotation_dir, "01_Primary_Annotations", basename(expected_annotation)
    )
  )
  message(
    "Genus pseudochromatograms: ", file.path(
      final_dir, "02_Primary_Figures"
    )
  )
  message(
    "Start here: ", file.path(output_root, "00_START_HERE.csv")
  )
  message("Raw-standard library QC: ", library_qc_dir)
}
