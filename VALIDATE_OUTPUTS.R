## ========================================================================== ##
## COLLABORATOR-HANDOFF OUTPUT VALIDATION
## ========================================================================== ##
## Run after the production raw-DAD pipeline:
##   Rscript VALIDATE_OUTPUTS.R Rhododendron_RawDAD_RetainedUniqueness
##   Rscript VALIDATE_OUTPUTS.R Rhododendron_RawDAD_ExpandedRescue
##
## With no arguments, every completed data/*_Agilent raw-DAD policy route is
## checked. The script validates input fingerprints, one-row-per-feature
## outputs, dictionary keys, the exact wide-to-long abundance conversion, the
## NA-free confident-abundance/uncertainty-flag contract, and policy propagation
## manifest through annotation.
## It also confirms that every CSV can be opened, every PDF is nonempty, and
## each workflow includes parallel class, chromophore, and class-first hybrid
## pseudochromatograms.
## Visual PDF and workbook review remains a separate release step because file
## structure alone cannot detect overlapping labels or clipped plot elements.
## ========================================================================== ##

source("00_library_inputs.R")
library_bundle <- load_library_inputs()

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

read_output_csv <- function(path) {
  if (!file.exists(path)) {
    candidates <- list.files(
      dirname(path), pattern = paste0("^", gsub(
        "([][{}()+*^$|\\?.])", "\\\\\\1", basename(path)
      ), "$"), recursive = TRUE, full.names = TRUE
    )
    assert_true(length(candidates) == 1L,
                paste("Missing or non-unique required output:", path))
    path <- candidates
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

resolve_output_paths <- function(paths) {
  vapply(paths, function(path) {
    if (file.exists(path)) return(path)
    candidates <- list.files(
      dirname(path), pattern = paste0("^", gsub(
        "([][{}()+*^$|\\?.])", "\\\\\\1", basename(path)
      ), "$"), recursive = TRUE, full.names = TRUE
    )
    if (length(candidates) == 1L) candidates else path
  }, character(1))
}

same_numeric_matrix <- function(left, right, tolerance = 1e-10) {
  left <- as.matrix(left)
  right <- as.matrix(right)
  storage.mode(left) <- "double"
  storage.mode(right) <- "double"
  identical(dim(left), dim(right)) &&
    isTRUE(all.equal(left, right, tolerance = tolerance,
                     check.attributes = FALSE))
}

validate_readable_outputs <- function(folder) {
  csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE,
                          recursive = TRUE)
  pdf_files <- list.files(folder, pattern = "\\.pdf$", full.names = TRUE,
                          recursive = TRUE)
  assert_true(length(csv_files) > 0L,
              paste("No CSV outputs were found in", folder))
  invisible(lapply(csv_files, read_output_csv))
  if (length(pdf_files)) {
    pdf_sizes <- file.info(pdf_files)$size
    assert_true(all(is.finite(pdf_sizes) & pdf_sizes > 1000),
                paste("An empty or implausibly small PDF was found in", folder))
  }
  c(csv_files = length(csv_files), pdf_files = length(pdf_files))
}

# New runs record every native-to-library wavelength conversion. Keep this
# check optional so output folders created by an earlier pipeline release can
# still be audited, while enforcing the conservative interpolation contract
# whenever the provenance table is present.
validate_wavelength_alignment <- function(folder, dataset) {
  path <- file.path(folder, "wavelength_grid_alignment.csv")
  path <- resolve_output_paths(path)
  if (!file.exists(path)) return(invisible(FALSE))
  alignment <- read_output_csv(path)
  required <- c(
    "input_id", "native_wavelength_min_nm", "native_wavelength_max_nm",
    "native_wavelength_count", "native_median_step_nm",
    "native_maximum_gap_nm", "analysis_wavelength_min_nm",
    "analysis_wavelength_max_nm", "analysis_wavelength_count",
    "analysis_median_step_nm", "interpolation_required",
    "lower_endpoint_carry_nm", "upper_endpoint_carry_nm",
    "alignment_method"
  )
  assert_true(nrow(alignment) > 0L && all(required %in% names(alignment)),
              paste(dataset, "has an incomplete wavelength-grid audit."))
  assert_true(
    all(alignment$native_wavelength_count >= 2L) &&
      all(alignment$analysis_wavelength_count >= 2L) &&
      all(alignment$native_median_step_nm > 0) &&
      all(alignment$analysis_median_step_nm > 0) &&
      all(alignment$native_maximum_gap_nm <=
            2 * alignment$analysis_median_step_nm + 1e-8),
    paste(dataset, "contains a wavelength grid too sparse for interpolation.")
  )
  edge_limit <- pmax(alignment$native_median_step_nm,
                     alignment$analysis_median_step_nm)
  assert_true(
    all(alignment$lower_endpoint_carry_nm >= 0) &&
      all(alignment$upper_endpoint_carry_nm >= 0) &&
      all(alignment$lower_endpoint_carry_nm <= edge_limit + 1e-8) &&
      all(alignment$upper_endpoint_carry_nm <= edge_limit + 1e-8),
    paste(dataset, "uses wavelength extrapolation beyond one grid interval.")
  )
  invisible(TRUE)
}

validate_base_dataset <- function(
    dataset,
    input_file = file.path("data", paste0(dataset, "_MaxAbsorbances.csv")),
    analysis_mode = "maximum_absorbance_only",
    output_folder = dataset) {
  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", dataset))
  assert_true(file.exists(input_file), paste("Missing sample input:", input_file))
  assert_true(dir.exists(output_folder), paste("Missing output folder:", output_folder))

  sample_input <- read.csv(
    input_file, stringsAsFactors = FALSE, check.names = FALSE
  )
  comprehensive <- read_output_csv(file.path(
    output_folder, paste0(slug, "_comprehensive_results.csv")
  ))
  metadata <- read_output_csv(file.path(
    output_folder, paste0(slug, "_peak_annotation_metadata.csv")
  ))
  annotated_input <- read_output_csv(file.path(
    output_folder, paste0(slug, "_original_data_with_peak_annotations.csv")
  ))

  input_ids <- as.character(sample_input[[1]])
  assert_true(!anyDuplicated(input_ids),
              paste(dataset, "input peak identifiers are not unique."))
  assert_true(identical(input_ids, as.character(comprehensive$peak_id)),
              paste(dataset, "comprehensive results do not preserve peak order."))
  assert_true(identical(input_ids, as.character(metadata$peak_id)),
              paste(dataset, "annotation metadata do not preserve peak order."))
  assert_true(identical(input_ids, as.character(annotated_input[[1]])),
              paste(dataset, "annotated original data do not preserve peak order."))
  assert_true(ncol(annotated_input) >= ncol(sample_input),
              paste(dataset, "annotated original data lost submitted columns."))
  assert_true(
    same_numeric_matrix(
      sample_input[, -1, drop = FALSE],
      annotated_input[, seq.int(2, ncol(sample_input)), drop = FALSE]
    ),
    paste(dataset, "annotated original data changed submitted spectral values.")
  )

  manifest <- read_output_csv(file.path(
    output_folder, "library_input_validation.csv"
  ))
  current_manifest <- library_bundle$input_manifest
  manifest_position <- match(current_manifest$input_role, manifest$input_role)
  assert_true(!anyNA(manifest_position) && all(
    current_manifest$md5 == manifest$md5[manifest_position]
  ), paste(dataset, "outputs do not match the current three-file library."))

  parameters <- read_output_csv(file.path(output_folder,
                                           "analysis_parameters.csv"))
  assert_true("sample_input_md5" %in% names(parameters) &&
                nrow(parameters) == 1L &&
                identical(as.character(parameters$sample_input_md5),
                          unname(tools::md5sum(input_file))),
              paste(dataset, "outputs do not match the current sample input."))
  validate_wavelength_alignment(output_folder, dataset)

  # Every compound-first run must retain the three parallel reporting views.
  # These differ only in how supported spectral evidence is displayed; they do
  # not refit the model or use retention time as a predictor.
  required_base_pseudochromatograms <- file.path(
    output_folder,
    paste0(
      slug,
      c("_class_pseudochromatogram.pdf",
        "_chromophore_pseudochromatogram.pdf",
        "_hybrid_class_chromophore_pseudochromatogram.pdf",
        "_class_signal_composition_mirror.pdf",
        "_chromophore_signal_composition_mirror.pdf",
        "_hybrid_signal_composition_mirror.pdf")
    )
  )
  required_base_pseudochromatograms <- resolve_output_paths(
    required_base_pseudochromatograms
  )
  assert_true(
    all(file.exists(required_base_pseudochromatograms)) &&
      all(file.info(required_base_pseudochromatograms)$size > 1000),
    paste(dataset, "is missing a filled or mirror class/chromophore figure.")
  )
  scaling_audit <- read_output_csv(file.path(
    output_folder, paste0(slug, "_pseudochromatogram_scaling_audit.csv")
  ))
  assert_true(
    nrow(scaling_audit) == 3L &&
      all(scaling_audit$maximum_reconstruction_error <= 1e-6) &&
      all(scaling_audit$black_sum_not_below_any_component),
    paste(dataset, "maximum-table pseudochromatogram scaling is inconsistent.")
  )

  file_counts <- validate_readable_outputs(output_folder)
  data.frame(
    dataset = dataset,
    analysis_mode = analysis_mode,
    input_peaks = nrow(sample_input),
    csv_outputs = unname(file_counts[["csv_files"]]),
    pdf_outputs = unname(file_counts[["pdf_files"]]),
    landscape_samples = NA_integer_,
    homologous_features = NA_integer_,
    abundance_cells = NA_integer_,
    stringsAsFactors = FALSE
  )
}

validate_landscape_dataset <- function(
    dataset, summary_row,
    landscape_input_folder = file.path("data", paste0(dataset, "_Landscape")),
    peak_input_file = file.path("data", paste0(dataset, "_MaxAbsorbances.csv"))) {
  if (!dir.exists(landscape_input_folder)) return(summary_row)

  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", dataset))
  landscape_folder <- paste0(dataset, "_Landscape_Analysis")
  variant_folder <- paste0(dataset, "_Spectral_Variants")
  assert_true(dir.exists(landscape_folder),
              paste("Missing landscape output folder:", landscape_folder))
  assert_true(dir.exists(variant_folder),
              paste("Missing variant output folder:", variant_folder))

  required_variant_pseudochromatograms <- file.path(
    variant_folder,
    c("spectral_variant_prevalence_pseudochromatogram.pdf",
      "spectral_variant_class_pseudochromatogram.pdf",
      "spectral_variant_chromophore_pseudochromatogram.pdf",
      "spectral_variant_hybrid_class_chromophore_pseudochromatogram.pdf",
      "spectral_variant_prevalence_signal_composition_mirror.pdf",
      "spectral_variant_class_signal_composition_mirror.pdf",
      "spectral_variant_chromophore_signal_composition_mirror.pdf",
      "spectral_variant_hybrid_signal_composition_mirror.pdf")
  )
  assert_true(
    all(file.exists(required_variant_pseudochromatograms)) &&
      all(file.info(required_variant_pseudochromatograms)$size > 1000),
    paste(dataset, "is missing a filled or mirror variant pseudochromatogram.")
  )
  hybrid_components <- read_output_csv(file.path(
    variant_folder, "variant_hybrid_pseudochromatogram_components.csv"
  ))
  variant_scaling_audit <- read_output_csv(file.path(
    variant_folder, "variant_pseudochromatogram_scaling_audit.csv"
  ))
  assert_true(
    nrow(variant_scaling_audit) == 4L &&
      all(variant_scaling_audit$maximum_component_fraction_sum_error <=
            1e-10) &&
      all(variant_scaling_audit$black_sum_not_below_any_peak_apex),
    paste(dataset, "variant pseudochromatogram scaling is inconsistent.")
  )
  assert_true(
    all(c("class_interpretation", "chromophore_interpretation",
          "hybrid_interpretation") %in% names(hybrid_components)),
    paste(dataset, "variant hybrid component audit is incomplete.")
  )

  input_peaks <- read.csv(
    peak_input_file,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  peak_ids <- as.character(input_peaks[[1]])
  landscape_samples <- tools::file_path_sans_ext(basename(list.files(
    landscape_input_folder, pattern = "\\.csv$", full.names = TRUE
  )))

  decisions <- read_output_csv(file.path(
    variant_folder, "retention_time_variant_decisions.csv"
  ))
  assert_true(identical(peak_ids, as.character(decisions$peak_id)),
              paste(dataset, "variant decisions do not preserve every peak."))

  abundance_wide <- read_output_csv(file.path(
    variant_folder, "sample_by_homologous_peak_summed_absorbance.csv"
  ))
  feature_ids <- names(abundance_wide)[-1]
  assert_true(identical(sort(landscape_samples),
                        sort(as.character(abundance_wide[[1]]))),
              paste(dataset, "abundance matrix does not contain every sample."))
  assert_true(!anyDuplicated(feature_ids),
              paste(dataset, "abundance feature identifiers are not unique."))

  dictionary <- read_output_csv(file.path(
    variant_folder, "homologous_peak_abundance_dictionary.csv"
  ))
  metadata <- read_output_csv(file.path(
    variant_folder, "homologous_peak_annotation_metadata.csv"
  ))
  assert_true(identical(feature_ids, as.character(dictionary$homologous_peak_id)),
              paste(dataset, "abundance dictionary order does not match the matrix."))
  assert_true(identical(feature_ids, as.character(metadata$homologous_peak_id)),
              paste(dataset, "annotation metadata order does not match the matrix."))

  abundance_long <- read_output_csv(file.path(
    variant_folder,
    "sample_by_homologous_peak_abundance_with_annotations_long.csv"
  ))
  expected_rows <- nrow(abundance_wide) * length(feature_ids)
  assert_true(nrow(abundance_long) == expected_rows,
              paste(dataset, "annotated long abundance table has the wrong size."))
  long_key <- paste(abundance_long$sample_id,
                    abundance_long$homologous_peak_id, sep = "||")
  assert_true(!anyDuplicated(long_key),
              paste(dataset, "annotated long abundance keys are not unique."))

  expected_values <- as.vector(t(as.matrix(
    abundance_wide[, -1, drop = FALSE]
  )))
  assert_true(isTRUE(all.equal(
    expected_values, abundance_long$summed_absorbance,
    tolerance = 1e-10, check.attributes = FALSE
  )), paste(dataset, "wide and long abundance values differ."))
  assert_true(all(is.na(abundance_long$summed_absorbance[
    abundance_long$abundance_status == "Detected but unresolved"
  ])), paste(dataset, "unresolved abundance cells must be NA."))
  assert_true(all(is.na(abundance_long$summed_absorbance[
    abundance_long$abundance_status ==
      "Not evaluable: outside retention-time coverage"
  ])), paste(dataset, "RT-unavailable abundance cells must be NA."))
  assert_true(all(abundance_long$summed_absorbance[
    abundance_long$abundance_status == "Not detected"
  ] == 0), paste(dataset, "not-detected abundance cells must be zero."))
  quantified <- abundance_long$summed_absorbance[
    abundance_long$abundance_status == "Quantified"
  ]
  assert_true(all(is.finite(quantified) & quantified > 0),
              paste(dataset, "quantified abundance cells must be positive."))

  assignment_audit <- read_output_csv(file.path(
    variant_folder, "sample_peak_homologous_assignment_and_abundance.csv"
  ))
  assert_true(nrow(assignment_audit) ==
                length(peak_ids) * length(landscape_samples),
              paste(dataset, "assignment audit is not peak-by-sample complete."))
  assignment_key <- paste(assignment_audit$peak_id,
                          assignment_audit$sample_id, sep = "||")
  assert_true(!anyDuplicated(assignment_key),
              paste(dataset, "assignment audit keys are not unique."))
  assert_true(
    all(c("variant_letter", "assignment_source",
          "retention_time_window_available",
          "distance_to_assigned_variant") %in% names(assignment_audit)),
    paste(dataset, "assignment audit lacks variant-attachment evidence.")
  )
  weakly_attached <- assignment_audit$assignment_source ==
    "SNR-weighted attachment"
  assert_true(
    all(is.finite(
      assignment_audit$distance_to_assigned_variant[weakly_attached]
    )) && all(
      assignment_audit$variant_letter[weakly_attached] != "Unresolved"
    ),
    paste(dataset, "contains a non-finite or unresolved weak attachment.")
  )
  explicitly_nonfinite <- grepl(
    "no finite variant attachment distance|insufficient finite variant distances",
    assignment_audit$assignment_source
  )
  assert_true(
    all(assignment_audit$variant_letter[explicitly_nonfinite] == "Unresolved"),
    paste(dataset, "assigned a spectrum lacking finite variant evidence.")
  )
  outside_rt_coverage <- !assignment_audit$retention_time_window_available
  assert_true(
    all(assignment_audit$variant_letter[outside_rt_coverage] == "Unresolved") &&
      all(grepl(
        "outside landscape coverage",
        assignment_audit$assignment_source[outside_rt_coverage]
      )),
    paste(dataset, "misclassified an RT-unavailable sample/peak pair.")
  )
  validate_wavelength_alignment(landscape_folder, dataset)

  sample_figure_folder <- file.path(
    variant_folder, "Sample_Pseudochromatograms"
  )
  assert_true(dir.exists(sample_figure_folder),
              paste(dataset, "is missing sample pseudochromatograms."))
  sample_figure_manifest <- read_output_csv(file.path(
    sample_figure_folder, "sample_pseudochromatogram_manifest.csv"
  ))
  sample_figure_audit <- read_output_csv(file.path(
    sample_figure_folder,
    "sample_pseudochromatogram_peak_annotations.csv"
  ))
  assert_true(identical(sort(landscape_samples),
                        sort(as.character(sample_figure_manifest$sample_id))),
              paste(dataset, "sample figure manifest does not cover every sample."))
  assert_true(all(sample_figure_manifest$pdf_pages == 9L),
              paste(dataset, "sample pseudochromatograms must have nine pages."))
  expected_sample_pdfs <- file.path(
    sample_figure_folder, sample_figure_manifest$pdf_file
  )
  assert_true(all(file.exists(expected_sample_pdfs)) &&
                all(file.info(expected_sample_pdfs)$size > 1000),
              paste(dataset, "a sample pseudochromatogram PDF is missing or empty."))
  assert_true(nrow(sample_figure_audit) ==
                length(peak_ids) * length(landscape_samples),
              paste(dataset, "sample pseudochromatogram audit is incomplete."))
  sample_figure_key <- paste(sample_figure_audit$peak_id,
                             sample_figure_audit$sample_id, sep = "||")
  assert_true(!anyDuplicated(sample_figure_key),
              paste(dataset, "sample pseudochromatogram audit keys are not unique."))
  assert_true(all(!sample_figure_audit$retention_time_used_for_library_matching),
              paste(dataset, "sample figures incorrectly report RT prediction use."))
  assert_true(
    all(c("class_interpretation", "chromophore_interpretation",
          "hybrid_interpretation", "plot_class_interpretation",
          "plot_chromophore_interpretation",
          "plot_hybrid_interpretation") %in% names(sample_figure_audit)),
    paste(dataset, "sample figure audit is missing a parallel interpretation view.")
  )
  retained_audit <- sample_figure_audit[
    sample_figure_audit$included_once_in_sample_pseudochromatogram, ]
  retained_rt_key <- paste(retained_audit$sample_id,
                           retained_audit$pseudo_apex_key, sep = "||")
  assert_true(!anyDuplicated(retained_rt_key),
              paste(dataset, "plots duplicate a selected sample apex."))
  positive_height_samples <- split(
    retained_audit$sample_relative_height_percent[
      retained_audit$sample_relative_height_percent > 0
    ],
    retained_audit$sample_id[
      retained_audit$sample_relative_height_percent > 0
    ]
  )
  assert_true(all(vapply(positive_height_samples, function(x) {
    isTRUE(all.equal(max(x), 100, tolerance = 1e-8))
  }, logical(1))), paste(dataset, "sample plot heights are not normalized to 100."))

  landscape_annotation <- read_output_csv(file.path(
    variant_folder,
    paste0(slug, "_original_data_with_landscape_annotations.csv")
  ))
  assert_true(identical(peak_ids, as.character(landscape_annotation[[1]])),
              paste(dataset, "landscape annotations do not preserve every peak."))

  landscape_counts <- validate_readable_outputs(landscape_folder)
  variant_counts <- validate_readable_outputs(variant_folder)
  sample_figure_counts <- validate_readable_outputs(sample_figure_folder)
  summary_row$csv_outputs <- summary_row$csv_outputs +
    unname(landscape_counts[["csv_files"]]) +
    unname(variant_counts[["csv_files"]]) +
    unname(sample_figure_counts[["csv_files"]])
  summary_row$pdf_outputs <- summary_row$pdf_outputs +
    unname(landscape_counts[["pdf_files"]]) +
    unname(variant_counts[["pdf_files"]]) +
    unname(sample_figure_counts[["pdf_files"]])
  summary_row$landscape_samples <- length(landscape_samples)
  summary_row$homologous_features <- length(feature_ids)
  summary_row$abundance_cells <- expected_rows
  if (identical(summary_row$analysis_mode, "maximum_absorbance_only")) {
    summary_row$analysis_mode <- "maximum_absorbance_plus_landscapes"
  }
  summary_row
}

validate_landscape_only_dataset <- function(source_dataset) {
  dataset <- paste0(source_dataset, "_LandscapeOnly")
  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", dataset))
  discovery_folder <- paste0(dataset, "_Discovery")
  landscape_input_folder <- file.path(
    "data", paste0(source_dataset, "_Landscape")
  )
  reference_file <- file.path(
    discovery_folder,
    paste0(slug, "_landscape_derived_reference_spectra.csv")
  )
  reference_file <- resolve_output_paths(reference_file)
  assert_true(dir.exists(discovery_folder),
              paste("Missing discovery output folder:", discovery_folder))
  assert_true(file.exists(reference_file),
              paste("Missing landscape-derived reference table:", reference_file))

  discovery <- read_output_csv(file.path(
    discovery_folder, "landscape_only_peak_discovery_summary.csv"
  ))
  reference <- read_output_csv(reference_file)
  detections <- read_output_csv(file.path(
    discovery_folder, "sample_landscape_peak_detections.csv"
  ))
  assert_true(identical(as.character(discovery$peak_id),
                        as.character(reference[[1]])),
              paste(dataset, "discovery and reference peak order differ."))
  assert_true(all(detections$peak_id %in% discovery$peak_id),
              paste(dataset, "a sample detection lacks a discovery peak."))
  assert_true(all(discovery$n_detected_samples >= 1L),
              paste(dataset, "contains an unsupported zero-observation peak."))
  assert_true(all(!discovery$retention_time_used_for_library_matching),
              paste(dataset, "incorrectly reports RT use in library matching."))
  validate_wavelength_alignment(discovery_folder, dataset)

  summary_row <- validate_base_dataset(
    dataset, reference_file, analysis_mode = "landscape_only"
  )
  parameters <- read_output_csv(file.path(dataset, "analysis_parameters.csv"))
  assert_true("analysis_input_mode" %in% names(parameters) &&
                identical(as.character(parameters$analysis_input_mode),
                          "landscape_only"),
              paste(dataset, "does not record landscape-only provenance."))
  summary_row <- validate_landscape_dataset(
    dataset, summary_row,
    landscape_input_folder = landscape_input_folder,
    peak_input_file = reference_file
  )
  discovery_counts <- validate_readable_outputs(discovery_folder)
  summary_row$csv_outputs <- summary_row$csv_outputs +
    unname(discovery_counts[["csv_files"]])
  summary_row$pdf_outputs <- summary_row$pdf_outputs +
    unname(discovery_counts[["pdf_files"]])
  summary_row
}

# Validate the raw Agilent route independently from the older landscape routes.
# The raw workflow creates its own peak list and quantitative matrix, then sends
# only representative spectral shapes to the library-only annotation stage.
validate_raw_dad_dataset <- function(source_dataset, alignment_policy) {
  policy_folder_label <- c(
    retained_uniqueness = "RetainedUniqueness",
    expanded_rescue = "ExpandedRescue"
  )[alignment_policy]
  assert_true(
    length(policy_folder_label) == 1L && !is.na(policy_folder_label),
    paste("Unknown raw DAD alignment policy:", alignment_policy)
  )
  raw_input_folder <- file.path("data", paste0(source_dataset, "_Agilent"))
  output_root <- paste0(
    source_dataset, "_RawDAD_", policy_folder_label
  )
  alignment_folder <- Sys.getenv(
    "UV_RAW_ALIGNMENT_OUTPUT_DIR",
    unset = file.path(output_root, "01_Alignment_QC")
  )
  annotation_dataset <- paste0(
    source_dataset, "_RawDAD_", policy_folder_label
  )
  annotation_folder <- Sys.getenv(
    "UV_OUTPUT_DIR",
    unset = file.path(output_root, "02_Library_Annotation")
  )
  library_qc_folder <- file.path(output_root, "00_Library_QC")
  source_slug <- tolower(gsub("[^A-Za-z0-9]+", "_", source_dataset))
  reference_file <- file.path(
    alignment_folder,
    paste0(source_slug, "_rawdad_aligned_reference_spectra.csv")
  )
  reference_file <- resolve_output_paths(reference_file)

  assert_true(dir.exists(raw_input_folder),
              paste("Missing raw Agilent input folder:", raw_input_folder))
  assert_true(dir.exists(alignment_folder),
              paste("Missing raw DAD alignment folder:", alignment_folder))
  required_library_qc <- file.path(library_qc_folder, c(
    "Absorbance_StandardLibrary_RawDAD.csv",
    "raw_dad_library_manifest.csv",
    "raw_dad_identity_corrections.csv",
    "standard_curve_manifest.csv",
    "raw_standard_compound_admission_decisions.csv",
    "raw_standard_replicate_extraction_qc.csv",
    "raw_standard_identity_conflict_audit.csv",
    "standard_curve_readiness.csv",
    "library_version_validation_summary.csv",
    "all_raw_standard_peak_extractions.pdf",
    "raw_vs_authenticated_standard_spectra_atlas.pdf",
    "raw_standard_identity_conflict_audit.pdf",
    "library_version_validation_comparison.pdf",
    "library_version_separation_comparison.pdf"
  ))
  assert_true(
    all(file.exists(required_library_qc)) &&
      all(file.info(required_library_qc)$size > 100),
    paste(source_dataset, "is missing its raw-standard library QC snapshot.")
  )
  identity_snapshot <- read_output_csv(file.path(
    library_qc_folder, "raw_standard_identity_conflict_audit.csv"
  ))
  assert_true(
    nrow(identity_snapshot) == 23L &&
      "identity_conflict_flag" %in% names(identity_snapshot) &&
      !any(identity_snapshot$identity_conflict_flag),
    paste(source_dataset,
          "contains an incomplete or unresolved raw-standard identity audit.")
  )
  standard_curve_snapshot <- read_output_csv(file.path(
    library_qc_folder, "standard_curve_readiness.csv"
  ))
  assert_true(
    nrow(standard_curve_snapshot) == 1L &&
      all(c("calibration_rows_available",
            "quantitative_calibration_enabled", "status") %in%
            names(standard_curve_snapshot)),
    paste(source_dataset, "has an invalid standard-curve readiness record.")
  )
  active_snapshot_md5 <- unname(tools::md5sum(required_library_qc[1]))
  active_current_md5 <- unname(tools::md5sum(
    "library/Absorbance_StandardLibrary_RawDAD.csv"
  ))
  assert_true(
    identical(active_snapshot_md5, active_current_md5),
    paste(source_dataset, "raw-standard library snapshot is not current.")
  )

  sample_directories <- list.dirs(raw_input_folder, recursive = FALSE,
                                   full.names = TRUE)
  sample_directories <- sample_directories[grepl(
    "[.]D$", sample_directories, ignore.case = TRUE
  )]
  raw_sample_ids <- tools::file_path_sans_ext(basename(sample_directories))
  assert_true(length(raw_sample_ids) > 0L,
              paste(source_dataset, "has no .D sample folders."))
  input_exclusions <- read_output_csv(file.path(
    alignment_folder, "raw_agilent_input_exclusions.csv"
  ))
  required_exclusion_fields <- c(
    "sample_id", "sample_directory", "dad_file", "file_size_bytes",
    "exclusion_reason", "exclusion_action"
  )
  assert_true(
    all(required_exclusion_fields %in% names(input_exclusions)) &&
      !anyDuplicated(input_exclusions$sample_id) &&
      all(input_exclusions$sample_id %in% raw_sample_ids) &&
      all(nzchar(input_exclusions$exclusion_reason)),
    paste(source_dataset, "has an invalid raw-input exclusion audit.")
  )
  sample_ids <- setdiff(raw_sample_ids, input_exclusions$sample_id)
  assert_true(length(sample_ids) > 0L,
              paste(source_dataset, "has no readable raw DAD samples."))

  # Reuse the compound-first validator for the representative spectra and
  # current-library fingerprint, then add raw alignment/matrix contracts.
  summary_row <- validate_base_dataset(
    annotation_dataset, reference_file,
    analysis_mode = "raw_agilent_dad_joint_alignment",
    output_folder = annotation_folder
  )

  policy_manifest <- read_output_csv(file.path(
    alignment_folder, "raw_alignment_policy_manifest.csv"
  ))
  required_policy_fields <- c(
    "source_dataset", "alignment_policy",
    "alignment_policy_display_name", "total_raw_sample_folders",
    "excluded_unreadable_inputs", "total_samples",
    "rescue_review_maximum_source_group_samples",
    "rescue_target_minimum_samples",
    "rescue_requires_target_larger_than_source",
    "recommended_primary_matrix", "downstream_annotation_dataset",
    "retention_time_used_for_library_annotation"
  )
  assert_true(
    nrow(policy_manifest) == 1L &&
      all(required_policy_fields %in% names(policy_manifest)) &&
      identical(as.character(policy_manifest$source_dataset), source_dataset) &&
      identical(as.character(policy_manifest$alignment_policy),
                alignment_policy) &&
      identical(as.character(policy_manifest$downstream_annotation_dataset),
                annotation_dataset),
    paste(source_dataset, alignment_policy,
          "has a missing or mismatched alignment-policy manifest.")
  )
  expected_review_maximum <- if (alignment_policy == "retained_uniqueness") {
    min(length(sample_ids), 3L)
  } else if (length(sample_ids) < 20L) {
    max(1L, ceiling(0.20 * length(sample_ids)))
  } else 20L
  assert_true(
    policy_manifest$total_raw_sample_folders == length(raw_sample_ids) &&
      policy_manifest$excluded_unreadable_inputs == nrow(input_exclusions) &&
      policy_manifest$total_samples == length(sample_ids) &&
      policy_manifest$rescue_review_maximum_source_group_samples ==
        expected_review_maximum &&
      policy_manifest$rescue_target_minimum_samples == 4L &&
      isTRUE(policy_manifest$rescue_requires_target_larger_than_source) &&
      isTRUE(policy_manifest$iterative_rescue_until_convergence) &&
      isTRUE(policy_manifest$rescue_targets_restricted_to_original_anchors) &&
      identical(
        as.character(policy_manifest$recommended_primary_matrix),
        "sample_by_aligned_feature_integrated_dad_area_strict_primary.csv"
      ) &&
      !isTRUE(policy_manifest$retention_time_used_for_library_annotation),
    paste(source_dataset, alignment_policy,
          "does not implement its documented alignment-policy preset.")
  )

  annotation_parameters <- read_output_csv(file.path(
    annotation_folder, "analysis_parameters.csv"
  ))
  assert_true(
    "source_alignment_policy" %in% names(annotation_parameters) &&
      identical(
        as.character(annotation_parameters$source_alignment_policy),
        alignment_policy
      ),
    paste(source_dataset, alignment_policy,
          "was not propagated into downstream library annotation.")
  )

  features <- read_output_csv(file.path(
    alignment_folder, "aligned_feature_summary.csv"
  ))
  annotated_features <- read_output_csv(file.path(
    alignment_folder, "aligned_feature_summary_with_library_annotations.csv"
  ))
  detections <- read_output_csv(file.path(
    alignment_folder, "sample_peak_detections_and_alignment.csv"
  ))
  ambiguous_detections <- read_output_csv(file.path(
    alignment_folder, "ambiguous_unforced_detections.csv"
  ))
  reconciliation <- read_output_csv(file.path(
    alignment_folder, "rare_fragment_reconciliation_audit.csv"
  ))
  reconciliation_iterations <- read_output_csv(file.path(
    alignment_folder, "rare_fragment_reconciliation_iteration_audit.csv"
  ))
  split_merges <- read_output_csv(file.path(
    alignment_folder, "probable_split_feature_merge_audit.csv"
  ))
  neighboring_pairs <- read_output_csv(file.path(
    alignment_folder, "neighboring_feature_pair_decisions.csv"
  ))
  competition_families <- read_output_csv(file.path(
    alignment_folder, "unresolved_competition_family_dictionary.csv"
  ))
  abundance_wide <- read_output_csv(file.path(
    alignment_folder, "sample_by_aligned_feature_integrated_dad_area.csv"
  ))
  primary_filtered_wide <- read_output_csv(file.path(
    alignment_folder,
    "sample_by_aligned_feature_integrated_dad_area_primary_filtered.csv"
  ))
  strict_primary_wide <- read_output_csv(file.path(
    alignment_folder,
    "sample_by_aligned_feature_integrated_dad_area_strict_primary.csv"
  ))
  uncertainty_flag_wide <- read_output_csv(file.path(
    alignment_folder,
    "sample_by_aligned_feature_assignment_uncertainty_flag.csv"
  ))
  strict_uncertainty_flag_wide <- read_output_csv(file.path(
    alignment_folder,
    paste0(
      "sample_by_aligned_feature_assignment_uncertainty_flag_",
      "strict_primary.csv"
    )
  ))
  upper_bound_wide <- read_output_csv(file.path(
    alignment_folder,
    "sample_by_aligned_feature_integrated_dad_area_upper_bound.csv"
  ))
  strict_upper_bound_wide <- read_output_csv(file.path(
    alignment_folder,
    paste0(
      "sample_by_aligned_feature_integrated_dad_area_upper_bound_",
      "strict_primary.csv"
    )
  ))
  family_sensitivity_wide <- read_output_csv(file.path(
    alignment_folder,
    "sample_by_feature_with_summed_unresolved_families_sensitivity.csv"
  ))
  abundance_long <- read_output_csv(file.path(
    alignment_folder, "sample_by_aligned_feature_integrated_dad_area_long.csv"
  ))
  peak_detection_summary <- read_output_csv(file.path(
    alignment_folder, "sample_peak_detection_summary.csv"
  ))
  peak_profile_manifest <- read_output_csv(file.path(
    alignment_folder, "peak_picking_profile_atlas_manifest.csv"
  ))
  peak_picker_sensitivity <- read_output_csv(file.path(
    alignment_folder, "peak_picking_parameter_sensitivity.csv"
  ))
  raw_alignment_parameters <- read_output_csv(file.path(
    alignment_folder, "raw_dad_alignment_parameters.csv"
  ))
  aligned_peak_map <- read_output_csv(file.path(
    alignment_folder, "sample_by_aligned_peak_map_data.csv"
  ))
  feature_color_dictionary <- read_output_csv(file.path(
    alignment_folder, "aligned_feature_figure_color_dictionary.csv"
  ))
  identity_atlas_manifest <- read_output_csv(file.path(
    alignment_folder,
    "sample_by_aligned_feature_identity_atlas_manifest.csv"
  ))
  feature_ids <- as.character(features$aligned_feature_id)
  matrix_feature_ids <- names(abundance_wide)[-1]

  assert_true(
    identical(as.character(peak_profile_manifest$sample_id), sort(sample_ids)) &&
      identical(peak_profile_manifest$sample_order,
                seq_along(sort(sample_ids))) &&
      all(peak_profile_manifest$atlas_page == floor(
        (seq_along(sort(sample_ids)) - 1) * ceiling(length(sample_ids) / 6) /
          length(sample_ids)
      ) + 1L) &&
      all(peak_profile_manifest$total_atlas_pages ==
            ceiling(length(sample_ids) / 6)),
    paste(source_dataset,
          "peak-picking profile atlas does not cover every sample once.")
  )

  assert_true(!anyDuplicated(feature_ids) &&
                identical(feature_ids, matrix_feature_ids),
              paste(source_dataset,
                    "raw DAD feature identifiers are duplicated or reordered."))
  assert_true(identical(feature_ids,
                        as.character(annotated_features$aligned_feature_id)),
              paste(source_dataset,
                    "annotated feature dictionary does not preserve feature order."))
  expected_conservative_features <- feature_ids[
    as.logical(features$include_in_primary_filtered_matrix)
  ]
  assert_true(identical(names(primary_filtered_wide)[-1],
                        expected_conservative_features) &&
                identical(as.character(primary_filtered_wide[[1]]),
                          as.character(abundance_wide[[1]])) &&
                same_numeric_matrix(
                  primary_filtered_wide[, -1, drop = FALSE],
                  abundance_wide[, expected_conservative_features,
                                 drop = FALSE]
                ),
              paste(source_dataset,
                    "primary-filtered raw DAD matrix does not match its feature flags."))
  expected_strict_features <- feature_ids[
    as.logical(features$include_in_strict_primary_matrix)
  ]
  assert_true(
    identical(names(strict_primary_wide)[-1], expected_strict_features) &&
      identical(as.character(strict_primary_wide[[1]]),
                as.character(abundance_wide[[1]])) &&
      same_numeric_matrix(
        strict_primary_wide[, -1, drop = FALSE],
        abundance_wide[, expected_strict_features, drop = FALSE]
      ),
    paste(source_dataset,
          "strict-primary matrix does not match its feature flags.")
  )
  assert_true(
    identical(names(uncertainty_flag_wide), names(abundance_wide)) &&
      identical(names(upper_bound_wide), names(abundance_wide)) &&
      identical(as.character(uncertainty_flag_wide[[1]]),
                as.character(abundance_wide[[1]])) &&
      identical(as.character(upper_bound_wide[[1]]),
                as.character(abundance_wide[[1]])) &&
      identical(names(strict_uncertainty_flag_wide),
                names(strict_primary_wide)) &&
      identical(names(strict_upper_bound_wide), names(strict_primary_wide)) &&
      same_numeric_matrix(
        strict_uncertainty_flag_wide[, -1, drop = FALSE],
        uncertainty_flag_wide[, expected_strict_features, drop = FALSE]
      ) &&
      same_numeric_matrix(
        strict_upper_bound_wide[, -1, drop = FALSE],
        upper_bound_wide[, expected_strict_features, drop = FALSE]
      ),
    paste(source_dataset,
          "uncertainty flags or upper bounds do not match the matrix subsets.")
  )

  family_members <- character()
  if (nrow(competition_families)) {
    family_members <- unique(unlist(strsplit(
      competition_families$member_features, ";", fixed = TRUE
    )))
  }
  assert_true(
    all(is.na(features$unresolved_competition_family[
      features$include_in_strict_primary_matrix
    ])) &&
      identical(
        sort(features$aligned_feature_id[
          !is.na(features$unresolved_competition_family)
        ]), sort(family_members)
      ),
    paste(source_dataset,
          "unresolved competition families do not match strict exclusions.")
  )
  expected_sensitivity_columns <- c(
    features$aligned_feature_id[is.na(features$unresolved_competition_family)],
    competition_families$competition_family_id
  )
  assert_true(
    setequal(names(family_sensitivity_wide)[-1],
             expected_sensitivity_columns) &&
      identical(as.character(family_sensitivity_wide[[1]]),
                as.character(abundance_wide[[1]])),
    paste(source_dataset,
          "family-level sensitivity matrix has incorrect columns or samples.")
  )
  if (nrow(competition_families)) {
    for (row in seq_len(nrow(competition_families))) {
      family_id <- competition_families$competition_family_id[row]
      members <- strsplit(
        competition_families$member_features[row], ";", fixed = TRUE
      )[[1]]
      member_values <- as.matrix(abundance_wide[, members, drop = FALSE])
      storage.mode(member_values) <- "double"
      expected_family <- rowSums(member_values)
      assert_true(isTRUE(all.equal(
        expected_family, as.numeric(family_sensitivity_wide[[family_id]]),
        tolerance = 1e-10, check.attributes = FALSE
      )), paste(source_dataset, family_id,
                "sensitivity-family sums are incorrect."))
    }
  }

  accepted_split_merges <- split_merges$automatic_merge %in% TRUE
  if (any(accepted_split_merges)) {
    assert_true(all(
      split_merges$same_sample_doublet_count[accepted_split_merges] == 0L &
        split_merges$consensus_rt_difference_minutes[accepted_split_merges] <=
          0.035 + 1e-10 &
        split_merges$consensus_spectral_distance[accepted_split_merges] <=
          0.12 + 1e-10 &
        split_merges$combined_rt_span_minutes[accepted_split_merges] <=
          0.10 + 1e-10 &
        split_merges$wavelength_bootstrap_equivalence_support[
          accepted_split_merges
        ] >= 0.85 - 1e-10
    ), paste(source_dataset,
             "contains an automatic split-feature merge that failed a gate."))
  }
  doublet_rows <- neighboring_pairs$pair_decision ==
    "Protected resolved double peak: same-sample co-occurrence"
  assert_true(
    !any(doublet_rows) || all(
      neighboring_pairs$same_sample_doublet_count[doublet_rows] > 0L
    ),
    paste(source_dataset,
          "contains a protected doublet without same-sample evidence.")
  )

  policy_tables <- list(
    aligned_feature_summary = features,
    annotated_feature_summary = annotated_features,
    detections = detections,
    ambiguous_detections = ambiguous_detections,
    reconciliation = reconciliation,
    reconciliation_iterations = reconciliation_iterations,
    split_merges = split_merges,
    neighboring_pairs = neighboring_pairs,
    competition_families = competition_families,
    abundance_long = abundance_long
  )
  for (table_name in names(policy_tables)) {
    table <- policy_tables[[table_name]]
    assert_true(
      "alignment_policy" %in% names(table) &&
        all(table$alignment_policy == alignment_policy),
      paste(source_dataset, table_name,
            "does not preserve the selected alignment policy.")
    )
  }
  assert_true(identical(sort(sample_ids),
                        sort(as.character(abundance_wide[[1]]))),
              paste(source_dataset,
                    "raw DAD abundance matrix does not contain every .D sample."))
  assert_true(nrow(abundance_long) == length(sample_ids) * length(feature_ids),
              paste(source_dataset, "raw DAD long matrix is incomplete."))
  assert_true(all(detections$sample_id %in% sample_ids),
              paste(source_dataset,
                    "raw DAD detection audit contains an unknown sample."))
  assert_true(
    identical(sort(as.character(peak_detection_summary$sample_id)),
              sort(sample_ids)) &&
      sum(peak_detection_summary$detected_peaks) == nrow(detections) &&
      sum(peak_detection_summary$core_detections) ==
        sum(detections$core_detection) &&
      sum(peak_detection_summary$assigned_peaks) ==
        sum(!is.na(detections$aligned_feature_id)),
    paste(source_dataset, "peak-picking sample summary is inconsistent.")
  )
  assert_true(
    sum(peak_picker_sensitivity$is_deployed_setting) == 1L &&
      peak_picker_sensitivity$total_candidate_peaks[
        peak_picker_sensitivity$is_deployed_setting
      ] >= nrow(detections),
    paste(source_dataset, "peak-picking sensitivity audit is invalid.")
  )
  parameter_value <- function(parameter_name) {
    values <- raw_alignment_parameters$value[
      raw_alignment_parameters$parameter == parameter_name
    ]
    if (length(values) != 1L) return(NA_real_)
    suppressWarnings(as.numeric(values))
  }
  deployed_picker <- peak_picker_sensitivity[
    peak_picker_sensitivity$is_deployed_setting, , drop = FALSE
  ]
  assert_true(
    isTRUE(all.equal(parameter_value("minimum_profile_snr"), 4)) &&
      isTRUE(all.equal(parameter_value("minimum_prominence_snr"), 2.5)) &&
      isTRUE(all.equal(deployed_picker$minimum_profile_snr, 4)) &&
      isTRUE(all.equal(deployed_picker$minimum_prominence_snr, 2.5)),
    paste(source_dataset,
          "does not preserve the deployed 4 / 2.5 peak-detection gates.")
  )
  assigned_detections <- detections[
    !is.na(detections$aligned_feature_id) &
      nzchar(detections$aligned_feature_id), , drop = FALSE
  ]
  assert_true(!anyDuplicated(paste(
    assigned_detections$sample_id,
    assigned_detections$aligned_feature_id, sep = "\r"
  )), paste(source_dataset,
            "raw DAD alignment assigns two peaks from one sample to one feature."))
  required_review_fields <- c(
    "rare_feature_review_status", "nearest_recurrent_feature",
    "nearest_recurrent_rt_difference_minutes",
    "nearest_recurrent_spectral_distance", "downstream_use_recommendation",
    "rescued_detection_count"
  )
  assert_true(all(required_review_fields %in% names(features)),
              paste(source_dataset,
                    "raw DAD rare-feature review metadata are incomplete."))
  rescued <- reconciliation$reconciliation_decision ==
    "Rescued into recurrent feature"
  if (any(rescued)) {
    assert_true(all(
      reconciliation$rt_difference_minutes[rescued] <= 0.12 + 1e-10 &
        reconciliation$spectral_distance[rescued] <= 0.07 + 1e-10 &
        reconciliation$joint_score[rescued] <= 1 + 1e-10 &
        reconciliation$wavelength_bootstrap_support[rescued] >= 0.85 - 1e-10
    ), paste(source_dataset,
             "contains a rare-fragment rescue that failed its evidence gates."))
  }
  required_iteration_fields <- c(
    "rescue_iteration", "detection_id", "source_group",
    "current_source_group", "proposed_target_group",
    "target_was_original_anchor", "reconciliation_decision"
  )
  assert_true(
    all(required_iteration_fields %in% names(reconciliation_iterations)),
    paste(source_dataset,
          "iterative rare-fragment audit is missing required fields.")
  )
  manifest_rescues <- as.integer(strsplit(
    as.character(policy_manifest$rescues_by_iteration), ";", fixed = TRUE
  )[[1]])
  iteration_numbers <- if (nrow(reconciliation_iterations)) {
    sort(unique(reconciliation_iterations$rescue_iteration))
  } else seq_len(policy_manifest$rescue_iterations_completed)
  assert_true(
    length(iteration_numbers) >= 1L &&
      identical(iteration_numbers, seq_len(max(iteration_numbers))) &&
      policy_manifest$rescue_iterations_completed == max(iteration_numbers),
    paste(source_dataset,
          "iterative rescue rounds are incomplete or disagree with the manifest.")
  )
  iteration_rescued <- reconciliation_iterations$reconciliation_decision ==
    "Rescued into recurrent feature"
  rescues_per_iteration <- if (nrow(reconciliation_iterations)) {
    vapply(
      iteration_numbers,
      function(iteration) sum(
        iteration_rescued &
          reconciliation_iterations$rescue_iteration == iteration
      ),
      integer(1)
    )
  } else rep(0L, length(iteration_numbers))
  assert_true(
    identical(rescues_per_iteration, manifest_rescues) &&
      tail(rescues_per_iteration, 1) == 0L &&
      (length(rescues_per_iteration) == 1L ||
         all(head(rescues_per_iteration, -1L) > 0L)),
    paste(source_dataset,
          "iterative rescue did not terminate with one complete zero-rescue round.")
  )
  assert_true(
    !anyDuplicated(
      reconciliation_iterations$detection_id[iteration_rescued]
    ),
    paste(source_dataset,
          "a detection was rescued more than once during iteration.")
  )
  if (any(iteration_rescued)) {
    assert_true(all(
      reconciliation_iterations$target_was_original_anchor[
        iteration_rescued
      ] %in% TRUE &
        reconciliation_iterations$rt_difference_minutes[
          iteration_rescued
        ] <= 0.12 + 1e-10 &
        reconciliation_iterations$spectral_distance[
          iteration_rescued
        ] <= 0.07 + 1e-10 &
        reconciliation_iterations$joint_score[iteration_rescued] <=
          1 + 1e-10 &
        reconciliation_iterations$wavelength_bootstrap_support[
          iteration_rescued
        ] >= 0.85 - 1e-10
    ), paste(source_dataset,
             "contains an iterative rescue that failed a safety gate."))
  }

  abundance_matrix <- as.matrix(abundance_wide[, -1, drop = FALSE])
  storage.mode(abundance_matrix) <- "double"
  rownames(abundance_matrix) <- abundance_wide[[1]]
  uncertainty_flag_matrix <- as.matrix(
    uncertainty_flag_wide[, -1, drop = FALSE]
  )
  storage.mode(uncertainty_flag_matrix) <- "integer"
  rownames(uncertainty_flag_matrix) <- uncertainty_flag_wide[[1]]
  upper_bound_matrix <- as.matrix(upper_bound_wide[, -1, drop = FALSE])
  storage.mode(upper_bound_matrix) <- "double"
  rownames(upper_bound_matrix) <- upper_bound_wide[[1]]
  assert_true(
    !anyNA(abundance_matrix) && all(abundance_matrix >= 0) &&
      !anyNA(uncertainty_flag_matrix) &&
      all(uncertainty_flag_matrix %in% c(0L, 1L)) &&
      !anyNA(upper_bound_matrix) &&
      all(upper_bound_matrix + 1e-10 >= abundance_matrix),
    paste(source_dataset,
          "confident abundance, uncertainty flags, or upper bounds are invalid.")
  )

  required_ambiguity_fields <- c(
    "detection_id", "sample_id", "integrated_dad_area",
    "candidate_aligned_features", "candidate_feature_count"
  )
  assert_true(
    all(required_ambiguity_fields %in% names(ambiguous_detections)),
    paste(source_dataset, "ambiguous-detection audit lacks quantitative fields.")
  )
  expected_uncertainty_flags <- matrix(
    0L, nrow = nrow(abundance_matrix), ncol = ncol(abundance_matrix),
    dimnames = dimnames(abundance_matrix)
  )
  expected_upper_bounds <- abundance_matrix
  if (nrow(ambiguous_detections)) {
    for (row in seq_len(nrow(ambiguous_detections))) {
      candidate_features <- strsplit(
        ambiguous_detections$candidate_aligned_features[row], ";", fixed = TRUE
      )[[1]]
      candidate_features <- candidate_features[
        nzchar(candidate_features) & candidate_features %in%
          colnames(abundance_matrix)
      ]
      sample_position <- match(
        ambiguous_detections$sample_id[row], rownames(abundance_matrix)
      )
      assert_true(
        length(candidate_features) ==
          ambiguous_detections$candidate_feature_count[row] &&
          is.finite(sample_position) &&
          is.finite(ambiguous_detections$integrated_dad_area[row]) &&
          ambiguous_detections$integrated_dad_area[row] >= 0,
        paste(source_dataset,
              "contains an invalid unresolved-detection candidate record.")
      )
      expected_uncertainty_flags[sample_position, candidate_features] <- 1L
      expected_upper_bounds[sample_position, candidate_features] <-
        expected_upper_bounds[sample_position, candidate_features, drop = TRUE] +
        ambiguous_detections$integrated_dad_area[row]
    }
  }
  assert_true(
    same_numeric_matrix(uncertainty_flag_matrix, expected_uncertainty_flags) &&
      same_numeric_matrix(upper_bound_matrix, expected_upper_bounds),
    paste(source_dataset,
          "uncertainty flags or upper bounds disagree with unresolved peaks.")
  )
  expected_values <- as.numeric(abundance_matrix[cbind(
    match(abundance_long$sample_id, rownames(abundance_matrix)),
    match(abundance_long$aligned_feature_id, colnames(abundance_matrix))
  )])
  expected_flags <- as.integer(uncertainty_flag_matrix[cbind(
    match(abundance_long$sample_id, rownames(uncertainty_flag_matrix)),
    match(abundance_long$aligned_feature_id,
          colnames(uncertainty_flag_matrix))
  )])
  expected_upper_values <- as.numeric(upper_bound_matrix[cbind(
    match(abundance_long$sample_id, rownames(upper_bound_matrix)),
    match(abundance_long$aligned_feature_id, colnames(upper_bound_matrix))
  )])
  assert_true(
    isTRUE(all.equal(
      expected_values, abundance_long$integrated_dad_area,
      tolerance = 1e-10, check.attributes = FALSE
    )) &&
      isTRUE(all.equal(
        expected_values, abundance_long$confident_integrated_dad_area,
        tolerance = 1e-10, check.attributes = FALSE
      )) &&
      identical(expected_flags,
                as.integer(abundance_long$unresolved_candidate_flag)) &&
      isTRUE(all.equal(
        expected_upper_values,
        abundance_long$featurewise_upper_bound_integrated_dad_area,
        tolerance = 1e-10, check.attributes = FALSE
      )),
    paste(source_dataset, "raw DAD wide and long uncertainty tables differ.")
  )
  expected_status <- ifelse(
    expected_flags > 0 & expected_values == 0,
    "Detected locally but alignment unresolved",
    ifelse(expected_values > 0, "Aligned peak detected",
           "Not detected above thresholds")
  )
  assert_true(identical(expected_status,
                        as.character(abundance_long$detection_status)),
              paste(source_dataset,
                    "raw DAD uncertainty status labels are inconsistent."))
  assert_true(all(features$detected_samples == colSums(
    abundance_matrix > 0
  )) && all(features$unresolved_candidate_samples == colSums(
    uncertainty_flag_matrix > 0
  )), paste(source_dataset,
            "feature prevalence or unresolved counts disagree with matrices."))

  # The simplified map must contain exactly one row for every assigned peak,
  # while the atlas manifest and color dictionary must cover every aligned
  # feature once and preserve chromatographic order. These checks prevent a
  # pretty diagnostic from silently omitting peaks or reusing a misleading key.
  expected_assigned_key <- paste(
    assigned_detections$sample_id, assigned_detections$detection_id,
    assigned_detections$aligned_feature_id, sep = "\r"
  )
  observed_map_key <- paste(
    aligned_peak_map$sample_id, aligned_peak_map$detection_id,
    aligned_peak_map$aligned_feature_id, sep = "\r"
  )
  assert_true(
    setequal(expected_assigned_key, observed_map_key) &&
      !anyDuplicated(observed_map_key) &&
      all(abs(
        aligned_peak_map$observed_minus_consensus_minutes -
          (aligned_peak_map$aligned_retention_time -
             aligned_peak_map$feature_consensus_retention_time)
      ) <= 1e-10),
    paste(source_dataset,
          "simplified aligned-peak map data are incomplete or inconsistent.")
  )
  assert_true(
    identical(as.character(feature_color_dictionary$aligned_feature_id),
              feature_ids) &&
      identical(as.character(identity_atlas_manifest$aligned_feature_id),
                feature_ids) &&
      !anyDuplicated(feature_color_dictionary$figure_color_hex) &&
      all(grepl("^#[0-9A-Fa-f]{6}$",
                feature_color_dictionary$figure_color_hex)) &&
      all(identity_atlas_manifest$figure_color_hex ==
            feature_color_dictionary$figure_color_hex) &&
      all(identity_atlas_manifest$feature_order == seq_along(feature_ids)) &&
      all(identity_atlas_manifest$atlas_total_pages ==
            ceiling(length(feature_ids) / 16)),
    paste(source_dataset,
          "same-feature identity atlas key is incomplete or inconsistent.")
  )

  required_figures <- file.path(alignment_folder, c(
    "genus_prevalence_alignment_pseudochromatogram.pdf",
    "genus_class_pseudochromatogram.pdf",
    "genus_chromophore_pseudochromatogram.pdf",
    "genus_hybrid_pseudochromatogram.pdf",
    "genus_class_signal_composition_mirror.pdf",
    "genus_chromophore_signal_composition_mirror.pdf",
    "genus_hybrid_signal_composition_mirror.pdf",
    "retention_time_correction_audit.pdf",
    "feature_prevalence_and_alignment_support.pdf",
    "peak_picking_parameter_sensitivity.pdf",
    "sample_peak_picking_audit.pdf",
    "representative_peak_picking_profiles.pdf",
    "sample_by_aligned_feature_detection_map.pdf",
    "sample_by_aligned_peak_map.pdf",
    "sample_by_aligned_feature_identity_atlas.pdf",
    "rare_fragment_reconciliation_diagnostics.pdf",
    "rare_fragment_rescue_convergence.pdf",
    "neighboring_feature_reconciliation_diagnostics.pdf",
    "all_aligned_feature_representative_spectra_with_uncertainty.pdf"
  ))
  required_figures <- resolve_output_paths(required_figures)
  assert_true(all(file.exists(required_figures)) &&
                all(file.info(required_figures)$size > 1000),
              paste(source_dataset,
                    "raw DAD route is missing a required diagnostic figure."))
  raw_scaling_audit <- read_output_csv(file.path(
    alignment_folder, "genus_pseudochromatogram_scaling_audit.csv"
  ))
  assert_true(
    nrow(raw_scaling_audit) == 3L &&
      all(raw_scaling_audit$maximum_reconstruction_error <= 1e-6) &&
      all(raw_scaling_audit$black_sum_not_below_any_component),
    paste(source_dataset, "raw-DAD pseudochromatogram scaling is inconsistent.")
  )
  final_folder <- Sys.getenv(
    "UV_RAW_FINAL_OUTPUT_DIR",
    unset = file.path(output_root, "03_Final_Results")
  )
  final_index <- read_output_csv(file.path(
    final_folder, "READ_FIRST_output_order.csv"
  ))
  final_strict_matrix <- read_output_csv(file.path(
    final_folder, "sample_by_peak_absorbance_strict.csv"
  ))
  final_complete_matrix <- read_output_csv(file.path(
    final_folder, "sample_by_peak_absorbance_complete.csv"
  ))
  final_strict_uncertainty_flags <- read_output_csv(file.path(
    final_folder, "sample_by_peak_assignment_uncertainty_flag_strict.csv"
  ))
  final_complete_uncertainty_flags <- read_output_csv(file.path(
    final_folder, "sample_by_peak_assignment_uncertainty_flag_complete.csv"
  ))
  final_strict_upper_bounds <- read_output_csv(file.path(
    final_folder, "sample_by_peak_absorbance_upper_bound_strict.csv"
  ))
  final_complete_upper_bounds <- read_output_csv(file.path(
    final_folder, "sample_by_peak_absorbance_upper_bound_complete.csv"
  ))
  final_strict_annotations <- read_output_csv(file.path(
    final_folder, "annotated_peak_table_with_sample_absorbance_strict.csv"
  ))
  final_complete_annotations <- read_output_csv(file.path(
    final_folder, "annotated_peak_table_with_sample_absorbance_complete.csv"
  ))
  expected_primary_annotation_fields <- c(
    "aligned_feature_id", "consensus_retention_time", "detected_samples",
    "total_samples", "detection_fraction", "unresolved_candidate_samples",
    "unresolved_candidate_fraction", "recurrence_status",
    "alignment_confidence", "practical_spectral_decision",
    "nearest_library_spectral_analogue", "nearest_analogue_library_class",
    "class_only_label", "chromophore_only_label",
    "hybrid_class_then_chromophore_label", "accepted_analogue_set",
    "ensemble_compound_probability", "ensemble_class_probability",
    "library_domain_p_value", "library_domain_accepted",
    "full_nearest_compound_bootstrap_support",
    "wavelength_block_compound_support", "possible_coelution",
    "recommended_followup"
  )
  assert_true(
    identical(names(final_strict_matrix), names(strict_primary_wide)) &&
      same_numeric_matrix(
        final_strict_matrix[, -1, drop = FALSE],
        strict_primary_wide[, -1, drop = FALSE]
      ) &&
      identical(names(final_complete_matrix), names(abundance_wide)) &&
      same_numeric_matrix(
        final_complete_matrix[, -1, drop = FALSE],
        abundance_wide[, -1, drop = FALSE]
      ) &&
      !anyNA(final_strict_matrix[, -1, drop = FALSE]) &&
      !anyNA(final_complete_matrix[, -1, drop = FALSE]) &&
      same_numeric_matrix(
        final_strict_uncertainty_flags[, -1, drop = FALSE],
        strict_uncertainty_flag_wide[, -1, drop = FALSE]
      ) &&
      same_numeric_matrix(
        final_complete_uncertainty_flags[, -1, drop = FALSE],
        uncertainty_flag_wide[, -1, drop = FALSE]
      ) &&
      same_numeric_matrix(
        final_strict_upper_bounds[, -1, drop = FALSE],
        strict_upper_bound_wide[, -1, drop = FALSE]
      ) &&
      same_numeric_matrix(
        final_complete_upper_bounds[, -1, drop = FALSE],
        upper_bound_wide[, -1, drop = FALSE]
      ) &&
      identical(as.character(final_complete_annotations$aligned_feature_id),
                feature_ids) &&
      identical(as.character(final_strict_annotations$aligned_feature_id),
                expected_strict_features) &&
      identical(
        names(final_strict_annotations),
        c(expected_primary_annotation_fields,
          as.character(strict_primary_wide[[1]]))
      ) &&
      !anyNA(final_strict_annotations) &&
      identical(final_index$interpretation_order, seq_len(nrow(final_index))),
    paste(source_dataset, "final collaborator tables are inconsistent.")
  )
  final_figures <- file.path(final_folder, c(
    "sample_by_peak_detection_matrix_strict.pdf",
    "sample_by_peak_detection_matrix_review_features.pdf",
    "genus_class_pseudochromatogram.pdf",
    "genus_chromophore_pseudochromatogram.pdf",
    "genus_hybrid_pseudochromatogram.pdf"
  ))
  final_figures <- resolve_output_paths(final_figures)
  assert_true(
    all(file.exists(final_figures)) && all(file.info(final_figures)$size > 1000),
    paste(source_dataset, "final collaborator figure package is incomplete.")
  )
  validate_wavelength_alignment(alignment_folder, source_dataset)
  alignment_counts <- validate_readable_outputs(alignment_folder)

  summary_row$dataset <- paste0(source_dataset, " [", alignment_policy, "]")
  summary_row$analysis_mode <- paste0(
    "raw_agilent_dad_joint_alignment:", alignment_policy
  )
  summary_row$input_peaks <- nrow(detections)
  summary_row$csv_outputs <- summary_row$csv_outputs +
    unname(alignment_counts[["csv_files"]])
  summary_row$pdf_outputs <- summary_row$pdf_outputs +
    unname(alignment_counts[["pdf_files"]])
  summary_row$landscape_samples <- length(sample_ids)
  summary_row$homologous_features <- length(feature_ids)
  summary_row$abundance_cells <- length(abundance_matrix)
  summary_row
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments)) {
  raw_pattern <- "_RawDAD_(RetainedUniqueness|ExpandedRescue)$"
  raw_arguments <- arguments[grepl(raw_pattern, arguments)]
  if (any(grepl("_RawDAD$", arguments))) {
    stop(
      "Raw DAD validation now requires the policy-qualified dataset name, ",
      "for example Genus_RawDAD_RetainedUniqueness or ",
      "Genus_RawDAD_ExpandedRescue."
    )
  }
  raw_sources <- sub(raw_pattern, "", raw_arguments)
  raw_policy_labels <- sub("^.*_RawDAD_", "", raw_arguments)
  raw_policies <- c(
    RetainedUniqueness = "retained_uniqueness",
    ExpandedRescue = "expanded_rescue"
  )[raw_policy_labels]
  non_raw_arguments <- arguments[!grepl(raw_pattern, arguments)]
  if (length(non_raw_arguments)) {
    stop(
      "Only policy-qualified raw-DAD datasets are supported, for example ",
      "Genus_RawDAD_RetainedUniqueness."
    )
  }
} else {
  completed_raw_roots <- list.dirs(
    ".", recursive = FALSE, full.names = FALSE
  )
  completed_raw_roots <- completed_raw_roots[grepl(
    "_RawDAD_(RetainedUniqueness|ExpandedRescue)$", completed_raw_roots
  )]
  raw_policy_labels <- sub(
    "^.*_RawDAD_([^_]+)$", "\\1", completed_raw_roots
  )
  raw_policies <- c(
    RetainedUniqueness = "retained_uniqueness",
    ExpandedRescue = "expanded_rescue"
  )[raw_policy_labels]
  raw_sources <- sub(
    "_RawDAD_(RetainedUniqueness|ExpandedRescue)$", "", completed_raw_roots
  )
  alignment_folders <- file.path(completed_raw_roots, "01_Alignment_QC")
  annotation_folders <- file.path(completed_raw_roots, "02_Library_Annotation")
  completed <- dir.exists(file.path("data", paste0(raw_sources, "_Agilent"))) &
    dir.exists(alignment_folders) & dir.exists(annotation_folders) &
    !is.na(raw_policies)
  raw_sources <- raw_sources[completed]
  raw_policies <- unname(raw_policies[completed])
}
assert_true(length(raw_sources) > 0L,
            "No completed policy-qualified raw-DAD datasets were available.")

raw_dad_summary <- Map(validate_raw_dad_dataset, raw_sources, raw_policies)
validation_summary <- do.call(rbind, raw_dad_summary)

cat("Collaborator-handoff validation passed.\n")
print(validation_summary, row.names = FALSE)
