## ========================================================================== ##
## ANNOTATED GENUS-LEVEL PSEUDOCHROMATOGRAMS FOR RAW DAD ALIGNMENT
## ========================================================================== ##
## Purpose
##   Join the fixed-library annotation results to the features produced by
##   01_raw_agilent_dad_alignment.R and draw genus-level pseudochromatograms
##   using actual across-sample DAD-area summaries. The alignment was created
##   without the library; this stage only adds interpretation and color.
##
## Entry point
##   Called automatically by RUN_AGILENT_DAD.R <Genus>
##   --alignment-policy=<retained_uniqueness|expanded_rescue>.
## ========================================================================== ##

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for genus pseudochromatograms.")
}
library(ggplot2)
source("00_library_inputs.R")

arguments <- commandArgs(trailingOnly = TRUE)
SOURCE_DATASET_NAME <- if (length(arguments)) arguments[1] else "Rhododendron"
requested_policy <- Sys.getenv("UV_RAW_ALIGNMENT_POLICY", unset = "")
if (!requested_policy %in% c("retained_uniqueness", "expanded_rescue")) {
  stop(
    "Raw DAD pseudochromatograms require UV_RAW_ALIGNMENT_POLICY. Run this ",
    "stage through RUN_AGILENT_DAD.R with an explicit --alignment-policy."
  )
}
policy_folder_label <- c(
  retained_uniqueness = "RetainedUniqueness",
  expanded_rescue = "ExpandedRescue"
)[requested_policy]
ANNOTATION_DATASET_NAME <- Sys.getenv(
  "UV_DATASET_NAME",
  unset = paste0(SOURCE_DATASET_NAME, "_RawDAD_", policy_folder_label)
)
ANNOTATION_SLUG <- tolower(gsub("[^A-Za-z0-9]+", "_", ANNOTATION_DATASET_NAME))
DISPLAY_NAME <- gsub(
  "([a-z])([A-Z])", "\\1 \\2", gsub("[_-]+", " ", SOURCE_DATASET_NAME)
)
ALIGNMENT_DIR <- Sys.getenv(
  "UV_RAW_ALIGNMENT_OUTPUT_DIR",
  unset = file.path(
    paste0(SOURCE_DATASET_NAME, "_RawDAD_", policy_folder_label),
    "01_Alignment_QC"
  )
)
ANNOTATION_DIR <- Sys.getenv(
  "UV_OUTPUT_DIR",
  unset = file.path(
    paste0(SOURCE_DATASET_NAME, "_RawDAD_", policy_folder_label),
    "02_Library_Annotation"
  )
)
FINAL_DIR <- Sys.getenv(
  "UV_RAW_FINAL_OUTPUT_DIR",
  unset = file.path(
    paste0(SOURCE_DATASET_NAME, "_RawDAD_", policy_folder_label),
    "03_Final_Results"
  )
)
dir.create(FINAL_DIR, recursive = TRUE, showWarnings = FALSE)

# The alignment manifest is the authoritative handoff contract. Refuse to join
# annotations or draw chromatograms if the requested policy and alignment
# folder disagree.
policy_manifest_file <- file.path(
  ALIGNMENT_DIR, "raw_alignment_policy_manifest.csv"
)
if (!file.exists(policy_manifest_file)) {
  stop("Missing raw alignment policy manifest: ", policy_manifest_file)
}
policy_manifest <- read.csv(
  policy_manifest_file, stringsAsFactors = FALSE, check.names = FALSE
)
if (nrow(policy_manifest) != 1L ||
    !policy_manifest$alignment_policy %in%
      c("retained_uniqueness", "expanded_rescue")) {
  stop("Raw alignment policy manifest is invalid: ", policy_manifest_file)
}
if (!identical(requested_policy, policy_manifest$alignment_policy)) {
  stop(
    "Requested alignment policy ('", requested_policy,
    "') does not match the alignment manifest ('",
    policy_manifest$alignment_policy, "')."
  )
}
ALIGNMENT_POLICY <- policy_manifest$alignment_policy
ALIGNMENT_POLICY_DISPLAY_NAME <-
  policy_manifest$alignment_policy_display_name
FIGURE_TITLE_PREFIX <- paste(
  DISPLAY_NAME, "-", ALIGNMENT_POLICY_DISPLAY_NAME
)

feature_file <- file.path(ALIGNMENT_DIR, "aligned_feature_summary.csv")
annotation_file <- file.path(
  ANNOTATION_DIR, paste0(ANNOTATION_SLUG, "_comprehensive_results.csv")
)
if (!file.exists(feature_file)) stop("Missing aligned feature summary: ", feature_file)
if (!file.exists(annotation_file)) stop("Missing annotation results: ", annotation_file)

features <- read.csv(feature_file, stringsAsFactors = FALSE,
                     check.names = FALSE)
annotations <- read.csv(annotation_file, stringsAsFactors = FALSE,
                        check.names = FALSE)
if (!all(features$alignment_policy == ALIGNMENT_POLICY)) {
  stop("Aligned feature table does not match its policy manifest.")
}
if (!"source_alignment_policy" %in% names(annotations) ||
    !all(annotations$source_alignment_policy == ALIGNMENT_POLICY)) {
  stop("Library annotation results do not match the selected alignment policy.")
}
annotation_position <- match(features$aligned_feature_id, annotations$peak_id)
if (anyNA(annotation_position) || anyDuplicated(annotations$peak_id)) {
  stop("Raw alignment features do not map one-to-one to annotation results.")
}
annotations <- annotations[annotation_position, , drop = FALSE]

# This is the principal feature dictionary for downstream ecological or
# comparative analyses. It retains every alignment field and every library
# interpretation in one row per aligned peak. Duplicate annotation column names
# receive an explicit prefix instead of being silently overwritten.
annotation_columns <- setdiff(names(annotations), "peak_id")
annotation_export <- annotations[annotation_columns]
overlapping_names <- names(annotation_export) %in% names(features)
names(annotation_export)[overlapping_names] <- paste0(
  "annotation_", names(annotation_export)[overlapping_names]
)
annotated_features <- cbind(features, annotation_export)
annotated_features <- annotated_features[
  order(annotated_features$consensus_retention_time,
        annotated_features$aligned_feature_id), , drop = FALSE
]
write.csv(
  annotated_features,
  file.path(ALIGNMENT_DIR,
            "aligned_feature_summary_with_library_annotations.csv"),
  row.names = FALSE
)

# ------------------------ FINAL COLLABORATOR TABLES -----------------------
# The main handoff contains the recommended strict feature set, ordered by
# consensus RT. A second complete table retains every discovery feature for
# audit. Both combine feature metadata and one absorbance column per sample so
# downstream users do not need to perform their first join correctly by hand.
complete_matrix <- read.csv(
  file.path(
    ALIGNMENT_DIR, "sample_by_aligned_feature_integrated_dad_area.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE
)
strict_matrix <- read.csv(
  file.path(
    ALIGNMENT_DIR,
    "sample_by_aligned_feature_integrated_dad_area_strict_primary.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE
)
complete_uncertainty_flags <- read.csv(
  file.path(
    ALIGNMENT_DIR,
    "sample_by_aligned_feature_assignment_uncertainty_flag.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE
)
strict_uncertainty_flags <- read.csv(
  file.path(
    ALIGNMENT_DIR,
    paste0(
      "sample_by_aligned_feature_assignment_uncertainty_flag_",
      "strict_primary.csv"
    )
  ), stringsAsFactors = FALSE, check.names = FALSE
)
complete_upper_bounds <- read.csv(
  file.path(
    ALIGNMENT_DIR,
    "sample_by_aligned_feature_integrated_dad_area_upper_bound.csv"
  ), stringsAsFactors = FALSE, check.names = FALSE
)
strict_upper_bounds <- read.csv(
  file.path(
    ALIGNMENT_DIR,
    paste0(
      "sample_by_aligned_feature_integrated_dad_area_upper_bound_",
      "strict_primary.csv"
    )
  ), stringsAsFactors = FALSE, check.names = FALSE
)
feature_order <- annotated_features$aligned_feature_id
strict_feature_order <- feature_order[
  annotated_features$include_in_strict_primary_matrix
]
if (!setequal(feature_order, names(complete_matrix)[-1]) ||
    !setequal(strict_feature_order, names(strict_matrix)[-1]) ||
    !setequal(feature_order, names(complete_uncertainty_flags)[-1]) ||
    !setequal(strict_feature_order, names(strict_uncertainty_flags)[-1]) ||
    !setequal(feature_order, names(complete_upper_bounds)[-1]) ||
    !setequal(strict_feature_order, names(strict_upper_bounds)[-1])) {
  stop("Final abundance matrices do not match the annotated feature dictionary.")
}
complete_matrix <- complete_matrix[, c("sample_id", feature_order), drop = FALSE]
strict_matrix <- strict_matrix[, c(
  "sample_id", strict_feature_order
), drop = FALSE]
complete_uncertainty_flags <- complete_uncertainty_flags[, c(
  "sample_id", feature_order
), drop = FALSE]
strict_uncertainty_flags <- strict_uncertainty_flags[, c(
  "sample_id", strict_feature_order
), drop = FALSE]
complete_upper_bounds <- complete_upper_bounds[, c(
  "sample_id", feature_order
), drop = FALSE]
strict_upper_bounds <- strict_upper_bounds[, c(
  "sample_id", strict_feature_order
), drop = FALSE]
write.csv(
  strict_matrix,
  file.path(FINAL_DIR, "sample_by_peak_absorbance_strict.csv"),
  row.names = FALSE
)
write.csv(
  complete_matrix,
  file.path(FINAL_DIR, "sample_by_peak_absorbance_complete.csv"),
  row.names = FALSE
)
write.csv(
  strict_uncertainty_flags,
  file.path(FINAL_DIR, "sample_by_peak_assignment_uncertainty_flag_strict.csv"),
  row.names = FALSE
)
write.csv(
  complete_uncertainty_flags,
  file.path(FINAL_DIR, "sample_by_peak_assignment_uncertainty_flag_complete.csv"),
  row.names = FALSE
)
write.csv(
  strict_upper_bounds,
  file.path(FINAL_DIR, "sample_by_peak_absorbance_upper_bound_strict.csv"),
  row.names = FALSE
)
write.csv(
  complete_upper_bounds,
  file.path(FINAL_DIR, "sample_by_peak_absorbance_upper_bound_complete.csv"),
  row.names = FALSE
)

metadata_first <- c(
  "alignment_policy", "alignment_policy_display_name", "aligned_feature_id",
  "consensus_retention_time", "detected_samples", "total_samples",
  "detection_fraction", "recurrence_status", "alignment_confidence",
  "rare_feature_review_status", "downstream_use_recommendation",
  "include_in_strict_primary_matrix", "practical_spectral_decision",
  "nearest_library_spectral_analogue", "class_only_label",
  "chromophore_only_label", "hybrid_class_then_chromophore_label"
)
metadata_first <- metadata_first[metadata_first %in% names(annotated_features)]
remaining_metadata <- setdiff(names(annotated_features), metadata_first)
annotated_features <- annotated_features[, c(
  metadata_first, remaining_metadata
), drop = FALSE]
sample_absorbance_by_feature <- as.data.frame(
  t(as.matrix(complete_matrix[, -1, drop = FALSE])),
  check.names = FALSE, stringsAsFactors = FALSE
)
names(sample_absorbance_by_feature) <- complete_matrix$sample_id
sample_absorbance_by_feature$aligned_feature_id <- rownames(
  sample_absorbance_by_feature
)
sample_absorbance_by_feature <- sample_absorbance_by_feature[, c(
  "aligned_feature_id", complete_matrix$sample_id
), drop = FALSE]
annotated_sample_table <- merge(
  annotated_features, sample_absorbance_by_feature,
  by = "aligned_feature_id", all.x = TRUE, sort = FALSE
)
annotated_sample_table <- annotated_sample_table[
  match(feature_order, annotated_sample_table$aligned_feature_id), , drop = FALSE
]
# Keep the primary annotated handoff concise and NA-free. Optional diagnostics
# such as nearest recurrent-feature distances and tie-specific wavelength
# statistics are meaningful only for some features; placing all of them in the
# primary table created a visually confusing field of NA values. They remain
# available in the complete audit table and the numbered QC folders.
primary_annotation_fields <- c(
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
primary_annotation_fields <- primary_annotation_fields[
  primary_annotation_fields %in% names(annotated_sample_table)
]
strict_annotated_sample_table <- annotated_sample_table[
  annotated_sample_table$include_in_strict_primary_matrix,
  c(primary_annotation_fields, complete_matrix$sample_id), drop = FALSE
]
if (anyNA(strict_annotated_sample_table)) {
  stop(
    "The concise primary annotated table contains an unexpected missing value; ",
    "inspect the annotation and abundance inputs before collaborator handoff."
  )
}
write.csv(
  strict_annotated_sample_table,
  file.path(FINAL_DIR, "annotated_peak_table_with_sample_absorbance_strict.csv"),
  row.names = FALSE
)
write.csv(
  annotated_sample_table,
  file.path(FINAL_DIR, "annotated_peak_table_with_sample_absorbance_complete.csv"),
  row.names = FALSE
)

# Broad chromophore families keep the legend readable. Fine-group assignments
# remain available in the comprehensive results and joined output table.
chromophore_label <- annotations$chromophore_only_label
supported_chromophore <- !chromophore_label %in%
  c("Ambiguous", "Outside library coverage") &
  nzchar(annotations$analogue_set_modal_chromophore_family)
chromophore_label[supported_chromophore] <-
  annotations$analogue_set_modal_chromophore_family[supported_chromophore]

hybrid_label <- annotations$hybrid_class_then_chromophore_label
hybrid_class_rows <- annotations$hybrid_label_type %in%
  c("Class from nearest analogue", "Class tendency")
hybrid_chromophore_rows <- annotations$hybrid_label_type %in% c(
  "Bootstrap chromophore tendency", "Fine chromophore fallback",
  "Broad chromophore fallback"
)
hybrid_label[hybrid_class_rows] <- paste0(
  "Class: ", hybrid_label[hybrid_class_rows]
)
hybrid_label[hybrid_chromophore_rows] <- paste0(
  "Chromophore: ",
  annotations$analogue_set_modal_chromophore_family[hybrid_chromophore_rows]
)

pseudo_data <- data.frame(
  alignment_policy = ALIGNMENT_POLICY,
  alignment_policy_display_name = ALIGNMENT_POLICY_DISPLAY_NAME,
  aligned_feature_id = features$aligned_feature_id,
  consensus_retention_time = features$consensus_retention_time,
  median_fwhm_minutes = features$median_fwhm_minutes,
  detected_samples = features$detected_samples,
  total_samples = features$total_samples,
  detection_fraction = features$detection_fraction,
  mean_integrated_dad_area_all_samples =
    features$mean_integrated_dad_area_all_samples,
  alignment_confidence = features$alignment_confidence,
  practical_spectral_decision = annotations$practical_spectral_decision,
  nearest_library_spectral_analogue =
    annotations$nearest_library_spectral_analogue,
  class_interpretation = annotations$class_only_label,
  chromophore_interpretation = chromophore_label,
  hybrid_interpretation = hybrid_label,
  stringsAsFactors = FALSE
)
write.csv(
  pseudo_data,
  file.path(ALIGNMENT_DIR, "genus_annotated_pseudochromatogram_data.csv"),
  row.names = FALSE
)

neutral_colors <- c(
  "Ambiguous" = "#A6A6A6",
  "Outside library coverage" = "#4D4D4D"
)
class_levels <- sort(unique(pseudo_data$class_interpretation))
class_non_neutral <- setdiff(class_levels, names(neutral_colors))
class_colors <- c(library_category_colors(class_non_neutral), neutral_colors)

chromophore_levels <- sort(unique(pseudo_data$chromophore_interpretation))
chromophore_non_neutral <- setdiff(chromophore_levels, names(neutral_colors))
chromophore_colors <- c(
  library_category_colors(
    chromophore_non_neutral,
    base_colors = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#FFD92F", "#A65628", "#F781BF"
    )
  ),
  neutral_colors
)

hybrid_levels <- sort(unique(pseudo_data$hybrid_interpretation))
hybrid_class_levels <- hybrid_levels[grepl("^Class: ", hybrid_levels)]
hybrid_chromophore_levels <- hybrid_levels[grepl("^Chromophore: ", hybrid_levels)]
hybrid_neutral_levels <- setdiff(
  hybrid_levels, c(hybrid_class_levels, hybrid_chromophore_levels)
)
hybrid_colors <- c(
  setNames(
    library_category_colors(sub("^Class: ", "", hybrid_class_levels)),
    hybrid_class_levels
  ),
  setNames(
    library_category_colors(
      sub("^Chromophore: ", "", hybrid_chromophore_levels),
      base_colors = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
        "#FF7F00", "#FFD92F", "#A65628", "#F781BF"
      )
    ),
    hybrid_chromophore_levels
  ),
  neutral_colors[intersect(names(neutral_colors), hybrid_neutral_levels)],
  setNames(rep("#A6A6A6", length(setdiff(
    hybrid_neutral_levels, names(neutral_colors)
  ))), setdiff(hybrid_neutral_levels, names(neutral_colors)))
)

alignment_colors <- c(
  "Supported alignment" = "#2C7FB8",
  "Neighboring-feature assignment unstable" = "#F28E2B",
  "RT-dispersed alignment" = "#756BB1",
  "Rare; limited recurrence support" = "#6A51A3",
  "Rare; possible alignment fragment" = "#C51B7D"
)

plot_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom", legend.box = "vertical",
    legend.text = element_text(size = 8.5),
    legend.key.height = grid::unit(0.42, "cm"),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.margin = margin(12, 14, 12, 12)
  )

make_curves <- function(height_values) {
  rows <- lapply(seq_len(nrow(pseudo_data)), function(i) {
    center <- pseudo_data$consensus_retention_time[i]
    fwhm <- max(pseudo_data$median_fwhm_minutes[i], 0.012, na.rm = TRUE)
    sigma <- max(fwhm / 2.355, 0.006)
    time <- seq(center - 4 * sigma, center + 4 * sigma, length.out = 101)
    data.frame(
      aligned_feature_id = pseudo_data$aligned_feature_id[i],
      retention_time = time,
      intensity = height_values[i] * exp(-0.5 * ((time - center) / sigma)^2),
      class_interpretation = pseudo_data$class_interpretation[i],
      chromophore_interpretation = pseudo_data$chromophore_interpretation[i],
      hybrid_interpretation = pseudo_data$hybrid_interpretation[i],
      alignment_confidence = pseudo_data$alignment_confidence[i],
      stringsAsFactors = FALSE
    )
  })
  curves <- do.call(rbind, rows)
  curves
}

# Square-root scaling is visualization-only and prevents one abundant feature
# from hiding the rest of a chemically diverse genus. Exact areas are retained
# in the CSV and abundance matrix.
abundance_curves <- make_curves(sqrt(pmax(
  pseudo_data$mean_integrated_dad_area_all_samples, 0
)))

# The prevalence figure is deliberately not drawn from Gaussian pseudopeaks.
# A mirror plot places every aligned feature at its exact consensus retention
# time: observed prevalence is above zero and relative mean DAD area is below.
# The lower half is square-root scaled only for display, as in the annotated
# abundance pseudochromatograms; exact values remain in the CSV and matrix.
make_prevalence_mirror_plot <- function(filename) {
  mirror <- pseudo_data
  mirror$prevalence_percent <- 100 * pmax(
    pmin(mirror$detection_fraction, 1), 0
  )
  abundance_display <- sqrt(pmax(
    mirror$mean_integrated_dad_area_all_samples, 0
  ))
  abundance_maximum <- max(abundance_display, na.rm = TRUE)
  mirror$relative_mean_area_percent <- if (
    is.finite(abundance_maximum) && abundance_maximum > 0
  ) {
    100 * abundance_display / abundance_maximum
  } else rep(0, nrow(mirror))
  alignment_display_labels <- c(
    "Supported alignment" = "Supported alignment",
    "Neighboring-feature assignment unstable" =
      "Neighboring-feature assignment unstable",
    "RT-dispersed alignment" = "Retention-time-dispersed alignment",
    "Rare; limited recurrence support" = "Rare; limited recurrence support",
    "Rare; possible alignment fragment" = "Rare; possible alignment fragment"
  )

  figure <- ggplot(mirror, aes(x = consensus_retention_time)) +
    # Black underlays keep the mirrored sticks crisp while the narrower color
    # overlay preserves the alignment-confidence interpretation.
    geom_segment(
      aes(xend = consensus_retention_time, y = 0,
          yend = prevalence_percent),
      color = "#000000", linewidth = 0.90, lineend = "butt"
    ) +
    geom_segment(
      aes(xend = consensus_retention_time, y = 0,
          yend = prevalence_percent, color = alignment_confidence),
      linewidth = 0.52, lineend = "butt"
    ) +
    geom_segment(
      aes(xend = consensus_retention_time, y = 0,
          yend = -relative_mean_area_percent),
      color = "#000000", linewidth = 0.90, lineend = "butt"
    ) +
    geom_segment(
      aes(xend = consensus_retention_time, y = 0,
          yend = -relative_mean_area_percent,
          color = alignment_confidence),
      linewidth = 0.52, lineend = "butt"
    ) +
    geom_hline(yintercept = 0, color = "#000000", linewidth = 0.55) +
    scale_color_manual(
      values = alignment_colors, labels = alignment_display_labels,
      name = "Alignment interpretation"
    ) +
    scale_y_continuous(
      breaks = seq(-100, 100, by = 25),
      labels = function(values) abs(values),
      limits = c(-105, 105),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    annotate("text", x = Inf, y = 98, label = "Detection prevalence",
             hjust = 1.05, vjust = 1, fontface = "bold", size = 3.4) +
    annotate("text", x = Inf, y = -98, label = "Relative mean DAD area",
             hjust = 1.05, vjust = 0, fontface = "bold", size = 3.4) +
    labs(
      title = paste(FIGURE_TITLE_PREFIX,
                    "genus-level prevalence/abundance mirror plot"),
      subtitle = paste(
        "Each vertical stick is one aligned feature; prevalence is above zero",
        "and mean DAD area is below zero"
      ),
      caption = paste(
        "Mean DAD area is square-root transformed and rescaled for display;",
        "exact values remain in the abundance matrix."
      ),
      x = "Aligned retention time (minutes)",
      y = "Mirror scale (%)"
    ) +
    guides(color = guide_legend(
      nrow = 2L, byrow = TRUE, title.position = "top", title.hjust = 0.5,
      override.aes = list(linewidth = 2.2)
    )) +
    plot_theme +
    theme(plot.caption = element_text(hjust = 0.5, margin = margin(t = 6)))
  ggsave(file.path(ALIGNMENT_DIR, filename), figure,
         width = 14, height = 8)
}

make_pseudochromatogram <- function(curves, group_column, colors, title,
                                    subtitle, legend_title, filename,
                                    mirror_filename,
                                    legend_rows = 3L, height = 9) {
  curves$display_group <- curves[[group_column]]
  # Sum after interpolation to one common time grid. Each feature curve has a
  # different native grid, so direct aggregation would draw an artificial
  # black zig-zag rather than the envelope of a chromatogram.
  feature_centers <- vapply(
    split(seq_len(nrow(curves)), curves$aligned_feature_id),
    function(rows) curves$retention_time[rows[which.max(curves$intensity[rows])]],
    numeric(1)
  )
  trace_time <- sort(unique(c(
    seq(min(curves$retention_time), max(curves$retention_time),
        length.out = 5000L),
    feature_centers
  )))
  trace_intensity <- Reduce(
    `+`,
    lapply(split(curves, curves$aligned_feature_id), function(feature_curve) {
      stats::approx(
        feature_curve$retention_time, feature_curve$intensity,
        xout = trace_time, yleft = 0, yright = 0
      )$y
    })
  )
  # Use one common display scale for both the component curves and their sum.
  # Consequently the black line is exactly the pointwise sum of the colored
  # pseudopeaks and can never fall below an individual colored component.
  common_scale <- 100 / max(trace_intensity, na.rm = TRUE)
  curves$intensity <- curves$intensity * common_scale
  trace <- data.frame(
    retention_time = trace_time,
    intensity = trace_intensity * common_scale
  )
  interpolated_components <- lapply(
    split(curves, curves$aligned_feature_id), function(component_curve) {
      stats::approx(
        component_curve$retention_time, component_curve$intensity,
        xout = trace$retention_time, yleft = 0, yright = 0
      )$y
    }
  )
  reconstructed_trace <- Reduce(`+`, interpolated_components)
  if (max(abs(reconstructed_trace - trace$intensity), na.rm = TRUE) > 1e-6) {
    stop("Raw-DAD pseudochromatogram components do not reconstruct the black sum.")
  }
  figure <- ggplot() +
    geom_area(
      data = curves,
      aes(retention_time, intensity, group = aligned_feature_id,
          fill = display_group),
      position = "identity", alpha = 0.94, color = NA
    ) +
    # Draw the summed chromatographic envelope last so its solid black outline
    # remains visible above the chemically informative feature fills.
    geom_line(data = trace, aes(retention_time, intensity),
              color = "#000000", linewidth = 0.55) +
    geom_hline(yintercept = 0, color = "#666666", linewidth = 0.28) +
    scale_fill_manual(
      values = colors, name = legend_title,
      labels = function(labels) vapply(
        labels, function(label) paste(strwrap(label, width = 34), collapse = "\n"),
        character(1), USE.NAMES = FALSE
      )
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, subtitle = subtitle,
         x = "Aligned retention time (minutes)",
         y = "Relative genus-level signal (%)") +
    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    plot_theme
  ggsave(file.path(ALIGNMENT_DIR, filename), figure,
         width = 14, height = height)

  # Companion mirror view: the upper black curve is the same pointwise summed
  # envelope. Each lower colored stick is the scaled apex contribution of one
  # aligned feature. Sticks describe feature composition; they are deliberately
  # not smoothed into a second chromatogram.
  stick_rows <- vapply(
    split(seq_len(nrow(curves)), curves$aligned_feature_id),
    function(rows) rows[which.max(curves$intensity[rows])], integer(1)
  )
  sticks <- curves[stick_rows, c(
    "aligned_feature_id", "retention_time", "intensity", "display_group"
  )]
  mirror <- ggplot() +
    geom_area(
      data = trace, aes(retention_time, intensity),
      fill = "#D9D9D9", alpha = 0.42, color = NA
    ) +
    geom_line(
      data = trace, aes(retention_time, intensity),
      color = "#000000", linewidth = 0.72
    ) +
    geom_segment(
      data = sticks,
      aes(x = retention_time, xend = retention_time, y = 0,
          yend = -intensity, color = display_group),
      linewidth = 0.78, lineend = "butt", alpha = 0.92
    ) +
    geom_hline(yintercept = 0, color = "#000000", linewidth = 0.45) +
    scale_color_manual(
      values = colors, name = legend_title,
      labels = function(labels) vapply(
        labels, function(label) paste(strwrap(label, width = 34), collapse = "\n"),
        character(1), USE.NAMES = FALSE
      )
    ) +
    scale_y_continuous(
      labels = function(values) abs(values),
      limits = c(-105, 105), breaks = seq(-100, 100, by = 25),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = paste0(title, " - signal/composition mirror"),
      subtitle = paste(
        "Top: summed black pseudochromatogram; bottom: individual feature",
        "apex contributions colored by interpretation"
      ),
      caption = paste(
        "Both halves use the same scale. The lower sticks are discrete feature",
        "contributions, not a second chromatographic trace."
      ),
      x = "Aligned retention time (minutes)",
      y = "Relative genus-level signal (%)"
    ) +
    guides(color = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    plot_theme +
    theme(plot.caption = element_text(hjust = 0, size = 8))
  ggsave(file.path(ALIGNMENT_DIR, mirror_filename), mirror,
         width = 14, height = height)

  data.frame(
    figure_file = filename,
    mirror_file = mirror_filename,
    shared_scale_factor = common_scale,
    black_trace_maximum = max(trace$intensity, na.rm = TRUE),
    largest_component_apex = max(sticks$intensity, na.rm = TRUE),
    maximum_reconstruction_error = max(
      abs(reconstructed_trace - trace$intensity), na.rm = TRUE
    ),
    black_sum_not_below_any_component = all(vapply(
      interpolated_components,
      function(values) all(values <= trace$intensity + 1e-8), logical(1)
    )),
    stringsAsFactors = FALSE
  )
}

make_prevalence_mirror_plot(
  "genus_prevalence_alignment_pseudochromatogram.pdf"
)
pseudochromatogram_scaling_audit <- list()
pseudochromatogram_scaling_audit[["class"]] <- make_pseudochromatogram(
  abundance_curves, "class_interpretation", class_colors,
  paste(FIGURE_TITLE_PREFIX, "genus-level class pseudochromatogram"),
  "Square-root-scaled mean DAD area across all samples; class is inherited from spectral analogue evidence",
  "Class interpretation", "genus_class_pseudochromatogram.pdf",
  "genus_class_signal_composition_mirror.pdf",
  legend_rows = 2L, height = 8.5
)
pseudochromatogram_scaling_audit[["chromophore"]] <- make_pseudochromatogram(
  abundance_curves, "chromophore_interpretation", chromophore_colors,
  paste(FIGURE_TITLE_PREFIX, "genus-level chromophore pseudochromatogram"),
  "Square-root-scaled mean DAD area; broad chromophore families keep the legend interpretable",
  "Chromophore interpretation", "genus_chromophore_pseudochromatogram.pdf",
  "genus_chromophore_signal_composition_mirror.pdf",
  legend_rows = 3L, height = 9
)
pseudochromatogram_scaling_audit[["hybrid"]] <- make_pseudochromatogram(
  abundance_curves, "hybrid_interpretation", hybrid_colors,
  paste(FIGURE_TITLE_PREFIX,
        "genus-level class-first hybrid pseudochromatogram"),
  "Class evidence is shown first; chromophore evidence is used when class support is unresolved",
  "Interpretation shown", "genus_hybrid_pseudochromatogram.pdf",
  "genus_hybrid_signal_composition_mirror.pdf",
  legend_rows = 4L, height = 10
)
write.csv(
  do.call(rbind, pseudochromatogram_scaling_audit),
  file.path(ALIGNMENT_DIR, "genus_pseudochromatogram_scaling_audit.csv"),
  row.names = FALSE
)

# ---------------------- FINAL SAMPLE-BY-PEAK FIGURES ----------------------
# The main matrix contains only the strict analysis set and is explicitly
# sorted by consensus retention time. Unlike the complete discovery atlas, it
# does not devote columns to provisional fragments or unresolved neighboring
# families. A companion rare-feature atlas preserves those observations for
# review instead of making the main figure look empty.
draw_sample_by_peak_matrix <- function(feature_ids, filename, title,
                                       subtitle, features_per_page = 24L) {
  if (!length(feature_ids)) return(invisible(NULL))
  plot_data <- merge(
    expand.grid(
      sample_id = complete_matrix$sample_id,
      aligned_feature_id = feature_ids,
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    ),
    annotated_features[, c(
      "aligned_feature_id", "consensus_retention_time"
    )],
    by = "aligned_feature_id", all.x = TRUE, sort = FALSE
  )
  complete_long <- do.call(rbind, lapply(feature_ids, function(feature_id) {
    data.frame(
      sample_id = complete_matrix$sample_id,
      aligned_feature_id = feature_id,
      integrated_dad_area = complete_matrix[[feature_id]],
      unresolved_candidate_flag = complete_uncertainty_flags[[feature_id]],
      stringsAsFactors = FALSE
    )
  }))
  abundance_position <- match(
    paste(plot_data$sample_id, plot_data$aligned_feature_id, sep = "\r"),
    paste(complete_long$sample_id, complete_long$aligned_feature_id,
          sep = "\r")
  )
  plot_data$integrated_dad_area <-
    complete_long$integrated_dad_area[abundance_position]
  plot_data$unresolved_candidate_flag <-
    complete_long$unresolved_candidate_flag[abundance_position]
  plot_data$cell_status <- ifelse(
    plot_data$unresolved_candidate_flag > 0 &
      plot_data$integrated_dad_area == 0,
    "Detected but unresolved",
    ifelse(plot_data$integrated_dad_area > 0, "Aligned peak", "Not detected")
  )
  page_count <- ceiling(length(feature_ids) / features_per_page)
  output_path <- file.path(FINAL_DIR, filename)
  grDevices::pdf(
    output_path, width = 14,
    height = max(10, min(18, 3 + 0.13 * length(complete_matrix$sample_id))),
    onefile = TRUE
  )
  for (page in seq_len(page_count)) {
    first <- (page - 1L) * features_per_page + 1L
    last <- min(page * features_per_page, length(feature_ids))
    page_features <- feature_ids[first:last]
    page_data <- plot_data[
      plot_data$aligned_feature_id %in% page_features, , drop = FALSE
    ]
    page_data$aligned_feature_id <- factor(
      page_data$aligned_feature_id, levels = page_features
    )
    page_data$sample_id <- factor(
      page_data$sample_id, levels = rev(complete_matrix$sample_id)
    )
    page_data$cell_status <- factor(
      page_data$cell_status,
      levels = c("Aligned peak", "Detected but unresolved", "Not detected")
    )
    page_metadata <- annotated_features[
      match(page_features, annotated_features$aligned_feature_id), , drop = FALSE
    ]
    axis_labels <- paste0(
      page_metadata$aligned_feature_id, "\n",
      formatC(page_metadata$consensus_retention_time,
              digits = 3, format = "f"), " min"
    )
    names(axis_labels) <- page_features
    figure <- ggplot(page_data, aes(
      aligned_feature_id, sample_id, shape = cell_status,
      color = cell_status
    )) +
      geom_point(size = 1.65, stroke = 0.6) +
      scale_shape_manual(values = c(
        "Aligned peak" = 16, "Detected but unresolved" = 4,
        "Not detected" = 1
      ), name = "Sample-feature result", drop = FALSE) +
      scale_color_manual(values = c(
        "Aligned peak" = "#2C7FB8",
        "Detected but unresolved" = "#E67E22",
        "Not detected" = "#D9D9D9"
      ), name = "Sample-feature result", drop = FALSE) +
      scale_x_discrete(labels = axis_labels, drop = FALSE) +
      labs(
        title = paste0(title, " - page ", page, " of ", page_count),
        subtitle = subtitle,
        x = "Aligned feature and consensus retention time",
        y = "Sample"
      ) +
      theme_minimal(base_size = 9) +
      theme(
        panel.grid.major = element_line(color = "#F0F0F0", linewidth = 0.20),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 48, hjust = 1, vjust = 1,
                                   size = 7),
        axis.text.y = element_text(size = 6.5),
        legend.position = "bottom", legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        plot.title.position = "plot"
      ) +
      guides(
        shape = guide_legend(
          nrow = 1L, order = 1L,
          override.aes = list(
            color = c("#2C7FB8", "#E67E22", "#D9D9D9")
          )
        ),
        color = "none"
      )
    print(figure)
  }
  grDevices::dev.off()
  invisible(output_path)
}

draw_sample_by_peak_matrix(
  strict_feature_order,
  "sample_by_peak_detection_matrix_strict.pdf",
  paste(FIGURE_TITLE_PREFIX, "strict sample-by-peak detection matrix"),
  paste(
    "Only analysis-ready features are shown; columns and pages progress from",
    "low to high consensus retention time"
  )
)
review_feature_ids <- setdiff(feature_order, strict_feature_order)
draw_sample_by_peak_matrix(
  review_feature_ids,
  "sample_by_peak_detection_matrix_review_features.pdf",
  paste(FIGURE_TITLE_PREFIX, "excluded/review feature detection matrix"),
  paste(
    "These provisional or unresolved features remain in the complete table",
    "but are excluded from the strict primary analysis"
  )
)

# A concise handoff index makes the interpretation order explicit and is also
# machine-readable for automated reports.
final_index <- data.frame(
  interpretation_order = 1:12,
  file = c(
    "sample_by_peak_detection_matrix_strict.pdf",
    "annotated_peak_table_with_sample_absorbance_strict.csv",
    "sample_by_peak_absorbance_strict.csv",
    "sample_by_peak_assignment_uncertainty_flag_strict.csv",
    "sample_by_peak_absorbance_upper_bound_strict.csv",
    "genus_hybrid_pseudochromatogram.pdf",
    "sample_by_peak_detection_matrix_review_features.pdf",
    "annotated_peak_table_with_sample_absorbance_complete.csv",
    "sample_by_peak_absorbance_complete.csv",
    "sample_by_peak_assignment_uncertainty_flag_complete.csv",
    "sample_by_peak_absorbance_upper_bound_complete.csv",
    "peak_picking_and_alignment_QC"
  ),
  purpose = c(
    "Check which analysis-ready peaks occur in each sample, in RT order.",
    "Interpret every strict peak and inspect each sample's DAD area.",
    "Use the NA-free confidently assigned abundances for primary analyses.",
    "Identify cells where detected signal had multiple plausible features.",
    paste(
      "Test each feature separately at its largest plausible abundance; do",
      "not sum upper-bound columns."
    ),
    "View the genus-level class-first chemical overview.",
    "Review provisional fragments and unresolved neighboring features.",
    "Retain all discovery features and annotations for audit.",
    "Retain all confidently assigned discovery-feature abundances for review.",
    "Locate unresolved candidate cells across all discovery features.",
    "Review feature-wise upper abundance limits; never sum these columns.",
    "If anything looks unusual, return to 01_Alignment_QC."
  ),
  stringsAsFactors = FALSE
)
write.csv(
  final_index, file.path(FINAL_DIR, "READ_FIRST_output_order.csv"),
  row.names = FALSE
)

# Copy only the three principal annotated chromatograms into the final package;
# detailed variants and scaling audits remain in the QC folder.
for (filename in c(
  "genus_class_pseudochromatogram.pdf",
  "genus_chromophore_pseudochromatogram.pdf",
  "genus_hybrid_pseudochromatogram.pdf"
)) {
  file.copy(
    file.path(ALIGNMENT_DIR, filename), file.path(FINAL_DIR, filename),
    overwrite = TRUE
  )
}

message("Annotated genus pseudochromatograms written to: ", ALIGNMENT_DIR)
message("Final collaborator results written to: ", FINAL_DIR)
