## ========================================================================== ##
## PRESENCE/ABSENCE-AWARE ALIGNMENT OF RAW AGILENT HPLC-DAD FILES
## ========================================================================== ##
## Purpose
##   Read a genus-level collection of Agilent `.D` folders directly, detect
##   chromatographic peaks within every sample, and align the resulting peaks
##   with both retention-time and full-spectrum evidence. The method is designed
##   for datasets that may contain different species and therefore may exhibit
##   substantial, biologically real peak presence/absence variation.
##
## Why this is separate from library annotation
##   Retention time is useful for deciding whether two observed peaks are
##   chromatographic homologues. It is not evidence that either peak is a
##   particular library compound. This script performs alignment without using
##   compound, chemical-class, or chromophore labels. The representative spectra
##   it writes can subsequently be passed to the fixed-library annotation model.
##
## Scientific safeguards
##   * No chromatogram is globally warped to look like one reference sample.
##   * A missing peak is allowed; nearby peaks are never forced into a group.
##   * Only high-quality RT/spectrum matches estimate conservative corrections.
##   * Corrections are limited to a small shift or a tightly constrained affine
##     adjustment; severe stretching is prohibited.
##   * Every aligned feature contains at most one peak from each sample.
##   * High-quality single-sample features are retained and labeled as rare.
##   * Low-quality ambiguous detections are reported rather than forced.
##
## Normal entry point
##   Rscript RUN_AGILENT_DAD.R Rhododendron \
##     --alignment-policy=retained_uniqueness
##   Rscript RUN_AGILENT_DAD.R Rhododendron \
##     --alignment-policy=expanded_rescue
##
## Direct use
##   Rscript 01_raw_agilent_dad_alignment.R Rhododendron retained_uniqueness
## ========================================================================== ##

if (!requireNamespace("chromConverter", quietly = TRUE)) {
  stop(
    "Package 'chromConverter' is required. Install it before running the ",
    "raw Agilent DAD route."
  )
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for alignment diagnostic figures.")
}

library(ggplot2)
source("00_library_inputs.R")

# ------------------------------- SETTINGS ---------------------------------
arguments <- commandArgs(trailingOnly = TRUE)
DATASET_NAME <- Sys.getenv(
  "UV_DATASET_NAME", unset = if (length(arguments)) arguments[1] else "Rhododendron"
)
DATASET_SLUG <- tolower(gsub("[^A-Za-z0-9]+", "_", DATASET_NAME))
DATASET_DISPLAY_NAME <- gsub(
  "([a-z])([A-Z])", "\\1 \\2", gsub("[_-]+", " ", DATASET_NAME)
)
RAW_DIR <- Sys.getenv(
  "UV_AGILENT_DIR", unset = file.path("data", paste0(DATASET_NAME, "_Agilent"))
)

# Require an explicit scientific policy. The retained-uniqueness preset reviews
# only 1-3-sample candidate groups for rescue. The expanded-rescue preset
# reviews groups up to 20 samples (or ceiling(20% of samples) for datasets with
# fewer than 20 samples). All spectral, RT, margin, bootstrap, and sample-
# collision safeguards are identical between presets.
ALIGNMENT_POLICY <- tolower(Sys.getenv(
  "UV_RAW_ALIGNMENT_POLICY",
  unset = if (length(arguments) >= 2L) arguments[2] else ""
))
VALID_ALIGNMENT_POLICIES <- c("retained_uniqueness", "expanded_rescue")
if (!ALIGNMENT_POLICY %in% VALID_ALIGNMENT_POLICIES) {
  stop(
    "Raw DAD alignment requires one explicit policy: retained_uniqueness or ",
    "expanded_rescue. Prefer RUN_AGILENT_DAD.R so policy-specific outputs ",
    "are routed through every downstream stage."
  )
}
ALIGNMENT_POLICY_FOLDER_LABEL <- c(
  retained_uniqueness = "RetainedUniqueness",
  expanded_rescue = "ExpandedRescue"
)[ALIGNMENT_POLICY]
ALIGNMENT_POLICY_DISPLAY_NAME <- Sys.getenv(
  "UV_RAW_ALIGNMENT_POLICY_DISPLAY_NAME",
  unset = c(
    retained_uniqueness = "Retained uniqueness",
    expanded_rescue = "Expanded rescue"
  )[ALIGNMENT_POLICY]
)
ALIGNMENT_FIGURE_TITLE_PREFIX <- paste(
  DATASET_DISPLAY_NAME, "-", ALIGNMENT_POLICY_DISPLAY_NAME
)
OUTPUT_DIR <- Sys.getenv(
  "UV_RAW_ALIGNMENT_OUTPUT_DIR",
  unset = file.path(
    paste0(DATASET_NAME, "_RawDAD_", ALIGNMENT_POLICY_FOLDER_LABEL),
    "01_Alignment_QC"
  )
)
ANNOTATION_DATASET <- Sys.getenv(
  "UV_RAW_ANNOTATION_DATASET",
  unset = paste0(DATASET_NAME, "_RawDAD_", ALIGNMENT_POLICY_FOLDER_LABEL)
)
ANNOTATION_OUTPUT_DIR <- Sys.getenv(
  "UV_RAW_ANNOTATION_OUTPUT_DIR",
  unset = file.path(
    paste0(DATASET_NAME, "_RawDAD_", ALIGNMENT_POLICY_FOLDER_LABEL),
    "02_Library_Annotation"
  )
)

# Peak detection is deliberately sensitive enough to retain rare chemistry but
# still requires both absolute profile signal and local chromatographic
# prominence. Spectral quality determines whether a detection may seed a group.
BASELINE_WINDOW_MINUTES <- 0.80
PROFILE_SMOOTHING_MINUTES <- 0.035
PROMINENCE_HALF_WINDOW_MINUTES <- 0.10
BACKGROUND_HALF_WINDOW_MINUTES <- 0.18
MAXIMUM_HALF_PEAK_WINDOW_MINUTES <- 0.12
MINIMUM_WITHIN_SAMPLE_PEAK_SEPARATION <- 0.020
# The permissive detection gate is intentionally below the core-landmark gate.
# A full 91-sample sensitivity audit showed that 4 / 2.5 retained every peak
# found at 5 / 3 while adding only 88 candidates (3.45%). This recovers weak,
# visually resolved peaks without allowing them to seed an alignment feature:
# the stricter 8 / 4 plus spectral-quality criteria below still define cores.
MINIMUM_PROFILE_SNR <- 4
MINIMUM_PROMINENCE_SNR <- 2.5
CORE_PROFILE_SNR <- 8
CORE_PROMINENCE_SNR <- 4
CORE_SPECTRAL_SNR <- 4
CORE_INFORMATIVE_FRACTION <- 0.025

# Peak-picking sensitivity is evaluated during every ordinary run. These
# neighboring settings are diagnostics, not competing hidden models: the
# deployed selection remains the explicitly named default above. A stable
# count plateau across these combinations is stronger evidence that the result
# is not an arbitrary threshold artifact than simply maximizing peak number.
PEAK_PICKING_SENSITIVITY_SMOOTHING_MINUTES <- c(0.020, 0.035, 0.050, 0.070)
PEAK_PICKING_SENSITIVITY_PROFILE_SNR <- c(3, 4, 5, 6, 8)
PEAK_PICKING_SENSITIVITY_PROMINENCE_SNR <- c(2, 2.5, 3, 4, 5)

# Coarse pairwise matches estimate only a sample-level shift. Final landmarks
# may support a small affine correction, but the slope is tightly constrained.
COARSE_RT_WINDOW <- 0.14
COARSE_SPECTRAL_DISTANCE <- 0.10
MINIMUM_PAIRWISE_LANDMARKS <- 3L
MAXIMUM_PAIRWISE_RT_MAD <- 0.030
MINIMUM_AFFINE_LANDMARKS <- 5L
MINIMUM_AFFINE_RT_SPAN <- 1.50
MINIMUM_ALLOWED_RT_SLOPE <- 0.97
MAXIMUM_ALLOWED_RT_SLOPE <- 1.03
MAXIMUM_ALLOWED_RT_INTERCEPT <- 0.12

# Complete-linkage grouping ensures every pair of core members satisfies both
# limits. This avoids transitive chains joining dissimilar peaks.
FINAL_CORE_RT_TOLERANCE <- 0.060
FINAL_CORE_SPECTRAL_TOLERANCE <- 0.15
# Complete linkage is intentionally strict, so a second agglomerative pass may
# reunite disjoint groups whose consensus spectra and RTs are effectively
# indistinguishable. The one-peak-per-sample rule is never relaxed.
CONSENSUS_MERGE_RT_TOLERANCE <- 0.035
CONSENSUS_MERGE_SPECTRAL_TOLERANCE <- 0.08
CONSENSUS_MERGE_MAXIMUM_RT_SPAN <- 0.090

# A strict first tree can leave a low-prevalence group just to the left or right
# of a better-supported homologue. The selected preset determines how broadly
# such groups are reviewed; it never changes the evidence required for rescue.
RARE_RECONCILIATION_ABSOLUTE_GROUP_CAP <- 20L
RARE_RECONCILIATION_SMALL_DATASET_FRACTION <- 0.20
RARE_RECONCILIATION_RETAINED_UNIQUENESS_CAP <- 3L
RARE_RECONCILIATION_MIN_ANCHOR_SAMPLES <- 4L
RARE_RECONCILIATION_RT_WINDOW <- 0.12
RARE_RECONCILIATION_SPECTRAL_TOLERANCE <- 0.07
RARE_RECONCILIATION_SCORE_MARGIN <- 0.12
RARE_RECONCILIATION_BOOTSTRAP_SUPPORT <- 0.85
RARE_REVIEW_RT_WINDOW <- 0.14
RARE_DISTINCT_SPECTRAL_DISTANCE <- 0.15
LOW_QUALITY_RT_TOLERANCE <- 0.065
LOW_QUALITY_SPECTRAL_TOLERANCE <- 0.21
MINIMUM_LOW_QUALITY_SCORE_MARGIN <- 0.12

# A final split-feature reconciliation addresses a different failure mode from
# rare-fragment rescue. Two groups are merged only when no sample contains both
# (mutual exclusion), their consensus RTs are extremely close, their spectra
# are very similar, the combined RT span remains compact, and wavelength-block
# resampling repeatedly supports spectral equivalence. Same-sample doublets are
# protected because co-occurrence is direct evidence for two resolved peaks.
SPLIT_FEATURE_MERGE_RT_TOLERANCE <- 0.035
SPLIT_FEATURE_MERGE_SPECTRAL_TOLERANCE <- 0.12
SPLIT_FEATURE_MERGE_MAXIMUM_RT_SPAN <- 0.10
SPLIT_FEATURE_MERGE_BOOTSTRAP_SUPPORT <- 0.85

# After strong merges, unstable neighboring pairs that remain mutually
# exclusive are treated as unresolved feature families. Their member columns
# stay separate in the complete matrix, both are excluded from the strict
# primary matrix, and a summed family column is supplied for sensitivity work.
COMPETITION_FAMILY_RT_WINDOW <- 0.12
COMPETITION_FAMILY_SPECTRAL_TOLERANCE <- 0.15
COMPETITION_ASSIGNMENT_SUPPORT_THRESHOLD <- 0.70

BOOTSTRAP_REPETITIONS <- as.integer(Sys.getenv(
  "UV_ALIGNMENT_BOOTSTRAPS", unset = "100"
))
BOOTSTRAP_BLOCK_WIDTH_NM <- 40
MAXIMUM_EXCLUDED_INPUT_FRACTION <- 0.05
RANDOM_SEED <- 20260811L
set.seed(RANDOM_SEED)

if (!dir.exists(RAW_DIR)) stop("Raw Agilent folder not found: ", RAW_DIR)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------- SMALL UTILITIES ------------------------------
odd_window_points <- function(minutes, retention_times, minimum = 3L) {
  step <- stats::median(diff(retention_times))
  width <- max(minimum, as.integer(round(minutes / step)))
  if (width %% 2L == 0L) width <- width + 1L
  maximum <- if (length(retention_times) %% 2L) {
    length(retention_times)
  } else length(retention_times) - 1L
  min(width, maximum)
}

centered_mean <- function(x, width) {
  width <- max(1L, as.integer(width))
  if (width %% 2L == 0L) width <- width + 1L
  filtered <- as.numeric(stats::filter(x, rep(1 / width, width), sides = 2))
  filtered[!is.finite(filtered)] <- x[!is.finite(filtered)]
  filtered
}

smooth_spectrum <- function(x, width = 5L) {
  if (length(x) < 3L) return(as.numeric(x))
  width <- min(as.integer(width), if (length(x) %% 2L) length(x) else length(x) - 1L)
  if (width < 3L) return(as.numeric(x))
  centered_mean(as.numeric(x), width)
}

normalize_spectrum <- function(x) {
  x <- smooth_spectrum(pmax(as.numeric(x), 0), 5L)
  span <- diff(range(x, na.rm = TRUE))
  if (!is.finite(span) || span <= sqrt(.Machine$double.eps)) {
    return(rep(NA_real_, length(x)))
  }
  (x - min(x, na.rm = TRUE)) / span
}

trapz_integral <- function(x, y) {
  if (length(x) < 2L) return(0)
  sum(diff(x) * (head(y, -1L) + tail(y, -1L)) / 2)
}

safe_mad <- function(x, fallback = 1e-8) {
  value <- stats::mad(x, center = stats::median(x, na.rm = TRUE),
                      constant = 1.4826, na.rm = TRUE)
  if (!is.finite(value) || value <= 0) fallback else value
}

nonmaximum_suppression <- function(indices, scores, retention_times,
                                   minimum_separation) {
  if (!length(indices)) return(integer())
  selected <- integer()
  for (index in indices[order(scores, decreasing = TRUE)]) {
    if (!length(selected) || all(abs(retention_times[index] -
                                     retention_times[selected]) >=
                                 minimum_separation)) {
      selected <- c(selected, index)
    }
  }
  sort(selected)
}

# Locate chromatographic maxima and calculate the two threshold statistics used
# by the picker. Keeping this logic in one helper guarantees that the deployed
# run and its parameter-sensitivity audit evaluate exactly the same candidates.
peak_candidate_metrics <- function(profile, retention_times) {
  profile_noise <- stats::mad(
    diff(profile), center = 0, constant = 1.4826, na.rm = TRUE
  ) / sqrt(2)
  if (!is.finite(profile_noise) || profile_noise <= 0) {
    positive <- profile[profile > 0]
    profile_noise <- if (length(positive)) {
      max(stats::quantile(positive, 0.10, na.rm = TRUE), 1e-8)
    } else 1e-8
  }
  maxima <- which(
    profile > c(-Inf, head(profile, -1L)) &
      profile >= c(tail(profile, -1L), -Inf)
  )
  edge_keep <- retention_times[maxima] >=
    min(retention_times) + PROMINENCE_HALF_WINDOW_MINUTES &
    retention_times[maxima] <=
    max(retention_times) - PROMINENCE_HALF_WINDOW_MINUTES
  maxima <- maxima[edge_keep]
  prominence <- vapply(maxima, function(index) {
    left <- which(
      retention_times >= retention_times[index] -
        PROMINENCE_HALF_WINDOW_MINUTES &
        retention_times < retention_times[index]
    )
    right <- which(
      retention_times > retention_times[index] &
        retention_times <= retention_times[index] +
          PROMINENCE_HALF_WINDOW_MINUTES
    )
    if (!length(left) || !length(right)) return(0)
    profile[index] - max(min(profile[left]), min(profile[right]))
  }, numeric(1))
  data.frame(
    apex_index = maxima,
    profile_snr = profile[maxima] / profile_noise,
    prominence = prominence,
    prominence_snr = prominence / profile_noise,
    profile_noise = rep(profile_noise, length(maxima)),
    stringsAsFactors = FALSE
  )
}

select_peak_candidates <- function(candidate_metrics, retention_times,
                                   minimum_profile_snr,
                                   minimum_prominence_snr) {
  keep <- candidate_metrics$profile_snr >= minimum_profile_snr &
    candidate_metrics$prominence_snr >= minimum_prominence_snr
  if (!any(keep)) return(integer())
  selected_rows <- which(keep)
  selected_apices <- nonmaximum_suppression(
    candidate_metrics$apex_index[selected_rows],
    candidate_metrics$prominence_snr[selected_rows],
    retention_times, MINIMUM_WITHIN_SAMPLE_PEAK_SEPARATION
  )
  selected_apices
}

format_feature_id <- function(retention_time, rank_at_time) {
  rt_text <- sub("0+$", "", sub("[.]0+$", "", sprintf("%.3f", retention_time)))
  paste0("X", rt_text, "_A", sprintf("%03d", rank_at_time))
}

matrix_to_long <- function(matrix_value, row_name = "sample_id",
                           column_name = "aligned_feature_id",
                           value_name = "value") {
  indices <- which(matrix(TRUE, nrow(matrix_value), ncol(matrix_value)),
                   arr.ind = TRUE)
  output <- data.frame(
    stringsAsFactors = FALSE,
    row = rownames(matrix_value)[indices[, 1]],
    column = colnames(matrix_value)[indices[, 2]],
    value = matrix_value[indices]
  )
  names(output) <- c(row_name, column_name, value_name)
  output
}

# Composite spectral distance uses concentration-invariant shape. Cosine,
# centered correlation, and first-derivative agreement capture complementary
# aspects of a DAD spectrum without fitting a supervised model.
normalize_rows_l2 <- function(x) {
  x <- as.matrix(x)
  norms <- sqrt(rowSums(x^2))
  norms[!is.finite(norms) | norms <= sqrt(.Machine$double.eps)] <- 1
  x / norms
}

spectral_distance_matrix <- function(spectra) {
  spectra <- as.matrix(spectra)
  l2 <- normalize_rows_l2(spectra)
  cosine <- 1 - tcrossprod(l2)
  cosine[cosine < 0] <- 0
  cosine[cosine > 1] <- 1

  centered <- spectra - rowMeans(spectra)
  centered <- normalize_rows_l2(centered)
  correlation <- (1 - tcrossprod(centered)) / 2
  correlation[correlation < 0] <- 0
  correlation[correlation > 1] <- 1

  derivative <- t(apply(spectra, 1, diff))
  derivative <- derivative - rowMeans(derivative)
  derivative <- normalize_rows_l2(derivative)
  derivative_distance <- (1 - tcrossprod(derivative)) / 2
  derivative_distance[derivative_distance < 0] <- 0
  derivative_distance[derivative_distance > 1] <- 1

  output <- 0.55 * cosine + 0.30 * correlation +
    0.15 * derivative_distance
  diag(output) <- 0
  output
}

spectral_distance_to_references <- function(query, references) {
  query <- matrix(query, nrow = 1)
  references <- as.matrix(references)

  query_l2 <- normalize_rows_l2(query)
  reference_l2 <- normalize_rows_l2(references)
  cosine <- pmax(0, pmin(1, 1 - as.numeric(query_l2 %*% t(reference_l2))))

  query_centered <- normalize_rows_l2(query - rowMeans(query))
  reference_centered <- normalize_rows_l2(references - rowMeans(references))
  correlation <- pmax(
    0, pmin(1, (1 - as.numeric(query_centered %*% t(reference_centered))) / 2)
  )

  query_derivative <- diff(as.numeric(query))
  query_derivative <- matrix(query_derivative - mean(query_derivative), nrow = 1)
  reference_derivative <- t(apply(references, 1, diff))
  reference_derivative <- reference_derivative - rowMeans(reference_derivative)
  derivative_distance <- pmax(
    0,
    pmin(1, (1 - as.numeric(normalize_rows_l2(query_derivative) %*%
                              t(normalize_rows_l2(reference_derivative)))) / 2)
  )

  0.55 * cosine + 0.30 * correlation + 0.15 * derivative_distance
}

connected_components <- function(number_of_nodes, edges) {
  neighbors <- vector("list", number_of_nodes)
  if (nrow(edges)) {
    for (row in seq_len(nrow(edges))) {
      first <- edges$sample_i[row]
      second <- edges$sample_j[row]
      neighbors[[first]] <- unique(c(neighbors[[first]], second))
      neighbors[[second]] <- unique(c(neighbors[[second]], first))
    }
  }
  component <- rep(NA_integer_, number_of_nodes)
  component_id <- 0L
  for (node in seq_len(number_of_nodes)) {
    if (!is.na(component[node])) next
    component_id <- component_id + 1L
    queue <- node
    component[node] <- component_id
    while (length(queue)) {
      current <- queue[1]
      queue <- queue[-1]
      unseen <- neighbors[[current]][is.na(component[neighbors[[current]]])]
      if (length(unseen)) {
        component[unseen] <- component_id
        queue <- c(queue, unseen)
      }
    }
  }
  component
}

# -------------------------- RAW FILE DISCOVERY -----------------------------
sample_directories <- sort(list.dirs(RAW_DIR, recursive = FALSE,
                                      full.names = TRUE))
sample_directories <- sample_directories[grepl("[.]D$", sample_directories,
                                                ignore.case = TRUE)]
if (!length(sample_directories)) {
  stop("No Agilent .D folders were found in: ", RAW_DIR)
}

find_dad_file <- function(directory) {
  candidates <- list.files(directory, pattern = "[.]uv$", full.names = TRUE,
                           ignore.case = TRUE)
  preferred <- candidates[toupper(basename(candidates)) == "DAD1.UV"]
  if (length(preferred)) return(preferred[1])
  if (length(candidates) == 1L) return(candidates)
  if (!length(candidates)) stop("No .uv DAD file found in: ", directory)
  stop("More than one .uv file was found and DAD1.uv was absent in: ", directory)
}

dad_files <- vapply(sample_directories, find_dad_file, character(1))
discovered_sample_ids <- tools::file_path_sans_ext(basename(sample_directories))
if (anyDuplicated(discovered_sample_ids)) {
  stop("Sample names derived from .D folders are not unique.")
}

# Validate every binary DAD cube before alignment. A truncated `.uv` file can
# fail inside the vendor decoder before ordinary coordinate checks are reached.
# Such an injection must be excluded, not represented as biological absence or
# an all-zero sample. The exclusion table remains part of the run manifest, and
# a dataset with more than 5% unreadable inputs stops for upstream review.
preflight_errors <- vapply(seq_along(dad_files), function(i) {
  tryCatch({
    test <- chromConverter::read_chemstation_uv(
      dad_files[i], format_out = "matrix", data_format = "wide",
      read_metadata = FALSE
    )
    test_rt <- suppressWarnings(as.numeric(rownames(test)))
    test_wavelengths <- suppressWarnings(as.numeric(colnames(test)))
    if (nrow(test) < 3L || ncol(test) < 3L ||
        any(!is.finite(test_rt)) || any(!is.finite(test_wavelengths)) ||
        any(!is.finite(test))) {
      stop("decoded cube has non-finite coordinates, absorbance, or dimensions")
    }
    ""
  }, error = function(e) paste(
    "Binary DAD decode failed (truncated file or malformed wavelength metadata):",
    conditionMessage(e)
  ))
}, character(1))
excluded_input <- nzchar(preflight_errors)
raw_input_exclusions <- data.frame(
  sample_id = discovered_sample_ids[excluded_input],
  sample_directory = normalizePath(
    sample_directories[excluded_input], mustWork = FALSE
  ),
  dad_file = normalizePath(dad_files[excluded_input], mustWork = FALSE),
  file_size_bytes = unname(file.info(dad_files[excluded_input])$size),
  exclusion_reason = preflight_errors[excluded_input],
  exclusion_action = rep(
    "Excluded from alignment and all sample matrices", sum(excluded_input)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  raw_input_exclusions,
  file.path(OUTPUT_DIR, "raw_agilent_input_exclusions.csv"),
  row.names = FALSE
)
if (mean(excluded_input) > MAXIMUM_EXCLUDED_INPUT_FRACTION) {
  stop(
    sum(excluded_input), " of ", length(excluded_input),
    " raw DAD inputs were unreadable (> ",
    100 * MAXIMUM_EXCLUDED_INPUT_FRACTION,
    "%). See raw_agilent_input_exclusions.csv."
  )
}
if (any(excluded_input)) {
  warning(
    "Excluded ", sum(excluded_input), " unreadable raw DAD input(s): ",
    paste(discovered_sample_ids[excluded_input], collapse = ", "),
    ". These samples are absent (not zero-filled) in final matrices."
  )
}
sample_directories <- sample_directories[!excluded_input]
dad_files <- dad_files[!excluded_input]
sample_ids <- discovered_sample_ids[!excluded_input]
if (!length(sample_ids)) stop("No readable raw DAD inputs remain after preflight.")
names(dad_files) <- sample_ids

# Compute the low-prevalence review boundary from the chosen policy. The
# boundary controls which groups receive extra scrutiny, not which detections
# are automatically moved or deleted.
RARE_RECONCILIATION_MAX_GROUP_SIZE <- if (
  ALIGNMENT_POLICY == "retained_uniqueness"
) {
  min(length(sample_ids), RARE_RECONCILIATION_RETAINED_UNIQUENESS_CAP)
} else if (length(sample_ids) < RARE_RECONCILIATION_ABSOLUTE_GROUP_CAP) {
  max(1L, as.integer(ceiling(
    RARE_RECONCILIATION_SMALL_DATASET_FRACTION * length(sample_ids)
  )))
} else RARE_RECONCILIATION_ABSOLUTE_GROUP_CAP
message(
  "Alignment policy: ", ALIGNMENT_POLICY_DISPLAY_NAME,
  ". Fragment review covers groups of 1-",
  RARE_RECONCILIATION_MAX_GROUP_SIZE, " samples (", length(sample_ids),
  " total samples)."
)

# The authenticated library supplies a stable comparison grid. This also makes
# future raw files sampled every 3 nm compatible with the current 5 nm library.
library_bundle <- load_library_inputs()
library_column_values <- suppressWarnings(as.numeric(sub(
  "(?i)nm$", "", names(library_bundle$spectra), perl = TRUE
)))
analysis_wavelengths <- library_column_values[
  is.finite(library_column_values) & library_column_values >= 190 &
    library_column_values <= 900
]
analysis_wavelengths <- sort(unique(analysis_wavelengths))
if (length(analysis_wavelengths) < 3L) {
  stop("The deployed spectral library does not provide a usable 190-900 nm grid.")
}

# ------------------- DETECT PEAKS IN EACH CHROMATOGRAM --------------------
all_detection_rows <- list()
all_shape_rows <- list()
all_integrated_rows <- list()
input_audit_rows <- list()
wavelength_audit_rows <- list()
peak_picking_sensitivity_rows <- list()
peak_profile_rows <- list()
detection_counter <- 1L

detect_sample_peaks <- function(sample_id, path) {
  raw <- chromConverter::read_chemstation_uv(
    path, format_out = "matrix", data_format = "wide", read_metadata = TRUE
  )
  retention_times <- suppressWarnings(as.numeric(rownames(raw)))
  native_wavelengths <- suppressWarnings(as.numeric(colnames(raw)))
  if (any(!is.finite(retention_times)) || any(!is.finite(native_wavelengths)) ||
      any(!is.finite(raw))) {
    stop("Non-finite coordinates or absorbance in: ", path)
  }
  if (is.unsorted(retention_times, strictly = TRUE)) {
    stop("Retention times are not strictly increasing in: ", path)
  }

  aligned_grid <- harmonize_wavelength_matrix(
    raw, native_wavelengths, analysis_wavelengths, input_id = sample_id
  )
  absorbance <- aligned_grid$values
  wavelengths <- analysis_wavelengths

  baseline_width <- odd_window_points(
    BASELINE_WINDOW_MINUTES, retention_times, minimum = 9L
  )
  baseline <- vapply(seq_len(ncol(absorbance)), function(column) {
    stats::runmed(absorbance[, column], baseline_width, endrule = "median")
  }, numeric(nrow(absorbance)))
  corrected <- pmax(absorbance - baseline, 0)

  raw_profile <- sqrt(rowMeans(corrected^2))
  smoothing_width <- odd_window_points(
    PROFILE_SMOOTHING_MINUTES, retention_times, minimum = 3L
  )
  profile <- centered_mean(raw_profile, smoothing_width)
  candidate_metrics <- peak_candidate_metrics(profile, retention_times)
  profile_noise <- unique(candidate_metrics$profile_noise)
  if (!length(profile_noise)) {
    profile_noise <- stats::mad(
      diff(profile), center = 0, constant = 1.4826, na.rm = TRUE
    ) / sqrt(2)
  }
  maxima <- select_peak_candidates(
    candidate_metrics, retention_times,
    MINIMUM_PROFILE_SNR, MINIMUM_PROMINENCE_SNR
  )

  # Evaluate small, scientifically plausible perturbations around the deployed
  # picker. This works on the same baseline-corrected chromatographic profile
  # and therefore adds negligible memory compared with rereading every `.D`
  # file in separate trial runs.
  sensitivity_grid <- expand.grid(
    smoothing_minutes = PEAK_PICKING_SENSITIVITY_SMOOTHING_MINUTES,
    minimum_profile_snr = PEAK_PICKING_SENSITIVITY_PROFILE_SNR,
    minimum_prominence_snr = PEAK_PICKING_SENSITIVITY_PROMINENCE_SNR,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  sensitivity_grid <- sensitivity_grid[
    abs(sensitivity_grid$minimum_profile_snr - MINIMUM_PROFILE_SNR) <= 4 &
      abs(sensitivity_grid$minimum_prominence_snr -
            MINIMUM_PROMINENCE_SNR) <= 2.5, , drop = FALSE
  ]
  sensitivity_grid$sample_id <- sample_id
  sensitivity_results <- lapply(seq_len(nrow(sensitivity_grid)), function(row) {
      width <- odd_window_points(
        sensitivity_grid$smoothing_minutes[row], retention_times, minimum = 3L
      )
      trial_profile <- centered_mean(raw_profile, width)
      trial_metrics <- peak_candidate_metrics(trial_profile, retention_times)
      trial_apices <- select_peak_candidates(
        trial_metrics, retention_times,
        sensitivity_grid$minimum_profile_snr[row],
        sensitivity_grid$minimum_prominence_snr[row]
      )
      match_window <- MINIMUM_WITHIN_SAMPLE_PEAK_SEPARATION / 2
      deployed_recovered <- if (length(maxima)) sum(vapply(
        maxima, function(apex) {
          any(abs(retention_times[apex] -
                    retention_times[trial_apices]) <= match_window)
        }, logical(1)
      )) else 0L
      trial_matched <- if (length(trial_apices)) sum(vapply(
        trial_apices, function(apex) {
          any(abs(retention_times[apex] -
                    retention_times[maxima]) <= match_window)
        }, logical(1)
      )) else 0L
      c(
        selected_peak_count = length(trial_apices),
        deployed_peaks_recovered = deployed_recovered,
        trial_peaks_matching_deployed = trial_matched,
        trial_peaks_not_matching_deployed = length(trial_apices) - trial_matched,
        deployed_peak_recovery_fraction = if (length(maxima)) {
          deployed_recovered / length(maxima)
        } else NA_real_,
        trial_match_fraction = if (length(trial_apices)) {
          trial_matched / length(trial_apices)
        } else NA_real_
      )
    })
  sensitivity_results <- as.data.frame(do.call(rbind, sensitivity_results))
  sensitivity_grid <- cbind(sensitivity_grid, sensitivity_results)
  sensitivity_grid$is_deployed_setting <- with(
    sensitivity_grid,
    abs(smoothing_minutes - PROFILE_SMOOTHING_MINUTES) < 1e-10 &
      minimum_profile_snr == MINIMUM_PROFILE_SNR &
      minimum_prominence_snr == MINIMUM_PROMINENCE_SNR
  )

  rows <- list()
  shapes <- list()
  integrated_spectra <- list()
  for (apex in maxima) {
    local <- which(retention_times >= retention_times[apex] -
                     PROMINENCE_HALF_WINDOW_MINUTES &
                     retention_times <= retention_times[apex] +
                     PROMINENCE_HALF_WINDOW_MINUTES)
    local_floor <- min(profile[local])
    local_prominence <- max(profile[apex] - local_floor, 0)
    half_height <- local_floor + 0.5 * local_prominence
    maximum_distance <- MAXIMUM_HALF_PEAK_WINDOW_MINUTES

    left <- apex
    while (left > 1L && profile[left - 1L] > half_height &&
           retention_times[apex] - retention_times[left - 1L] <=
             maximum_distance) left <- left - 1L
    right <- apex
    while (right < length(profile) && profile[right + 1L] > half_height &&
           retention_times[right + 1L] - retention_times[apex] <=
             maximum_distance) right <- right + 1L
    peak_indices <- seq.int(left, right)

    background_indices <- which(
      retention_times >= retention_times[apex] - BACKGROUND_HALF_WINDOW_MINUTES &
        retention_times <= retention_times[apex] + BACKGROUND_HALF_WINDOW_MINUTES &
        (retention_times < retention_times[left] - 0.01 |
           retention_times > retention_times[right] + 0.01)
    )
    if (length(background_indices) < 6L) {
      background_indices <- setdiff(local, peak_indices)
    }
    if (length(background_indices) < 3L) next

    local_background <- apply(absorbance[background_indices, , drop = FALSE],
                              2, stats::median)
    background_noise <- apply(
      absorbance[background_indices, , drop = FALSE], 2, safe_mad
    )
    positive_peak <- pmax(
      absorbance[peak_indices, , drop = FALSE] -
        matrix(local_background, nrow = length(peak_indices),
               ncol = length(local_background), byrow = TRUE),
      0
    )
    apex_spectrum <- pmax(absorbance[apex, ] - local_background, 0)
    integrated_spectrum <- vapply(seq_len(ncol(positive_peak)), function(j) {
      trapz_integral(retention_times[peak_indices], positive_peak[, j])
    }, numeric(1))
    normalized_shape <- normalize_spectrum(integrated_spectrum)
    if (any(!is.finite(normalized_shape))) next

    spectral_ratios <- apex_spectrum / pmax(background_noise, 1e-8)
    spectral_snr <- as.numeric(stats::quantile(
      spectral_ratios, 0.90, names = FALSE, na.rm = TRUE
    ))
    informative_fraction <- mean(spectral_ratios >= 3)
    wavelength_step <- stats::median(diff(wavelengths))
    integrated_dad_area <- sum(integrated_spectrum) * wavelength_step
    apex_summed_absorbance <- sum(apex_spectrum) * wavelength_step
    fwhm <- retention_times[right] - retention_times[left]
    local_profile_snr <- profile[apex] / profile_noise
    local_prominence_snr <- local_prominence / profile_noise
    quality_score <- sqrt(log1p(local_profile_snr) *
                            log1p(local_prominence_snr)) *
      log1p(max(spectral_snr, 0)) * (0.5 + informative_fraction)

    rows[[length(rows) + 1L]] <- data.frame(
      sample_id = sample_id,
      raw_retention_time = retention_times[apex],
      apex_index = apex,
      left_retention_time = retention_times[left],
      right_retention_time = retention_times[right],
      fwhm_minutes = fwhm,
      profile_height = profile[apex],
      profile_snr = local_profile_snr,
      prominence = local_prominence,
      prominence_snr = local_prominence_snr,
      spectral_shape_snr = spectral_snr,
      informative_wavelength_fraction = informative_fraction,
      integrated_dad_area = integrated_dad_area,
      apex_summed_absorbance = apex_summed_absorbance,
      maximum_spectral_wavelength_nm = wavelengths[which.max(normalized_shape)],
      quality_score = quality_score,
      stringsAsFactors = FALSE
    )
    shapes[[length(shapes) + 1L]] <- normalized_shape
    integrated_spectra[[length(integrated_spectra) + 1L]] <- integrated_spectrum
  }

  detection_table <- if (length(rows)) do.call(rbind, rows) else data.frame()
  shape_matrix <- if (length(shapes)) do.call(rbind, shapes) else
    matrix(numeric(), nrow = 0, ncol = length(wavelengths))
  integrated_matrix <- if (length(integrated_spectra)) {
    do.call(rbind, integrated_spectra)
  } else matrix(numeric(), nrow = 0, ncol = length(wavelengths))
  colnames(shape_matrix) <- colnames(integrated_matrix) <- wavelengths

  input_audit <- data.frame(
    sample_id = sample_id,
    dad_file = normalizePath(path),
    file_version = as.character(attr(raw, "file_version")),
    acquisition_method = as.character(attr(raw, "method")),
    native_time_points = nrow(raw),
    native_retention_time_minimum = min(retention_times),
    native_retention_time_maximum = max(retention_times),
    native_median_time_step_seconds = 60 * stats::median(diff(retention_times)),
    native_wavelength_count = length(native_wavelengths),
    native_wavelength_minimum_nm = min(native_wavelengths),
    native_wavelength_maximum_nm = max(native_wavelengths),
    detected_peaks = nrow(detection_table),
    stringsAsFactors = FALSE
  )
  peak_profile <- data.frame(
    sample_id = sample_id,
    retention_time = retention_times,
    baseline_corrected_rms_profile = raw_profile,
    smoothed_rms_profile = profile,
    profile_noise = profile_noise,
    profile_detection_height = MINIMUM_PROFILE_SNR * profile_noise,
    stringsAsFactors = FALSE
  )
  list(
    detections = detection_table,
    shapes = shape_matrix,
    integrated = integrated_matrix,
    input_audit = input_audit,
    wavelength_audit = aligned_grid$audit,
    peak_picking_sensitivity = sensitivity_grid,
    peak_profile = peak_profile
  )
}

message("Reading raw Agilent DAD files and detecting peaks: ", DATASET_NAME)
for (sample_index in seq_along(dad_files)) {
  sample_id <- names(dad_files)[sample_index]
  result <- detect_sample_peaks(sample_id, dad_files[sample_index])
  if (nrow(result$detections)) {
    identifiers <- paste0("D", sprintf("%05d", detection_counter +
                                         seq_len(nrow(result$detections)) - 1L))
    result$detections$detection_id <- identifiers
    rownames(result$shapes) <- rownames(result$integrated) <- identifiers
    all_detection_rows[[length(all_detection_rows) + 1L]] <- result$detections
    all_shape_rows[[length(all_shape_rows) + 1L]] <- result$shapes
    all_integrated_rows[[length(all_integrated_rows) + 1L]] <- result$integrated
    detection_counter <- detection_counter + nrow(result$detections)
  }
  input_audit_rows[[length(input_audit_rows) + 1L]] <- result$input_audit
  wavelength_audit_rows[[length(wavelength_audit_rows) + 1L]] <-
    result$wavelength_audit
  peak_picking_sensitivity_rows[[
    length(peak_picking_sensitivity_rows) + 1L
  ]] <- result$peak_picking_sensitivity
  peak_profile_rows[[length(peak_profile_rows) + 1L]] <- result$peak_profile
  # A raw DAD cube is much larger than its extracted peaks. Explicit collection
  # prevents dozens of already-processed binary matrices from remaining in R's
  # allocation pool during genus-scale runs on memory-constrained computers.
  rm(result)
  if (sample_index %% 5L == 0L) invisible(gc(verbose = FALSE))
  if (sample_index %% 10L == 0L || sample_index == length(dad_files)) {
    message("  processed ", sample_index, " / ", length(dad_files), " samples")
  }
}

if (!length(all_detection_rows)) stop("No peaks passed the detection thresholds.")
detections <- do.call(rbind, all_detection_rows)
rownames(detections) <- NULL
shape_spectra <- do.call(rbind, all_shape_rows)
integrated_spectra <- do.call(rbind, all_integrated_rows)
input_audit <- do.call(rbind, input_audit_rows)
wavelength_audit <- do.call(rbind, wavelength_audit_rows)
peak_picking_sensitivity <- do.call(rbind, peak_picking_sensitivity_rows)
peak_profiles <- do.call(rbind, peak_profile_rows)

order_index <- match(detections$detection_id, rownames(shape_spectra))
if (anyNA(order_index)) stop("Internal detection-spectrum ordering failure.")
shape_spectra <- shape_spectra[order_index, , drop = FALSE]
integrated_spectra <- integrated_spectra[order_index, , drop = FALSE]
rownames(shape_spectra) <- rownames(integrated_spectra) <- detections$detection_id

detections$core_detection <- with(
  detections,
  profile_snr >= CORE_PROFILE_SNR &
    prominence_snr >= CORE_PROMINENCE_SNR &
    spectral_shape_snr >= CORE_SPECTRAL_SNR &
    informative_wavelength_fraction >= CORE_INFORMATIVE_FRACTION
)

# -------------------- PAIRWISE CONSERVATIVE RT SHIFTS ---------------------
message("Estimating conservative pairwise retention-time shifts")
core_indices <- which(detections$core_detection)
message(
  "  ", nrow(detections), " peaks passed detection; ", length(core_indices),
  " high-quality peaks may seed alignment groups"
)
if (!length(core_indices)) {
  stop("No high-quality core peaks were available for RT-spectral alignment.")
}
# Only core detections require an all-by-all distance matrix. Computing a dense
# matrix for every weak candidate can exhaust memory in chemically diverse
# genera and is unnecessary because weak detections are compared with compact
# group consensuses later.
core_spectral_distances <- spectral_distance_matrix(
  shape_spectra[core_indices, , drop = FALSE]
)
rownames(core_spectral_distances) <- colnames(core_spectral_distances) <-
  as.character(core_indices)
sample_index_by_detection <- match(detections$sample_id, sample_ids)
pairwise_rows <- list()
pair_counter <- 1L

for (first_sample in seq_len(length(sample_ids) - 1L)) {
  first_indices <- which(sample_index_by_detection == first_sample &
                           detections$core_detection)
  if (!length(first_indices)) next
  for (second_sample in seq.int(first_sample + 1L, length(sample_ids))) {
    second_indices <- which(sample_index_by_detection == second_sample &
                              detections$core_detection)
    if (!length(second_indices)) next
    rt_difference <- outer(
      detections$raw_retention_time[first_indices],
      detections$raw_retention_time[second_indices], "-"
    )
    spectral <- core_spectral_distances[
      as.character(first_indices), as.character(second_indices), drop = FALSE
    ]
    eligible <- abs(rt_difference) <= COARSE_RT_WINDOW &
      spectral <= COARSE_SPECTRAL_DISTANCE
    candidates <- which(eligible, arr.ind = TRUE)
    if (!nrow(candidates)) next
    costs <- abs(rt_difference[candidates]) / COARSE_RT_WINDOW +
      spectral[candidates] / COARSE_SPECTRAL_DISTANCE
    candidates <- candidates[order(costs), , drop = FALSE]
    used_first <- used_second <- integer()
    selected <- matrix(integer(), ncol = 2)
    for (candidate_row in seq_len(nrow(candidates))) {
      current_first <- candidates[candidate_row, 1]
      current_second <- candidates[candidate_row, 2]
      if (current_first %in% used_first || current_second %in% used_second) next
      selected <- rbind(selected, c(current_first, current_second))
      used_first <- c(used_first, current_first)
      used_second <- c(used_second, current_second)
    }
    if (nrow(selected) < MINIMUM_PAIRWISE_LANDMARKS) next
    shifts <- rt_difference[selected]
    median_shift <- stats::median(shifts)
    residual_mad <- stats::mad(shifts, center = median_shift,
                               constant = 1.4826)
    if (!is.finite(residual_mad)) residual_mad <- Inf
    if (residual_mad > MAXIMUM_PAIRWISE_RT_MAD) next
    pairwise_rows[[pair_counter]] <- data.frame(
      sample_i = first_sample,
      sample_j = second_sample,
      sample_i_id = sample_ids[first_sample],
      sample_j_id = sample_ids[second_sample],
      matched_landmarks = nrow(selected),
      median_rt_difference = median_shift,
      rt_difference_mad = residual_mad,
      edge_weight = nrow(selected) / (residual_mad + 0.005)^2,
      stringsAsFactors = FALSE
    )
    pair_counter <- pair_counter + 1L
  }
}

pairwise_edges <- if (length(pairwise_rows)) {
  do.call(rbind, pairwise_rows)
} else data.frame(
  sample_i = integer(), sample_j = integer(), sample_i_id = character(),
  sample_j_id = character(), matched_landmarks = integer(),
  median_rt_difference = numeric(), rt_difference_mad = numeric(),
  edge_weight = numeric()
)

components <- connected_components(length(sample_ids), pairwise_edges)
initial_shift <- numeric(length(sample_ids))
for (component_id in sort(unique(components))) {
  members <- which(components == component_id)
  component_edges <- pairwise_edges[
    pairwise_edges$sample_i %in% members &
      pairwise_edges$sample_j %in% members, , drop = FALSE
  ]
  if (!nrow(component_edges) || length(members) == 1L) next
  incidence <- matrix(0, nrow(component_edges), length(members))
  incidence[cbind(seq_len(nrow(component_edges)),
                  match(component_edges$sample_i, members))] <- 1
  incidence[cbind(seq_len(nrow(component_edges)),
                  match(component_edges$sample_j, members))] <- -1
  sqrt_weight <- sqrt(component_edges$edge_weight /
                        stats::median(component_edges$edge_weight))
  weighted_incidence <- incidence * sqrt_weight
  weighted_difference <- component_edges$median_rt_difference * sqrt_weight
  # The final row centers the component at zero rather than privileging one
  # species or chromatogram as the reference.
  solution <- qr.solve(
    rbind(weighted_incidence, rep(1, length(members))),
    c(weighted_difference, 0), tol = 1e-10
  )
  initial_shift[members] <- solution
}
names(initial_shift) <- sample_ids
detections$initial_rt_shift <- initial_shift[detections$sample_id]
detections$initial_aligned_retention_time <-
  detections$raw_retention_time - detections$initial_rt_shift

# --------------------- PRELIMINARY LANDMARK GROUPS ------------------------
complete_linkage_groups <- function(indices, aligned_times, spectral_distances,
                                    rt_tolerance, spectral_tolerance,
                                    sample_labels) {
  if (!length(indices)) return(integer())
  if (length(indices) == 1L) return(setNames(1L, indices))
  rt_distance <- abs(outer(aligned_times[indices], aligned_times[indices], "-")) /
    rt_tolerance
  spectral_distance <- spectral_distances[
    as.character(indices), as.character(indices), drop = FALSE
  ] /
    spectral_tolerance
  joint <- pmax(rt_distance, spectral_distance)
  same_sample <- outer(sample_labels[indices], sample_labels[indices], "==")
  joint[same_sample & row(joint) != col(joint)] <- 1e6
  diag(joint) <- 0
  tree <- stats::hclust(stats::as.dist(joint), method = "complete")
  groups <- stats::cutree(tree, h = 1)
  setNames(as.integer(groups), indices)
}

preliminary_groups <- complete_linkage_groups(
  core_indices, detections$initial_aligned_retention_time,
  core_spectral_distances, FINAL_CORE_RT_TOLERANCE + 0.015,
  FINAL_CORE_SPECTRAL_TOLERANCE, detections$sample_id
)
detections$preliminary_group <- NA_integer_
detections$preliminary_group[as.integer(names(preliminary_groups))] <-
  preliminary_groups

preliminary_summary <- do.call(rbind, lapply(
  sort(unique(preliminary_groups)), function(group_id) {
    members <- as.integer(names(preliminary_groups)[preliminary_groups == group_id])
    data.frame(
      preliminary_group = group_id,
      n_samples = length(unique(detections$sample_id[members])),
      consensus_rt = stats::median(
        detections$initial_aligned_retention_time[members]
      ),
      rt_span = diff(range(detections$initial_aligned_retention_time[members])),
      stringsAsFactors = FALSE
    )
  }
))
landmark_groups <- preliminary_summary$preliminary_group[
  preliminary_summary$n_samples >= max(4L, ceiling(0.05 * length(sample_ids))) &
    preliminary_summary$rt_span <= 0.10
]
landmark_reference <- setNames(
  preliminary_summary$consensus_rt, preliminary_summary$preliminary_group
)

# -------------------- CONSTRAINED SAMPLE-SPECIFIC RT MODEL ----------------
rt_model_rows <- list()
for (sample_id in sample_ids) {
  indices <- which(detections$sample_id == sample_id &
                     detections$preliminary_group %in% landmark_groups)
  if (length(indices)) {
    # Complete linkage prevents duplicate sample members within a group, but the
    # explicit check makes this contract auditable if settings change later.
    indices <- indices[!duplicated(detections$preliminary_group[indices])]
  }
  reference_rt <- landmark_reference[
    as.character(detections$preliminary_group[indices])
  ]
  observed_rt <- detections$raw_retention_time[indices]
  starting_shift <- initial_shift[sample_id]
  intercept <- starting_shift
  slope <- 1
  model_type <- if (length(indices)) "Pairwise shift" else "No supported correction"
  rejected_affine <- FALSE

  if (length(indices) >= MINIMUM_AFFINE_LANDMARKS &&
      diff(range(reference_rt)) >= MINIMUM_AFFINE_RT_SPAN) {
    keep <- rep(TRUE, length(indices))
    fitted <- NULL
    for (iteration in seq_len(3L)) {
      # Fit with ordinary column names so predict() can safely evaluate all
      # landmarks, including those temporarily excluded as residual outliers.
      fit_data <- data.frame(
        observed = observed_rt[keep],
        reference = reference_rt[keep]
      )
      fitted <- stats::lm(observed ~ reference, data = fit_data)
      residuals <- observed_rt - stats::predict(
        fitted, newdata = data.frame(reference = reference_rt)
      )
      scale <- safe_mad(residuals[keep], fallback = 0.01)
      new_keep <- abs(residuals) <= max(0.025, 3 * scale)
      if (sum(new_keep) < MINIMUM_AFFINE_LANDMARKS || identical(new_keep, keep)) break
      keep <- new_keep
    }
    coefficients <- stats::coef(fitted)
    proposed_intercept <- unname(coefficients[1])
    proposed_slope <- unname(coefficients[2])
    proposed_rmse <- sqrt(mean((observed_rt -
      (proposed_intercept + proposed_slope * reference_rt))^2))
    shift_rmse <- sqrt(mean((observed_rt -
      (starting_shift + reference_rt))^2))
    allowed <- is.finite(proposed_intercept) && is.finite(proposed_slope) &&
      abs(proposed_intercept) <= MAXIMUM_ALLOWED_RT_INTERCEPT &&
      proposed_slope >= MINIMUM_ALLOWED_RT_SLOPE &&
      proposed_slope <= MAXIMUM_ALLOWED_RT_SLOPE &&
      proposed_rmse <= shift_rmse
    if (allowed) {
      intercept <- proposed_intercept
      slope <- proposed_slope
      model_type <- "Constrained affine landmark correction"
    } else {
      rejected_affine <- TRUE
    }
  }

  predicted_reference <- (observed_rt - intercept) / slope
  residual_rmse <- if (length(indices)) {
    sqrt(mean((predicted_reference - reference_rt)^2))
  } else NA_real_
  rt_model_rows[[length(rt_model_rows) + 1L]] <- data.frame(
    sample_id = sample_id,
    connected_component = components[match(sample_id, sample_ids)],
    pairwise_initial_shift_minutes = starting_shift,
    landmark_count = length(indices),
    landmark_rt_span_minutes = if (length(indices) > 1L) diff(range(reference_rt)) else 0,
    correction_model = model_type,
    intercept_minutes = intercept,
    slope = slope,
    landmark_rmse_minutes = residual_rmse,
    rejected_affine_candidate = rejected_affine,
    stringsAsFactors = FALSE
  )
}
rt_models <- do.call(rbind, rt_model_rows)
model_position <- match(detections$sample_id, rt_models$sample_id)
detections$rt_model_intercept <- rt_models$intercept_minutes[model_position]
detections$rt_model_slope <- rt_models$slope[model_position]
detections$aligned_retention_time <-
  (detections$raw_retention_time - detections$rt_model_intercept) /
  detections$rt_model_slope
detections$rt_correction_minutes <-
  detections$aligned_retention_time - detections$raw_retention_time

# --------------------------- FINAL CORE GROUPS -----------------------------
message("Building complete-linkage RT-spectral peak groups")
final_core_groups <- complete_linkage_groups(
  core_indices, detections$aligned_retention_time,
  core_spectral_distances, FINAL_CORE_RT_TOLERANCE,
  FINAL_CORE_SPECTRAL_TOLERANCE, detections$sample_id
)
detections$alignment_group <- NA_integer_
detections$alignment_group[as.integer(names(final_core_groups))] <-
  final_core_groups
detections$assignment_status <- ifelse(
  detections$core_detection, "Core complete-linkage assignment", "Not yet assigned"
)
detections$assignment_score <- ifelse(detections$core_detection, 0, NA_real_)
detections$assignment_margin <- ifelse(detections$core_detection, Inf, NA_real_)
detections$competing_groups <- ""

group_consensus <- function(group_ids) {
  groups <- sort(unique(group_ids[is.finite(group_ids)]))
  spectra <- t(vapply(groups, function(group_id) {
    members <- which(detections$alignment_group == group_id)
    normalize_spectrum(apply(shape_spectra[members, , drop = FALSE], 2,
                             stats::median))
  }, numeric(ncol(shape_spectra))))
  rownames(spectra) <- groups
  list(
    groups = groups,
    retention_times = setNames(vapply(groups, function(group_id) {
      members <- which(detections$alignment_group == group_id)
      stats::median(detections$aligned_retention_time[
        members
      ])
    }, numeric(1)), groups),
    spectra = spectra,
    samples = setNames(lapply(groups, function(group_id) {
      members <- which(detections$alignment_group == group_id)
      detections$sample_id[members]
    }), groups)
  )
}

# Strict complete linkage can fragment a real feature when one noisy member is
# near the tolerance boundary. Agglomerate only groups with disjoint sample
# membership, nearly identical consensus spectra, and a still-compact combined
# RT span. If two groups occur in the same sample, they cannot be merged.
number_of_consensus_merges <- 0L
repeat {
  current_consensus <- group_consensus(detections$alignment_group)
  if (length(current_consensus$groups) < 2L) break
  consensus_distance <- spectral_distance_matrix(current_consensus$spectra)
  rt_distance <- abs(outer(
    current_consensus$retention_times,
    current_consensus$retention_times, "-"
  ))
  candidate_pairs <- which(
    upper.tri(consensus_distance) &
      rt_distance <= CONSENSUS_MERGE_RT_TOLERANCE &
      consensus_distance <= CONSENSUS_MERGE_SPECTRAL_TOLERANCE,
    arr.ind = TRUE
  )
  if (!nrow(candidate_pairs)) break
  eligible_rows <- vapply(seq_len(nrow(candidate_pairs)), function(row) {
    first_position <- candidate_pairs[row, 1]
    second_position <- candidate_pairs[row, 2]
    first_group <- current_consensus$groups[first_position]
    second_group <- current_consensus$groups[second_position]
    disjoint <- !length(intersect(
      current_consensus$samples[[first_position]],
      current_consensus$samples[[second_position]]
    ))
    members <- which(detections$alignment_group %in%
                       c(first_group, second_group))
    combined_span <- diff(range(detections$aligned_retention_time[members]))
    disjoint && combined_span <= CONSENSUS_MERGE_MAXIMUM_RT_SPAN
  }, logical(1))
  candidate_pairs <- candidate_pairs[eligible_rows, , drop = FALSE]
  if (!nrow(candidate_pairs)) break
  candidate_score <- pmax(
    rt_distance[candidate_pairs] / CONSENSUS_MERGE_RT_TOLERANCE,
    consensus_distance[candidate_pairs] /
      CONSENSUS_MERGE_SPECTRAL_TOLERANCE
  )
  selected <- candidate_pairs[which.min(candidate_score), ]
  first_group <- current_consensus$groups[selected[1]]
  second_group <- current_consensus$groups[selected[2]]
  first_size <- sum(detections$alignment_group == first_group, na.rm = TRUE)
  second_size <- sum(detections$alignment_group == second_group, na.rm = TRUE)
  retained_group <- if (first_size >= second_size) first_group else second_group
  removed_group <- if (first_size >= second_size) second_group else first_group
  detections$alignment_group[
    !is.na(detections$alignment_group) &
      detections$alignment_group == removed_group
  ] <- retained_group
  number_of_consensus_merges <- number_of_consensus_merges + 1L
}
message("  agglomerated ", number_of_consensus_merges,
        " disjoint, consensus-equivalent core groups")

# ---------------- ITERATIVE RARE-FRAGMENT RECONCILIATION ------------------
# Rebuild the target consensuses after every accepted round and continue until
# a complete round accepts no additional rescue. Three safeguards prevent
# snowball merging: only groups that were recurrent anchors before iteration
# may receive detections, every detection may move at most once, and every move
# must pass all evidence gates using the consensus from that exact round.
detections$pre_reconciliation_group <- detections$alignment_group
initial_reconciliation_consensus <- group_consensus(detections$alignment_group)
initial_reconciliation_group_sizes <- setNames(vapply(
  initial_reconciliation_consensus$groups, function(group_id) {
    members <- which(
      is.finite(detections$alignment_group) &
        detections$alignment_group == group_id
    )
    length(unique(detections$sample_id[members]))
  }, integer(1)), initial_reconciliation_consensus$groups
)
original_anchor_groups <- as.integer(names(initial_reconciliation_group_sizes)[
  initial_reconciliation_group_sizes >= RARE_RECONCILIATION_MIN_ANCHOR_SAMPLES
])
fragment_groups <- as.integer(names(initial_reconciliation_group_sizes)[
  initial_reconciliation_group_sizes <= RARE_RECONCILIATION_MAX_GROUP_SIZE
])
fragment_indices <- which(
  detections$core_detection & detections$alignment_group %in% fragment_groups
)
reconciliation_blocks <- split(
  seq_along(analysis_wavelengths),
  floor((analysis_wavelengths - min(analysis_wavelengths)) /
          BOOTSTRAP_BLOCK_WIDTH_NM)
)
iterative_reconciliation_rows <- list()
moved_in_iteration <- setNames(rep(NA_integer_, length(fragment_indices)),
                               fragment_indices)
reconciliation_iteration <- 0L
rescues_by_iteration <- integer()
maximum_convergence_iterations <- length(fragment_indices) + 1L

repeat {
  reconciliation_iteration <- reconciliation_iteration + 1L
  if (reconciliation_iteration > maximum_convergence_iterations) {
    stop("Iterative rare-fragment rescue exceeded its convergence bound.")
  }
  reconciliation_consensus <- group_consensus(detections$alignment_group)
  reconciliation_group_sizes <- setNames(vapply(
    reconciliation_consensus$groups, function(group_id) {
      members <- which(
        is.finite(detections$alignment_group) &
          detections$alignment_group == group_id
      )
      length(unique(detections$sample_id[members]))
    }, integer(1)), reconciliation_consensus$groups
  )
  anchor_positions <- which(
    reconciliation_consensus$groups %in% original_anchor_groups &
      reconciliation_group_sizes[
        as.character(reconciliation_consensus$groups)
      ] >= RARE_RECONCILIATION_MIN_ANCHOR_SAMPLES
  )
  remaining_indices <- fragment_indices[
    is.na(moved_in_iteration[as.character(fragment_indices)])
  ]
  round_rows <- list()

  for (index in remaining_indices) {
    current_source_group <- detections$alignment_group[index]
    source_group_size <- reconciliation_group_sizes[
      as.character(current_source_group)
    ]
    original_source_group_size <- initial_reconciliation_group_sizes[
      as.character(detections$pre_reconciliation_group[index])
    ]
    eligible_positions <- anchor_positions[vapply(
      anchor_positions, function(position) {
        target_group <- reconciliation_consensus$groups[position]
        sample_available <- !(detections$sample_id[index] %in%
          reconciliation_consensus$samples[[position]])
        target_is_different <- target_group != current_source_group
        target_is_better_supported <- reconciliation_group_sizes[
          as.character(target_group)
        ] > source_group_size
        rt_close <- abs(
          detections$aligned_retention_time[index] -
            reconciliation_consensus$retention_times[position]
        ) <= RARE_RECONCILIATION_RT_WINDOW
        sample_available && target_is_different &&
          target_is_better_supported && rt_close
      }, logical(1)
    )]

    if (!length(eligible_positions)) {
      round_rows[[length(round_rows) + 1L]] <- data.frame(
        rescue_iteration = reconciliation_iteration,
        detection_id = detections$detection_id[index],
        sample_id = detections$sample_id[index],
        source_group = detections$pre_reconciliation_group[index],
        current_source_group = current_source_group,
        source_group_samples = original_source_group_size,
        current_source_group_samples = source_group_size,
        proposed_target_group = NA_integer_,
        target_group_samples = NA_integer_,
        target_was_original_anchor = NA,
        rt_difference_minutes = NA_real_, spectral_distance = NA_real_,
        joint_score = NA_real_, best_vs_second_margin = NA_real_,
        wavelength_bootstrap_support = NA_real_,
        reconciliation_decision =
          "Retained rare: no eligible original anchor in expanded RT window",
        stringsAsFactors = FALSE
      )
      next
    }

    candidate_groups <- reconciliation_consensus$groups[eligible_positions]
    candidate_spectra <- reconciliation_consensus$spectra[
      eligible_positions, , drop = FALSE
    ]
    spectral <- spectral_distance_to_references(
      shape_spectra[index, ], candidate_spectra
    )
    rt_difference <- abs(
      detections$aligned_retention_time[index] -
        reconciliation_consensus$retention_times[eligible_positions]
    )
    score <- 0.85 * spectral / RARE_RECONCILIATION_SPECTRAL_TOLERANCE +
      0.15 * rt_difference / RARE_RECONCILIATION_RT_WINDOW
    score_order <- order(score)
    best <- score_order[1]
    margin <- if (length(score_order) > 1L) {
      score[score_order[2]] - score[best]
    } else Inf

    bootstrap_wins <- logical(BOOTSTRAP_REPETITIONS)
    for (bootstrap_iteration in seq_len(BOOTSTRAP_REPETITIONS)) {
      sampled_blocks <- sample(
        names(reconciliation_blocks), length(reconciliation_blocks),
        replace = TRUE
      )
      columns <- unlist(reconciliation_blocks[sampled_blocks], use.names = FALSE)
      query <- shape_spectra[index, columns]
      references <- candidate_spectra[, columns, drop = FALSE]
      query_l2 <- query / max(sqrt(sum(query^2)), sqrt(.Machine$double.eps))
      reference_l2 <- references /
        pmax(sqrt(rowSums(references^2)), sqrt(.Machine$double.eps))
      resampled_spectral <- pmax(
        0, 1 - as.numeric(reference_l2 %*% query_l2)
      )
      resampled_score <-
        0.85 * resampled_spectral / RARE_RECONCILIATION_SPECTRAL_TOLERANCE +
        0.15 * rt_difference / RARE_RECONCILIATION_RT_WINDOW
      bootstrap_wins[bootstrap_iteration] <- which.min(resampled_score) == best
    }
    bootstrap_support <- mean(bootstrap_wins)
    full_evidence_passes <-
      spectral[best] <= RARE_RECONCILIATION_SPECTRAL_TOLERANCE &&
      score[best] <= 1 &&
      (!is.finite(margin) || margin >= RARE_RECONCILIATION_SCORE_MARGIN) &&
      bootstrap_support >= RARE_RECONCILIATION_BOOTSTRAP_SUPPORT
    decision <- if (full_evidence_passes) {
      "Eligible for recurrent-feature rescue"
    } else if (spectral[best] > RARE_RECONCILIATION_SPECTRAL_TOLERANCE) {
      "Retained rare: spectrum differs from recurrent feature"
    } else if (is.finite(margin) &&
               margin < RARE_RECONCILIATION_SCORE_MARGIN) {
      "Retained rare: competing recurrent features"
    } else if (bootstrap_support < RARE_RECONCILIATION_BOOTSTRAP_SUPPORT) {
      "Retained rare: wavelength bootstrap unstable"
    } else {
      "Retained rare: joint RT-spectral score outside rescue limit"
    }
    round_rows[[length(round_rows) + 1L]] <- data.frame(
      rescue_iteration = reconciliation_iteration,
      detection_id = detections$detection_id[index],
      sample_id = detections$sample_id[index],
      source_group = detections$pre_reconciliation_group[index],
      current_source_group = current_source_group,
      source_group_samples = original_source_group_size,
      current_source_group_samples = source_group_size,
      proposed_target_group = candidate_groups[best],
      target_group_samples = reconciliation_group_sizes[
        as.character(candidate_groups[best])
      ],
      target_was_original_anchor = candidate_groups[best] %in%
        original_anchor_groups,
      rt_difference_minutes = rt_difference[best],
      spectral_distance = spectral[best], joint_score = score[best],
      best_vs_second_margin = margin,
      wavelength_bootstrap_support = bootstrap_support,
      reconciliation_decision = decision,
      stringsAsFactors = FALSE
    )
  }

  round_audit <- if (length(round_rows)) do.call(rbind, round_rows) else NULL
  accepted_this_round <- 0L
  if (!is.null(round_audit)) {
    proposals <- which(round_audit$reconciliation_decision ==
                         "Eligible for recurrent-feature rescue")
    proposals <- proposals[order(
      round_audit$joint_score[proposals],
      -round_audit$wavelength_bootstrap_support[proposals]
    )]
    for (audit_row in proposals) {
      index <- match(round_audit$detection_id[audit_row],
                     detections$detection_id)
      target_group <- round_audit$proposed_target_group[audit_row]
      target_samples <- detections$sample_id[
        is.finite(detections$alignment_group) &
          detections$alignment_group == target_group
      ]
      if (detections$sample_id[index] %in% target_samples) {
        round_audit$reconciliation_decision[audit_row] <-
          "Retained rare: target already contains this sample"
        next
      }
      detections$alignment_group[index] <- target_group
      detections$assignment_status[index] <-
        "Iterative rare-fragment RT-spectral rescue"
      detections$assignment_score[index] <- round_audit$joint_score[audit_row]
      detections$assignment_margin[index] <-
        round_audit$best_vs_second_margin[audit_row]
      moved_in_iteration[as.character(index)] <- reconciliation_iteration
      round_audit$reconciliation_decision[audit_row] <-
        "Rescued into recurrent feature"
      accepted_this_round <- accepted_this_round + 1L
    }
    iterative_reconciliation_rows[[length(iterative_reconciliation_rows) + 1L]] <-
      round_audit
  }
  rescues_by_iteration <- c(rescues_by_iteration, accepted_this_round)
  if (accepted_this_round == 0L) break
}

rare_fragment_reconciliation_iteration_audit <- if (
  length(iterative_reconciliation_rows)
) {
  do.call(rbind, iterative_reconciliation_rows)
} else data.frame(
  rescue_iteration = integer(), detection_id = character(),
  sample_id = character(), source_group = integer(),
  current_source_group = integer(), source_group_samples = integer(),
  current_source_group_samples = integer(),
  proposed_target_group = integer(), target_group_samples = integer(),
  target_was_original_anchor = logical(), rt_difference_minutes = numeric(),
  spectral_distance = numeric(), joint_score = numeric(),
  best_vs_second_margin = numeric(), wavelength_bootstrap_support = numeric(),
  reconciliation_decision = character()
)

# Preserve the original one-row-per-detection contract for routine use. A
# rescued row is authoritative; otherwise the final convergence-round result is
# retained. The full history remains in the iteration audit.
final_reconciliation_rows <- lapply(fragment_indices, function(index) {
  attempts <- rare_fragment_reconciliation_iteration_audit[
    rare_fragment_reconciliation_iteration_audit$detection_id ==
      detections$detection_id[index], , drop = FALSE
  ]
  rescued_row <- which(attempts$reconciliation_decision ==
                         "Rescued into recurrent feature")
  selected <- if (length(rescued_row)) {
    attempts[rescued_row[1], , drop = FALSE]
  } else attempts[nrow(attempts), , drop = FALSE]
  selected$reconciliation_iterations_evaluated <- nrow(attempts)
  selected$rescue_iteration_accepted <- if (length(rescued_row)) {
    selected$rescue_iteration
  } else NA_integer_
  selected$convergence_iteration <- reconciliation_iteration
  selected$converged_without_additional_rescue <-
    tail(rescues_by_iteration, 1) == 0L
  selected$target_group_samples_after_convergence <- if (
    is.finite(selected$proposed_target_group)
  ) {
    length(unique(detections$sample_id[
      is.finite(detections$alignment_group) &
        detections$alignment_group == selected$proposed_target_group
    ]))
  } else NA_integer_
  selected
})
fragment_reconciliation <- if (length(final_reconciliation_rows)) {
  do.call(rbind, final_reconciliation_rows)
} else rare_fragment_reconciliation_iteration_audit

rescued_detections <- sum(is.finite(moved_in_iteration))
remaining_small_groups <- sum(vapply(
  sort(unique(detections$alignment_group[is.finite(detections$alignment_group)])),
  function(group_id) {
    members <- which(
      is.finite(detections$alignment_group) &
        detections$alignment_group == group_id
    )
    length(unique(detections$sample_id[members])) <=
      RARE_RECONCILIATION_MAX_GROUP_SIZE
  }, logical(1)
))
message("  iterative reconciliation converged after ",
        reconciliation_iteration, " rounds (rescues by round: ",
        paste(rescues_by_iteration, collapse = ", "), "); ",
        rescued_detections,
        " displaced detections entered original recurrent anchors; ",
        remaining_small_groups,
        " candidate groups remain for explicit low-prevalence review")

# Lower-quality peaks may join a supported core only if the best joint score is
# sufficiently separated from the alternative. Unmatched peaks are not turned
# into low-quality singleton features.
noncore_indices <- which(!detections$core_detection)
ambiguous_rows <- list()
# Weak detections are scored against the fixed, high-quality core consensuses.
# Rebuilding every median spectrum after each weak assignment is both
# statistically undesirable (weak spectra could move the target) and expensive
# for genus-scale datasets. Only the per-group sample occupancy is updated.
assignment_consensus <- group_consensus(detections$alignment_group)
for (index in noncore_indices[order(detections$quality_score[noncore_indices],
                                    decreasing = TRUE)]) {
  if (!length(assignment_consensus$groups)) break
  eligible_group <- vapply(seq_along(assignment_consensus$groups),
                           function(group_position) {
    group_rt <- unname(assignment_consensus$retention_times[group_position])
    if (length(group_rt) != 1L || !is.finite(group_rt)) return(FALSE)
    sample_available <- !(detections$sample_id[index] %in%
                            assignment_consensus$samples[[group_position]])
    isTRUE(sample_available &&
             abs(detections$aligned_retention_time[index] - group_rt) <=
               LOW_QUALITY_RT_TOLERANCE)
  }, logical(1))
  if (!any(eligible_group)) {
    detections$assignment_status[index] <- "Low-quality unmatched detection"
    next
  }
  candidate_groups <- assignment_consensus$groups[eligible_group]
  candidate_spectra <- assignment_consensus$spectra[
    match(candidate_groups, assignment_consensus$groups), , drop = FALSE
  ]
  spectral_distance <- spectral_distance_to_references(
    shape_spectra[index, ], candidate_spectra
  )
  candidate_positions <- match(candidate_groups, assignment_consensus$groups)
  rt_distance <- abs(detections$aligned_retention_time[index] -
    assignment_consensus$retention_times[candidate_positions])
  score <- 0.65 * spectral_distance / LOW_QUALITY_SPECTRAL_TOLERANCE +
    0.35 * rt_distance / LOW_QUALITY_RT_TOLERANCE
  acceptable <- spectral_distance <= LOW_QUALITY_SPECTRAL_TOLERANCE & score <= 1
  if (!any(acceptable)) {
    detections$assignment_status[index] <- "Low-quality spectrally unmatched detection"
    next
  }
  candidate_groups <- candidate_groups[acceptable]
  score <- score[acceptable]
  spectral_distance <- spectral_distance[acceptable]
  order_score <- order(score)
  best <- order_score[1]
  margin <- if (length(order_score) > 1L) {
    score[order_score[2]] - score[best]
  } else Inf
  detections$assignment_score[index] <- score[best]
  detections$assignment_margin[index] <- margin
  detections$competing_groups[index] <- paste(candidate_groups, collapse = ";")
  if (is.finite(margin) && margin < MINIMUM_LOW_QUALITY_SCORE_MARGIN) {
    detections$assignment_status[index] <- "Ambiguous low-quality detection"
    ambiguous_rows[[length(ambiguous_rows) + 1L]] <- data.frame(
      detection_id = detections$detection_id[index],
      sample_id = detections$sample_id[index],
      aligned_retention_time = detections$aligned_retention_time[index],
      candidate_groups = paste(candidate_groups, collapse = ";"),
      best_score = score[best],
      score_margin = margin,
      stringsAsFactors = FALSE
    )
  } else {
    detections$alignment_group[index] <- candidate_groups[best]
    detections$assignment_status[index] <- "Quality-weighted consensus assignment"
    group_position <- match(candidate_groups[best], assignment_consensus$groups)
    assignment_consensus$samples[[group_position]] <- c(
      assignment_consensus$samples[[group_position]], detections$sample_id[index]
    )
  }
}

ambiguous_detections <- if (length(ambiguous_rows)) {
  do.call(rbind, ambiguous_rows)
} else data.frame(
  detection_id = character(), sample_id = character(),
  aligned_retention_time = numeric(), candidate_groups = character(),
  best_score = numeric(), score_margin = numeric()
)

# ------------------- STRONG SPLIT-FEATURE RECONCILIATION ------------------
# The earlier agglomeration requires exceptionally tight consensus agreement.
# Here we evaluate every remaining mutually exclusive neighboring pair with the
# same composite spectral distance used elsewhere plus wavelength-block
# support. This is still conservative: evidence of a same-sample double peak
# makes a pair categorically ineligible for automatic merging.
split_reconciliation_blocks <- split(
  seq_along(analysis_wavelengths),
  floor((analysis_wavelengths - min(analysis_wavelengths)) /
          BOOTSTRAP_BLOCK_WIDTH_NM)
)
split_feature_merge_rows <- list()
number_of_split_feature_merges <- 0L
repeat {
  split_consensus <- group_consensus(detections$alignment_group)
  if (length(split_consensus$groups) < 2L) break
  split_spectral <- spectral_distance_matrix(split_consensus$spectra)
  split_rt <- abs(outer(
    split_consensus$retention_times,
    split_consensus$retention_times, "-"
  ))
  candidate_pairs <- which(
    upper.tri(split_spectral) &
      split_rt <= COMPETITION_FAMILY_RT_WINDOW &
      split_spectral <= COMPETITION_FAMILY_SPECTRAL_TOLERANCE,
    arr.ind = TRUE
  )
  if (!nrow(candidate_pairs)) break

  pair_rows <- lapply(seq_len(nrow(candidate_pairs)), function(row) {
    first_position <- candidate_pairs[row, 1]
    second_position <- candidate_pairs[row, 2]
    first_group <- split_consensus$groups[first_position]
    second_group <- split_consensus$groups[second_position]
    first_samples <- split_consensus$samples[[first_position]]
    second_samples <- split_consensus$samples[[second_position]]
    overlap <- intersect(first_samples, second_samples)
    members <- which(detections$alignment_group %in%
                       c(first_group, second_group))
    combined_span <- diff(range(detections$aligned_retention_time[members]))

    bootstrap_equivalent <- logical(BOOTSTRAP_REPETITIONS)
    for (iteration in seq_len(BOOTSTRAP_REPETITIONS)) {
      sampled_blocks <- sample(
        names(split_reconciliation_blocks),
        length(split_reconciliation_blocks), replace = TRUE
      )
      columns <- unlist(
        split_reconciliation_blocks[sampled_blocks], use.names = FALSE
      )
      bootstrap_distance <- spectral_distance_to_references(
        split_consensus$spectra[first_position, columns],
        matrix(split_consensus$spectra[second_position, columns], nrow = 1)
      )
      bootstrap_equivalent[iteration] <-
        bootstrap_distance <= SPLIT_FEATURE_MERGE_SPECTRAL_TOLERANCE
    }
    bootstrap_support <- mean(bootstrap_equivalent)
    automatic_merge <- !length(overlap) &&
      split_rt[first_position, second_position] <=
        SPLIT_FEATURE_MERGE_RT_TOLERANCE &&
      split_spectral[first_position, second_position] <=
        SPLIT_FEATURE_MERGE_SPECTRAL_TOLERANCE &&
      combined_span <= SPLIT_FEATURE_MERGE_MAXIMUM_RT_SPAN &&
      bootstrap_support >= SPLIT_FEATURE_MERGE_BOOTSTRAP_SUPPORT
    decision <- if (length(overlap)) {
      "Retained separately: same-sample double peak observed"
    } else if (automatic_merge) {
      "Merged: mutually exclusive probable split homologue"
    } else {
      "Retained separately: evidence inconclusive"
    }
    data.frame(
      first_group = first_group,
      second_group = second_group,
      first_consensus_retention_time =
        split_consensus$retention_times[first_position],
      second_consensus_retention_time =
        split_consensus$retention_times[second_position],
      first_group_samples = length(unique(first_samples)),
      second_group_samples = length(unique(second_samples)),
      same_sample_doublet_count = length(overlap),
      same_sample_doublet_samples = paste(overlap, collapse = ";"),
      consensus_rt_difference_minutes =
        split_rt[first_position, second_position],
      combined_rt_span_minutes = combined_span,
      consensus_spectral_distance =
        split_spectral[first_position, second_position],
      wavelength_bootstrap_equivalence_support = bootstrap_support,
      pair_decision = decision,
      automatic_merge = automatic_merge,
      retained_group = NA_integer_,
      merged_group = NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  pair_audit <- do.call(rbind, pair_rows)
  merge_candidates <- which(pair_audit$automatic_merge)
  if (!length(merge_candidates)) {
    split_feature_merge_rows[[length(split_feature_merge_rows) + 1L]] <-
      pair_audit
    break
  }

  merge_score <- pmax(
    pair_audit$consensus_rt_difference_minutes[merge_candidates] /
      SPLIT_FEATURE_MERGE_RT_TOLERANCE,
    pair_audit$consensus_spectral_distance[merge_candidates] /
      SPLIT_FEATURE_MERGE_SPECTRAL_TOLERANCE,
    pair_audit$combined_rt_span_minutes[merge_candidates] /
      SPLIT_FEATURE_MERGE_MAXIMUM_RT_SPAN
  )
  selected_row <- merge_candidates[which.min(merge_score)]
  selected <- pair_audit[selected_row, , drop = FALSE]
  first_size <- selected$first_group_samples
  second_size <- selected$second_group_samples
  retained_group <- if (first_size >= second_size) {
    selected$first_group
  } else selected$second_group
  removed_group <- if (first_size >= second_size) {
    selected$second_group
  } else selected$first_group
  selected$retained_group <- retained_group
  selected$merged_group <- removed_group
  split_feature_merge_rows[[length(split_feature_merge_rows) + 1L]] <- selected
  moved_indices <- which(
    is.finite(detections$alignment_group) &
      detections$alignment_group == removed_group
  )
  detections$alignment_group[moved_indices] <- retained_group
  detections$assignment_status[moved_indices] <-
    "Mutually exclusive split-feature merge"
  number_of_split_feature_merges <- number_of_split_feature_merges + 1L
}

split_feature_merge_audit <- if (length(split_feature_merge_rows)) {
  output <- do.call(rbind, split_feature_merge_rows)
  output[!duplicated(paste(
    pmin(output$first_group, output$second_group),
    pmax(output$first_group, output$second_group), sep = "\r"
  )), , drop = FALSE]
} else data.frame(
  first_group = integer(), second_group = integer(),
  first_consensus_retention_time = numeric(),
  second_consensus_retention_time = numeric(),
  first_group_samples = integer(), second_group_samples = integer(),
  same_sample_doublet_count = integer(),
  same_sample_doublet_samples = character(),
  consensus_rt_difference_minutes = numeric(),
  combined_rt_span_minutes = numeric(),
  consensus_spectral_distance = numeric(),
  wavelength_bootstrap_equivalence_support = numeric(),
  pair_decision = character(), automatic_merge = logical(),
  retained_group = integer(), merged_group = integer()
)
message(
  "  merged ", number_of_split_feature_merges,
  " mutually exclusive, RT-close, spectrally equivalent split-feature pairs"
)

# Redirect every pre-merge group reference so later ambiguity and rescue audits
# still point to a current feature. This matters when several accepted merges
# form a short chain.
redirect_split_group <- function(group_id) {
  current <- as.integer(group_id)
  visited <- integer()
  repeat {
    row <- which(
      is.finite(split_feature_merge_audit$merged_group) &
        split_feature_merge_audit$merged_group == current
    )
    if (!length(row) || current %in% visited) return(current)
    visited <- c(visited, current)
    current <- split_feature_merge_audit$retained_group[row[1]]
  }
}
if (nrow(ambiguous_detections)) {
  ambiguous_detections$candidate_groups <- vapply(
    strsplit(ambiguous_detections$candidate_groups, ";", fixed = TRUE),
    function(groups) paste(unique(vapply(
      as.integer(groups), redirect_split_group, integer(1)
    )), collapse = ";"),
    character(1)
  )
}
if (nrow(fragment_reconciliation)) {
  finite_target <- is.finite(fragment_reconciliation$proposed_target_group)
  fragment_reconciliation$proposed_target_group[finite_target] <- vapply(
    fragment_reconciliation$proposed_target_group[finite_target],
    redirect_split_group, integer(1)
  )
}
if (nrow(rare_fragment_reconciliation_iteration_audit)) {
  finite_target <- is.finite(
    rare_fragment_reconciliation_iteration_audit$proposed_target_group
  )
  rare_fragment_reconciliation_iteration_audit$proposed_target_group[
    finite_target
  ] <- vapply(
    rare_fragment_reconciliation_iteration_audit$proposed_target_group[
      finite_target
    ],
    redirect_split_group, integer(1)
  )
}

# --------------------- STABLE FEATURE IDENTIFIERS -------------------------
assigned_groups <- sort(unique(detections$alignment_group[
  is.finite(detections$alignment_group)
]))
group_rt <- vapply(assigned_groups, function(group_id) {
  members <- which(detections$alignment_group == group_id)
  stats::median(detections$aligned_retention_time[
    members
  ])
}, numeric(1))
group_order <- order(group_rt)
assigned_groups <- assigned_groups[group_order]
group_rt <- group_rt[group_order]
rounded_rt <- sprintf("%.3f", group_rt)
within_rt_rank <- ave(seq_along(group_rt), rounded_rt, FUN = seq_along)
feature_ids <- mapply(format_feature_id, group_rt, within_rt_rank,
                      USE.NAMES = FALSE)
feature_by_group <- setNames(feature_ids, assigned_groups)
detections$aligned_feature_id <- feature_by_group[
  as.character(detections$alignment_group)
]

# ------------------- ASSIGNMENT BOOTSTRAP SUPPORT -------------------------
message("Bootstrapping wavelength-block assignment support")
consensus <- group_consensus(detections$alignment_group)
consensus_feature_ids <- feature_by_group[as.character(consensus$groups)]
rownames(consensus$spectra) <- consensus_feature_ids
names(consensus$retention_times) <- consensus_feature_ids

# Final rare-feature review answers a different question from rescue: after all
# accepted moves, does a remaining dataset-eligible feature look distinct from
# every larger, better-supported feature in an expanded RT neighborhood, or
# should it be treated as a possible unresolved alignment fragment? The
# complete matrix retains both, while this label supports conservative filtering.
final_group_sizes <- setNames(vapply(consensus$groups, function(group_id) {
  members <- which(
    is.finite(detections$alignment_group) &
      detections$alignment_group == group_id
  )
  length(unique(detections$sample_id[members]))
}, integer(1)), consensus$groups)
group_review_rows <- lapply(seq_along(consensus$groups), function(position) {
  group_id <- consensus$groups[position]
  group_size <- final_group_sizes[as.character(group_id)]
  if (group_size > RARE_RECONCILIATION_MAX_GROUP_SIZE) {
    return(data.frame(
      alignment_group = group_id,
      rare_feature_review_status = "Recurrent feature",
      nearest_recurrent_feature = NA_character_,
      nearest_recurrent_rt_difference_minutes = NA_real_,
      nearest_recurrent_spectral_distance = NA_real_,
      downstream_use_recommendation = "Include in primary abundance analysis",
      stringsAsFactors = FALSE
    ))
  }
  better_supported_positions <- which(
    final_group_sizes[as.character(consensus$groups)] >=
      max(RARE_RECONCILIATION_MIN_ANCHOR_SAMPLES, group_size + 1L)
  )
  candidate_positions <- better_supported_positions[vapply(
    better_supported_positions, function(candidate_position) {
      no_sample_collision <- !length(intersect(
        consensus$samples[[position]],
        consensus$samples[[candidate_position]]
      ))
      rt_close <- abs(
        consensus$retention_times[position] -
          consensus$retention_times[candidate_position]
      ) <= RARE_REVIEW_RT_WINDOW
      no_sample_collision && rt_close
    }, logical(1)
  )]
  if (!length(candidate_positions)) {
    return(data.frame(
      alignment_group = group_id,
      rare_feature_review_status =
        "Rare; no recurrent feature in expanded RT window",
      nearest_recurrent_feature = NA_character_,
      nearest_recurrent_rt_difference_minutes = NA_real_,
      nearest_recurrent_spectral_distance = NA_real_,
      downstream_use_recommendation =
        "Retain as rare candidate; evaluate in sensitivity analysis",
      stringsAsFactors = FALSE
    ))
  }
  distances <- spectral_distance_to_references(
    consensus$spectra[position, ],
    consensus$spectra[candidate_positions, , drop = FALSE]
  )
  best <- which.min(distances)
  best_position <- candidate_positions[best]
  best_distance <- distances[best]
  status <- if (best_distance >= RARE_DISTINCT_SPECTRAL_DISTANCE) {
    "Rare; spectrally distinct from nearby recurrent features"
  } else {
    "Provisional rare; possible unresolved alignment fragment"
  }
  recommendation <- if (best_distance >= RARE_DISTINCT_SPECTRAL_DISTANCE) {
    "Retain as rare candidate; evaluate in sensitivity analysis"
  } else {
    "Exclude from primary models; retain in complete matrix for review"
  }
  data.frame(
    alignment_group = group_id,
    rare_feature_review_status = status,
    nearest_recurrent_feature = consensus_feature_ids[best_position],
    nearest_recurrent_rt_difference_minutes = abs(
      consensus$retention_times[position] -
        consensus$retention_times[best_position]
    ),
    nearest_recurrent_spectral_distance = best_distance,
    downstream_use_recommendation = recommendation,
    stringsAsFactors = FALSE
  )
})
group_review <- do.call(rbind, group_review_rows)

if (nrow(fragment_reconciliation)) {
  fragment_reconciliation$proposed_target_feature <- feature_by_group[
    as.character(fragment_reconciliation$proposed_target_group)
  ]
  fragment_reconciliation$final_aligned_feature <- detections$aligned_feature_id[
    match(fragment_reconciliation$detection_id, detections$detection_id)
  ]
  fragment_reconciliation$source_group_retained_as_feature <-
    fragment_reconciliation$source_group %in% assigned_groups
}
if (nrow(rare_fragment_reconciliation_iteration_audit)) {
  rare_fragment_reconciliation_iteration_audit$proposed_target_feature <-
    feature_by_group[as.character(
      rare_fragment_reconciliation_iteration_audit$proposed_target_group
    )]
  rare_fragment_reconciliation_iteration_audit$final_aligned_feature <-
    detections$aligned_feature_id[match(
      rare_fragment_reconciliation_iteration_audit$detection_id,
      detections$detection_id
    )]
  rare_fragment_reconciliation_iteration_audit$
    source_group_retained_as_feature <-
      rare_fragment_reconciliation_iteration_audit$source_group %in%
        assigned_groups
}

bootstrap_blocks <- split(
  seq_along(analysis_wavelengths),
  floor((analysis_wavelengths - min(analysis_wavelengths)) /
          BOOTSTRAP_BLOCK_WIDTH_NM)
)
detections$assignment_bootstrap_support <- NA_real_
detections$nearest_competing_feature <- NA_character_
detections$full_data_assignment_margin <- NA_real_

assigned_indices <- which(is.finite(detections$alignment_group))
for (index in assigned_indices) {
  current_feature <- detections$aligned_feature_id[index]
  nearby_features <- names(consensus$retention_times)[
    abs(consensus$retention_times - detections$aligned_retention_time[index]) <=
      max(0.10, RARE_RECONCILIATION_RT_WINDOW)
  ]
  # A successfully reconciled detection may sit just outside the ordinary
  # 0.10-minute bootstrap neighborhood. Always retain its assigned feature in
  # the candidate set so support is evaluated rather than producing an empty
  # current-feature score.
  candidate_features <- unique(c(current_feature, nearby_features))
  candidate_positions <- match(candidate_features, rownames(consensus$spectra))
  candidate_features <- candidate_features[!is.na(candidate_positions)]
  candidate_positions <- candidate_positions[!is.na(candidate_positions)]
  if (length(candidate_features) < 2L) next
  candidate_spectra <- consensus$spectra[candidate_positions, , drop = FALSE]
  rownames(candidate_spectra) <- candidate_features
  full_spectral <- spectral_distance_to_references(
    shape_spectra[index, ], candidate_spectra
  )
  full_score <- 0.70 * full_spectral / LOW_QUALITY_SPECTRAL_TOLERANCE +
    0.30 * abs(consensus$retention_times[candidate_features] -
                 detections$aligned_retention_time[index]) /
    LOW_QUALITY_RT_TOLERANCE
  score_order <- order(full_score)
  alternatives <- candidate_features[candidate_features != current_feature]
  if (length(alternatives)) {
    alternative_scores <- full_score[candidate_features != current_feature]
    nearest_position <- which.min(alternative_scores)
    detections$nearest_competing_feature[index] <- alternatives[nearest_position]
    detections$full_data_assignment_margin[index] <-
      min(alternative_scores) - full_score[candidate_features == current_feature]
  }

  wins <- logical(BOOTSTRAP_REPETITIONS)
  for (iteration in seq_len(BOOTSTRAP_REPETITIONS)) {
    sampled_blocks <- sample(names(bootstrap_blocks), length(bootstrap_blocks),
                             replace = TRUE)
    columns <- unlist(bootstrap_blocks[sampled_blocks], use.names = FALSE)
    query <- shape_spectra[index, columns]
    references <- candidate_spectra[, columns, drop = FALSE]
    query_l2 <- query / max(sqrt(sum(query^2)), sqrt(.Machine$double.eps))
    reference_l2 <- references /
      pmax(sqrt(rowSums(references^2)), sqrt(.Machine$double.eps))
    spectral <- pmax(0, 1 - as.numeric(reference_l2 %*% query_l2))
    score <- 0.70 * spectral / LOW_QUALITY_SPECTRAL_TOLERANCE +
      0.30 * abs(consensus$retention_times[candidate_features] -
                   detections$aligned_retention_time[index]) /
      LOW_QUALITY_RT_TOLERANCE
    wins[iteration] <- candidate_features[which.min(score)] == current_feature
  }
  detections$assignment_bootstrap_support[index] <- mean(wins)
}

# ---------------- UNRESOLVED NEIGHBORING-FEATURE FAMILIES -----------------
# Summarize instability at the feature level, then distinguish three cases:
#   1. same-sample co-occurrence -> protect as a resolved double peak;
#   2. mutually exclusive and close/similar -> unresolved sensitivity family;
#   3. no qualifying neighbor -> retain separately with its individual warning.
feature_bootstrap_support <- setNames(vapply(
  feature_ids, function(feature_id) {
    values <- detections$assignment_bootstrap_support[
      detections$aligned_feature_id == feature_id
    ]
    values <- values[is.finite(values)]
    if (length(values)) stats::median(values) else NA_real_
  }, numeric(1)
), feature_ids)

final_pair_rows <- list()
if (length(feature_ids) >= 2L) {
  final_spectral_distances <- spectral_distance_matrix(consensus$spectra)
  rownames(final_spectral_distances) <- feature_ids
  colnames(final_spectral_distances) <- feature_ids
  for (first_position in seq_len(length(feature_ids) - 1L)) {
    for (second_position in seq.int(first_position + 1L, length(feature_ids))) {
      first_feature <- feature_ids[first_position]
      second_feature <- feature_ids[second_position]
      rt_difference <- abs(
        consensus$retention_times[first_feature] -
          consensus$retention_times[second_feature]
      )
      if (rt_difference > COMPETITION_FAMILY_RT_WINDOW) next
      spectral_distance <- final_spectral_distances[
        first_feature, second_feature
      ]
      if (spectral_distance > COMPETITION_FAMILY_SPECTRAL_TOLERANCE) next
      first_samples <- unique(detections$sample_id[
        which(detections$aligned_feature_id == first_feature)
      ])
      second_samples <- unique(detections$sample_id[
        which(detections$aligned_feature_id == second_feature)
      ])
      overlap <- intersect(first_samples, second_samples)
      first_unstable <- is.finite(feature_bootstrap_support[first_feature]) &&
        feature_bootstrap_support[first_feature] <
          COMPETITION_ASSIGNMENT_SUPPORT_THRESHOLD
      second_unstable <- is.finite(feature_bootstrap_support[second_feature]) &&
        feature_bootstrap_support[second_feature] <
          COMPETITION_ASSIGNMENT_SUPPORT_THRESHOLD
      decision <- if (length(overlap)) {
        "Protected resolved double peak: same-sample co-occurrence"
      } else if (first_unstable || second_unstable) {
        "Unresolved neighboring-feature family"
      } else {
        "Retained separately: assignments stable"
      }
      final_pair_rows[[length(final_pair_rows) + 1L]] <- data.frame(
        first_feature = first_feature,
        second_feature = second_feature,
        first_feature_samples = length(first_samples),
        second_feature_samples = length(second_samples),
        same_sample_doublet_count = length(overlap),
        same_sample_doublet_samples = paste(overlap, collapse = ";"),
        consensus_rt_difference_minutes = rt_difference,
        consensus_spectral_distance = spectral_distance,
        first_feature_assignment_support =
          feature_bootstrap_support[first_feature],
        second_feature_assignment_support =
          feature_bootstrap_support[second_feature],
        pair_decision = decision,
        stringsAsFactors = FALSE
      )
    }
  }
}
neighboring_feature_pair_audit <- if (length(final_pair_rows)) {
  do.call(rbind, final_pair_rows)
} else data.frame(
  first_feature = character(), second_feature = character(),
  first_feature_samples = integer(), second_feature_samples = integer(),
  same_sample_doublet_count = integer(),
  same_sample_doublet_samples = character(),
  consensus_rt_difference_minutes = numeric(),
  consensus_spectral_distance = numeric(),
  first_feature_assignment_support = numeric(),
  second_feature_assignment_support = numeric(), pair_decision = character()
)

# Build complete-linkage families: every pair admitted to one family must be
# mutually exclusive, RT-close, spectrally similar, and connected to observed
# assignment instability. This prevents a transitive chain from collapsing a
# genuine same-sample doublet at its ends.
family_candidate_features <- unique(c(
  neighboring_feature_pair_audit$first_feature[
    neighboring_feature_pair_audit$pair_decision ==
      "Unresolved neighboring-feature family"
  ],
  neighboring_feature_pair_audit$second_feature[
    neighboring_feature_pair_audit$pair_decision ==
      "Unresolved neighboring-feature family"
  ]
))
competition_family_members <- list()
if (length(family_candidate_features)) {
  family_candidate_features <- family_candidate_features[order(
    consensus$retention_times[family_candidate_features]
  )]
  eligible_pair_key <- paste(
    pmin(
      neighboring_feature_pair_audit$first_feature,
      neighboring_feature_pair_audit$second_feature
    ),
    pmax(
      neighboring_feature_pair_audit$first_feature,
      neighboring_feature_pair_audit$second_feature
    ), sep = "\r"
  )[neighboring_feature_pair_audit$pair_decision ==
      "Unresolved neighboring-feature family"]
  for (feature_id in family_candidate_features) {
    placed <- FALSE
    if (length(competition_family_members)) {
      for (family_position in seq_along(competition_family_members)) {
        member_ids <- competition_family_members[[family_position]]
        keys <- paste(
          pmin(feature_id, member_ids), pmax(feature_id, member_ids), sep = "\r"
        )
        if (all(keys %in% eligible_pair_key)) {
          competition_family_members[[family_position]] <- c(
            member_ids, feature_id
          )
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) competition_family_members[[length(competition_family_members) + 1L]] <-
      feature_id
  }
  competition_family_members <- competition_family_members[
    lengths(competition_family_members) >= 2L
  ]
}

competition_family_by_feature <- setNames(rep(NA_character_, length(feature_ids)),
                                          feature_ids)
competition_family_summary_rows <- list()
if (length(competition_family_members)) {
  family_order <- order(vapply(competition_family_members, function(members) {
    stats::median(consensus$retention_times[members])
  }, numeric(1)))
  competition_family_members <- competition_family_members[family_order]
  for (family_position in seq_along(competition_family_members)) {
    members <- competition_family_members[[family_position]]
    family_id <- paste0("UF", sprintf("%03d", family_position))
    competition_family_by_feature[members] <- family_id
    competition_family_summary_rows[[family_position]] <- data.frame(
      competition_family_id = family_id,
      member_features = paste(members, collapse = ";"),
      number_of_member_features = length(members),
      family_consensus_retention_time = stats::median(
        consensus$retention_times[members]
      ),
      minimum_member_assignment_support = min(
        feature_bootstrap_support[members], na.rm = TRUE
      ),
      maximum_pair_rt_difference_minutes = max(abs(outer(
        consensus$retention_times[members],
        consensus$retention_times[members], "-"
      ))),
      family_interpretation = paste(
        "Members remain separate in the complete matrix, are excluded from",
        "the strict primary matrix, and are summed only for sensitivity",
        "analysis because their sample assignments are not reliably separable."
      ),
      stringsAsFactors = FALSE
    )
  }
}
competition_family_summary <- if (length(competition_family_summary_rows)) {
  do.call(rbind, competition_family_summary_rows)
} else data.frame(
  competition_family_id = character(), member_features = character(),
  number_of_member_features = integer(),
  family_consensus_retention_time = numeric(),
  minimum_member_assignment_support = numeric(),
  maximum_pair_rt_difference_minutes = numeric(),
  family_interpretation = character()
)

# ---------------- FEATURE SUMMARIES AND REPRESENTATIVE SPECTRA ------------
feature_summary_rows <- list()
representative_rows <- list()
uncertainty_rows <- list()
for (feature_position in seq_along(feature_ids)) {
  feature_id <- feature_ids[feature_position]
  group_id <- assigned_groups[feature_position]
  members <- which(detections$alignment_group == group_id)
  member_shapes <- shape_spectra[members, , drop = FALSE]
  member_integrated <- integrated_spectra[members, , drop = FALSE]
  representative_shape <- normalize_spectrum(apply(
    member_shapes, 2, stats::median
  ))
  # Scale the consensus shape to the strongest well-observed member. The fixed
  # library model is shape-normalized later; retaining an observed magnitude
  # keeps the compatibility pseudochromatogram useful without inventing values.
  strongest_member <- members[which.max(detections$quality_score[members])]
  reference_scale <- max(integrated_spectra[strongest_member, ], na.rm = TRUE)
  representative_spectrum <- representative_shape * reference_scale
  number_detected <- length(unique(detections$sample_id[members]))
  prevalence <- number_detected / length(sample_ids)
  bootstrap_values <- detections$assignment_bootstrap_support[members]
  bootstrap_values <- bootstrap_values[is.finite(bootstrap_values)]
  support <- if (length(bootstrap_values)) stats::median(bootstrap_values) else NA_real_
  review <- group_review[
    match(group_id, group_review$alignment_group), , drop = FALSE
  ]
  rescued_member_count <- sum(
    grepl(
      "rare-fragment RT-spectral rescue",
      detections$assignment_status[members], ignore.case = TRUE
    )
  )
  recurrence <- if (number_detected == 1L) {
    "High-quality single-sample feature"
  } else if (prevalence < 0.10) {
    "Rare recurrent feature"
  } else if (prevalence < 0.25) {
    "Moderately recurrent feature"
  } else "Widely recurrent feature"
  alignment_confidence <- if (
    grepl("^Provisional rare", review$rare_feature_review_status)
  ) {
    "Rare; possible alignment fragment"
  } else if (review$rare_feature_review_status != "Recurrent feature") {
    "Rare; limited recurrence support"
  } else if (is.finite(support) && support < 0.70) {
    "Neighboring-feature assignment unstable"
  } else if (diff(range(detections$aligned_retention_time[members])) > 0.08) {
    "RT-dispersed alignment"
  } else "Supported alignment"

  observed_areas <- detections$integrated_dad_area[members]
  feature_summary_rows[[feature_position]] <- data.frame(
    aligned_feature_id = feature_id,
    alignment_group = group_id,
    consensus_retention_time = stats::median(
      detections$aligned_retention_time[members]
    ),
    raw_retention_time_median = stats::median(
      detections$raw_retention_time[members]
    ),
    aligned_retention_time_iqr = stats::IQR(
      detections$aligned_retention_time[members]
    ),
    aligned_retention_time_minimum = min(
      detections$aligned_retention_time[members]
    ),
    aligned_retention_time_maximum = max(
      detections$aligned_retention_time[members]
    ),
    detected_samples = number_detected,
    total_samples = length(sample_ids),
    detection_fraction = prevalence,
    recurrence_status = recurrence,
    alignment_confidence = alignment_confidence,
    rare_feature_review_status = review$rare_feature_review_status,
    nearest_recurrent_feature = review$nearest_recurrent_feature,
    nearest_recurrent_rt_difference_minutes =
      review$nearest_recurrent_rt_difference_minutes,
    nearest_recurrent_spectral_distance =
      review$nearest_recurrent_spectral_distance,
    downstream_use_recommendation = review$downstream_use_recommendation,
    unresolved_competition_family = unname(
      competition_family_by_feature[feature_id]
    ),
    neighboring_feature_decision = if (
      is.finite(match(feature_id, family_candidate_features)) &&
        !is.na(competition_family_by_feature[feature_id])
    ) {
      "Exclude from strict primary matrix; use summed family in sensitivity analysis"
    } else if (is.finite(support) &&
               support < COMPETITION_ASSIGNMENT_SUPPORT_THRESHOLD) {
      "Retain separately; no mutually exclusive close-similar family established"
    } else "No unresolved neighboring-feature family",
    rescued_detection_count = rescued_member_count,
    median_assignment_bootstrap_support = support,
    median_full_data_assignment_margin = stats::median(
      detections$full_data_assignment_margin[members], na.rm = TRUE
    ),
    median_integrated_dad_area_detected = stats::median(observed_areas),
    mean_integrated_dad_area_all_samples = sum(observed_areas) / length(sample_ids),
    integrated_dad_area_iqr_detected = stats::IQR(observed_areas),
    median_fwhm_minutes = stats::median(detections$fwhm_minutes[members]),
    representative_source_sample = detections$sample_id[strongest_member],
    representative_source_detection = detections$detection_id[strongest_member],
    representative_quality_score = detections$quality_score[strongest_member],
    stringsAsFactors = FALSE
  )
  representative_rows[[feature_position]] <- data.frame(
    aligned_feature_id = feature_id,
    t(representative_spectrum), check.names = FALSE
  )
  names(representative_rows[[feature_position]])[-1] <- analysis_wavelengths
  uncertainty_rows[[feature_position]] <- data.frame(
    aligned_feature_id = feature_id,
    wavelength_nm = analysis_wavelengths,
    median_normalized_absorbance = apply(member_shapes, 2, stats::median),
    q10_normalized_absorbance = apply(member_shapes, 2, stats::quantile,
                                      probs = 0.10),
    q25_normalized_absorbance = apply(member_shapes, 2, stats::quantile,
                                      probs = 0.25),
    q75_normalized_absorbance = apply(member_shapes, 2, stats::quantile,
                                      probs = 0.75),
    q90_normalized_absorbance = apply(member_shapes, 2, stats::quantile,
                                      probs = 0.90),
    number_of_spectra = nrow(member_shapes),
    stringsAsFactors = FALSE
  )
}

feature_summary <- do.call(rbind, feature_summary_rows)
feature_summary$median_full_data_assignment_margin[
  !is.finite(feature_summary$median_full_data_assignment_margin)
] <- NA_real_
feature_summary$include_in_conservative_alignment_matrix <- !grepl(
  "^Provisional rare", feature_summary$rare_feature_review_status
)
# The clearer name is used in new outputs. Keep the historical column as an
# exact alias so older analysis code remains reproducible.
feature_summary$include_in_primary_filtered_matrix <-
  feature_summary$include_in_conservative_alignment_matrix
feature_summary$include_in_strict_primary_matrix <-
  feature_summary$include_in_primary_filtered_matrix &
  is.na(feature_summary$unresolved_competition_family)
reference_spectra <- do.call(rbind, representative_rows)
representative_uncertainty <- do.call(rbind, uncertainty_rows)

# ---------------------------- ABUNDANCE MATRIX -----------------------------
# The primary matrix contains confidently assigned signal only and is always
# complete: zero means that no signal was confidently assigned to that feature
# in that sample. An unresolved local detection is not hidden inside an NA.
# Instead, two companion matrices retain its meaning explicitly:
#   * uncertainty_flag_matrix: 1 for every plausible destination feature;
#   * upper_bound_matrix: the largest per-feature value if all unresolved
#     candidate signal for that cell belonged to that feature.
# Upper bounds are feature-wise sensitivity limits. They must not be summed
# across candidate features because one unresolved peak appears in every
# plausible destination column.
abundance_matrix <- matrix(
  0, nrow = length(sample_ids), ncol = length(feature_ids),
  dimnames = list(sample_ids, feature_ids)
)
for (index in assigned_indices) {
  abundance_matrix[detections$sample_id[index],
                   detections$aligned_feature_id[index]] <-
    detections$integrated_dad_area[index]
}

uncertainty_flag_matrix <- matrix(
  0L, nrow = length(sample_ids), ncol = length(feature_ids),
  dimnames = list(sample_ids, feature_ids)
)
upper_bound_matrix <- abundance_matrix

# Record every plausible feature explicitly for each ambiguous detection. The
# confident matrix stays numeric and NA-free, while the flag and upper-bound
# matrices preserve the evidence needed for review and sensitivity analyses.
if (nrow(ambiguous_detections)) {
  ambiguous_detections$integrated_dad_area <- detections$integrated_dad_area[
    match(ambiguous_detections$detection_id, detections$detection_id)
  ]
  ambiguous_detections$candidate_aligned_features <- ""
  ambiguous_detections$candidate_feature_count <- 0L
  for (row in seq_len(nrow(ambiguous_detections))) {
    groups <- as.integer(strsplit(ambiguous_detections$candidate_groups[row],
                                  ";", fixed = TRUE)[[1]])
    features <- feature_by_group[as.character(groups)]
    features <- features[is.finite(match(features, colnames(abundance_matrix)))]
    features <- unique(as.character(features))
    ambiguous_detections$candidate_aligned_features[row] <- paste(
      features, collapse = ";"
    )
    ambiguous_detections$candidate_feature_count[row] <- length(features)
    sample_index <- match(
      ambiguous_detections$sample_id[row], rownames(abundance_matrix)
    )
    feature_indices <- match(features, colnames(abundance_matrix))
    feature_indices <- feature_indices[is.finite(feature_indices)]
    unresolved_area <- ambiguous_detections$integrated_dad_area[row]
    if (length(feature_indices) && is.finite(sample_index) &&
        is.finite(unresolved_area) && unresolved_area >= 0) {
      uncertainty_flag_matrix[sample_index, feature_indices] <- 1L
      upper_bound_matrix[sample_index, feature_indices] <-
        upper_bound_matrix[sample_index, feature_indices, drop = TRUE] +
        unresolved_area
    }
  }
} else {
  ambiguous_detections$integrated_dad_area <- numeric()
  ambiguous_detections$candidate_aligned_features <- character()
  ambiguous_detections$candidate_feature_count <- integer()
}
abundance_output <- data.frame(
  sample_id = rownames(abundance_matrix), abundance_matrix,
  check.names = FALSE, row.names = NULL
)
uncertainty_flag_output <- data.frame(
  sample_id = rownames(uncertainty_flag_matrix), uncertainty_flag_matrix,
  check.names = FALSE, row.names = NULL
)
upper_bound_output <- data.frame(
  sample_id = rownames(upper_bound_matrix), upper_bound_matrix,
  check.names = FALSE, row.names = NULL
)
conservative_feature_ids <- feature_summary$aligned_feature_id[
  feature_summary$include_in_conservative_alignment_matrix
]
conservative_abundance_output <- abundance_output[, c(
  "sample_id", conservative_feature_ids
), drop = FALSE]
strict_primary_feature_ids <- feature_summary$aligned_feature_id[
  feature_summary$include_in_strict_primary_matrix
]
strict_primary_abundance_output <- abundance_output[, c(
  "sample_id", strict_primary_feature_ids
), drop = FALSE]
strict_primary_uncertainty_flag_output <- uncertainty_flag_output[, c(
  "sample_id", strict_primary_feature_ids
), drop = FALSE]
strict_primary_upper_bound_output <- upper_bound_output[, c(
  "sample_id", strict_primary_feature_ids
), drop = FALSE]

# Feature-level summaries tell users how often an otherwise analysis-ready
# feature was a plausible destination for signal that could not be assigned.
feature_summary$unresolved_candidate_samples <- colSums(
  uncertainty_flag_matrix > 0
)[feature_summary$aligned_feature_id]
feature_summary$unresolved_candidate_fraction <-
  feature_summary$unresolved_candidate_samples / length(sample_ids)

# Sensitivity matrix: stable individual features remain individual columns,
# while each unresolved competition family becomes one summed column. This
# matrix also contains confident assignments only and therefore has no NAs;
# unresolved candidate signals remain in the separate flag/upper-bound files.
sensitivity_stable_feature_ids <- feature_summary$aligned_feature_id[
  is.na(feature_summary$unresolved_competition_family)
]
sensitivity_abundance_matrix <- abundance_matrix[
  , sensitivity_stable_feature_ids, drop = FALSE
]
if (nrow(competition_family_summary)) {
  for (family_position in seq_len(nrow(competition_family_summary))) {
    family_id <- competition_family_summary$competition_family_id[family_position]
    members <- strsplit(
      competition_family_summary$member_features[family_position],
      ";", fixed = TRUE
    )[[1]]
    member_values <- abundance_matrix[, members, drop = FALSE]
    family_values <- rowSums(member_values)
    sensitivity_abundance_matrix <- cbind(
      sensitivity_abundance_matrix,
      setNames(data.frame(family_values), family_id)
    )
  }
}
sensitivity_column_rt <- c(
  setNames(
    feature_summary$consensus_retention_time[
      match(sensitivity_stable_feature_ids, feature_summary$aligned_feature_id)
    ], sensitivity_stable_feature_ids
  ),
  setNames(
    competition_family_summary$family_consensus_retention_time,
    competition_family_summary$competition_family_id
  )
)
sensitivity_abundance_matrix <- sensitivity_abundance_matrix[
  , names(sort(sensitivity_column_rt)), drop = FALSE
]
sensitivity_abundance_output <- data.frame(
  sample_id = rownames(sensitivity_abundance_matrix),
  sensitivity_abundance_matrix,
  check.names = FALSE, row.names = NULL
)
abundance_long <- matrix_to_long(
  abundance_matrix, value_name = "integrated_dad_area"
)
long_cell_index <- cbind(
  match(abundance_long$sample_id, rownames(abundance_matrix)),
  match(abundance_long$aligned_feature_id, colnames(abundance_matrix))
)
abundance_long$confident_integrated_dad_area <-
  abundance_long$integrated_dad_area
abundance_long$unresolved_candidate_flag <- as.integer(
  uncertainty_flag_matrix[long_cell_index]
)
abundance_long$featurewise_upper_bound_integrated_dad_area <-
  as.numeric(upper_bound_matrix[long_cell_index])
abundance_long$detection_status <- ifelse(
  abundance_long$unresolved_candidate_flag > 0 &
    abundance_long$integrated_dad_area == 0,
  "Detected locally but alignment unresolved",
  ifelse(abundance_long$integrated_dad_area > 0,
         "Aligned peak detected", "Not detected above thresholds")
)
feature_metadata_position <- match(
  abundance_long$aligned_feature_id, feature_summary$aligned_feature_id
)
abundance_long$alignment_confidence <-
  feature_summary$alignment_confidence[feature_metadata_position]
abundance_long$rare_feature_review_status <-
  feature_summary$rare_feature_review_status[feature_metadata_position]
abundance_long$downstream_use_recommendation <-
  feature_summary$downstream_use_recommendation[feature_metadata_position]
abundance_long$include_in_conservative_alignment_matrix <-
  feature_summary$include_in_conservative_alignment_matrix[
    feature_metadata_position
  ]
abundance_long$include_in_primary_filtered_matrix <-
  feature_summary$include_in_primary_filtered_matrix[feature_metadata_position]
abundance_long$include_in_strict_primary_matrix <-
  feature_summary$include_in_strict_primary_matrix[feature_metadata_position]
abundance_long$unresolved_competition_family <-
  feature_summary$unresolved_competition_family[feature_metadata_position]

# Peak-picker audit tables are assembled only after alignment so they can
# distinguish detected, core, assigned, unresolved, and unmatched outcomes.
sample_peak_detection_summary <- do.call(rbind, lapply(sample_ids, function(id) {
  rows <- detections$sample_id == id
  assigned <- rows & !is.na(detections$aligned_feature_id)
  unresolved <- rows & detections$assignment_status ==
    "Ambiguous low-quality detection"
  data.frame(
    sample_id = id,
    detected_peaks = sum(rows),
    core_detections = sum(rows & detections$core_detection),
    assigned_peaks = sum(assigned),
    unresolved_local_detections = sum(unresolved),
    unmatched_low_quality_detections = sum(
      rows & grepl("unmatched detection$", detections$assignment_status,
                   ignore.case = TRUE)
    ),
    assigned_fraction = if (sum(rows)) sum(assigned) / sum(rows) else NA_real_,
    median_profile_snr = stats::median(detections$profile_snr[rows]),
    median_prominence_snr = stats::median(detections$prominence_snr[rows]),
    stringsAsFactors = FALSE
  )
}))

sensitivity_setting_key <- with(
  peak_picking_sensitivity,
  paste(smoothing_minutes, minimum_profile_snr,
        minimum_prominence_snr, sep = "\r")
)
peak_picking_sensitivity_summary <- do.call(rbind, lapply(
  split(seq_len(nrow(peak_picking_sensitivity)), sensitivity_setting_key),
  function(rows) {
    values <- peak_picking_sensitivity$selected_peak_count[rows]
    data.frame(
      smoothing_minutes = peak_picking_sensitivity$smoothing_minutes[rows[1]],
      minimum_profile_snr =
        peak_picking_sensitivity$minimum_profile_snr[rows[1]],
      minimum_prominence_snr =
        peak_picking_sensitivity$minimum_prominence_snr[rows[1]],
      is_deployed_setting =
        peak_picking_sensitivity$is_deployed_setting[rows[1]],
      total_candidate_peaks = sum(values),
      total_deployed_peaks_recovered = sum(
        peak_picking_sensitivity$deployed_peaks_recovered[rows]
      ),
      total_trial_peaks_not_matching_deployed = sum(
        peak_picking_sensitivity$trial_peaks_not_matching_deployed[rows]
      ),
      median_deployed_peak_recovery_percent = 100 * stats::median(
        peak_picking_sensitivity$deployed_peak_recovery_fraction[rows],
        na.rm = TRUE
      ),
      median_trial_match_percent = 100 * stats::median(
        peak_picking_sensitivity$trial_match_fraction[rows], na.rm = TRUE
      ),
      minimum_peaks_per_sample = min(values),
      q10_peaks_per_sample = unname(stats::quantile(values, 0.10)),
      median_peaks_per_sample = stats::median(values),
      q90_peaks_per_sample = unname(stats::quantile(values, 0.90)),
      maximum_peaks_per_sample = max(values),
      stringsAsFactors = FALSE
    )
  }
))
deployed_count <- peak_picking_sensitivity_summary$total_candidate_peaks[
  peak_picking_sensitivity_summary$is_deployed_setting
]
if (length(deployed_count) != 1L) {
  stop("Peak-picking sensitivity audit did not identify one deployed setting.")
}
peak_picking_sensitivity_summary$change_from_deployed_percent <-
  100 * (peak_picking_sensitivity_summary$total_candidate_peaks -
           deployed_count) / deployed_count
peak_picking_sensitivity_summary <- peak_picking_sensitivity_summary[order(
  peak_picking_sensitivity_summary$smoothing_minutes,
  peak_picking_sensitivity_summary$minimum_profile_snr,
  peak_picking_sensitivity_summary$minimum_prominence_snr
), , drop = FALSE]

# Put policy provenance directly in every row-oriented alignment table. Wide
# quantitative matrices remain strictly sample ID plus numeric feature columns;
# their policy is recorded in the adjacent manifest instead.
add_alignment_policy <- function(table) {
  data.frame(
    alignment_policy = rep(ALIGNMENT_POLICY, nrow(table)),
    alignment_policy_display_name = rep(
      ALIGNMENT_POLICY_DISPLAY_NAME, nrow(table)
    ),
    table,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
input_audit <- add_alignment_policy(input_audit)
wavelength_audit <- add_alignment_policy(wavelength_audit)
pairwise_edges <- add_alignment_policy(pairwise_edges)
rt_models <- add_alignment_policy(rt_models)
detections <- add_alignment_policy(detections)
ambiguous_detections <- add_alignment_policy(ambiguous_detections)
fragment_reconciliation <- add_alignment_policy(fragment_reconciliation)
rare_fragment_reconciliation_iteration_audit <- add_alignment_policy(
  rare_fragment_reconciliation_iteration_audit
)
split_feature_merge_audit <- add_alignment_policy(split_feature_merge_audit)
neighboring_feature_pair_audit <- add_alignment_policy(
  neighboring_feature_pair_audit
)
competition_family_summary <- add_alignment_policy(
  competition_family_summary
)
feature_summary <- add_alignment_policy(feature_summary)
representative_uncertainty <- add_alignment_policy(representative_uncertainty)
abundance_long <- add_alignment_policy(abundance_long)
sample_peak_detection_summary <- add_alignment_policy(
  sample_peak_detection_summary
)
peak_picking_sensitivity_summary <- add_alignment_policy(
  peak_picking_sensitivity_summary
)

# ------------------------------- OUTPUTS ----------------------------------
write.csv(input_audit, file.path(OUTPUT_DIR, "raw_agilent_input_audit.csv"),
          row.names = FALSE)
write.csv(wavelength_audit,
          file.path(OUTPUT_DIR, "wavelength_grid_alignment.csv"),
          row.names = FALSE)
write.csv(
  sample_peak_detection_summary,
  file.path(OUTPUT_DIR, "sample_peak_detection_summary.csv"),
  row.names = FALSE
)
write.csv(
  peak_picking_sensitivity_summary,
  file.path(OUTPUT_DIR, "peak_picking_parameter_sensitivity.csv"),
  row.names = FALSE
)
write.csv(pairwise_edges,
          file.path(OUTPUT_DIR, "pairwise_rt_landmark_edges.csv"),
          row.names = FALSE)
write.csv(rt_models,
          file.path(OUTPUT_DIR, "sample_retention_time_corrections.csv"),
          row.names = FALSE)
write.csv(detections,
          file.path(OUTPUT_DIR, "sample_peak_detections_and_alignment.csv"),
          row.names = FALSE)
write.csv(ambiguous_detections,
          file.path(OUTPUT_DIR, "ambiguous_unforced_detections.csv"),
          row.names = FALSE)
write.csv(fragment_reconciliation,
          file.path(OUTPUT_DIR, "rare_fragment_reconciliation_audit.csv"),
          row.names = FALSE)
write.csv(
  rare_fragment_reconciliation_iteration_audit,
  file.path(
    OUTPUT_DIR, "rare_fragment_reconciliation_iteration_audit.csv"
  ),
  row.names = FALSE
)
write.csv(
  split_feature_merge_audit,
  file.path(OUTPUT_DIR, "probable_split_feature_merge_audit.csv"),
  row.names = FALSE
)
write.csv(
  neighboring_feature_pair_audit,
  file.path(OUTPUT_DIR, "neighboring_feature_pair_decisions.csv"),
  row.names = FALSE
)
write.csv(
  competition_family_summary,
  file.path(OUTPUT_DIR, "unresolved_competition_family_dictionary.csv"),
  row.names = FALSE
)
write.csv(feature_summary,
          file.path(OUTPUT_DIR, "aligned_feature_summary.csv"),
          row.names = FALSE)
write.csv(
  feature_summary[
    !feature_summary$include_in_conservative_alignment_matrix, , drop = FALSE
  ],
  file.path(OUTPUT_DIR, "provisional_rare_features_for_manual_review.csv"),
  row.names = FALSE
)
write.csv(reference_spectra, file.path(
  OUTPUT_DIR, paste0(DATASET_SLUG, "_rawdad_aligned_reference_spectra.csv")
), row.names = FALSE)
write.csv(representative_uncertainty, file.path(
  OUTPUT_DIR, "aligned_feature_representative_spectra_with_uncertainty.csv"
), row.names = FALSE)
write.csv(abundance_output, file.path(
  OUTPUT_DIR, "sample_by_aligned_feature_integrated_dad_area.csv"
), row.names = FALSE)
write.csv(uncertainty_flag_output, file.path(
  OUTPUT_DIR, "sample_by_aligned_feature_assignment_uncertainty_flag.csv"
), row.names = FALSE)
write.csv(upper_bound_output, file.path(
  OUTPUT_DIR, "sample_by_aligned_feature_integrated_dad_area_upper_bound.csv"
), row.names = FALSE)
write.csv(strict_primary_abundance_output, file.path(
  OUTPUT_DIR,
  "sample_by_aligned_feature_integrated_dad_area_strict_primary.csv"
), row.names = FALSE)
write.csv(strict_primary_uncertainty_flag_output, file.path(
  OUTPUT_DIR,
  "sample_by_aligned_feature_assignment_uncertainty_flag_strict_primary.csv"
), row.names = FALSE)
write.csv(strict_primary_upper_bound_output, file.path(
  OUTPUT_DIR,
  paste0(
    "sample_by_aligned_feature_integrated_dad_area_upper_bound_",
    "strict_primary.csv"
  )
), row.names = FALSE)
write.csv(sensitivity_abundance_output, file.path(
  OUTPUT_DIR,
  "sample_by_feature_with_summed_unresolved_families_sensitivity.csv"
), row.names = FALSE)
# Intermediate fragment-filtered matrix. The strict matrix remains the primary
# output because it also removes unresolved neighboring-feature families.
write.csv(conservative_abundance_output, file.path(
  OUTPUT_DIR,
  "sample_by_aligned_feature_integrated_dad_area_primary_filtered.csv"
), row.names = FALSE)
write.csv(abundance_long, file.path(
  OUTPUT_DIR, "sample_by_aligned_feature_integrated_dad_area_long.csv"
), row.names = FALSE)

policy_description <- if (ALIGNMENT_POLICY == "retained_uniqueness") {
  paste(
    "Reviews only candidate groups observed in 1-3 samples; preserves more",
    "low-prevalence groups as distinct unless the strict rescue evidence is",
    "met. Recommended when retaining potential biological uniqueness is the",
    "primary concern."
  )
} else {
  paste(
    "Reviews candidate groups observed in up to 20 samples, or ceiling(20%",
    "of samples) when fewer than 20 samples are available; applies the same",
    "strict rescue evidence gates. Recommended when minimizing displaced",
    "alignment fragments is the primary concern."
  )
}
policy_manifest <- data.frame(
  source_dataset = DATASET_NAME,
  alignment_policy = ALIGNMENT_POLICY,
  alignment_policy_display_name = ALIGNMENT_POLICY_DISPLAY_NAME,
  alignment_policy_description = policy_description,
  total_raw_sample_folders = length(discovered_sample_ids),
  excluded_unreadable_inputs = nrow(raw_input_exclusions),
  total_samples = length(sample_ids),
  rescue_review_maximum_source_group_samples =
    RARE_RECONCILIATION_MAX_GROUP_SIZE,
  rescue_target_minimum_samples = RARE_RECONCILIATION_MIN_ANCHOR_SAMPLES,
  rescue_requires_target_larger_than_source = TRUE,
  iterative_rescue_until_convergence = TRUE,
  rescue_targets_restricted_to_original_anchors = TRUE,
  rescue_iterations_completed = reconciliation_iteration,
  rescues_by_iteration = paste(rescues_by_iteration, collapse = ";"),
  complete_matrix =
    "sample_by_aligned_feature_integrated_dad_area.csv",
  recommended_primary_matrix =
    "sample_by_aligned_feature_integrated_dad_area_strict_primary.csv",
  recommended_primary_uncertainty_flags =
    "sample_by_aligned_feature_assignment_uncertainty_flag_strict_primary.csv",
  primary_featurewise_upper_bounds = paste0(
    "sample_by_aligned_feature_integrated_dad_area_upper_bound_",
    "strict_primary.csv"
  ),
  primary_fragment_filtered_matrix =
    "sample_by_aligned_feature_integrated_dad_area_primary_filtered.csv",
  family_sensitivity_matrix =
    "sample_by_feature_with_summed_unresolved_families_sensitivity.csv",
  representative_spectra = paste0(
    DATASET_SLUG, "_rawdad_aligned_reference_spectra.csv"
  ),
  downstream_annotation_dataset = ANNOTATION_DATASET,
  downstream_annotation_folder = ANNOTATION_OUTPUT_DIR,
  retention_time_used_for_library_annotation = FALSE,
  stringsAsFactors = FALSE
)
write.csv(
  policy_manifest,
  file.path(OUTPUT_DIR, "raw_alignment_policy_manifest.csv"),
  row.names = FALSE
)
write.csv(data.frame(
  parameter = c(
    "analysis_mode", "alignment_policy", "alignment_policy_display_name",
    "raw_input_directory", "number_of_raw_sample_folders",
    "number_of_excluded_unreadable_inputs", "number_of_samples",
    "number_of_detected_peaks", "number_of_core_detections",
    "number_of_aligned_features", "baseline_window_minutes",
    "minimum_profile_snr", "minimum_prominence_snr",
    "core_rt_tolerance_minutes", "core_spectral_distance_tolerance",
    "maximum_allowed_rt_slope", "minimum_allowed_rt_slope",
    "rare_reconciliation_max_group_size",
    "retained_uniqueness_max_group_size",
    "rare_reconciliation_absolute_group_cap",
    "rare_reconciliation_small_dataset_fraction",
    "rare_reconciliation_min_anchor_samples",
    "rare_reconciliation_rt_window_minutes",
    "rare_reconciliation_spectral_tolerance",
    "rare_reconciliation_score_margin",
    "rare_reconciliation_bootstrap_support",
    "split_feature_merge_rt_tolerance_minutes",
    "split_feature_merge_spectral_tolerance",
    "split_feature_merge_maximum_rt_span_minutes",
    "split_feature_merge_bootstrap_support",
    "competition_family_rt_window_minutes",
    "competition_family_spectral_tolerance",
    "competition_assignment_support_threshold",
    "number_of_probable_split_feature_merges",
    "number_of_unresolved_competition_families",
    "number_of_rescued_displaced_detections",
    "rare_reconciliation_iterative_until_convergence",
    "rare_reconciliation_original_anchors_only",
    "rare_reconciliation_iterations_completed",
    "rare_reconciliation_rescues_by_iteration",
    "bootstrap_repetitions", "retention_time_used_for_library_annotation"
  ),
  value = c(
    "raw_agilent_presence_absence_aware_joint_alignment",
    ALIGNMENT_POLICY, ALIGNMENT_POLICY_DISPLAY_NAME,
    normalizePath(RAW_DIR), length(discovered_sample_ids),
    nrow(raw_input_exclusions), length(sample_ids), nrow(detections),
    sum(detections$core_detection), nrow(feature_summary),
    BASELINE_WINDOW_MINUTES, MINIMUM_PROFILE_SNR, MINIMUM_PROMINENCE_SNR,
    FINAL_CORE_RT_TOLERANCE, FINAL_CORE_SPECTRAL_TOLERANCE,
    MAXIMUM_ALLOWED_RT_SLOPE, MINIMUM_ALLOWED_RT_SLOPE,
    RARE_RECONCILIATION_MAX_GROUP_SIZE,
    RARE_RECONCILIATION_RETAINED_UNIQUENESS_CAP,
    RARE_RECONCILIATION_ABSOLUTE_GROUP_CAP,
    RARE_RECONCILIATION_SMALL_DATASET_FRACTION,
    RARE_RECONCILIATION_MIN_ANCHOR_SAMPLES,
    RARE_RECONCILIATION_RT_WINDOW,
    RARE_RECONCILIATION_SPECTRAL_TOLERANCE,
    RARE_RECONCILIATION_SCORE_MARGIN,
    RARE_RECONCILIATION_BOOTSTRAP_SUPPORT,
    SPLIT_FEATURE_MERGE_RT_TOLERANCE,
    SPLIT_FEATURE_MERGE_SPECTRAL_TOLERANCE,
    SPLIT_FEATURE_MERGE_MAXIMUM_RT_SPAN,
    SPLIT_FEATURE_MERGE_BOOTSTRAP_SUPPORT,
    COMPETITION_FAMILY_RT_WINDOW,
    COMPETITION_FAMILY_SPECTRAL_TOLERANCE,
    COMPETITION_ASSIGNMENT_SUPPORT_THRESHOLD,
    number_of_split_feature_merges,
    nrow(competition_family_summary),
    rescued_detections,
    TRUE,
    TRUE,
    reconciliation_iteration,
    paste(rescues_by_iteration, collapse = ";"),
    BOOTSTRAP_REPETITIONS, FALSE
  ), stringsAsFactors = FALSE
), file.path(OUTPUT_DIR, "raw_dad_alignment_parameters.csv"), row.names = FALSE)

saveRDS(list(
  alignment_policy = ALIGNMENT_POLICY,
  alignment_policy_display_name = ALIGNMENT_POLICY_DISPLAY_NAME,
  alignment_policy_manifest = policy_manifest,
  detections = detections,
  shape_spectra = shape_spectra,
  integrated_spectra = integrated_spectra,
  feature_summary = feature_summary,
  abundance_matrix = abundance_matrix,
  primary_filtered_abundance_matrix = abundance_matrix[
    , conservative_feature_ids, drop = FALSE
  ],
  strict_primary_abundance_matrix = abundance_matrix[
    , strict_primary_feature_ids, drop = FALSE
  ],
  family_sensitivity_abundance_matrix = sensitivity_abundance_matrix,
  conservative_abundance_matrix = abundance_matrix[
    , conservative_feature_ids, drop = FALSE
  ],
  fragment_reconciliation = fragment_reconciliation,
  rare_fragment_reconciliation_iteration_audit =
    rare_fragment_reconciliation_iteration_audit,
  rare_fragment_reconciliation_iterations_completed =
    reconciliation_iteration,
  rare_fragment_reconciliation_rescues_by_iteration =
    rescues_by_iteration,
  split_feature_merge_audit = split_feature_merge_audit,
  neighboring_feature_pair_audit = neighboring_feature_pair_audit,
  competition_family_summary = competition_family_summary,
  sample_peak_detection_summary = sample_peak_detection_summary,
  peak_picking_parameter_sensitivity = peak_picking_sensitivity_summary,
  wavelengths = analysis_wavelengths
), file.path(OUTPUT_DIR, "raw_dad_alignment_objects.rds"))

# ------------------------------- FIGURES ----------------------------------
alignment_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom", legend.box = "vertical",
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot", panel.grid.major.y = element_line(
      color = "#E6E6E6", linewidth = 0.25
    )
  )

# 1. Before/after RT correction. Points are retained detections and thin lines
# connect the same peak before and after correction.
rt_shift_plot <- rbind(
  data.frame(
    detection_id = detections$detection_id,
    sample_id = detections$sample_id,
    coordinate = "Raw retention time",
    retention_time = detections$raw_retention_time,
    core_detection = detections$core_detection
  ),
  data.frame(
    detection_id = detections$detection_id,
    sample_id = detections$sample_id,
    coordinate = "Conservatively corrected",
    retention_time = detections$aligned_retention_time,
    core_detection = detections$core_detection
  )
)
rt_shift_plot$coordinate <- factor(
  rt_shift_plot$coordinate,
  levels = c("Raw retention time", "Conservatively corrected")
)
figure <- ggplot(
  rt_shift_plot,
  aes(coordinate, retention_time, group = detection_id,
      color = core_detection)
) +
  geom_line(alpha = 0.10, linewidth = 0.25) +
  geom_point(alpha = 0.35, size = 0.7) +
  scale_color_manual(
    values = c(`TRUE` = "#1B6CA8", `FALSE` = "#999999"),
    labels = c(`TRUE` = "Landmark-eligible detection",
               `FALSE` = "Lower-quality detection"),
    name = NULL
  ) +
  labs(
    title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX, "RT corrections"),
    subtitle = "Large global stretching is prohibited; lines show each detected peak",
    x = NULL, y = "Retention time (minutes)"
  ) + alignment_theme
ggsave(file.path(OUTPUT_DIR, "retention_time_correction_audit.pdf"),
       figure, width = 10, height = 7)

# 2. Feature prevalence and assignment support. Lack of a nearby competitor is
# shown separately from a numerically high bootstrap value.
support_plot <- feature_summary[
  is.finite(feature_summary$median_assignment_bootstrap_support), , drop = FALSE
]
support_plot$bootstrap_display <- support_plot$median_assignment_bootstrap_support
figure <- ggplot(
  support_plot,
  aes(detection_fraction, bootstrap_display,
      size = median_integrated_dad_area_detected,
      color = alignment_confidence)
) +
  geom_point(alpha = 0.80) +
  scale_size_continuous(range = c(1.8, 8), name = "Median detected\nDAD area") +
  scale_color_manual(values = c(
    "Supported alignment" = "#1B9E77",
    "Neighboring-feature assignment unstable" = "#D95F02",
    "RT-dispersed alignment" = "#7570B3",
    "Rare; limited recurrence support" = "#4E79A7",
    "Rare; possible alignment fragment" = "#B07AA1"
  ), labels = c(
    "Supported alignment" = "Supported alignment",
    "Neighboring-feature assignment unstable" =
      "Neighboring-feature assignment unstable",
    "RT-dispersed alignment" = "Retention-time-dispersed alignment",
    "Rare; limited recurrence support" = "Rare; limited recurrence support",
    "Rare; possible alignment fragment" = "Rare; possible alignment fragment"
  ), name = "Alignment interpretation") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                  "feature recurrence and assignment support"),
    subtitle = "Missing y-values indicate no nearby competing feature, not failed support",
    x = "Samples containing the aligned feature",
    y = "Median wavelength-block assignment support"
  ) +
  guides(
    size = guide_legend(order = 1L, title.position = "top", title.hjust = 0.5),
    color = guide_legend(
      order = 2L, nrow = 2L, byrow = TRUE,
      title.position = "top", title.hjust = 0.5,
      override.aes = list(size = 3)
    )
  ) + alignment_theme
ggsave(file.path(OUTPUT_DIR, "feature_prevalence_and_alignment_support.pdf"),
       figure, width = 11, height = 8.5)

# 3. Peak-picker optimization audit. Facets show how many maxima survive small
# threshold changes at each smoothing width. The deployed setting is marked in
# red; a broad, smooth response around it indicates a stable operating region.
picker_plot <- peak_picking_sensitivity_summary
picker_plot$smoothing_label <- factor(
  paste0("Smoothing = ", picker_plot$smoothing_minutes, " min"),
  levels = paste0(
    "Smoothing = ", sort(unique(picker_plot$smoothing_minutes)), " min"
  )
)
picker_plot$profile_snr_label <- factor(
  picker_plot$minimum_profile_snr,
  levels = sort(unique(picker_plot$minimum_profile_snr))
)
picker_plot$prominence_snr_label <- factor(
  picker_plot$minimum_prominence_snr,
  levels = sort(unique(picker_plot$minimum_prominence_snr))
)
figure <- ggplot(
  picker_plot,
  aes(profile_snr_label, prominence_snr_label,
      fill = total_candidate_peaks)
) +
  geom_tile(color = "#FFFFFF", linewidth = 0.7) +
  geom_tile(
    data = picker_plot[picker_plot$is_deployed_setting, , drop = FALSE],
    fill = NA, color = "#D7301F", linewidth = 1.25
  ) +
  geom_text(
    aes(label = total_candidate_peaks), size = 3.1, fontface = "bold"
  ) +
  facet_wrap(~ smoothing_label, ncol = 2) +
  scale_fill_viridis_c(option = "C", direction = -1, guide = "none") +
  labs(
    title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                  "peak-picking parameter sensitivity"),
    subtitle = paste(
      "Numbers are total candidate maxima across all samples; red border marks",
      "the deployed setting"
    ),
    caption = paste(
      "This audit changes one peak-picking operating point at a time; it does",
      "not select thresholds by maximizing the number of peaks."
    ),
    x = "Minimum chromatographic profile SNR",
    y = "Minimum local-prominence SNR"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot",
    plot.caption = element_text(hjust = 0)
  )
ggsave(
  file.path(OUTPUT_DIR, "peak_picking_parameter_sensitivity.pdf"),
  figure, width = 11, height = 9
)

# Sample-level picker QC: samples are ordered by total detections and the
# stacked bars show whether those maxima became landmark-eligible detections,
# were assigned conservatively, or remained unresolved/unmatched.
picker_sample_plot <- sample_peak_detection_summary
picker_sample_plot$noncore_assigned <- pmax(
  picker_sample_plot$assigned_peaks - picker_sample_plot$core_detections, 0
)
picker_sample_plot$other_unassigned <- pmax(
  picker_sample_plot$detected_peaks - picker_sample_plot$core_detections -
    picker_sample_plot$noncore_assigned -
    picker_sample_plot$unresolved_local_detections -
    picker_sample_plot$unmatched_low_quality_detections,
  0
)
sample_status_long <- rbind(
  data.frame(sample_id = picker_sample_plot$sample_id,
             status = "Landmark-eligible detections",
             count = picker_sample_plot$core_detections),
  data.frame(sample_id = picker_sample_plot$sample_id,
             status = "Additional assigned detections",
             count = picker_sample_plot$noncore_assigned),
  data.frame(sample_id = picker_sample_plot$sample_id,
             status = "Unresolved local detections",
             count = picker_sample_plot$unresolved_local_detections),
  data.frame(sample_id = picker_sample_plot$sample_id,
             status = "Unmatched low-quality detections",
             count = picker_sample_plot$unmatched_low_quality_detections),
  data.frame(sample_id = picker_sample_plot$sample_id,
             status = "Other unassigned detections",
             count = picker_sample_plot$other_unassigned)
)
sample_order <- picker_sample_plot$sample_id[order(
  picker_sample_plot$detected_peaks, picker_sample_plot$sample_id
)]
sample_status_long$sample_id <- factor(
  sample_status_long$sample_id, levels = sample_order
)
sample_status_long$status <- factor(sample_status_long$status, levels = c(
  "Landmark-eligible detections", "Additional assigned detections",
  "Unresolved local detections", "Unmatched low-quality detections",
  "Other unassigned detections"
))
figure <- ggplot(sample_status_long, aes(sample_id, count, fill = status)) +
  geom_col(width = 0.86) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Landmark-eligible detections" = "#1B6CA8",
    "Additional assigned detections" = "#67A9CF",
    "Unresolved local detections" = "#F28E2B",
    "Unmatched low-quality detections" = "#B07AA1",
    "Other unassigned detections" = "#BDBDBD"
  ), name = "Detection outcome") +
  labs(
    title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                  "sample-level peak-picking audit"),
    subtitle = "Every detected maximum is counted once and separated by downstream alignment outcome",
    x = "Sample (ordered by detected peak count)",
    y = "Detected chromatographic peaks"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom", legend.box = "vertical",
    axis.text.y = element_text(size = 6.5),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot"
  ) +
  guides(fill = guide_legend(nrow = 2L, byrow = TRUE))
ggsave(
  file.path(OUTPUT_DIR, "sample_peak_picking_audit.pdf"),
  figure, width = 12,
  height = max(10, min(18, 3 + 0.13 * length(sample_ids)))
)

# The peak-profile atlas checks the deployed detector in every sample, rather
# than assuming that three representative chromatograms capture all unusual
# baselines or weak runs. Six samples per page preserve readable peak shapes,
# sample labels, thresholds, and marker types without creating a 91-page PDF.
profile_samples_per_page <- 6L
profile_sample_order <- sort(unique(as.character(peak_profiles$sample_id)))
profile_page_count <- ceiling(
  length(profile_sample_order) / profile_samples_per_page
)
# Divide sequential samples as evenly as possible across the required number
# of pages. This avoids a nearly empty last page when the sample count is one
# more than a multiple of the nominal six-panel capacity.
profile_page_assignment <- floor(
  (seq_along(profile_sample_order) - 1) * profile_page_count /
    length(profile_sample_order)
) + 1L
profile_manifest <- merge(
  data.frame(
    sample_id = profile_sample_order,
    sample_order = seq_along(profile_sample_order),
    atlas_page = profile_page_assignment,
    stringsAsFactors = FALSE
  ),
  sample_peak_detection_summary,
  by = "sample_id", all.x = TRUE, sort = FALSE
)
profile_manifest <- profile_manifest[match(
  profile_sample_order, profile_manifest$sample_id
), , drop = FALSE]
profile_manifest$total_atlas_pages <- profile_page_count
write.csv(
  profile_manifest,
  file.path(OUTPUT_DIR, "peak_picking_profile_atlas_manifest.csv"),
  row.names = FALSE
)

# A log1p display makes the height gate and small peaks visible without changing
# which maxima were detected. The dashed line is only the global profile-height
# component of the rule. A maximum must also have local-prominence SNR >= 2.5;
# consequently, a shoulder may lie above the line but correctly lack a marker.
# Exact untransformed signals remain in the raw `.D` files.
profile_pdf <- file.path(
  OUTPUT_DIR, "representative_peak_picking_profiles.pdf"
)
grDevices::pdf(profile_pdf, width = 12, height = 10, onefile = TRUE)
for (profile_page in seq_len(profile_page_count)) {
  page_samples <- profile_sample_order[
    profile_page_assignment == profile_page
  ]
  profile_plot <- peak_profiles[
    peak_profiles$sample_id %in% page_samples, , drop = FALSE
  ]
  profile_detection_points <- detections[
    detections$sample_id %in% page_samples, c(
      "sample_id", "raw_retention_time", "profile_height", "core_detection"
    ), drop = FALSE
  ]
  page_summary <- sample_peak_detection_summary[match(
    page_samples, sample_peak_detection_summary$sample_id
  ), , drop = FALSE]
  page_labels <- paste0(
    page_samples, " | detected ", page_summary$detected_peaks,
    " | landmark-eligible ", page_summary$core_detections
  )
  names(page_labels) <- page_samples
  profile_plot$sample_panel <- factor(
    page_labels[profile_plot$sample_id], levels = page_labels
  )
  profile_detection_points$sample_panel <- factor(
    page_labels[profile_detection_points$sample_id], levels = page_labels
  )
  profile_plot$display_signal <- log1p(profile_plot$smoothed_rms_profile)
  profile_plot$display_threshold <- log1p(
    profile_plot$profile_detection_height
  )
  profile_detection_points$display_height <- log1p(
    profile_detection_points$profile_height
  )
  figure <- ggplot(
    profile_plot, aes(retention_time, display_signal)
  ) +
    geom_line(color = "#111111", linewidth = 0.55) +
    geom_line(aes(y = display_threshold), color = "#D7301F",
              linetype = "dashed", linewidth = 0.45) +
    geom_point(
      data = profile_detection_points,
      aes(raw_retention_time, display_height, shape = core_detection),
      color = "#2C7FB8", size = 1.45
    ) +
    facet_wrap(~ sample_panel, ncol = 2, scales = "free_y") +
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      labels = c(`TRUE` = "Landmark-eligible detection",
                 `FALSE` = "Additional detection"),
      name = "Detected maximum"
    ) +
    labs(
      title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                    "all-sample peak-picking profiles"),
      subtitle = paste0(
        "Page ", profile_page, " of ", profile_page_count,
        " — red line: profile-height gate (SNR 4) — ",
        "y-axis: log1p display scaling"
      ),
      x = "Retention time (minutes)",
      y = "log1p baseline-corrected RMS wavelength signal",
      caption = paste(
        "Dots pass both profile SNR >= 4 and local-prominence SNR >= 2.5.",
        "A peak-like shoulder above the red line can lack a dot if it fails",
        "local prominence or minimum-separation filtering."
      )
    ) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 8.5),
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot",
      plot.caption = element_text(hjust = 0, size = 7.5),
      panel.spacing = grid::unit(0.75, "lines")
    )
  print(figure)
}
grDevices::dev.off()

# 4. Sample-by-feature detection-status map. This is a quality-control matrix:
# its colors describe the evidence/status of each cell, not feature identity or
# chemical class. Two blue cells in different columns are therefore not the
# same feature; cells in the same column refer to the same aligned feature.
detection_map <- abundance_long
detection_map$display_status <- detection_map$detection_status
is_detected <- detection_map$detection_status == "Aligned peak detected"
detection_map$display_status[is_detected] <- "Aligned recurrent detection"
detection_map$display_status[
  is_detected & detection_map$alignment_confidence ==
    "Rare; limited recurrence support"
] <- "Reviewed rare detection"
detection_map$display_status[
  is_detected & detection_map$alignment_confidence ==
    "Rare; possible alignment fragment"
] <- "Provisional rare detection"
detection_map$status <- factor(
  detection_map$display_status,
  levels = c(
    "Not detected above thresholds", "Aligned recurrent detection",
    "Reviewed rare detection", "Provisional rare detection",
    "Detected locally but alignment unresolved"
  )
)
detection_map$aligned_feature_id <- factor(
  detection_map$aligned_feature_id, levels = feature_summary$aligned_feature_id
)
detection_map$sample_id <- factor(detection_map$sample_id,
                                  levels = rev(sample_ids))
figure <- ggplot(detection_map,
                 aes(aligned_feature_id, sample_id, fill = status)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = c(
    "Not detected above thresholds" = "#F3F3F3",
    "Aligned recurrent detection" = "#2C7FB8",
    "Reviewed rare detection" = "#6A51A3",
    "Provisional rare detection" = "#C51B7D",
    "Detected locally but alignment unresolved" = "#F28E2B"
  ), name = "Cell meaning") +
  labs(
    title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                  "aligned-feature detection and review-status matrix"),
    subtitle = paste(
      "Rows are samples and columns are aligned features; color describes",
      "cell status, not feature identity or chemical class"
    ),
    x = "Aligned feature (retention-time order)", y = "Sample"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(), legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    plot.title = element_text(face = "bold"), plot.title.position = "plot"
  ) +
  guides(fill = guide_legend(nrow = 2L, byrow = TRUE))
ggsave(file.path(OUTPUT_DIR, "sample_by_aligned_feature_detection_map.pdf"),
       figure, width = 14,
       height = max(10, min(18, 3 + 0.13 * length(sample_ids))))

# 5. A deliberately simple aligned-peak map. Every assigned detection is moved
# to its feature's consensus RT for display, so observations of the same
# homologous feature form an exact vertical column. Color is a secondary
# feature-identity cue. Unlike the status matrix above, this figure contains
# assigned peaks only; unresolved local detections remain available in the
# status matrix and the detailed detection audit.
feature_count <- nrow(feature_summary)
feature_hues <- (seq_len(feature_count) * 137.508) %% 360
feature_luminance <- rep(c(45, 62, 52, 69), length.out = feature_count)
feature_chroma <- rep(c(72, 62, 68), length.out = feature_count)
feature_identity_colors <- grDevices::hcl(
  h = feature_hues, c = feature_chroma, l = feature_luminance,
  fixup = TRUE
)
names(feature_identity_colors) <- feature_summary$aligned_feature_id

feature_color_dictionary <- data.frame(
  alignment_policy = ALIGNMENT_POLICY,
  alignment_policy_display_name = ALIGNMENT_POLICY_DISPLAY_NAME,
  aligned_feature_id = feature_summary$aligned_feature_id,
  consensus_retention_time = feature_summary$consensus_retention_time,
  figure_color_hex = unname(feature_identity_colors),
  stringsAsFactors = FALSE
)
write.csv(
  feature_color_dictionary,
  file.path(OUTPUT_DIR, "aligned_feature_figure_color_dictionary.csv"),
  row.names = FALSE
)

assigned_peak_map <- detections[
  !is.na(detections$aligned_feature_id) &
    nzchar(detections$aligned_feature_id), , drop = FALSE
]
assigned_feature_position <- match(
  assigned_peak_map$aligned_feature_id, feature_summary$aligned_feature_id
)
assigned_peak_map$feature_consensus_retention_time <-
  feature_summary$consensus_retention_time[assigned_feature_position]
assigned_peak_map$observed_minus_consensus_minutes <-
  assigned_peak_map$aligned_retention_time -
  assigned_peak_map$feature_consensus_retention_time
assigned_peak_map$feature_alignment_confidence <-
  feature_summary$alignment_confidence[assigned_feature_position]
assigned_peak_map$feature_figure_color_hex <-
  feature_identity_colors[assigned_peak_map$aligned_feature_id]
assigned_peak_map <- assigned_peak_map[, c(
  "alignment_policy", "alignment_policy_display_name", "sample_id",
  "detection_id", "aligned_feature_id", "raw_retention_time",
  "aligned_retention_time", "feature_consensus_retention_time",
  "observed_minus_consensus_minutes", "assignment_status",
  "assignment_bootstrap_support", "feature_alignment_confidence",
  "feature_figure_color_hex"
)]
write.csv(
  assigned_peak_map,
  file.path(OUTPUT_DIR, "sample_by_aligned_peak_map_data.csv"),
  row.names = FALSE
)

assigned_peak_map$sample_id <- factor(
  assigned_peak_map$sample_id, levels = rev(sample_ids)
)
assigned_peak_map$aligned_feature_id <- factor(
  assigned_peak_map$aligned_feature_id,
  levels = feature_summary$aligned_feature_id
)
feature_guides <- feature_summary
figure <- ggplot() +
  geom_vline(
    data = feature_guides,
    aes(xintercept = consensus_retention_time),
    color = "#E2E2E2", linewidth = 0.20
  ) +
  geom_point(
    data = assigned_peak_map,
    aes(
      x = feature_consensus_retention_time, y = sample_id,
      color = aligned_feature_id
    ),
    size = 1.55, alpha = 0.92
  ) +
  scale_color_manual(values = feature_identity_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.012, 0.02))) +
  labs(
    title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                  "aligned peaks by sample"),
    subtitle = paste(
      "Each point is one assigned peak; the same aligned feature shares a",
      "color and forms one vertical consensus-RT column"
    ),
    x = "Aligned-feature consensus retention time (minutes)",
    y = "Sample"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#F0F0F0", linewidth = 0.20),
    axis.text.y = element_text(size = 6.5),
    plot.title = element_text(face = "bold"),
    plot.title.position = "plot"
  )
ggsave(
  file.path(OUTPUT_DIR, "sample_by_aligned_peak_map.pdf"),
  figure, width = 14,
  height = max(10, min(18, 3 + 0.13 * length(sample_ids)))
)

# 6. Paginated same-feature identity atlas. The overview above preserves the
# real RT axis, whereas this atlas allocates equal width to a small number of
# features per page. This makes rare and tightly spaced features legible. A
# colored point means that the sample contains that aligned feature; a pale
# open point means no qualifying detection; and an orange cross means a local
# detection was deliberately not forced into that feature. Feature colors are
# stable across the entire run and are recorded in the color dictionary.
FEATURES_PER_IDENTITY_PAGE <- 16L
identity_page_count <- ceiling(feature_count / FEATURES_PER_IDENTITY_PAGE)
identity_manifest <- vector("list", identity_page_count)
identity_atlas_path <- file.path(
  OUTPUT_DIR, "sample_by_aligned_feature_identity_atlas.pdf"
)
grDevices::pdf(
  identity_atlas_path, width = 14,
  height = max(10, min(18, 3 + 0.13 * length(sample_ids))),
  onefile = TRUE
)
for (page_number in seq_len(identity_page_count)) {
  first_feature <- (page_number - 1L) * FEATURES_PER_IDENTITY_PAGE + 1L
  last_feature <- min(
    page_number * FEATURES_PER_IDENTITY_PAGE, feature_count
  )
  page_features <- feature_summary$aligned_feature_id[
    first_feature:last_feature
  ]
  page_summary <- feature_summary[first_feature:last_feature, , drop = FALSE]
  page_data <- abundance_long[
    abundance_long$aligned_feature_id %in% page_features, , drop = FALSE
  ]
  page_data$aligned_feature_id <- factor(
    page_data$aligned_feature_id, levels = page_features
  )
  page_data$sample_id <- factor(page_data$sample_id, levels = rev(sample_ids))
  page_data$atlas_status <- ifelse(
    page_data$unresolved_candidate_flag > 0 &
      page_data$integrated_dad_area == 0,
    "Detected but unresolved",
    ifelse(page_data$integrated_dad_area > 0, "Aligned peak", "Not detected")
  )
  page_data$atlas_status <- factor(
    page_data$atlas_status,
    levels = c("Aligned peak", "Detected but unresolved", "Not detected")
  )
  detected_page <- page_data[page_data$atlas_status == "Aligned peak", ]
  unresolved_page <- page_data[
    page_data$atlas_status == "Detected but unresolved", , drop = FALSE
  ]
  absent_page <- page_data[page_data$atlas_status == "Not detected", ]
  feature_axis_labels <- paste0(
    page_summary$aligned_feature_id, "\n",
    formatC(page_summary$consensus_retention_time,
            digits = 3, format = "f"), " min"
  )
  names(feature_axis_labels) <- page_summary$aligned_feature_id

  atlas_figure <- ggplot() +
    geom_point(
      data = absent_page,
      aes(aligned_feature_id, sample_id, shape = atlas_status),
      color = "#D9D9D9", size = 0.85, stroke = 0.35
    ) +
    geom_point(
      data = detected_page,
      aes(
        aligned_feature_id, sample_id,
        color = aligned_feature_id, shape = atlas_status
      ),
      size = 2.05, alpha = 0.96
    ) +
    geom_point(
      data = unresolved_page,
      aes(aligned_feature_id, sample_id, shape = atlas_status),
      color = "#E67E22", size = 2.25, stroke = 0.75
    ) +
    scale_color_manual(
      values = feature_identity_colors[page_features], guide = "none"
    ) +
    scale_shape_manual(
      values = c(
        "Aligned peak" = 16, "Detected but unresolved" = 4,
        "Not detected" = 1
      ),
      name = "Sample-feature result",
      drop = FALSE
    ) +
    scale_x_discrete(labels = feature_axis_labels, drop = FALSE) +
    labs(
      title = paste(
        ALIGNMENT_FIGURE_TITLE_PREFIX,
        "same-feature identity atlas - page",
        page_number, "of", identity_page_count
      ),
      subtitle = paste(
        "Within and across pages, each feature keeps its recorded color;",
        "matching column and color means the same aligned feature"
      ),
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
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.title.position = "plot"
    ) +
    guides(
      shape = guide_legend(
        nrow = 1L,
        override.aes = list(
          color = c("#2C7FB8", "#E67E22", "#D9D9D9"),
          size = c(2.6, 2.8, 2.1)
        )
      )
    )
  print(atlas_figure)

  identity_manifest[[page_number]] <- data.frame(
    alignment_policy = ALIGNMENT_POLICY,
    alignment_policy_display_name = ALIGNMENT_POLICY_DISPLAY_NAME,
    atlas_page = page_number,
    atlas_total_pages = identity_page_count,
    feature_order = first_feature:last_feature,
    aligned_feature_id = page_summary$aligned_feature_id,
    consensus_retention_time = page_summary$consensus_retention_time,
    figure_color_hex = unname(
      feature_identity_colors[page_summary$aligned_feature_id]
    ),
    stringsAsFactors = FALSE
  )
}
grDevices::dev.off()
identity_manifest <- do.call(rbind, identity_manifest)
write.csv(
  identity_manifest,
  file.path(
    OUTPUT_DIR, "sample_by_aligned_feature_identity_atlas_manifest.csv"
  ),
  row.names = FALSE
)

# 7. Rare-fragment reconciliation audit. The shaded lower-left region is only
# the first two gates; accepted rescues must also satisfy margin, bootstrap, and
# one-peak-per-sample requirements recorded in the CSV audit.
reconciliation_plot <- fragment_reconciliation[
  is.finite(fragment_reconciliation$rt_difference_minutes) &
    is.finite(fragment_reconciliation$spectral_distance), , drop = FALSE
]
if (nrow(reconciliation_plot)) {
  reconciliation_plot$display_decision <- ifelse(
    reconciliation_plot$reconciliation_decision ==
      "Rescued into recurrent feature",
    "Reconciled with recurrent feature",
    ifelse(grepl("competing|bootstrap|target already",
                 reconciliation_plot$reconciliation_decision,
                 ignore.case = TRUE),
           "Not reconciled: ambiguous or unstable",
           "Retained as rare evidence")
  )
  reconciliation_plot$rescue_timing <- ifelse(
    reconciliation_plot$reconciliation_decision !=
      "Rescued into recurrent feature",
    "Not reconciled",
    ifelse(
      reconciliation_plot$rescue_iteration_accepted == 1L,
      "Reconciled in first round", "Reconciled after consensus update"
    )
  )
  figure <- ggplot(
    reconciliation_plot,
    aes(rt_difference_minutes, spectral_distance,
        color = display_decision,
        shape = rescue_timing,
        size = source_group_samples)
  ) +
    annotate(
      "rect", xmin = 0, xmax = RARE_RECONCILIATION_RT_WINDOW,
      ymin = 0, ymax = RARE_RECONCILIATION_SPECTRAL_TOLERANCE,
      fill = "#D9F0D3", alpha = 0.45
    ) +
    geom_vline(xintercept = RARE_RECONCILIATION_RT_WINDOW,
               color = "#777777", linetype = 2, linewidth = 0.35) +
    geom_hline(yintercept = RARE_RECONCILIATION_SPECTRAL_TOLERANCE,
               color = "#777777", linetype = 2, linewidth = 0.35) +
    geom_point(alpha = 0.88, stroke = 0.8) +
    scale_color_manual(values = c(
      "Reconciled with recurrent feature" = "#1B9E77",
      "Not reconciled: ambiguous or unstable" = "#D95F02",
      "Retained as rare evidence" = "#6A51A3"
    ), labels = c(
      "Reconciled with recurrent feature" =
        "Reconciled with better-supported feature",
      "Not reconciled: ambiguous or unstable" =
        "Not reconciled: ambiguous or unstable",
      "Retained as rare evidence" = "Retained as distinct/rare evidence"
    ), name = "Final decision") +
    scale_shape_manual(values = c(
      "Reconciled in first round" = 16,
      "Reconciled after consensus update" = 17,
      "Not reconciled" = 1
    ), name = "Iteration outcome") +
    scale_size_continuous(
      range = c(1.8, 6.5),
      breaks = unique(pmax(1, pretty(
        reconciliation_plot$source_group_samples, n = 4
      ))),
      name = "Original candidate-group\nsample count"
    ) +
    labs(
      title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                    "low-prevalence feature reconciliation audit"),
      subtitle = paste(
        "Consensus is rebuilt between rounds; the green region is necessary",
        "but margin, bootstrap support, and sample occupancy must also pass"
      ),
      x = "Difference from proposed recurrent feature (minutes)",
      y = "Composite spectral distance"
    ) +
    guides(
      color = guide_legend(
        order = 1L, nrow = 2L, byrow = TRUE,
        title.position = "top", title.hjust = 0.5,
        override.aes = list(size = 3.2)
      ),
      shape = guide_legend(
        order = 2L, nrow = 1L,
        title.position = "top", title.hjust = 0.5,
        override.aes = list(size = 3.2)
      ),
      size = guide_legend(
        order = 3L, nrow = 1L,
        title.position = "top", title.hjust = 0.5
      )
    ) + alignment_theme
  ggsave(file.path(OUTPUT_DIR,
                   "rare_fragment_reconciliation_diagnostics.pdf"),
         figure, width = 11, height = 9)
}

# 8. Iterative rescue convergence. The final gray, zero-height round is an
# explicit stopping test: a complete pass found no additional peak that met
# every rescue gate. Newly created groups cannot become targets because target
# eligibility is frozen to the original recurrent anchors.
convergence_plot <- data.frame(
  rescue_iteration = seq_along(rescues_by_iteration),
  rescued_detections = rescues_by_iteration,
  round_status = ifelse(
    rescues_by_iteration > 0L, "Reconciliation round", "Convergence round"
  )
)
figure <- ggplot(
  convergence_plot,
  aes(rescue_iteration, rescued_detections, fill = round_status)
) +
  geom_col(width = 0.68, color = "#222222", linewidth = 0.35) +
  geom_text(
    aes(label = rescued_detections), vjust = -0.45,
    fontface = "bold", size = 3.8
  ) +
  scale_fill_manual(
    values = c(
      "Reconciliation round" = "#2C7FB8",
      "Convergence round" = "#BDBDBD"
    ),
    name = NULL
  ) +
  scale_x_continuous(
    breaks = convergence_plot$rescue_iteration,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(n = 5),
    expand = expansion(mult = c(0, 0.14))
  ) +
  labs(
    title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                  "iterative reconciliation convergence"),
    subtitle = paste(
      "Original recurrent anchors only; each detection may move once, and",
      "the final zero-reconciliation round establishes convergence"
    ),
    x = "Reconciliation round",
    y = "Detections reconciled in that round"
  ) +
  guides(fill = guide_legend(
    nrow = 1L, title.position = "top", title.hjust = 0.5
  )) + alignment_theme
ggsave(
  file.path(OUTPUT_DIR, "rare_fragment_rescue_convergence.pdf"),
  figure, width = 10, height = 7
)

# 9. Neighboring-feature decisions. This makes the final rule auditable: the
# lower-left region is the RT/spectral neighborhood considered, while color and
# shape show whether a pair merged, became a sensitivity family, remained
# separate, or was protected by a same-sample double peak.
pair_plot <- neighboring_feature_pair_audit
if (nrow(split_feature_merge_audit)) {
  merged_pairs <- split_feature_merge_audit[
    split_feature_merge_audit$automatic_merge, , drop = FALSE
  ]
  if (nrow(merged_pairs)) {
    merged_plot <- data.frame(
      consensus_rt_difference_minutes =
        merged_pairs$consensus_rt_difference_minutes,
      consensus_spectral_distance =
        merged_pairs$consensus_spectral_distance,
      same_sample_doublet_count = merged_pairs$same_sample_doublet_count,
      pair_decision =
        "Merged probable split homologue",
      stringsAsFactors = FALSE
    )
    pair_plot <- rbind(
      pair_plot[, names(merged_plot), drop = FALSE], merged_plot
    )
  }
}
if (nrow(pair_plot)) {
  pair_plot$display_decision <- pair_plot$pair_decision
  pair_plot$display_decision[
    pair_plot$pair_decision == "Unresolved neighboring-feature family"
  ] <- "Unresolved family: summed only for sensitivity"
  pair_plot$display_decision[
    pair_plot$pair_decision ==
      "Protected resolved double peak: same-sample co-occurrence"
  ] <- "Protected double peak: observed in same sample"
  pair_plot$display_decision[
    pair_plot$pair_decision == "Retained separately: assignments stable"
  ] <- "Retained separately: assignments stable"
  pair_plot$display_decision <- factor(
    pair_plot$display_decision,
    levels = c(
      "Merged probable split homologue",
      "Unresolved family: summed only for sensitivity",
      "Protected double peak: observed in same sample",
      "Retained separately: assignments stable"
    )
  )
  figure <- ggplot(
    pair_plot,
    aes(consensus_rt_difference_minutes,
        consensus_spectral_distance,
        color = display_decision, shape = display_decision)
  ) +
    annotate(
      "rect", xmin = 0, xmax = SPLIT_FEATURE_MERGE_RT_TOLERANCE,
      ymin = 0, ymax = SPLIT_FEATURE_MERGE_SPECTRAL_TOLERANCE,
      fill = "#D9F0D3", alpha = 0.45
    ) +
    geom_vline(xintercept = SPLIT_FEATURE_MERGE_RT_TOLERANCE,
               color = "#777777", linetype = 2, linewidth = 0.35) +
    geom_hline(yintercept = SPLIT_FEATURE_MERGE_SPECTRAL_TOLERANCE,
               color = "#777777", linetype = 2, linewidth = 0.35) +
    geom_point(size = 3, alpha = 0.88) +
    scale_color_manual(values = c(
      "Merged probable split homologue" = "#1B9E77",
      "Unresolved family: summed only for sensitivity" = "#D95F02",
      "Protected double peak: observed in same sample" = "#7570B3",
      "Retained separately: assignments stable" = "#666666"
    ), name = "Pair decision", drop = FALSE) +
    scale_shape_manual(values = c(
      "Merged probable split homologue" = 16,
      "Unresolved family: summed only for sensitivity" = 17,
      "Protected double peak: observed in same sample" = 15,
      "Retained separately: assignments stable" = 1
    ), name = "Pair decision", drop = FALSE) +
    labs(
      title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                    "neighboring-feature decision audit"),
      subtitle = paste(
        "Automatic merging additionally requires mutual exclusion, compact",
        "combined RT span, and at least 85% wavelength-block support"
      ),
      x = "Consensus retention-time difference (minutes)",
      y = "Consensus spectral distance"
    ) +
    guides(
      color = guide_legend(
        nrow = 2L, byrow = TRUE, title.position = "top", title.hjust = 0.5
      ),
      shape = "none"
    ) + alignment_theme
  ggsave(file.path(OUTPUT_DIR,
                   "neighboring_feature_reconciliation_diagnostics.pdf"),
         figure, width = 11, height = 8.5)
} else {
  figure <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label = paste(
        "No neighboring feature pairs met the RT and spectral review window.",
        "No split-feature merges or unresolved sensitivity families were needed.",
        sep = "\n"
      ),
      size = 4.2, lineheight = 1.3
    ) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(
      title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                    "neighboring-feature decision audit"),
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.title.position = "plot")
  ggsave(file.path(OUTPUT_DIR,
                   "neighboring_feature_reconciliation_diagnostics.pdf"),
         figure, width = 11, height = 8.5)
}

# 10. Representative spectra. Every feature receives a line and uncertainty
# ribbon; pages keep labels readable even when a genus contains many features.
representative_uncertainty$aligned_feature_id <- factor(
  representative_uncertainty$aligned_feature_id,
  levels = feature_summary$aligned_feature_id
)
pdf(file.path(OUTPUT_DIR,
              "all_aligned_feature_representative_spectra_with_uncertainty.pdf"),
    width = 12, height = 9, onefile = TRUE)
features_per_page <- 12L
for (page_start in seq(1L, nrow(feature_summary), by = features_per_page)) {
  page_features <- feature_summary$aligned_feature_id[
    page_start:min(page_start + features_per_page - 1L, nrow(feature_summary))
  ]
  page_summary <- feature_summary[match(
    page_features, feature_summary$aligned_feature_id
  ), , drop = FALSE]
  # Prevalence is based on aligned detections across every input sample, not on
  # the number of spectra remaining after uncertainty summarization. Showing it
  # in the panel header prevents a precise-looking ribbon from obscuring that a
  # feature may occur in only a small fraction of the dataset.
  page_panel_labels <- setNames(
    paste0(
      page_summary$aligned_feature_id, " | ",
      round(100 * page_summary$detection_fraction), "% of samples"
    ),
    page_summary$aligned_feature_id
  )
  page_data <- representative_uncertainty[
    representative_uncertainty$aligned_feature_id %in% page_features, , drop = FALSE
  ]
  page_data$feature_panel <- factor(
    page_panel_labels[as.character(page_data$aligned_feature_id)],
    levels = unname(page_panel_labels[page_features])
  )
  figure <- ggplot(
    page_data,
    aes(wavelength_nm, median_normalized_absorbance)
  ) +
    geom_ribbon(aes(ymin = q10_normalized_absorbance,
                    ymax = q90_normalized_absorbance),
                fill = "#9ECAE1", alpha = 0.30) +
    geom_ribbon(aes(ymin = q25_normalized_absorbance,
                    ymax = q75_normalized_absorbance),
                fill = "#3182BD", alpha = 0.28) +
    geom_line(color = "#08519C", linewidth = 0.65) +
    facet_wrap(~ feature_panel, ncol = 3, scales = "free_y") +
    labs(
      title = paste(ALIGNMENT_FIGURE_TITLE_PREFIX,
                    "aligned-feature representative spectra"),
      subtitle = "Median line; dark ribbon = IQR; light ribbon = 10th-90th percentile",
      x = "Wavelength (nm)", y = "Normalized absorbance"
    ) + alignment_theme +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  print(figure)
}
grDevices::dev.off()

message(
  "Raw DAD alignment complete: ", nrow(detections), " detected peaks -> ",
  nrow(feature_summary), " aligned features across ", length(sample_ids),
  " samples."
)
message("Output folder: ", OUTPUT_DIR)
