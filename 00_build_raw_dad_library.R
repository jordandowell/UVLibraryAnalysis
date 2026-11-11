## ========================================================================== ##
## BUILD THE ACTIVE SPECTRAL LIBRARY FROM AUTHENTICATED RAW STANDARD RUNS
## ========================================================================== ##
## Purpose
##   Convert one Agilent DAD chromatogram per authenticated standard replicate
##   into one concentration-invariant auxiliary spectrum. The existing CSV
##   library is retained as the identity anchor and active point estimate
##   because leave-one-replicate-out testing showed clearer spectral separation.
##   Raw spectra add QC and within-standard uncertainty evidence but are not
##   appended as extra training rows, so the same injections are never counted
##   twice.
##
## Important separation of evidence
##   Standard files may have been collected with different HPLC methods.
##   Consequently, retention time and peak width are used only to find and
##   integrate a peak inside its own raw file. Neither enters sample annotation.
##   The authenticated CSV retention time defines a tight primary search zone
##   and a wider recovery zone. Spectral agreement with both the compound
##   median and the matching explicit replicate then chooses the peak. The
##   full chromatogram is only a last-resort diagnostic fallback.
##
## Compound-level fallback
##   All three raw replicates for a compound must pass extraction and replicate
##   concordance checks before they are retained as auxiliary evidence. Failed
##   sets are quarantined, while the original CSV rows remain active point
##   estimates for every compound.
##
## Outputs
##   library/Absorbance_StandardLibrary_RawDAD.csv        active model input
##   library/raw_dad_derived/*                            QC tables and PDFs
## ========================================================================== ##

if (!requireNamespace("chromConverter", quietly = TRUE)) {
  stop("Package 'chromConverter' is required to build the raw DAD library.")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package 'ggplot2' is required for raw-library QC figures.")
}

suppressPackageStartupMessages(library(ggplot2))

ANCHOR_FILE <- file.path("library", "Absorbance_StandardLibrary_Test.csv")
MANIFEST_FILE <- file.path("library", "raw_dad_library_manifest.csv")
CLASS_FILE <- file.path("library", "compound_class_dictionary.csv")
STANDARD_CURVE_FILE <- file.path("library", "standard_curve_manifest.csv")
OUTPUT_DIR <- file.path("library", "raw_dad_derived")
ACTIVE_FILE <- file.path("library", "Absorbance_StandardLibrary_RawDAD.csv")
ANALYSIS_MIN_NM <- 190
ANALYSIS_MAX_NM <- 500
BASELINE_WINDOWS_MINUTES <- c(0.30, 0.60, 1.20, 2.40)
MINIMUM_PROFILE_SNR <- 6
MAXIMUM_ANCHOR_DISTANCE <- 0.16
MAXIMUM_EXPECTED_REGION_ANCHOR_DISTANCE <- 0.45
DOMINANT_RECOVERY_SIGNAL_RATIO <- 25
MINIMUM_WITHIN_PEAK_CORRELATION <- 0.85
# This threshold assesses reproducibility among raw injections acquired with
# the same raw-standard method. It does not compare their retention times with
# the curated CSV reference, which may have been acquired by another method.
MAXIMUM_RAW_REPLICATE_DISTANCE <- 0.12
MAXIMUM_RAW_APEX_RT_RANGE_MINUTES <- 0.30
PRIMARY_EXPECTED_RT_SEARCH_HALF_WINDOW_MINUTES <- 0.45
EXPECTED_RT_SEARCH_HALF_WINDOW_MINUTES <- 1.30
# Identity conflicts require all three conditions. This combination avoids
# flagging benign spectral-equivalence ties: the assigned anchor must disagree
# substantially, another authenticated anchor must be a close match, and the
# alternative must improve the distance by at least threefold.
IDENTITY_LARGE_ASSIGNED_DISTANCE <- 0.15
IDENTITY_STRONG_ALTERNATIVE_DISTANCE <- 0.08
IDENTITY_MINIMUM_IMPROVEMENT_RATIO <- 3

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

required_files <- c(ANCHOR_FILE, MANIFEST_FILE, CLASS_FILE)
if (any(!file.exists(required_files))) {
  stop("Raw-library builder is missing: ",
       paste(required_files[!file.exists(required_files)], collapse = ", "))
}

anchor <- read.csv(ANCHOR_FILE, check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.csv(MANIFEST_FILE, check.names = FALSE,
                     stringsAsFactors = FALSE)
class_dictionary <- read.csv(CLASS_FILE, stringsAsFactors = FALSE)

# The standard-curve manifest is deliberately allowed to contain only its
# header until concentration-series data become available. It establishes a
# stable future input contract without enabling quantitative calibration or
# implying that the current library supports concentration estimates.
standard_curve_required_columns <- c(
  "compound", "raw_d_folder", "replicate_folder", "dad1_uv",
  "concentration_value", "concentration_unit", "injection_volume_uL",
  "role", "notes"
)
standard_curve_manifest <- if (file.exists(STANDARD_CURVE_FILE)) {
  read.csv(STANDARD_CURVE_FILE, stringsAsFactors = FALSE,
           check.names = FALSE)
} else {
  data.frame(matrix(nrow = 0L, ncol = length(standard_curve_required_columns),
                    dimnames = list(NULL, standard_curve_required_columns)))
}
missing_standard_curve_columns <- setdiff(
  standard_curve_required_columns, names(standard_curve_manifest)
)
if (length(missing_standard_curve_columns)) {
  stop(
    "Standard-curve manifest is missing required columns: ",
    paste(missing_standard_curve_columns, collapse = ", ")
  )
}
required_manifest <- c("compound", "replicate_folder", "role", "dad1_uv")
if (length(setdiff(required_manifest, names(manifest)))) {
  stop("Raw DAD manifest is missing required columns.")
}
standard_manifest <- manifest[manifest$role == "standard_replicate", , drop = FALSE]
if (!nrow(standard_manifest) || any(!file.exists(standard_manifest$dad1_uv))) {
  stop("The raw DAD manifest contains missing standard replicate files.")
}

compound_key <- function(x) tolower(gsub("[[:space:]]+", " ", trimws(x)))
anchor_wavelengths <- suppressWarnings(as.numeric(sub(
  "nm$", "", names(anchor), ignore.case = TRUE
)))
wavelength_columns <- which(is.finite(anchor_wavelengths))
target_wavelengths <- anchor_wavelengths[wavelength_columns]
model_columns <- which(target_wavelengths >= ANALYSIS_MIN_NM &
                         target_wavelengths <= ANALYSIS_MAX_NM)
if (length(model_columns) < 20L) stop("Anchor library lacks a usable UV grid.")

normalize_spectrum <- function(x) {
  x <- pmax(as.numeric(x), 0)
  span <- diff(range(x, na.rm = TRUE))
  if (!is.finite(span) || span <= sqrt(.Machine$double.eps)) {
    return(rep(NA_real_, length(x)))
  }
  (x - min(x, na.rm = TRUE)) / span
}

safe_correlation <- function(x, y) {
  value <- suppressWarnings(stats::cor(x, y))
  if (is.finite(value)) value else -1
}

# This lightweight distance mirrors the concentration-invariant components
# used during raw-sample alignment. The full annotation model later recomputes
# its six-metric ensemble and reliability weights from the active library.
shape_distance <- function(x, y) {
  x <- normalize_spectrum(x)
  y <- normalize_spectrum(y)
  if (anyNA(x) || anyNA(y)) return(Inf)
  cosine <- 1 - sum(x * y) / sqrt(sum(x^2) * sum(y^2))
  correlation <- (1 - safe_correlation(x, y)) / 2
  derivative <- (1 - safe_correlation(diff(x), diff(y))) / 2
  max(0, 0.55 * cosine + 0.30 * correlation + 0.15 * derivative)
}

shape_distance_rows <- function(values, reference) {
  vapply(seq_len(nrow(values)), function(i) {
    shape_distance(values[i, ], reference)
  }, numeric(1))
}

odd_window_points <- function(minutes, retention_times, minimum = 9L) {
  points <- max(minimum, as.integer(round(minutes /
                                            stats::median(diff(retention_times)))))
  if (points %% 2L == 0L) points <- points + 1L
  maximum <- if (length(retention_times) %% 2L) length(retention_times) else
    length(retention_times) - 1L
  min(points, maximum)
}

local_maxima <- function(x) {
  which(x > c(-Inf, head(x, -1L)) & x >= c(tail(x, -1L), -Inf))
}

profile_noise <- function(profile) {
  value <- stats::mad(diff(profile), center = 0, constant = 1.4826,
                      na.rm = TRUE) / sqrt(2)
  if (!is.finite(value) || value <= 0) {
    value <- max(stats::quantile(profile[profile > 0], 0.10, na.rm = TRUE), 1e-8)
  }
  value
}

interpolate_grid <- function(values, native_wavelengths) {
  t(vapply(seq_len(nrow(values)), function(i) {
    stats::approx(native_wavelengths, values[i, ], xout = target_wavelengths,
                  method = "linear", rule = 2, ties = "ordered")$y
  }, numeric(length(target_wavelengths))))
}

anchor_medians <- do.call(rbind, lapply(unique(anchor$Compound), function(cp) {
  rows <- compound_key(anchor$Compound) == compound_key(cp)
  apply(anchor[rows, wavelength_columns, drop = FALSE], 2, median)
}))
rownames(anchor_medians) <- unique(anchor$Compound)

# Extract one raw replicate. Multiple baseline windows allow narrow and broad
# standard peaks to compete fairly. Candidate peaks must be chromatographic
# maxima with measurable signal; spectral agreement selects among them.
extract_raw_replicate <- function(manifest_row, replicate_number) {
  compound <- manifest_row$compound
  raw <- chromConverter::read_chemstation_uv(
    manifest_row$dad1_uv, format_out = "matrix", data_format = "wide",
    read_metadata = TRUE
  )
  retention_times <- suppressWarnings(as.numeric(rownames(raw)))
  native_wavelengths <- suppressWarnings(as.numeric(colnames(raw)))
  if (any(!is.finite(retention_times)) || any(!is.finite(native_wavelengths)) ||
      any(!is.finite(raw))) stop("Invalid raw standard file: ", manifest_row$dad1_uv)
  values <- interpolate_grid(raw, native_wavelengths)
  reference <- anchor_medians[
    match(compound_key(compound), compound_key(rownames(anchor_medians))),
    model_columns
  ]
  anchor_rows <- which(compound_key(anchor$Compound) == compound_key(compound))
  replicate_match <- anchor_rows[match(
    as.character(replicate_number), as.character(anchor$Rep[anchor_rows])
  )]
  if (is.na(replicate_match)) replicate_match <- anchor_rows[replicate_number]
  explicit_reference <- as.numeric(anchor[
    replicate_match, wavelength_columns[model_columns], drop = TRUE
  ])
  explicit_retention_time <- as.numeric(anchor$RT[replicate_match])
  if (!is.finite(explicit_retention_time)) {
    stop("Explicit anchor retention time is missing for ", compound,
         " replicate ", replicate_number, ".")
  }

  candidates <- list()
  raw_profile <- sqrt(rowMeans(pmax(values, 0)^2))
  candidate_counter <- 1L
  for (window in BASELINE_WINDOWS_MINUTES) {
    width <- odd_window_points(window, retention_times)
    profile_baseline <- stats::runmed(raw_profile, width, endrule = "median")
    profile <- pmax(raw_profile - profile_baseline, 0)
    noise <- profile_noise(profile)
    maxima <- local_maxima(profile)
    maxima <- maxima[profile[maxima] / noise >= MINIMUM_PROFILE_SNR]
    if (!length(maxima)) maxima <- order(profile, decreasing = TRUE)[seq_len(5)]
    candidate_spectra <- t(vapply(maxima, function(apex) {
      # Estimate a constant spectrum immediately outside the putative peak.
      # This is fast, preserves broad bands, and avoids fitting a global
      # chromatographic baseline to every wavelength independently.
      flank <- which(
        abs(retention_times - retention_times[apex]) >= 0.55 * window &
          abs(retention_times - retention_times[apex]) <= 0.90 * window
      )
      if (length(flank) < 4L) flank <- c(max(1L, apex - width),
                                         min(nrow(values), apex + width))
      background <- apply(values[flank, , drop = FALSE], 2, median)
      pmax(values[apex, ] - background, 0)
    }, numeric(ncol(values))))
    median_distances <- shape_distance_rows(
      candidate_spectra[, model_columns, drop = FALSE], reference
    )
    explicit_distances <- shape_distance_rows(
      candidate_spectra[, model_columns, drop = FALSE], explicit_reference
    )
    candidates[[candidate_counter]] <- data.frame(
      baseline_window_minutes = window,
      apex_index = maxima,
      candidate_retention_time_minutes = retention_times[maxima],
      candidate_origin = "Chromatographic local maximum",
      profile_snr = profile[maxima] / noise,
      candidate_corrected_rms = sqrt(rowMeans(
        candidate_spectra[, model_columns, drop = FALSE]^2
      )),
      compound_median_anchor_distance = median_distances,
      explicit_replicate_anchor_distance = explicit_distances,
      anchor_distance = 0.5 * median_distances + 0.5 * explicit_distances,
      within_primary_explicit_rt_search_region = abs(
        retention_times[maxima] - explicit_retention_time
      ) <= PRIMARY_EXPECTED_RT_SEARCH_HALF_WINDOW_MINUTES,
      within_explicit_rt_search_region = abs(
        retention_times[maxima] - explicit_retention_time
      ) <= EXPECTED_RT_SEARCH_HALF_WINDOW_MINUTES,
      stringsAsFactors = FALSE
    )
    candidate_counter <- candidate_counter + 1L
  }
  candidate_table <- do.call(rbind, candidates)
  candidate_table$selection_score <- candidate_table$anchor_distance +
    0.01 / log1p(pmax(candidate_table$profile_snr, 1))
  # The curated replicate spectrum and its stated RT jointly authenticate the
  # intended standard peak. First use a tight zone around the explicit RT. If
  # that zone genuinely lacks a credible peak, expand the search to accommodate
  # a method-dependent RT shift (important for broad or shifted standards such
  # as luteolin). Only then permit a full-chromatogram fallback. This hierarchy
  # prevents a spectrally convenient neighboring impurity from overriding the
  # experimentally identified standard peak.
  primary_candidates <- candidate_table[
    candidate_table$within_primary_explicit_rt_search_region &
      candidate_table$profile_snr >= MINIMUM_PROFILE_SNR, ,
    drop = FALSE
  ]
  # Inside the primary zone, RT proximity is real authentication evidence and
  # receives a modest penalty term. This prevents a neighboring impurity at the
  # edge of the zone from winning solely because its spectrum is cleaner.
  if (nrow(primary_candidates)) {
    primary_candidates$selection_score <-
      primary_candidates$selection_score + 0.08 * abs(
        primary_candidates$candidate_retention_time_minutes -
          explicit_retention_time
      ) / PRIMARY_EXPECTED_RT_SEARCH_HALF_WINDOW_MINUTES
  }
  recovery_candidates <- candidate_table[
    candidate_table$within_explicit_rt_search_region &
      candidate_table$profile_snr >= MINIMUM_PROFILE_SNR &
      candidate_table$anchor_distance <=
        MAXIMUM_EXPECTED_REGION_ANCHOR_DISTANCE, , drop = FALSE
  ]
  primary_pick <- if (nrow(primary_candidates)) {
    primary_candidates[which.min(primary_candidates$selection_score), ,
                       drop = FALSE]
  } else primary_candidates
  strongest_recovery_signal <- if (nrow(recovery_candidates)) {
    max(recovery_candidates$candidate_corrected_rms, na.rm = TRUE)
  } else -Inf
  primary_signal <- if (nrow(primary_pick)) {
    primary_pick$candidate_corrected_rms
  } else -Inf
  dominant_recovery_exists <- nrow(primary_pick) && nrow(recovery_candidates) &&
    is.finite(strongest_recovery_signal) && is.finite(primary_signal) &&
    strongest_recovery_signal >=
      DOMINANT_RECOVERY_SIGNAL_RATIO * max(primary_signal, 1e-8)
  selection_tier <- if (nrow(primary_candidates) && !dominant_recovery_exists) {
    "Primary explicit-RT region"
  } else if (nrow(recovery_candidates)) {
    "Expanded explicit-RT recovery region"
  } else "Full-chromatogram fallback"
  selection_pool <- switch(
    selection_tier,
    "Primary explicit-RT region" = primary_pick,
    "Expanded explicit-RT recovery region" = recovery_candidates,
    candidate_table
  )
  used_primary_region <- selection_tier == "Primary explicit-RT region"
  used_explicit_region <- selection_tier != "Full-chromatogram fallback"
  selected <- selection_pool[which.min(selection_pool$selection_score), ]
  apex <- selected$apex_index
  window <- selected$baseline_window_minutes
  flank <- which(
    abs(retention_times - retention_times[apex]) >= 0.55 * window &
      abs(retention_times - retention_times[apex]) <= 0.90 * window
  )
  if (length(flank) < 4L) {
    width <- odd_window_points(window, retention_times)
    flank <- c(max(1L, apex - width), min(nrow(values), apex + width))
  }
  background <- apply(values[flank, , drop = FALSE], 2, median)
  corrected <- pmax(sweep(values, 2, background, "-"), 0)
  profile <- pmax(
    raw_profile - stats::runmed(
      raw_profile, odd_window_points(window, retention_times),
      endrule = "median"
    ), 0
  )

  # Integrate the contiguous high-signal portion of the selected peak. The
  # region is bounded at 1.5 baseline windows so a broad standard may remain
  # broad without absorbing unrelated chromatographic peaks.
  maximum_half_width <- 1.5 * window
  allowed <- which(abs(retention_times - retention_times[apex]) <=
                     maximum_half_width)
  threshold <- max(0.15 * profile[apex], 3 * profile_noise(profile))
  high <- allowed[profile[allowed] >= threshold]
  left <- apex
  while (left > 1L && (left - 1L) %in% high) left <- left - 1L
  right <- apex
  while (right < length(retention_times) && (right + 1L) %in% high) {
    right <- right + 1L
  }
  region <- seq.int(left, right)
  if (length(region) < 3L) {
    region <- which(abs(retention_times - retention_times[apex]) <= 0.02)
  }
  high_region <- region[profile[region] >= 0.40 * max(profile[region])]
  if (!length(high_region)) high_region <- apex

  # Normalizing every high-signal scan before averaging makes the consensus
  # insensitive to concentration and gives broad peaks a stable shape without
  # allowing the apex alone to dominate.
  normalized_scans <- t(apply(
    corrected[high_region, , drop = FALSE], 1, normalize_spectrum
  ))
  consensus <- apply(normalized_scans, 2, median, na.rm = TRUE)
  apex_spectrum <- normalize_spectrum(corrected[apex, ])
  consensus_distance <- 0.5 *
    shape_distance(consensus[model_columns], reference) + 0.5 *
    shape_distance(consensus[model_columns], explicit_reference)
  apex_distance <- 0.5 *
    shape_distance(apex_spectrum[model_columns], reference) + 0.5 *
    shape_distance(apex_spectrum[model_columns], explicit_reference)
  selected_spectrum <- if (is.finite(consensus_distance) &&
                           consensus_distance <= apex_distance + 0.03) {
    consensus
  } else apex_spectrum
  extraction_method <- if (identical(selected_spectrum, consensus)) {
    "Median normalized scans across high-signal peak region"
  } else "Baseline-corrected apex spectrum"

  scan_correlations <- vapply(seq_len(nrow(normalized_scans)), function(i) {
    safe_correlation(normalized_scans[i, model_columns],
                     selected_spectrum[model_columns])
  }, numeric(1))
  half_height <- region[profile[region] >= 0.50 * profile[apex]]
  fwhm <- if (length(half_height) >= 2L) {
    diff(range(retention_times[half_height]))
  } else 0
  all_anchor_distances <- vapply(seq_len(nrow(anchor_medians)), function(i) {
    shape_distance(selected_spectrum[model_columns],
                   anchor_medians[i, model_columns])
  }, numeric(1))
  distance_order <- order(all_anchor_distances)
  own_index <- match(compound_key(compound), compound_key(rownames(anchor_medians)))
  own_rank <- match(own_index, distance_order)

  qc <- data.frame(
    Compound = compound,
    Rep = replicate_number,
    raw_file = manifest_row$replicate_folder,
    raw_path = manifest_row$dad1_uv,
    explicit_anchor_retention_time_minutes = explicit_retention_time,
    primary_explicit_rt_search_half_window_minutes =
      PRIMARY_EXPECTED_RT_SEARCH_HALF_WINDOW_MINUTES,
    explicit_rt_search_half_window_minutes =
      EXPECTED_RT_SEARCH_HALF_WINDOW_MINUTES,
    peak_selection_tier = selection_tier,
    selected_within_primary_explicit_rt_search_region = used_primary_region,
    selected_within_explicit_rt_search_region = used_explicit_region,
    selected_rt_minus_explicit_anchor_minutes =
      retention_times[apex] - explicit_retention_time,
    raw_apex_retention_time_minutes = retention_times[apex],
    raw_peak_start_minutes = retention_times[min(region)],
    raw_peak_end_minutes = retention_times[max(region)],
    raw_peak_fwhm_minutes = fwhm,
    selected_baseline_window_minutes = selected$baseline_window_minutes,
    profile_snr = selected$profile_snr,
    scans_in_consensus = length(high_region),
    median_within_peak_correlation = median(scan_correlations, na.rm = TRUE),
    compound_median_anchor_distance =
      shape_distance(selected_spectrum[model_columns], reference),
    explicit_replicate_anchor_distance =
      shape_distance(selected_spectrum[model_columns], explicit_reference),
    anchor_distance = 0.5 *
      shape_distance(selected_spectrum[model_columns], reference) + 0.5 *
      shape_distance(selected_spectrum[model_columns], explicit_reference),
    nearest_anchor_compound = rownames(anchor_medians)[distance_order[1]],
    nearest_anchor_distance = all_anchor_distances[distance_order[1]],
    own_anchor_rank = own_rank,
    nearest_to_second_distance_margin =
      all_anchor_distances[distance_order[2]] - all_anchor_distances[distance_order[1]],
    extraction_method = extraction_method,
    replicate_qc_pass = selected$profile_snr >= MINIMUM_PROFILE_SNR &
      median(scan_correlations, na.rm = TRUE) >= MINIMUM_WITHIN_PEAK_CORRELATION &
      (if (used_explicit_region) {
        0.5 * shape_distance(selected_spectrum[model_columns], reference) +
          0.5 * shape_distance(
            selected_spectrum[model_columns], explicit_reference
          ) <= MAXIMUM_EXPECTED_REGION_ANCHOR_DISTANCE
      } else {
        all_anchor_distances[own_index] <= MAXIMUM_ANCHOR_DISTANCE
      }),
    stringsAsFactors = FALSE
  )
  profile_output <- data.frame(
    Compound = compound, Rep = replicate_number,
    raw_file = manifest_row$replicate_folder,
    retention_time_minutes = retention_times,
    explicit_anchor_retention_time_minutes = explicit_retention_time,
    baseline_corrected_rms = profile,
    selected_region = seq_along(retention_times) %in% region,
    selected_apex = seq_along(retention_times) == apex,
    stringsAsFactors = FALSE
  )
  candidate_table$Compound <- compound
  candidate_table$Rep <- replicate_number
  candidate_table$explicit_anchor_retention_time_minutes <-
    explicit_retention_time
  candidate_table$selected_candidate <-
    candidate_table$apex_index == apex &
    candidate_table$baseline_window_minutes == window
  candidate_table$peak_selection_tier <- selection_tier
  list(spectrum = selected_spectrum, qc = qc, profile = profile_output,
       candidates = candidate_table)
}

extractions <- list()
counter <- 1L
for (compound in unique(standard_manifest$compound)) {
  local <- standard_manifest[
    compound_key(standard_manifest$compound) == compound_key(compound), ,
    drop = FALSE
  ]
  local <- local[order(local$replicate_folder), , drop = FALSE]
  if (nrow(local) != 3L) {
    stop("Expected exactly three raw standard replicates for: ", compound)
  }
  message("Extracting raw standard: ", compound)
  for (replicate_number in seq_len(nrow(local))) {
    extractions[[counter]] <- extract_raw_replicate(
      local[replicate_number, , drop = FALSE], replicate_number
    )
    counter <- counter + 1L
  }
}

qc <- do.call(rbind, lapply(extractions, `[[`, "qc"))
raw_spectra <- do.call(rbind, lapply(extractions, `[[`, "spectrum"))
profile_data <- do.call(rbind, lapply(extractions, `[[`, "profile"))
candidate_audit <- do.call(rbind, lapply(extractions, `[[`, "candidates"))
colnames(raw_spectra) <- paste0(target_wavelengths, "nm")

# ----------------------- RAW IDENTITY-CONFLICT AUDIT -----------------------
# Compare the median raw spectrum assigned to each manifest label against every
# authenticated CSV anchor. A strong off-diagonal match is not automatically a
# label error because some standards are spectral equivalents. We flag only a
# large assigned-anchor disagreement combined with a close, at least threefold
# better alternative. Reciprocal A-to-B and B-to-A conflicts receive an
# explicit possible-label-swap warning.
raw_compounds <- unique(qc$Compound)
raw_identity_consensus <- do.call(rbind, lapply(raw_compounds, function(cp) {
  rows <- compound_key(qc$Compound) == compound_key(cp)
  apply(raw_spectra[rows, model_columns, drop = FALSE], 2, median)
}))
rownames(raw_identity_consensus) <- raw_compounds

identity_distance_matrix <- t(vapply(seq_along(raw_compounds), function(i) {
  vapply(seq_len(nrow(anchor_medians)), function(j) {
    shape_distance(raw_identity_consensus[i, ],
                   anchor_medians[j, model_columns])
  }, numeric(1))
}, numeric(nrow(anchor_medians))))
rownames(identity_distance_matrix) <- raw_compounds
colnames(identity_distance_matrix) <- rownames(anchor_medians)

identity_audit <- do.call(rbind, lapply(seq_along(raw_compounds), function(i) {
  cp <- raw_compounds[i]
  assigned_index <- match(
    compound_key(cp), compound_key(colnames(identity_distance_matrix))
  )
  if (is.na(assigned_index)) {
    stop("Raw standard has no authenticated CSV anchor: ", cp)
  }
  distances <- identity_distance_matrix[i, ]
  other_indices <- setdiff(seq_along(distances), assigned_index)
  nearest_other_index <- other_indices[which.min(distances[other_indices])]
  assigned_distance <- distances[assigned_index]
  alternative_distance <- distances[nearest_other_index]
  improvement_ratio <- assigned_distance /
    max(alternative_distance, sqrt(.Machine$double.eps))
  conflict <- assigned_distance >= IDENTITY_LARGE_ASSIGNED_DISTANCE &&
    alternative_distance <= IDENTITY_STRONG_ALTERNATIVE_DISTANCE &&
    improvement_ratio >= IDENTITY_MINIMUM_IMPROVEMENT_RATIO
  data.frame(
    Compound = cp,
    assigned_authenticated_anchor = colnames(identity_distance_matrix)[
      assigned_index
    ],
    assigned_anchor_distance = assigned_distance,
    nearest_alternative_anchor = colnames(identity_distance_matrix)[
      nearest_other_index
    ],
    nearest_alternative_distance = alternative_distance,
    assigned_to_alternative_distance_ratio = improvement_ratio,
    large_assigned_mismatch_threshold = IDENTITY_LARGE_ASSIGNED_DISTANCE,
    strong_alternative_match_threshold =
      IDENTITY_STRONG_ALTERNATIVE_DISTANCE,
    minimum_improvement_ratio = IDENTITY_MINIMUM_IMPROVEMENT_RATIO,
    identity_conflict_flag = conflict,
    stringsAsFactors = FALSE
  )
}))

identity_audit$possible_reciprocal_label_swap <- FALSE
identity_audit$suggested_swap_partner <- NA_character_
for (i in seq_len(nrow(identity_audit))) {
  if (!identity_audit$identity_conflict_flag[i]) next
  partner_index <- match(
    compound_key(identity_audit$nearest_alternative_anchor[i]),
    compound_key(identity_audit$Compound)
  )
  reciprocal <- !is.na(partner_index) &&
    identity_audit$identity_conflict_flag[partner_index] &&
    compound_key(identity_audit$nearest_alternative_anchor[partner_index]) ==
      compound_key(identity_audit$Compound[i])
  identity_audit$possible_reciprocal_label_swap[i] <- reciprocal
  if (reciprocal) {
    identity_audit$suggested_swap_partner[i] <-
      identity_audit$Compound[partner_index]
  }
}
identity_audit$recommended_action <- ifelse(
  identity_audit$possible_reciprocal_label_swap,
  paste0(
    "Review a possible reciprocal naming switch with ",
    identity_audit$suggested_swap_partner,
    "; raw auxiliary evidence is withheld until resolved."
  ),
  ifelse(
    identity_audit$identity_conflict_flag,
    paste0(
      "Review raw-file identity; the raw consensus matches ",
      identity_audit$nearest_alternative_anchor,
      " substantially better than its assigned anchor."
    ),
    "No strong raw-versus-authenticated identity conflict detected."
  )
)

if (any(identity_audit$identity_conflict_flag)) {
  warning(paste(
    "Raw DAD identity conflict(s):",
    paste(
      paste0(
        identity_audit$Compound[identity_audit$identity_conflict_flag],
        " -> ",
        identity_audit$nearest_alternative_anchor[
          identity_audit$identity_conflict_flag
        ]
      ),
      collapse = "; "
    ),
    "Review raw_standard_identity_conflict_audit.csv. Conflicting raw",
    "evidence will not be admitted."
  ), call. = FALSE)
}

# Assess raw replicate concordance before deciding at compound level. Spectral
# identity ties are allowed and reported but not used as a rejection criterion:
# several standards deliberately form spectral-equivalence groups, so requiring
# rank 1 would discard genuine raw spectra precisely where downstream tie logic
# is designed to preserve ambiguity.
compound_decisions <- do.call(rbind, lapply(unique(qc$Compound), function(cp) {
  rows <- compound_key(qc$Compound) == compound_key(cp)
  local <- raw_spectra[rows, model_columns, drop = FALSE]
  pairwise <- c(
    shape_distance(local[1, ], local[2, ]),
    shape_distance(local[1, ], local[3, ]),
    shape_distance(local[2, ], local[3, ])
  )
  all_replicates_pass <- all(qc$replicate_qc_pass[rows])
  replicate_set_pass <- max(pairwise) <= MAXIMUM_RAW_REPLICATE_DISTANCE
  identity_row <- match(compound_key(cp), compound_key(identity_audit$Compound))
  identity_conflict <- identity_audit$identity_conflict_flag[identity_row]
  # Admission deliberately ignores the absolute offset from the curated CSV
  # retention time. Those reference spectra and the raw .D standards may come
  # from different HPLC methods. Retention time is relevant here only for
  # checking that the three raw injections locate a reproducible peak under
  # their shared raw-standard method.
  admitted <- all_replicates_pass && replicate_set_pass && !identity_conflict &&
    diff(range(qc$raw_apex_retention_time_minutes[rows])) <=
      MAXIMUM_RAW_APEX_RT_RANGE_MINUTES
  data.frame(
    Compound = cp,
    raw_replicates_available = sum(rows),
    raw_replicates_passing_qc = sum(qc$replicate_qc_pass[rows]),
    median_raw_to_anchor_distance = median(qc$anchor_distance[rows]),
    maximum_raw_replicate_pair_distance = max(pairwise),
    raw_apex_rt_range_minutes = diff(range(
      qc$raw_apex_retention_time_minutes[rows]
    )),
    median_own_anchor_rank = median(qc$own_anchor_rank[rows]),
    maximum_own_anchor_rank = max(qc$own_anchor_rank[rows]),
    all_raw_replicates_pass = all_replicates_pass,
    raw_replicate_set_concordant = replicate_set_pass,
    maximum_allowed_raw_replicate_distance =
      MAXIMUM_RAW_REPLICATE_DISTANCE,
    raw_peak_location_concordant = diff(range(
      qc$raw_apex_retention_time_minutes[rows]
    )) <= MAXIMUM_RAW_APEX_RT_RANGE_MINUTES,
    curated_rt_offset_used_for_admission = FALSE,
    # Diagnostic provenance only; this field is not an admission criterion.
    all_selected_within_explicit_rt_region = all(
      qc$selected_within_explicit_rt_search_region[rows]
    ),
    spectral_equivalence_tie_present = any(qc$own_anchor_rank[rows] > 1L),
    identity_conflict_flag = identity_conflict,
    possible_reciprocal_label_swap =
      identity_audit$possible_reciprocal_label_swap[identity_row],
    suggested_swap_partner =
      identity_audit$suggested_swap_partner[identity_row],
    raw_auxiliary_evidence_status = if (admitted) {
      "Admitted raw DAD auxiliary evidence"
    } else if (identity_conflict) {
      "Raw DAD set quarantined - identity conflict"
    } else "Raw DAD set quarantined",
    active_point_estimate = "Authenticated CSV spectrum",
    # Retained for compatibility with previously distributed QC tables. It now
    # names raw-evidence admission, not the active training-row source.
    active_source = if (admitted) {
      "Raw DAD auxiliary evidence admitted"
    } else "Raw DAD auxiliary evidence quarantined",
    stringsAsFactors = FALSE
  )
}))

qc$compound_active_source <- compound_decisions$active_source[
  match(compound_key(qc$Compound), compound_key(compound_decisions$Compound))
]
qc$used_in_active_library <- qc$compound_active_source ==
  "Raw DAD auxiliary evidence admitted"

# Preserve exactly three independent replicate rows per compound. Validation
# showed that the curated spectra provide slightly clearer nearest-versus-
# second-nearest separation than raw peak consensus. Therefore they remain the
# active point estimates. Raw spectra are retained as admitted auxiliary
# within-standard evidence for QC, uncertainty envelopes, and future model
# comparisons; they are not appended or averaged into the training rows.
active <- anchor
raw_candidate <- anchor
active$Spectrum_Source <- "Authenticated CSV point estimate"
active$Raw_File <- NA_character_
active$Raw_Apex_RT_Minutes <- NA_real_
active$Raw_Peak_FWHM_Minutes <- NA_real_
for (i in seq_len(nrow(active))) {
  cp <- active$Compound[i]
  decision <- compound_decisions$active_source[
    compound_key(compound_decisions$Compound) == compound_key(cp)
  ]
  if (identical(decision, "Raw DAD auxiliary evidence admitted")) {
    rep_number <- as.integer(active$Rep[i])
    raw_row <- which(compound_key(qc$Compound) == compound_key(cp) &
                       qc$Rep == rep_number)
    raw_candidate[i, wavelength_columns] <- raw_spectra[raw_row, ]
    active$Spectrum_Source[i] <-
      "Authenticated CSV point estimate + admitted raw-DAD auxiliary evidence"
    active$Raw_File[i] <- qc$raw_file[raw_row]
    active$Raw_Apex_RT_Minutes[i] <- qc$raw_apex_retention_time_minutes[raw_row]
    active$Raw_Peak_FWHM_Minutes[i] <- qc$raw_peak_fwhm_minutes[raw_row]
  }
}

# Retention-time columns remain useful provenance but are explicitly named as
# method-specific metadata. The model discovers wavelengths by their `nm`
# suffix and therefore cannot accidentally use these fields as predictors.
write.csv(active, ACTIVE_FILE, row.names = FALSE, na = "")
write.csv(qc, file.path(OUTPUT_DIR, "raw_standard_replicate_extraction_qc.csv"),
          row.names = FALSE, na = "")
write.csv(candidate_audit,
          file.path(OUTPUT_DIR, "raw_standard_peak_candidate_audit.csv"),
          row.names = FALSE, na = "")
write.csv(compound_decisions,
          file.path(OUTPUT_DIR, "raw_standard_compound_admission_decisions.csv"),
          row.names = FALSE, na = "")
write.csv(identity_audit,
          file.path(OUTPUT_DIR, "raw_standard_identity_conflict_audit.csv"),
          row.names = FALSE, na = "")

standard_curve_status <- data.frame(
  manifest_file = STANDARD_CURVE_FILE,
  calibration_rows_available = nrow(standard_curve_manifest),
  quantitative_calibration_enabled = FALSE,
  status = if (nrow(standard_curve_manifest)) {
    paste(
      "Calibration rows supplied, but quantitative calibration is not yet",
      "enabled; validate curve design and implement the calibration module."
    )
  } else {
    paste(
      "No concentration-series rows supplied; current outputs are qualitative",
      "spectral analogues and relative DAD-area summaries only."
    )
  },
  stringsAsFactors = FALSE
)
write.csv(standard_curve_status,
          file.path(OUTPUT_DIR, "standard_curve_readiness.csv"),
          row.names = FALSE)
if (nrow(standard_curve_manifest)) {
  warning(
    "Standard-curve rows are present but quantitative calibration is not enabled.",
    call. = FALSE
  )
}

raw_export <- data.frame(
  Compound = qc$Compound, Rep = qc$Rep, raw_file = qc$raw_file,
  raw_spectra, check.names = FALSE, stringsAsFactors = FALSE
)
write.csv(raw_export,
          file.path(OUTPUT_DIR, "all_extracted_raw_standard_spectra.csv"),
          row.names = FALSE)

# -------------------------- VALIDATION COMPARISON --------------------------
# Leave one replicate out, rebuild the compound medians, and ask whether the
# held spectrum returns to its own compound. The evaluation is performed for
# both the original and active libraries with the same distance, so adoption of
# raw spectra cannot be justified merely by adding observations.
validate_library <- function(data, label) {
  spectra <- as.matrix(data[, wavelength_columns[model_columns], drop = FALSE])
  storage.mode(spectra) <- "double"
  rows <- lapply(seq_len(nrow(data)), function(i) {
    training <- seq_len(nrow(data)) != i
    representatives <- do.call(rbind, lapply(unique(data$Compound), function(cp) {
      local <- training & compound_key(data$Compound) == compound_key(cp)
      apply(spectra[local, , drop = FALSE], 2, median)
    }))
    rownames(representatives) <- unique(data$Compound)
    distances <- vapply(seq_len(nrow(representatives)), function(j) {
      shape_distance(spectra[i, ], representatives[j, ])
    }, numeric(1))
    ordered <- order(distances)
    predicted <- rownames(representatives)[ordered[1]]
    true_class <- class_dictionary$class_name[
      match(compound_key(data$Compound[i]),
            compound_key(class_dictionary$Compound))
    ]
    predicted_class <- class_dictionary$class_name[
      match(compound_key(predicted), compound_key(class_dictionary$Compound))
    ]
    data.frame(
      library_version = label,
      Compound = data$Compound[i], Rep = data$Rep[i],
      nearest_compound = predicted,
      exact_compound_correct = compound_key(predicted) ==
        compound_key(data$Compound[i]),
      nearest_class = predicted_class,
      class_correct = predicted_class == true_class,
      minimum_distance = distances[ordered[1]],
      minimum_to_second_ratio = distances[ordered[1]] /
        max(distances[ordered[2]], sqrt(.Machine$double.eps)),
      minimum_to_mean_other_ratio = distances[ordered[1]] /
        mean(distances[ordered[-1]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

validation <- rbind(
  validate_library(active, "Authenticated CSV active"),
  validate_library(raw_candidate, "Raw-DAD candidate")
)
write.csv(validation,
          file.path(OUTPUT_DIR, "leave_one_replicate_out_library_comparison.csv"),
          row.names = FALSE)
validation_summary <- do.call(rbind, lapply(
  split(validation, validation$library_version), function(local) {
    data.frame(
      library_version = local$library_version[1],
      n_replicates = nrow(local),
      exact_compound_accuracy = mean(local$exact_compound_correct),
      class_accuracy = mean(local$class_correct),
      median_minimum_distance = median(local$minimum_distance),
      median_minimum_to_second_ratio = median(local$minimum_to_second_ratio),
      median_minimum_to_mean_other_ratio =
        median(local$minimum_to_mean_other_ratio),
      stringsAsFactors = FALSE
    )
  }
))
write.csv(validation_summary,
          file.path(OUTPUT_DIR, "library_version_validation_summary.csv"),
          row.names = FALSE)

# ------------------------------- FIGURES -----------------------------------
theme_library <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey30"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Identity audit: points above the diagonal are closer to their assigned CSV
# anchor than to any alternative. Points well below it are candidates for a
# manifest error; reciprocal conflicts are specifically reported as possible
# two-compound naming switches. Luteolin and naringenin remain labeled after
# correction so the resolved pair is easy to verify visually.
identity_audit$audit_label <- ifelse(
  identity_audit$identity_conflict_flag |
    compound_key(identity_audit$Compound) %in% c("luteolin", "naringenin"),
  identity_audit$Compound, ""
)
identity_audit_limit <- 1.08 * max(
  identity_audit$assigned_anchor_distance,
  identity_audit$nearest_alternative_distance,
  IDENTITY_LARGE_ASSIGNED_DISTANCE,
  na.rm = TRUE
)
identity_figure <- ggplot(
  identity_audit,
  aes(assigned_anchor_distance, nearest_alternative_distance,
      color = identity_conflict_flag)
) +
  geom_abline(slope = 1, intercept = 0, color = "#666666",
              linetype = 2, linewidth = 0.45) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_text(
    data = identity_audit[nzchar(identity_audit$audit_label), , drop = FALSE],
    aes(assigned_anchor_distance, nearest_alternative_distance,
        label = audit_label),
    inherit.aes = FALSE, check_overlap = TRUE,
    nudge_x = 0.006, hjust = 0, size = 3.1
  ) +
  scale_color_manual(
    values = c("FALSE" = "#0072B2", "TRUE" = "#D55E00"),
    labels = c("FALSE" = "No strong conflict", "TRUE" = "Identity conflict"),
    name = "Raw-versus-anchor audit"
  ) +
  scale_x_continuous(
    limits = c(0, identity_audit_limit),
    breaks = scales::breaks_pretty(n = 5)
  ) +
  scale_y_continuous(
    limits = c(0, identity_audit_limit),
    breaks = scales::breaks_pretty(n = 5)
  ) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Raw-standard identity conflict audit",
    subtitle = paste(
      "Above diagonal = assigned anchor is better; below diagonal = another",
      "authenticated anchor is better"
    ),
    caption = paste(
      "Strong conflicts require a large assigned mismatch, a close alternative,",
      "and at least a threefold distance improvement."
    ),
    x = "Distance to assigned authenticated anchor",
    y = "Distance to nearest alternative anchor"
  ) +
  theme_library +
  theme(plot.caption = element_text(hjust = 0, size = 8.5))
ggsave(
  file.path(OUTPUT_DIR, "raw_standard_identity_conflict_audit.pdf"),
  identity_figure, width = 8.5, height = 7
)

# All raw chromatographic profiles are included, six per page. A blue band is
# the integrated region, the red line is a detected chromatographic apex, and
# the purple line is the curated reference RT. Keeping the last two separate
# makes a standard/metadata conflict visible instead of drawing an expected RT
# as though it were an observed peak.
profile_data$panel <- paste(profile_data$Compound,
                            paste0("replicate ", profile_data$Rep), sep = " - ")
panels <- unique(profile_data$panel)
grDevices::pdf(file.path(OUTPUT_DIR, "all_raw_standard_peak_extractions.pdf"),
               width = 11, height = 8.5, onefile = TRUE)
for (start in seq(1L, length(panels), by = 6L)) {
  keep_panels <- panels[start:min(start + 5L, length(panels))]
  local <- profile_data[profile_data$panel %in% keep_panels, , drop = FALSE]
  regions <- do.call(rbind, lapply(split(local, local$panel), function(z) {
    selected <- z[z$selected_region, , drop = FALSE]
    data.frame(panel = z$panel[1], xmin = min(selected$retention_time_minutes),
               xmax = max(selected$retention_time_minutes))
  }))
  apices <- local[local$selected_apex, c("panel", "retention_time_minutes")]
  explicit_anchors <- unique(local[, c(
    "panel", "explicit_anchor_retention_time_minutes"
  )])
  print(
    ggplot(local, aes(retention_time_minutes,
                      log1p(baseline_corrected_rms))) +
      geom_rect(data = regions, aes(xmin = xmin, xmax = xmax,
                                    ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = "#56B4E9", alpha = 0.18) +
      geom_line(linewidth = 0.45, color = "black") +
      geom_vline(data = apices, aes(xintercept = retention_time_minutes),
                 color = "#D55E00", linewidth = 0.55, linetype = 2) +
      geom_vline(
        data = explicit_anchors,
        aes(xintercept = explicit_anchor_retention_time_minutes),
        color = "#6A3D9A", linewidth = 0.55, linetype = 3
      ) +
      facet_wrap(~panel, ncol = 2, scales = "free_y") +
      labs(
        title = "Authenticated raw-standard peak extractions",
        subtitle = paste(
          "Blue = integrated region; red dashed = detected candidate apex;",
          "purple dotted = curated reference RT"
        ),
        x = "Retention time (minutes)",
        y = "log(1 + baseline-corrected RMS absorbance)"
      ) + theme_library
  )
}
grDevices::dev.off()

# One page per compound shows the raw-derived median and IQR against the
# authenticated CSV median and IQR. This makes broad-peak extraction effects
# directly reviewable without implying that RT is transferable across methods.
anchor_long <- list()
raw_long <- list()
for (cp in unique(anchor$Compound)) {
  arows <- compound_key(anchor$Compound) == compound_key(cp)
  rrows <- compound_key(qc$Compound) == compound_key(cp)
  am <- as.matrix(anchor[arows, wavelength_columns, drop = FALSE])
  rm <- raw_spectra[rrows, , drop = FALSE]
  normalized_anchor <- t(apply(am, 1, normalize_spectrum))
  anchor_long[[cp]] <- data.frame(
    Compound = cp, wavelength_nm = target_wavelengths,
    anchor_median = apply(normalized_anchor, 2, median),
    anchor_q25 = apply(normalized_anchor, 2, quantile, 0.25),
    anchor_q75 = apply(normalized_anchor, 2, quantile, 0.75),
    raw_median = apply(rm, 2, median),
    raw_q25 = apply(rm, 2, quantile, 0.25),
    raw_q75 = apply(rm, 2, quantile, 0.75),
    active_source = compound_decisions$active_source[
      compound_key(compound_decisions$Compound) == compound_key(cp)
    ], stringsAsFactors = FALSE
  )
}
spectra_plot_data <- do.call(rbind, anchor_long)
grDevices::pdf(file.path(OUTPUT_DIR,
                         "raw_vs_authenticated_standard_spectra_atlas.pdf"),
               width = 9, height = 6.5, onefile = TRUE)
for (cp in unique(spectra_plot_data$Compound)) {
  local <- spectra_plot_data[spectra_plot_data$Compound == cp, ]
  print(
    ggplot(local, aes(wavelength_nm)) +
      geom_ribbon(aes(ymin = raw_q25, ymax = raw_q75),
                  fill = "#0072B2", alpha = 0.18) +
      geom_line(aes(y = raw_median, color = "Raw DAD peak consensus"),
                linewidth = 0.85) +
      geom_ribbon(aes(ymin = anchor_q25, ymax = anchor_q75),
                  fill = "#D55E00", alpha = 0.12) +
      geom_line(aes(y = anchor_median, color = "Authenticated CSV anchor"),
                linewidth = 0.75, linetype = 2) +
      scale_color_manual(values = c(
      "Raw DAD peak consensus" = "#0072B2",
        "Authenticated CSV anchor" = "#D55E00"
      ), name = NULL) +
      coord_cartesian(xlim = c(190, 500)) +
      labs(
        title = cp,
        subtitle = paste(
          "Raw evidence:", local$active_source[1],
          "| curated line remains active | shaded bands = replicate IQR"
        ),
        x = "Wavelength (nm)", y = "Normalized absorbance"
      ) + theme_library
  )
}
grDevices::dev.off()

metric_long <- rbind(
  data.frame(library_version = validation_summary$library_version,
             metric = "Exact compound", value =
               validation_summary$exact_compound_accuracy),
  data.frame(library_version = validation_summary$library_version,
             metric = "Chemical class", value = validation_summary$class_accuracy)
)
# Add the separation statistic because accuracy alone is saturated for this
# deliberately replicated library. Lower nearest/second ratios indicate a
# clearer return to the held-out compound.
ratio_long <- data.frame(
  library_version = validation_summary$library_version,
  metric = "Median nearest / second distance",
  value = validation_summary$median_minimum_to_second_ratio,
  stringsAsFactors = FALSE
)
comparison_figure <- ggplot(metric_long,
                           aes(library_version, value, fill = library_version)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = scales::percent(value, accuracy = 1)),
            vjust = -0.35, fontface = "bold") +
  facet_wrap(~metric) +
  scale_fill_manual(values = c("Authenticated CSV active" = "#999999",
                               "Raw-DAD candidate" = "#0072B2")) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format()) +
  labs(
    title = "Leave-one-replicate-out validation of the active library",
    subtitle = paste0(
      "Raw consensus is evaluated after peak-purity, anchor-agreement, and ",
      "replicate-concordance checks.\nCurated spectra remain active."
    ),
    x = NULL, y = "Held-out replicate accuracy"
  ) + theme_library +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "library_version_validation_comparison.pdf"),
       comparison_figure, width = 9, height = 5.8)

separation_figure <- ggplot(
  ratio_long, aes(library_version, value, fill = library_version)
) +
  geom_col(width = 0.58, show.legend = FALSE) +
  geom_text(aes(label = format(value, digits = 3, scientific = TRUE)),
            vjust = -0.35, fontface = "bold") +
  scale_fill_manual(values = c("Authenticated CSV active" = "#999999",
                               "Raw-DAD candidate" = "#0072B2")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "Replicate separation determines the active point estimate",
    subtitle = paste0(
      "Lower = clearer own-compound separation. Curated spectra remain the ",
      "active point estimates."
    ),
    x = NULL, y = "Median minimum / second-nearest distance"
  ) + theme_library +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "library_version_separation_comparison.pdf"),
       separation_figure, width = 8, height = 5.4)

cat("Raw DAD library build complete.\n")
cat("Compounds with admitted raw auxiliary evidence:",
    sum(compound_decisions$raw_auxiliary_evidence_status ==
          "Admitted raw DAD auxiliary evidence"), "of",
    nrow(compound_decisions), "\n")
print(validation_summary, row.names = FALSE)
