## ========================================================================== ##
## COMPREHENSIVE COMPOUND-FIRST, OPEN-SET UV SPECTRAL ANALYSIS
## ========================================================================== ##
## Production entry point: RUN_AGILENT_DAD.R
##
## Sequence
##   1. Fit preprocessing and references using library standards only.
##   2. Compare sample peaks using six complementary spectral distances.
##   3. Calibrate acceptance from held-out library information.
##   4. Measure method, replicate, wavelength, and cluster stability.
##   5. Assign Strong, Moderate, Ambiguous, or Outside library coverage support.
##   6. Recluster ambiguous peaks into discovery-oriented spectral families.
##   7. Write documented CSV tables and consistently styled PDF figures.
##
## "Nearest compound" means closest library analogue, not confirmed identity.
## Retention time is metadata used for the final display only; it is never used
## by a distance, probability, clustering, or acceptance calculation.
##
## Pipeline role and data contract
##   RUN_AGILENT_DAD.R supplies one representative spectrum per aligned feature
##   and invokes this engine once per genus and alignment policy. The engine
##   reads the active raw-DAD library and the validated class/chromophore
##   dictionaries from library/. This stage is invoked by RUN_AGILENT_DAD.R
##   and does not provide a separate historical input route.

source("02_nonlinear_compound_similarity.R")

# The production runner supplies one aligned-feature spectrum table and its
# dataset name. Requiring both values prevents accidental execution through a
# retired input route.
DATASET_FILE <- Sys.getenv("UV_DATASET_FILE", unset = "")
DATASET_NAME <- Sys.getenv("UV_DATASET_NAME", unset = "")
if (!nzchar(DATASET_FILE) || !file.exists(DATASET_FILE)) {
  stop("UV_DATASET_FILE must identify the aligned-feature spectrum table created by RUN_AGILENT_DAD.R.")
}
if (!nzchar(DATASET_NAME)) stop("Could not derive a dataset name from the CSV file.")
DATASET_SLUG <- tolower(gsub("[^A-Za-z0-9]+", "_", DATASET_NAME))
DATASET_DISPLAY_NAME <- gsub(
  "([a-z])([A-Z])", "\\1 \\2", gsub("[_-]+", " ", DATASET_NAME)
)
OUTPUT_DIR_3 <- Sys.getenv("UV_OUTPUT_DIR", unset = DATASET_NAME)
ANALYSIS_INPUT_MODE <- Sys.getenv(
  "UV_ANALYSIS_INPUT_MODE", unset = "submitted_maximum_absorbance_table"
)
# Raw-DAD runs pass the alignment policy through unchanged. Historical routes
# used `not_applicable`. This provenance never enters a
# spectral distance or probability calculation.
SOURCE_ALIGNMENT_POLICY <- Sys.getenv(
  "UV_RAW_ALIGNMENT_POLICY", unset = "not_applicable"
)
SOURCE_ALIGNMENT_POLICY_DISPLAY_NAME <- Sys.getenv(
  "UV_RAW_ALIGNMENT_POLICY_DISPLAY_NAME", unset = "Not applicable"
)
dataset_filename <- function(suffix) paste0(DATASET_SLUG, "_", suffix)

# Consistent display capitalization for library compounds. Only the first
# alphabetic character is changed, preserving the remainder of the supplied
# chemical name. The positional prefix in p-Cresol is handled explicitly.
format_compound_name <- function(x) {
  vapply(x, function(value) {
    value <- trimws(value)
    if (tolower(value) == "p-cresol") return("p-Cresol")
    first_letter <- regexpr("[[:alpha:]]", value)[1]
    if (first_letter > 0) {
      substr(value, first_letter, first_letter) <-
        toupper(substr(value, first_letter, first_letter))
    }
    value
  }, character(1), USE.NAMES = FALSE)
}

# Standardize class labels for every table and legend. Preserve the supplied
# abbreviation NFP while using title case for ordinary multiword names.
format_class_name <- function(x) {
  vapply(x, function(value) {
    value <- trimws(value)
    if (toupper(value) == "NFP") "NFP" else tools::toTitleCase(tolower(value))
  }, character(1), USE.NAMES = FALSE)
}
BOOTSTRAP_REPETITIONS <- 1000
SEED_SENSITIVITY_REPETITIONS <- 250
SEED_SENSITIVITY_SEEDS <- 20260806 + 0:3
CONSENSUS_REPETITIONS <- 200
CONFORMAL_ALPHA <- 0.10
COMPOUND_CONFORMAL_ALPHA <- 0.25
set.seed(20260806)


# -------------------------- DISTANCE DEFINITIONS ---------------------------
# Each distance captures a different aspect of UV shape. All matrices use
# spectra as rows and matched wavelength/features as columns. Metric scales are
# learned from library-to-library comparisons before any sample is evaluated.

# Normalize each row to unit Euclidean length, protecting zero-length rows.
normalize_rows_l2 <- function(x) {
  norms <- sqrt(rowSums(x^2))
  norms[norms <= sqrt(.Machine$double.eps)] <- 1
  x / norms
}


# Angular shape dissimilarity, bounded to [0, 2] after numerical cleanup.
cosine_distances <- function(query, reference) {
  query <- normalize_rows_l2(query)
  reference <- normalize_rows_l2(reference)
  distances <- 1 - tcrossprod(query, reference)
  distances[distances < 0] <- 0
  distances[distances > 2] <- 2
  distances
}


# Correlation distance implemented as cosine distance after row centering.
correlation_distances <- function(query, reference) {
  query_centered <- query - rowMeans(query)
  reference_centered <- reference - rowMeans(reference)
  cosine_distances(query_centered, reference_centered)
}


# Minimum correlation distance across small opposing wavelength shifts. A shift
# of two columns accommodates minor band displacement without free peak warping.
shifted_correlation_distances <- function(query, reference, max_shift = 2) {
  output <- matrix(Inf, nrow(query), nrow(reference))
  for (shift in seq(-max_shift, max_shift)) {
    if (shift < 0) {
      query_columns <- seq_len(ncol(query) + shift)
      reference_columns <- (1 - shift):ncol(reference)
    } else if (shift > 0) {
      query_columns <- (1 + shift):ncol(query)
      reference_columns <- seq_len(ncol(reference) - shift)
    } else {
      query_columns <- reference_columns <- seq_len(ncol(query))
    }
    candidate <- correlation_distances(
      query[, query_columns, drop = FALSE],
      reference[, reference_columns, drop = FALSE]
    )
    output <- pmin(output, candidate)
  }
  output
}


# Extract the wavelength positions of the strongest local maxima. Spectra with
# fewer maxima repeat their final position so every row has fixed length.
top_peak_positions <- function(spectra, wavelengths, number = 3) {
  t(apply(spectra, 1, function(x) {
    local <- which(x > c(-Inf, x[-length(x)]) &
                   x >= c(x[-1], -Inf))
    if (!length(local)) local <- which.max(x)
    selected <- local[order(x[local], decreasing = TRUE)]
    selected <- head(selected, number)
    positions <- sort(wavelengths[selected])
    if (length(positions) < number) {
      positions <- c(positions, rep(positions[length(positions)],
                                    number - length(positions)))
    }
    positions
  }))
}


# Compare corresponding ordered maxima and scale displacement by the full
# wavelength range, producing a dimensionless peak-position distance.
peak_position_distances <- function(query_positions, reference_positions,
                                    wavelength_range) {
  output <- matrix(0, nrow(query_positions), nrow(reference_positions))
  for (i in seq_len(nrow(query_positions))) {
    output[i, ] <- rowMeans(abs(sweep(reference_positions, 2,
                                      query_positions[i, ], "-"))) /
      wavelength_range
  }
  output
}


# ---------------------------- ENSEMBLE MODEL --------------------------------
# Compound medians, metric normalizers, and moderate reliability weights are
# learned from library spectra. Weight estimation uses replicate compactness and
# never optimizes the broad class labels that performed poorly previously.

# Collapse replicate rows independently for each compound using pointwise medians.
median_compound_representation <- function(values, compounds) {
  compound_levels <- sort(unique(compounds))
  output <- t(vapply(compound_levels, function(compound) {
    apply(values[compounds == compound, , drop = FALSE], 2, median)
  }, numeric(ncol(values))))
  rownames(output) <- compound_levels
  output
}


# Estimate a robust scale from median nearest-neighbor separation in a library
# distance matrix; this puts heterogeneous metrics on comparable magnitudes.
metric_scale <- function(distance_matrix) {
  diagonal <- diag(distance_matrix)
  diag(distance_matrix) <- Inf
  value <- median(apply(distance_matrix, 1, min))
  diag(distance_matrix) <- diagonal
  max(value, sqrt(.Machine$double.eps))
}


# Derive unsupervised metric weights from replicate reproducibility. A metric
# receives more weight when each library spectrum lies substantially closer to
# its own compound median than to the nearest different compound. No class label
# or sample spectrum is used. Scores are capped before normalization so one
# nearly perfect metric cannot dominate the ensemble.
metric_reliability_weights <- function(metric_distances, compounds, scales) {
  scores <- vapply(names(metric_distances), function(metric_name) {
    distance <- metric_distances[[metric_name]] / scales[metric_name]
    pair_index <- which(upper.tri(distance), arr.ind = TRUE)
    pair_distance <- distance[pair_index]
    within <- compounds[pair_index[, 1]] == compounds[pair_index[, 2]]
    if (!any(within) || !any(!within)) return(1)
    # A 0.10 normalized-distance noise floor prevents nearly identical technical
    # replicates from giving every metric an effectively infinite ratio. Log
    # compression then keeps weights moderate and interpretable.
    ratio <- median(pair_distance[!within], na.rm = TRUE) /
      (median(pair_distance[within], na.rm = TRUE) + 0.10)
    score <- log1p(ratio)
    if (!is.finite(score) || score <= 0) 1 else score
  }, numeric(1))
  scores / sum(scores)
}


# Combine already normalized metric matrices with optional metric omission.
# This helper is reused by ablation diagnostics so those calculations exactly
# match the primary weighted ensemble apart from the intentionally removed view.
combine_normalized_metrics <- function(normalized, weights,
                                       omit_metric = character()) {
  keep <- setdiff(names(normalized), omit_metric)
  local_weights <- weights[keep]
  if (is.null(local_weights) || any(!is.finite(local_weights)) ||
      sum(local_weights) <= 0) {
    local_weights <- setNames(rep(1 / length(keep), length(keep)), keep)
  } else local_weights <- local_weights / sum(local_weights)
  Reduce(`+`, Map(`*`, normalized[keep], local_weights))
}


# Fit all library-only representatives and six metric normalizers. The returned
# list is the frozen ensemble model subsequently applied to sample spectra.
build_ensemble_model <- function(spectra, wavelet_features, compounds,
                                 wavelengths) {
  derivatives <- t(apply(spectra, 1, diff))
  peaks <- top_peak_positions(spectra, wavelengths)
  reps <- list(
    shape = median_compound_representation(spectra, compounds),
    derivative = median_compound_representation(derivatives, compounds),
    wavelet = median_compound_representation(wavelet_features, compounds),
    peaks = median_compound_representation(peaks, compounds)
  )
  wavelength_range <- diff(range(wavelengths))
  within <- list(
    cosine = cosine_distances(reps$shape, reps$shape),
    correlation = correlation_distances(reps$shape, reps$shape),
    shifted_correlation = shifted_correlation_distances(reps$shape, reps$shape),
    derivative = correlation_distances(reps$derivative, reps$derivative),
    wavelet = sqrt(squared_distances(reps$wavelet, reps$wavelet)),
    peak_position = peak_position_distances(reps$peaks, reps$peaks,
                                            wavelength_range)
  )
  scales <- vapply(within, metric_scale, numeric(1))
  replicate_metrics <- list(
    cosine = cosine_distances(spectra, spectra),
    correlation = correlation_distances(spectra, spectra),
    shifted_correlation = shifted_correlation_distances(spectra, spectra),
    derivative = correlation_distances(derivatives, derivatives),
    wavelet = sqrt(squared_distances(wavelet_features, wavelet_features)),
    peak_position = peak_position_distances(peaks, peaks, wavelength_range)
  )
  weights <- metric_reliability_weights(replicate_metrics, compounds, scales)
  list(representatives = reps, scales = scales, weights = weights,
       wavelengths = wavelengths, wavelength_range = wavelength_range)
}


# Calculate each raw and normalized query-to-compound distance. Normalized
# values are capped at 10 before equal-weight averaging to limit one outlier.
ensemble_distances <- function(spectra, wavelet_features, model) {
  derivatives <- t(apply(spectra, 1, diff))
  peaks <- top_peak_positions(spectra, model$wavelengths)
  metrics <- list(
    cosine = cosine_distances(spectra, model$representatives$shape),
    correlation = correlation_distances(spectra, model$representatives$shape),
    shifted_correlation = shifted_correlation_distances(
      spectra, model$representatives$shape
    ),
    derivative = correlation_distances(
      derivatives, model$representatives$derivative
    ),
    wavelet = sqrt(squared_distances(
      wavelet_features, model$representatives$wavelet
    )),
    peak_position = peak_position_distances(
      peaks, model$representatives$peaks, model$wavelength_range
    )
  )
  normalized <- Map(function(distance, scale) {
    pmin(distance / scale, 10)
  }, metrics, model$scales[names(metrics)])
  ensemble <- combine_normalized_metrics(normalized, model$weights)
  colnames(ensemble) <- rownames(model$representatives$shape)
  list(ensemble = ensemble, metrics = metrics, normalized = normalized)
}


# Convert ensemble distances to relative Gaussian-kernel support. This is a
# ranking aid over the current library, not a calibrated identity posterior.
distances_to_probabilities <- function(distance_matrix, bandwidth) {
  bandwidth <- max(bandwidth, sqrt(.Machine$double.eps))
  log_weight <- -(distance_matrix^2) / (2 * bandwidth^2)
  log_weight <- log_weight - apply(log_weight, 1, max)
  weight <- exp(log_weight)
  weight / rowSums(weight)
}


# -------------------------- CONFORMAL CALIBRATION ---------------------------
# Distance is treated as a nonconformity score. Compound calibration leaves out
# individual replicates; class calibration leaves out entire compounds so class
# acceptance is not merely memorization of a compound represented in training.

# Build compound- and class-specific nonconformity distributions. Compound
# representatives are rebuilt without the scored replicate; class scores cannot
# use the replicate's own compound as the nearest reference.
calibrate_conformal <- function(library_spectra, library_features,
                                compounds, classes, ensemble_model) {
  compound_levels <- rownames(ensemble_model$representatives$shape)
  class_levels <- sort(unique(classes))
  compound_scores <- setNames(vector("list", length(compound_levels)),
                              compound_levels)
  class_scores <- setNames(vector("list", length(class_levels)), class_levels)
  rows <- list()
  result_row <- 1L

  for (i in seq_len(nrow(library_spectra))) {
    held_compound <- compounds[i]
    comparison <- ensemble_distances(
      library_spectra[i, , drop = FALSE],
      library_features[i, , drop = FALSE], ensemble_model
    )$ensemble[1, ]

    # Compound conformity: representative rebuilt without this replicate.
    sibling_rows <- which(compounds == held_compound &
                          seq_len(length(compounds)) != i)
    sibling_shape <- matrix(apply(library_spectra[sibling_rows, , drop = FALSE],
                                  2, median), nrow = 1)
    sibling_wavelet <- matrix(apply(library_features[sibling_rows, , drop = FALSE],
                                    2, median), nrow = 1)
    local_model <- build_ensemble_model(
      rbind(sibling_shape, sibling_shape),
      rbind(sibling_wavelet, sibling_wavelet),
      c("target", "target"), ensemble_model$wavelengths
    )
    # Reuse global scales to avoid zero two-row scale estimates.
    local_model$scales <- ensemble_model$scales
    local_model$weights <- ensemble_model$weights
    compound_score <- ensemble_distances(
      library_spectra[i, , drop = FALSE],
      library_features[i, , drop = FALSE], local_model
    )$ensemble[1, 1]
    compound_scores[[held_compound]] <- c(
      compound_scores[[held_compound]], compound_score
    )

    comparison[held_compound] <- Inf
    observed_class <- classes[i]
    rep_classes <- classes[match(compound_levels, compounds)]
    within_class <- comparison[rep_classes == observed_class]
    class_score <- min(within_class)
    class_scores[[observed_class]] <- c(class_scores[[observed_class]], class_score)
    rows[[result_row]] <- data.frame(
      library_row = i,
      held_out_compound = held_compound,
      observed_class = observed_class,
      compound_nonconformity = compound_score,
      class_nonconformity = class_score
    )
    result_row <- result_row + 1L
  }
  list(compound_scores = compound_scores, class_scores = class_scores,
       validation = do.call(rbind, rows))
}


# Finite-sample smoothed conformal p-values. Larger values mean that a test score
# is more compatible with the corresponding calibration distribution.
conformal_p_values <- function(test_scores, calibration_scores) {
  output <- matrix(NA_real_, nrow(test_scores), ncol(test_scores),
                   dimnames = dimnames(test_scores))
  for (column in seq_len(ncol(test_scores))) {
    calibration <- calibration_scores[[colnames(test_scores)[column]]]
    output[, column] <- vapply(test_scores[, column], function(score) {
      (1 + sum(calibration >= score)) / (length(calibration) + 1)
    }, numeric(1))
  }
  output
}


# Partially pooled compound conformity addresses the coarse four-level p-values
# produced by only three replicates per compound. Scores are first divided by
# their compound-specific median calibration score, then compared with the
# pooled scaled calibration distribution. This retains compound scale while
# borrowing resolution across the fixed library. It remains a diagnostic and
# does not replace the grouped open-set decision.
pooled_scaled_compound_p_values <- function(test_scores, compound_scores) {
  scales <- vapply(compound_scores, function(scores) {
    max(median(scores), sqrt(.Machine$double.eps))
  }, numeric(1))
  pooled_calibration <- unlist(Map(`/`, compound_scores, scales),
                               use.names = FALSE)
  scaled_test <- sweep(test_scores, 2,
                       scales[colnames(test_scores)], "/")
  output <- apply(scaled_test, c(1, 2), function(score) {
    (1 + sum(pooled_calibration >= score)) /
      (length(pooled_calibration) + 1)
  })
  dim(output) <- dim(test_scores)
  dimnames(output) <- dimnames(test_scores)
  output
}


# Withhold complete compounds to mimic unseen natural products. Separation and
# agreement statistics from these folds calibrate the open-set decision.
calibrate_library_domain <- function(library_spectra, library_features,
                                     compounds, wavelengths, global_scales,
                                     global_weights = NULL,
                                     omit_metric = character()) {
  # Each complete compound is treated as an unseen natural product. Its
  # distance to the remaining standards defines the expected range for a new
  # spectrum that is still covered by the spectral library.
  compound_levels <- sort(unique(compounds))
  rows <- lapply(compound_levels, function(held_compound) {
    training <- compounds != held_compound
    held <- compounds == held_compound
    model <- build_ensemble_model(
      library_spectra[training, , drop = FALSE],
      library_features[training, , drop = FALSE],
      compounds[training], wavelengths
    )
    model$scales <- global_scales
    if (!is.null(global_weights)) model$weights <- global_weights
    distance_result <- ensemble_distances(
      library_spectra[held, , drop = FALSE],
      library_features[held, , drop = FALSE], model
    )
    local_ensemble <- if (length(omit_metric)) {
      combine_normalized_metrics(distance_result$normalized,
                                 model$weights, omit_metric)
    } else distance_result$ensemble
    colnames(local_ensemble) <- colnames(distance_result$ensemble)
    sorted <- t(apply(local_ensemble, 1, sort))
    method_neighbors <- sapply(distance_result$metrics, function(distance) {
      colnames(distance_result$ensemble)[max.col(-distance, ties.method = "first")]
    })
    agreement <- apply(method_neighbors, 1, function(x) max(table(x)) / length(x))
    data.frame(
      held_out_compound = held_compound,
      median_minimum_distance = median(sorted[, 1]),
      median_min_to_second_ratio = median(sorted[, 1] / sorted[, 2]),
      median_min_to_mean_other_ratio = median(
        sorted[, 1] / rowMeans(sorted[, -1, drop = FALSE])
      ),
      median_method_compound_agreement = median(agreement)
    )
  })
  do.call(rbind, rows)
}


# Empirical upper-tail compatibility where larger statistics are less library-like.
empirical_upper_tail_p <- function(values, calibration) {
  # Larger distances are less compatible with library coverage.
  vapply(values, function(value) {
    (1 + sum(calibration >= value)) / (length(calibration) + 1)
  }, numeric(1))
}


# Test whether supplied classes occupy tighter spectral neighborhoods than
# permuted labels before considering any class-supervised interpretation.
class_structure_diagnostics <- function(distance_matrix, class_labels,
                                        permutations = 999,
                                        seed = 20260806) {
  # Quantifies whether predefined chemical classes occupy tighter spectral
  # neighborhoods than expected after random relabeling.
  pair_index <- which(upper.tri(distance_matrix), arr.ind = TRUE)
  pair_distance <- distance_matrix[pair_index]
  same_class <- class_labels[pair_index[, 1]] == class_labels[pair_index[, 2]]
  observed_effect <- median(pair_distance[!same_class]) -
    median(pair_distance[same_class])
  set.seed(seed)
  null_effect <- replicate(permutations, {
    shuffled <- sample(class_labels)
    shuffled_same <- shuffled[pair_index[, 1]] == shuffled[pair_index[, 2]]
    median(pair_distance[!shuffled_same]) - median(pair_distance[shuffled_same])
  })
  nearest <- max.col(-distance_matrix, ties.method = "first")
  data.frame(
    statistic = c("median_within_distance", "median_between_distance",
                  "between_minus_within", "within_to_between_ratio",
                  "nearest_neighbor_same_class_fraction",
                  "permutation_p_value"),
    value = c(median(pair_distance[same_class]),
              median(pair_distance[!same_class]), observed_effect,
              median(pair_distance[same_class]) /
                median(pair_distance[!same_class]),
              mean(class_labels[nearest] == class_labels),
              (1 + sum(null_effect >= observed_effect)) /
                (permutations + 1))
  )
}


# ------------------------ STABILITY AND CLUSTERING --------------------------
# Contiguous wavelength omission asks whether conclusions depend on one spectral
# region. Consensus clustering resamples features and records co-clustering
# frequencies rather than relying on a single dendrogram cut.

# Omit one contiguous wavelength region and sample one replicate per compound in
# every iteration. Frequencies quantify neighbor stability and reveal matches
# dominated by one region of the spectrum.
wavelength_block_bootstrap <- function(query_spectra, library_spectra,
                                       compounds, classes, wavelengths,
                                       repetitions = 250, seed = 20260806) {
  set.seed(seed)
  breaks <- c(190, 230, 270, 320, 400, 501)
  block <- cut(wavelengths, breaks = breaks, include.lowest = TRUE,
               right = FALSE)
  block_levels <- levels(block)
  compound_levels <- sort(unique(compounds))
  compound_classes <- classes[match(compound_levels, compounds)]
  class_levels <- sort(unique(classes))
  compound_counts <- matrix(0L, nrow(query_spectra), length(compound_levels),
                            dimnames = list(NULL, compound_levels))
  class_counts <- matrix(0L, nrow(query_spectra), length(class_levels),
                         dimnames = list(NULL, class_levels))
  block_selected <- matrix(0L, repetitions, length(block_levels),
                           dimnames = list(NULL, block_levels))
  omission_class_counts <- array(
    0L,
    dim = c(nrow(query_spectra), length(class_levels), length(block_levels)),
    dimnames = list(NULL, class_levels, block_levels)
  )
  omission_counts <- setNames(integer(length(block_levels)), block_levels)

  for (iteration in seq_len(repetitions)) {
    selected_blocks <- sample(block_levels, length(block_levels) - 1)
    omitted_block <- setdiff(block_levels, selected_blocks)
    omission_counts[omitted_block] <- omission_counts[omitted_block] + 1L
    selected_columns <- which(block %in% selected_blocks)
    block_selected[iteration, selected_blocks] <- 1L
    selected_rows <- vapply(compound_levels, function(compound) {
      sample(which(compounds == compound), 1)
    }, integer(1))
    reps <- library_spectra[selected_rows, selected_columns, drop = FALSE]
    rownames(reps) <- compound_levels
    distances <- correlation_distances(
      query_spectra[, selected_columns, drop = FALSE], reps
    )
    nearest <- max.col(-distances, ties.method = "first")
    for (i in seq_len(nrow(query_spectra))) {
      compound_counts[i, nearest[i]] <- compound_counts[i, nearest[i]] + 1L
      nearest_class <- compound_classes[nearest[i]]
      class_counts[i, nearest_class] <- class_counts[i, nearest_class] + 1L
      omission_class_counts[i, nearest_class, omitted_block] <-
        omission_class_counts[i, nearest_class, omitted_block] + 1L
    }
  }
  omission_class_support <- omission_class_counts
  for (block_name in block_levels) {
    omission_class_support[, , block_name] <-
      omission_class_counts[, , block_name] / omission_counts[block_name]
  }
  list(compound_support = compound_counts / repetitions,
       class_support = class_counts / repetitions,
       block_selection = colMeans(block_selected),
       omission_class_support = omission_class_support)
}


# Wilson score interval for a binomial bootstrap frequency. Unlike a symmetric
# normal interval, it remains within [0, 1] near zero or one.
wilson_interval <- function(successes, trials, confidence = 0.95) {
  z <- qnorm(1 - (1 - confidence) / 2)
  estimate <- successes / trials
  denominator <- 1 + z^2 / trials
  center <- (estimate + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(
    estimate * (1 - estimate) / trials + z^2 / (4 * trials^2)
  ) / denominator
  cbind(lower = pmax(0, center - half_width),
        upper = pmin(1, center + half_width))
}


# Repeat shorter bootstrap analyses under independent seeds. This separates
# Monte Carlo instability from genuine spectral ambiguity and records whether
# the modal compound changes across random-number streams.
bootstrap_seed_sensitivity <- function(query_spectra, library_spectra,
                                       compounds, classes, wavelengths,
                                       repetitions, seeds) {
  runs <- lapply(seeds, function(seed) wavelength_block_bootstrap(
    query_spectra, library_spectra, compounds, classes, wavelengths,
    repetitions = repetitions, seed = seed
  ))
  compound_levels <- colnames(runs[[1]]$compound_support)
  long <- do.call(rbind, lapply(seq_along(runs), function(run_index) {
    support <- runs[[run_index]]$compound_support
    modal <- max.col(support, ties.method = "first")
    data.frame(
      peak_row = seq_len(nrow(support)),
      seed = seeds[run_index],
      modal_compound = compound_levels[modal],
      modal_support = support[cbind(seq_len(nrow(support)), modal)],
      stringsAsFactors = FALSE
    )
  }))
  summary <- do.call(rbind, lapply(split(long, long$peak_row), function(group) {
    votes <- sort(table(group$modal_compound), decreasing = TRUE)
    data.frame(
      peak_row = group$peak_row[1],
      seed_modal_compound = names(votes)[1],
      seed_modal_compound_fraction = max(votes) / sum(votes),
      seed_modal_support_mean = mean(group$modal_support),
      seed_modal_support_sd = sd(group$modal_support),
      seed_modal_compound_changed = length(votes) > 1,
      stringsAsFactors = FALSE
    )
  }))
  list(long = long, summary = summary)
}


# Nested grouped validation of the practical open-set rule. The outer fold
# withholds every replicate of one compound. All model scales, weights, domain
# distributions, and thresholds are then re-estimated inside the remaining
# library before the outer spectra are evaluated. Because the true compound is
# absent, this validates library coverage and analogue-family transfer rather
# than exact-compound recovery.
nested_leave_one_compound_validation <- function(
    library_spectra, library_features, compounds, classes, wavelengths,
    bootstrap_repetitions = 250, seed = 20260806) {
  compound_levels <- sort(unique(compounds))
  rows <- vector("list", length(compound_levels))
  for (outer_index in seq_along(compound_levels)) {
    held_compound <- compound_levels[outer_index]
    training <- compounds != held_compound
    held <- !training
    model <- build_ensemble_model(
      library_spectra[training, , drop = FALSE],
      library_features[training, , drop = FALSE],
      compounds[training], wavelengths
    )
    inner_calibration <- calibrate_library_domain(
      library_spectra[training, , drop = FALSE],
      library_features[training, , drop = FALSE],
      compounds[training], wavelengths, model$scales, model$weights
    )
    held_result <- ensemble_distances(
      library_spectra[held, , drop = FALSE],
      library_features[held, , drop = FALSE], model
    )
    sorted_distance <- t(apply(held_result$ensemble, 1, sort))
    minimum_distance <- sorted_distance[, 1]
    ratio <- sorted_distance[, 1] / sorted_distance[, 2]
    method_neighbors <- sapply(held_result$metrics, function(distance) {
      colnames(held_result$ensemble)[max.col(-distance, ties.method = "first")]
    })
    agreement <- apply(method_neighbors, 1,
                       function(x) max(table(x)) / length(x))
    local_bootstrap <- wavelength_block_bootstrap(
      library_spectra[held, , drop = FALSE],
      library_spectra[training, , drop = FALSE],
      compounds[training], classes[training], wavelengths,
      repetitions = bootstrap_repetitions,
      seed = seed + outer_index
    )
    modal_support <- apply(local_bootstrap$compound_support, 1, max)
    nearest_index <- max.col(-held_result$ensemble, ties.method = "first")
    nearest_compound <- colnames(held_result$ensemble)[nearest_index]
    nearest_class <- classes[match(nearest_compound, compounds)]
    held_class <- unique(classes[held])
    domain_p <- empirical_upper_tail_p(
      minimum_distance, inner_calibration$median_minimum_distance
    )
    domain_accepted <- domain_p > CONFORMAL_ALPHA
    strong_ratio <- median(inner_calibration$median_min_to_second_ratio)
    moderate_ratio <- unname(quantile(
      inner_calibration$median_min_to_second_ratio, 0.75
    ))
    agreement_threshold <- median(
      inner_calibration$median_method_compound_agreement
    )
    decision <- ifelse(
      !domain_accepted, "Outside library coverage",
      ifelse(ratio <= strong_ratio & agreement >= agreement_threshold &
               modal_support >= 0.80, "Strong analogue support",
             ifelse(ratio <= moderate_ratio & agreement >= (1 / 3) &
                      modal_support >= 0.60, "Moderate analogue support",
                    "Ambiguous"))
    )
    rows[[outer_index]] <- data.frame(
      held_out_compound = held_compound,
      held_out_class = held_class,
      replicate = seq_len(sum(held)),
      nearest_remaining_compound = nearest_compound,
      nearest_remaining_class = nearest_class,
      nearest_class_matches_held_class = nearest_class == held_class,
      minimum_distance = minimum_distance,
      minimum_to_second_ratio = ratio,
      method_agreement = agreement,
      bootstrap_modal_support = modal_support,
      nested_library_domain_p_value = domain_p,
      nested_decision = decision,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}


# Select k by silhouette width, then resample nonlinear feature columns to
# estimate pairwise co-clustering. The historical function name is retained
# internally, but the implementation is applied to every genus dataset.
betula_consensus_clustering <- function(features, repetitions = 200,
                                        feature_fraction = 0.8,
                                        minimum_cluster_size = 1,
                                        seed = 20260806) {
  original_distance <- dist(features)
  original_tree <- hclust(original_distance, method = "average")
  candidates <- 2:min(20, nrow(features) - 1)
  silhouettes <- vapply(candidates, function(k) {
    mean(cluster::silhouette(cutree(original_tree, k), original_distance)[, 3])
  }, numeric(1))
  eligible <- vapply(candidates, function(k) {
    min(table(cutree(original_tree, k))) >= minimum_cluster_size
  }, logical(1))
  if (!any(eligible)) eligible[1] <- TRUE
  selected_k <- candidates[eligible][which.max(silhouettes[eligible])]

  set.seed(seed)
  co_cluster <- matrix(0, nrow(features), nrow(features))
  for (iteration in seq_len(repetitions)) {
    columns <- sample(seq_len(ncol(features)),
                      ceiling(feature_fraction * ncol(features)))
    groups <- cutree(hclust(dist(features[, columns, drop = FALSE]),
                           method = "average"), k = selected_k)
    co_cluster <- co_cluster + outer(groups, groups, `==`)
  }
  co_cluster <- co_cluster / repetitions
  consensus_tree <- hclust(as.dist(1 - co_cluster), method = "average")
  groups <- cutree(consensus_tree, k = selected_k)
  while (selected_k > 2 && min(table(groups)) < minimum_cluster_size) {
    selected_k <- selected_k - 1L
    groups <- cutree(consensus_tree, k = selected_k)
  }
  cluster_stability <- vapply(sort(unique(groups)), function(group) {
    members <- which(groups == group)
    if (length(members) == 1) return(1)
    mean(co_cluster[members, members][upper.tri(co_cluster[members, members])])
  }, numeric(1))
  list(k = selected_k, silhouette = data.frame(k = candidates,
                                                average_silhouette = silhouettes),
       co_cluster = co_cluster, tree = consensus_tree, groups = groups,
       cluster_stability = cluster_stability)
}


# ---------------------- MIXTURE AND QUALITY CHECKS --------------------------
# Two-reference reconstruction is a diagnostic for possible co-elution, not a
# definitive deconvolution. Quality metrics separate noisy measurements from
# coherent spectra that are genuinely unsupported by the current library.

# Compare each query with all scaled single references and all nonnegative
# two-reference least-squares combinations. Thresholding is performed later
# against improvements observed for known library spectra.
fit_best_single_and_pair <- function(query_spectra, representatives) {
  compounds <- rownames(representatives)
  pairs <- combn(seq_len(nrow(representatives)), 2)
  results <- lapply(seq_len(nrow(query_spectra)), function(i) {
    y <- query_spectra[i, ]
    single_scale <- as.numeric(representatives %*% y) /
      rowSums(representatives^2)
    single_scale <- pmax(single_scale, 0)
    single_residual <- vapply(seq_len(nrow(representatives)), function(j) {
      sqrt(mean((y - single_scale[j] * representatives[j, ])^2))
    }, numeric(1))
    best_single <- which.min(single_residual)

    best_pair_residual <- Inf
    best_pair <- c(NA_integer_, NA_integer_)
    best_coefficients <- c(NA_real_, NA_real_)
    for (pair_index in seq_len(ncol(pairs))) {
      pair <- pairs[, pair_index]
      design <- t(representatives[pair, , drop = FALSE])
      coefficient <- tryCatch(
        as.numeric(solve(crossprod(design), crossprod(design, y))),
        error = function(e) c(0, 0)
      )
      coefficient <- pmax(coefficient, 0)
      if (all(coefficient == 0)) next
      residual <- sqrt(mean((y - design %*% coefficient)^2))
      if (residual < best_pair_residual) {
        best_pair_residual <- residual
        best_pair <- pair
        best_coefficients <- coefficient
      }
    }
    coefficient_fraction <- min(best_coefficients) / sum(best_coefficients)
    improvement <- if (single_residual[best_single] > sqrt(.Machine$double.eps)) {
      (single_residual[best_single] - best_pair_residual) /
        single_residual[best_single]
    } else {
      0
    }
    data.frame(
      best_single_compound = compounds[best_single],
      best_single_rmse = single_residual[best_single],
      mixture_compound_1 = compounds[best_pair[1]],
      mixture_compound_2 = compounds[best_pair[2]],
      mixture_coefficient_1 = best_coefficients[1],
      mixture_coefficient_2 = best_coefficients[2],
      mixture_minor_fraction = coefficient_fraction,
      best_pair_rmse = best_pair_residual,
      mixture_relative_improvement = improvement,
      possible_coelution = improvement >= 0.20 & coefficient_fraction >= 0.15
    )
  })
  do.call(rbind, results)
}


# Quantify signal magnitude, roughness, smoothing discrepancy, local maxima, and
# long-wavelength contribution without using any library annotation.
spectral_quality_metrics <- function(raw_spectra, smoothed_spectra,
                                     wavelengths) {
  normalized_raw <- t(apply(raw_spectra, 1, function(x) {
    span <- max(x) - min(x)
    (x - min(x)) / span
  }))
  data.frame(
    maximum_raw_absorbance = apply(raw_spectra, 1, max),
    raw_dynamic_range = apply(raw_spectra, 1, function(x) max(x) - min(x)),
    raw_roughness = apply(normalized_raw, 1,
                          function(x) sd(diff(x, differences = 2))),
    smoothing_rmse = sqrt(rowMeans((normalized_raw - smoothed_spectra)^2)),
    number_local_maxima = apply(smoothed_spectra, 1, function(x) {
      sum(x > c(-Inf, x[-length(x)]) & x >= c(x[-1], -Inf) & x >= 0.05)
    }),
    fraction_signal_above_400nm = if (any(wavelengths >= 400)) {
      rowSums(smoothed_spectra[, wavelengths >= 400, drop = FALSE]) /
        rowSums(smoothed_spectra)
    } else 0
  )
}


dir.create(OUTPUT_DIR_3, showWarnings = FALSE, recursive = TRUE)

# --------------------------- LIBRARY-ONLY FIT -------------------------------
# Everything through calibration is learned from authenticated standards. No
# sample spectrum contributes to representatives, feature scaling, bandwidths,
# conformal distributions, or open-set thresholds.
library_bundle <- load_library_inputs()
library_data <- library_bundle$spectra
library_data$Compound <- format_compound_name(library_data$Compound)
library_data$Short_Class <- format_class_name(library_data$Short_Class)
library_data$Known_Class <- format_class_name(library_data$Known_Class)
library_spec <- spectral_columns(library_data, WAVELENGTH_MIN, WAVELENGTH_MAX)
library_raw <- as.matrix(library_data[, library_spec$columns, drop = FALSE])
storage.mode(library_raw) <- "double"
library_spectra <- preprocess_spectra(library_data, library_spec$columns,
                                     SMOOTHING_WINDOW, SMOOTHING_DEGREE)
library_wavelet_raw <- extract_haar_features(library_spectra)
scaler <- fit_robust_scaler(library_wavelet_raw)
library_wavelet <- apply_robust_scaler(library_wavelet_raw, scaler)
compounds <- library_data$Compound
classes <- library_data$Short_Class
ensemble_model <- build_ensemble_model(
  library_spectra, library_wavelet, compounds, library_spec$wavelengths
)
conformal_model <- calibrate_conformal(
  library_spectra, library_wavelet, compounds, classes, ensemble_model
)
library_domain_calibration <- calibrate_library_domain(
  library_spectra, library_wavelet, compounds, library_spec$wavelengths,
  ensemble_model$scales, ensemble_model$weights
)

# --------------------------- SAMPLE DATASET ---------------------------------
# `betula_*` object names are retained from the original implementation for
# compatibility, but they contain the aligned-feature dataset selected by the
# raw-DAD runner. The
# first column supplies peak identifiers; parsed retention time remains metadata.
betula_data <- read_spectral_csv(DATASET_FILE)
# The library grid remains fixed. A denser or offset sample grid is aligned in
# raw absorbance space before normalization, smoothing, or any distance model.
# This supports instrument exports such as 3 nm Acer data against the current
# 5 nm library without changing the authenticated standards or their learned
# calibration distributions.
betula_native_spec <- spectral_columns(betula_data, -Inf, Inf)
betula_grid_alignment <- harmonize_wavelength_matrix(
  betula_data[, betula_native_spec$columns, drop = FALSE],
  betula_native_spec$wavelengths,
  library_spec$wavelengths,
  input_id = basename(DATASET_FILE)
)
betula_raw <- betula_grid_alignment$values
betula_spec <- list(
  columns = seq_along(library_spec$wavelengths),
  wavelengths = library_spec$wavelengths
)
betula_spectra <- preprocess_spectral_matrix(
  betula_raw, SMOOTHING_WINDOW, SMOOTHING_DEGREE
)
betula_wavelet <- apply_robust_scaler(extract_haar_features(betula_spectra), scaler)
peak_id <- betula_data[[1]]

ensemble_result <- ensemble_distances(betula_spectra, betula_wavelet, ensemble_model)
representative_ensemble_distance <- ensemble_distances(
  ensemble_model$representatives$shape,
  ensemble_model$representatives$wavelet,
  ensemble_model
)$ensemble
diag(representative_ensemble_distance) <- Inf
class_structure <- class_structure_diagnostics(
  representative_ensemble_distance,
  classes[match(rownames(ensemble_model$representatives$shape), compounds)]
)
ensemble_bandwidth <- median(apply(representative_ensemble_distance, 1, min))
ensemble_probability <- distances_to_probabilities(
  ensemble_result$ensemble, ensemble_bandwidth
)
colnames(ensemble_probability) <- colnames(ensemble_result$ensemble)
compound_levels <- colnames(ensemble_probability)
compound_classes <- classes[match(compound_levels, compounds)]
class_probability <- aggregate_class_probabilities(
  ensemble_probability, compound_classes
)

# Leave-one-metric-out sensitivity. Each ablated ensemble is renormalized over
# the remaining reliability weights and receives its own leave-one-compound-out
# domain calibration, rather than borrowing the full-model distance threshold.
metric_ablation_results <- do.call(rbind, lapply(
  names(ensemble_result$normalized), function(omitted_metric) {
    ablated_distance <- combine_normalized_metrics(
      ensemble_result$normalized, ensemble_model$weights, omitted_metric
    )
    colnames(ablated_distance) <- compound_levels
    ablated_nearest_index <- max.col(-ablated_distance, ties.method = "first")
    ablated_calibration <- calibrate_library_domain(
      library_spectra, library_wavelet, compounds, library_spec$wavelengths,
      ensemble_model$scales, ensemble_model$weights, omitted_metric
    )
    ablated_minimum <- apply(ablated_distance, 1, min)
    data.frame(
      peak_row = seq_len(nrow(ablated_distance)),
      omitted_metric = omitted_metric,
      ablated_nearest_compound = compound_levels[ablated_nearest_index],
      ablated_minimum_distance = ablated_minimum,
      ablated_library_domain_p_value = empirical_upper_tail_p(
        ablated_minimum, ablated_calibration$median_minimum_distance
      ),
      stringsAsFactors = FALSE
    )
  }
))
metric_ablation_results$peak_id <- peak_id[metric_ablation_results$peak_row]

# Method agreement records whether conclusions persist when spectral shape is
# represented in six different ways; agreement is evidence of robustness, not
# an additional supervised vote.
method_neighbor <- sapply(ensemble_result$metrics, function(distance) {
  colnames(ensemble_result$ensemble)[max.col(-distance, ties.method = "first")]
})
method_class <- apply(method_neighbor, 2, function(neighbor) {
  compound_classes[match(neighbor, compound_levels)]
})
modal_compound <- apply(method_neighbor, 1, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})
modal_class <- apply(method_class, 1, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})
method_compound_agreement <- apply(method_neighbor, 1, function(x) max(table(x)) / length(x))
method_class_agreement <- apply(method_class, 1, function(x) max(table(x)) / length(x))

# Conformal prediction sets retain every compound/class whose nonconformity is
# compatible with its library calibration distribution at the chosen cutoff.
compound_p <- conformal_p_values(ensemble_result$ensemble,
                                 conformal_model$compound_scores)
pooled_compound_p <- pooled_scaled_compound_p_values(
  ensemble_result$ensemble, conformal_model$compound_scores
)
class_distance <- sapply(sort(unique(classes)), function(class_name) {
  apply(ensemble_result$ensemble[, compound_classes == class_name, drop = FALSE],
        1, min)
})
class_p <- conformal_p_values(class_distance, conformal_model$class_scores)
compound_sets <- apply(compound_p, 1, function(p) {
  selected <- names(p)[p > COMPOUND_CONFORMAL_ALPHA]
  if (length(selected)) paste(selected, collapse = "; ") else "Unsupported"
})
class_sets <- apply(class_p, 1, function(p) {
  selected <- names(p)[p > CONFORMAL_ALPHA]
  if (length(selected)) paste(selected, collapse = "; ") else "Unsupported"
})

# Wavelength-region support and fully unsupervised sample consensus clustering.
block_support <- wavelength_block_bootstrap(
  betula_spectra, library_spectra, compounds, classes,
  library_spec$wavelengths, BOOTSTRAP_REPETITIONS
)
seed_sensitivity <- bootstrap_seed_sensitivity(
  betula_spectra, library_spectra, compounds, classes,
  library_spec$wavelengths, SEED_SENSITIVITY_REPETITIONS,
  SEED_SENSITIVITY_SEEDS
)
seed_sensitivity$long$peak_id <- peak_id[seed_sensitivity$long$peak_row]
seed_sensitivity$summary$peak_id <- peak_id[seed_sensitivity$summary$peak_row]
nested_validation <- nested_leave_one_compound_validation(
  library_spectra, library_wavelet, compounds, classes,
  library_spec$wavelengths,
  bootstrap_repetitions = SEED_SENSITIVITY_REPETITIONS
)
consensus <- betula_consensus_clustering(
  betula_wavelet, CONSENSUS_REPETITIONS
)
dimnames(consensus$co_cluster) <- list(peak_id, peak_id)
fine_candidates <- consensus$silhouette[
  consensus$silhouette$k >= 4 & consensus$silhouette$k <= 12, , drop = FALSE
]
fine_k <- fine_candidates$k[which.max(fine_candidates$average_silhouette)]
fine_groups <- cutree(consensus$tree, k = fine_k)
fine_stability <- vapply(sort(unique(fine_groups)), function(group) {
  members <- which(fine_groups == group)
  if (length(members) == 1) return(1)
  mean(consensus$co_cluster[members, members][
    upper.tri(consensus$co_cluster[members, members])
  ])
}, numeric(1))

# Mixture/co-elution and spectral-quality diagnostics.
shape_representatives <- ensemble_model$representatives$shape
mixture_results <- fit_best_single_and_pair(betula_spectra, shape_representatives)
library_mixture_reference <- fit_best_single_and_pair(
  library_spectra, shape_representatives
)
mixture_improvement_threshold <- max(
  0.20,
  unname(quantile(library_mixture_reference$mixture_relative_improvement, 0.95,
                  na.rm = TRUE))
)
mixture_results$possible_coelution <-
  mixture_results$mixture_relative_improvement >= mixture_improvement_threshold &
  mixture_results$mixture_minor_fraction >= 0.15
quality <- spectral_quality_metrics(betula_raw, betula_spectra,
                                    betula_spec$wavelengths)

best_compound_index <- max.col(ensemble_probability, ties.method = "first")
best_class_index <- max.col(class_probability, ties.method = "first")
best_compound_probability <- ensemble_probability[
  cbind(seq_len(nrow(ensemble_probability)), best_compound_index)
]
bootstrap_modal_compound_index <- max.col(
  block_support$compound_support, ties.method = "first"
)
bootstrap_modal_compound <- compound_levels[bootstrap_modal_compound_index]
bootstrap_modal_compound_support <- block_support$compound_support[
  cbind(seq_len(nrow(block_support$compound_support)),
        bootstrap_modal_compound_index)
]
full_nearest_compound_support <- block_support$compound_support[
  cbind(seq_len(nrow(block_support$compound_support)), best_compound_index)
]
full_nearest_support_interval <- wilson_interval(
  round(full_nearest_compound_support * BOOTSTRAP_REPETITIONS),
  BOOTSTRAP_REPETITIONS
)
bootstrap_modal_support_interval <- wilson_interval(
  round(bootstrap_modal_compound_support * BOOTSTRAP_REPETITIONS),
  BOOTSTRAP_REPETITIONS
)
full_nearest_matches_bootstrap_modal <-
  compound_levels[best_compound_index] == bootstrap_modal_compound
best_class_probability <- class_probability[
  cbind(seq_len(nrow(class_probability)), best_class_index)
]
sorted_compound_distances <- t(apply(ensemble_result$ensemble, 1, sort))
minimum_compound_distance <- sorted_compound_distances[, 1]
second_compound_distance <- sorted_compound_distances[, 2]
mean_other_compound_distance <-
  (rowSums(ensemble_result$ensemble) - minimum_compound_distance) /
  (ncol(ensemble_result$ensemble) - 1)
minimum_to_mean_other_ratio <- minimum_compound_distance /
  mean_other_compound_distance
minimum_to_second_ratio <- minimum_compound_distance / second_compound_distance
sorted_class_probabilities <- t(apply(class_probability, 1, sort,
                                      decreasing = TRUE))
class_probability_margin <- sorted_class_probabilities[, 1] -
  sorted_class_probabilities[, 2]
best_class_conformal_p <- class_p[
  cbind(seq_len(nrow(class_p)), best_class_index)
]
best_class_conformally_accepted <- best_class_conformal_p > CONFORMAL_ALPHA

# Primary decisions now concern spectral-library coverage and analogue
# stability, not recovery of broad chemical classes that do not separate well
# in the library dendrogram.
library_domain_p <- empirical_upper_tail_p(
  minimum_compound_distance,
  library_domain_calibration$median_minimum_distance
)
library_domain_accepted <- library_domain_p > CONFORMAL_ALPHA
metric_ablation_results$full_nearest_compound <- rep(
  compound_levels[best_compound_index],
  times = length(names(ensemble_result$normalized))
)
metric_ablation_results$nearest_compound_changed <-
  metric_ablation_results$ablated_nearest_compound !=
  metric_ablation_results$full_nearest_compound
metric_ablation_results$full_library_domain_accepted <- rep(
  library_domain_accepted,
  times = length(names(ensemble_result$normalized))
)
metric_ablation_results$ablated_library_domain_accepted <-
  metric_ablation_results$ablated_library_domain_p_value > CONFORMAL_ALPHA
metric_ablation_results$library_domain_decision_changed <-
  metric_ablation_results$ablated_library_domain_accepted !=
  metric_ablation_results$full_library_domain_accepted
strong_ratio_threshold <- median(
  library_domain_calibration$median_min_to_second_ratio
)
moderate_ratio_threshold <- unname(quantile(
  library_domain_calibration$median_min_to_second_ratio, 0.75
))
agreement_threshold <- median(
  library_domain_calibration$median_method_compound_agreement
)

strong_analogue_support <-
  library_domain_accepted &
  minimum_to_second_ratio <= strong_ratio_threshold &
  method_compound_agreement >= agreement_threshold &
  full_nearest_matches_bootstrap_modal &
  full_nearest_compound_support >= 0.80
moderate_analogue_support <-
  !strong_analogue_support &
  library_domain_accepted &
  minimum_to_second_ratio <= moderate_ratio_threshold &
  method_compound_agreement >= (1 / 3) &
  full_nearest_matches_bootstrap_modal &
  full_nearest_compound_support >= 0.60
bootstrap_support_threshold_uncertain <-
  (full_nearest_support_interval[, "lower"] < 0.80 &
     full_nearest_support_interval[, "upper"] >= 0.80) |
  (full_nearest_support_interval[, "lower"] < 0.60 &
     full_nearest_support_interval[, "upper"] >= 0.60)
practical_spectral_decision <- ifelse(
  !library_domain_accepted, "Outside library coverage",
  ifelse(strong_analogue_support, "Strong analogue support",
         ifelse(moderate_analogue_support, "Moderate analogue support",
                "Ambiguous"))
)

# Compound-first interpretation: class is inherited from the selected analogue,
# not predicted independently. For ambiguous peaks, the top-three analogue set
# is inspected to determine whether it at least points toward one shared class.
ordered_analogue_indices <- t(apply(ensemble_result$ensemble, 1, order))
top_analogue_indices <- ordered_analogue_indices[, 1:3, drop = FALSE]
top_analogue_names <- matrix(compound_levels[top_analogue_indices],
                             nrow = nrow(top_analogue_indices))
top_analogue_classes <- matrix(compound_classes[top_analogue_indices],
                               nrow = nrow(top_analogue_indices))
top_analogue_class_consistency <- apply(top_analogue_classes, 1, function(x) {
  max(table(x)) / length(x)
})
top_analogue_modal_class <- apply(top_analogue_classes, 1, function(x) {
  names(sort(table(x), decreasing = TRUE))[1]
})
nearest_analogue_class <- compound_classes[best_compound_index]
ambiguous_class_tendency <-
  top_analogue_class_consistency >= (2 / 3) &
  method_class_agreement >= 0.50 &
  apply(block_support$class_support, 1, max) >= 0.70
analogue_class_interpretation <- ifelse(
  practical_spectral_decision == "Outside library coverage",
  "No class inference: spectrum is outside calibrated library coverage",
  ifelse(
    practical_spectral_decision == "Ambiguous",
    ifelse(ambiguous_class_tendency,
           paste0("Ambiguous compound analogue; analogue set tends toward ",
                  top_analogue_modal_class),
           "Ambiguous compound analogue; class unresolved"),
    paste0("Class of nearest library analogue: ", nearest_analogue_class)
  )
)

# ------------------ SPECTRAL EQUIVALENCE AND ANALOGUE SETS -----------------
# Closely spaced standards are grouped at the median nearest-neighbor distance
# of the library. Complete linkage prevents long chains of only locally similar
# compounds from merging into one broad group.
equivalence_distance_threshold <- median(
  apply(representative_ensemble_distance, 1, min)
)
equivalence_distance <- representative_ensemble_distance
diag(equivalence_distance) <- 0
equivalence_tree <- hclust(as.dist(equivalence_distance), method = "complete")
# A fully Bayesian agglomerative model is not used here because the six-part
# nonlinear ensemble distance has no justified parametric likelihood. Treating
# exp(-distance) as a likelihood would produce posterior-looking values without
# valid calibration. Replicate resampling, wavelength-block bootstrapping, and
# consensus co-clustering instead provide empirical support for branches while
# retaining the useful agglomerative hierarchy. A Bayesian model becomes useful
# only after specifying and validating a generative model for spectral curves.
equivalence_groups <- cutree(equivalence_tree,
                             h = equivalence_distance_threshold)
equivalence_group_by_compound <- setNames(equivalence_groups,
                                          compound_levels)
equivalence_members <- vapply(equivalence_groups, function(group) {
  paste(compound_levels[equivalence_groups == group], collapse = "; ")
}, character(1))
equivalence_class_composition <- vapply(equivalence_groups, function(group) {
  member_classes <- compound_classes[equivalence_groups == group]
  paste(names(sort(table(member_classes), decreasing = TRUE)), collapse = "; ")
}, character(1))

# Build every compound lookup and grouped display dictionary from the deployed
# chromophore input. Category counts may be smaller or larger than the LSU
# library and no compound name is hard-coded into the active analysis.
chromophore_mapping <- library_bundle$chromophore_mapping
chromophore_mapping$Compound <- format_compound_name(
  chromophore_mapping$Compound
)
chromophore_lookup <- data.frame(
  compound = chromophore_mapping$Compound,
  chromophore_code = chromophore_mapping$broad_chromophore_code,
  stringsAsFactors = FALSE
)
fine_chromophore_lookup <- data.frame(
  compound = chromophore_mapping$Compound,
  fine_chromophore_code = chromophore_mapping$fine_chromophore_code,
  stringsAsFactors = FALSE
)
chromophore_dictionary <- library_bundle$broad_chromophore_dictionary
fine_chromophore_dictionary <- library_bundle$fine_chromophore_dictionary

# Expanded, human-readable reference vocabulary for future plant datasets.
# Only fine groups supplied for authenticated standards in the deployed
# chromophore dictionary may be assigned automatically. The additional groups
# below are deliberately
# marked "Reference-only": they support literature searches and future library
# expansion, but the present model must not claim them from spectral shape alone.
additional_chromophore_dictionary <- data.frame(
  chromophore_code = c(
    "MET", "ARY", "BNQ", "COU", "CHR", "CHA", "AUR", "STI", "LIG",
    "TAN", "ANT", "ANQ", "BEN", "NAP", "SIN", "BENF", "PYN", "PYE",
    "PUR", "NUC", "POR", "CAR", "BET", "AZO", "ENR", "ENO", "POLY",
    "FUR", "THI", "NIT", "IMI", "PHEA"
  ),
  chromophore_group = c(
    "Methoxy-substituted aromatic", "General substituted aromatic",
    "Benzoquinone or quinone-like conjugation", "Coumarin benzopyrone",
    "Chromone benzopyranone", "Chalcone open-chain conjugated aromatic",
    "Aurone benzofuranone", "Stilbene diarylethene",
    "Lignan or neolignan coupled aromatic", "Condensed or hydrolysable tannin",
    "Anthocyanin or flavylium extended aromatic", "Anthraquinone fused quinone",
    "Benzophenone diaryl carbonyl", "Naphthalene fused aromatic",
    "Indole simple fused heteroaromatic", "Benzofuran fused oxygen heteroaromatic",
    "Pyridine simple heteroaromatic", "Pyrimidine heteroaromatic",
    "Purine fused nitrogen heteroaromatic", "Nucleobase-like heteroaromatic",
    "Porphyrin or chlorophyll tetrapyrrole", "Carotenoid long-chain polyene",
    "Betalain conjugated nitrogen pigment", "Azo or azomethine conjugated system",
    "Alpha,beta-unsaturated carbonyl", "Conjugated enone or dienone",
    "General extended polyene", "Furan simple oxygen heteroaromatic",
    "Thiophene sulfur heteroaromatic", "Nitroaromatic charge-transfer system",
    "Imidazole-like nitrogen heteroaromatic",
    "Phenylacetic or weakly conjugated aromatic acid"
  ),
  broad_family = c(
    "Simple aromatic", "Simple aromatic", "Quinone", "Benzopyrone",
    "Benzopyrone", "Extended phenylpropanoid", "Extended phenylpropanoid",
    "Extended phenylpropanoid", "Coupled aromatic", "Polyphenolic",
    "Flavonoid pigment", "Quinone", "Aromatic carbonyl", "Fused aromatic",
    "Nitrogen heteroaromatic", "Oxygen heteroaromatic",
    "Nitrogen heteroaromatic", "Nitrogen heteroaromatic",
    "Nitrogen heteroaromatic", "Nitrogen heteroaromatic", "Tetrapyrrole pigment",
    "Polyene pigment", "Nitrogen pigment", "Conjugated nitrogen system",
    "Conjugated carbonyl", "Conjugated carbonyl", "Polyene", "Oxygen heteroaromatic",
    "Sulfur heteroaromatic", "Substituted aromatic",
    "Nitrogen heteroaromatic", "Simple aromatic acid"
  ),
  qualitative_spectral_cue = c(
    "Aromatic bands shifted by one or more methoxy substituents",
    "Dominant aromatic absorption with limited diagnostic long-wave structure",
    "Conjugated carbonyl system with broader near-UV contribution",
    "Fused aromatic lactone with structured near-UV absorption",
    "Fused aromatic ketone with extended near-UV absorption",
    "Broad extended-conjugation band distinct from closed flavonoid rings",
    "Extended benzofuranone conjugation with a visible-edge tendency",
    "Two aromatic systems coupled through an alkene",
    "Multiple coupled aromatic bands with substitution-dependent structure",
    "Broad overlapping polyphenolic absorption with substantial heterogeneity",
    "Strong extended absorption that may reach the visible region",
    "Fused aromatic quinone with multiple structured bands",
    "Two aromatic rings coupled through a carbonyl",
    "Fused-ring aromatic pattern with structured short- and near-UV bands",
    "Fused benzopyrrole pattern; substitution strongly shifts band shape",
    "Fused benzene-furan pattern with structured aromatic absorption",
    "Predominantly short-UV heteroaromatic absorption",
    "Short-UV nitrogen heteroaromatic absorption",
    "Fused nitrogen heteroaromatic bands with carbonyl-dependent shifts",
    "Compact heteroaromatic pattern dominated by short wavelengths",
    "Multiple intense bands spanning UV and visible wavelengths",
    "Very broad long-wavelength absorption from extensive conjugation",
    "Broad UV-visible absorption from a conjugated nitrogen-containing pigment",
    "Extended conjugation with potentially strong visible-edge absorption",
    "Carbonyl-conjugated alkene producing a broader near-UV shoulder or band",
    "Longer conjugated carbonyl system than a single enone",
    "Broad absorption that shifts to longer wavelengths as conjugation increases",
    "Primarily short-UV heteroaromatic absorption",
    "Short-UV heteroaromatic absorption with sulfur-dependent shifts",
    "Substituted aromatic spectrum with strong electron-withdrawing effects",
    "Compact nitrogen heteroaromatic absorption dominated by short wavelengths",
    "Aromatic absorption with weak separation from other simple aromatic acids"
  ),
  dictionary_status = "Reference-only; requires a matching future standard",
  library_compounds = "None in current library",
  n_library_compounds = 0L,
  stringsAsFactors = FALSE
)
expanded_chromophore_dictionary <- rbind(
  data.frame(
    chromophore_code = fine_chromophore_dictionary$fine_chromophore_code,
    chromophore_group = fine_chromophore_dictionary$fine_chromophore_group,
    broad_family = unname(setNames(
      chromophore_dictionary$chromophore_family,
      chromophore_dictionary$chromophore_code
    )[fine_chromophore_dictionary$broad_chromophore_code]),
    qualitative_spectral_cue = fine_chromophore_dictionary$structural_rationale,
    dictionary_status = "Library-anchored; eligible for automatic assignment",
    library_compounds = fine_chromophore_dictionary$library_compounds,
    n_library_compounds = fine_chromophore_dictionary$n_library_compounds,
    stringsAsFactors = FALSE
  ),
  additional_chromophore_dictionary
)
compound_chromophore_code <- setNames(
  chromophore_lookup$chromophore_code, chromophore_lookup$compound
)[compound_levels]
if (anyNA(compound_chromophore_code)) {
  stop("Every library compound must have a curated chromophore annotation.")
}
compound_fine_chromophore_code <- setNames(
  fine_chromophore_lookup$fine_chromophore_code,
  fine_chromophore_lookup$compound
)[compound_levels]
if (anyNA(compound_fine_chromophore_code)) {
  stop("Every library compound must have a fine chromophore annotation.")
}
# Translate each compound's fine-group code to its descriptive group name.
# Reapply the compound names explicitly: indexing the code lookup otherwise
# leaves the vector named by group code, which breaks later compound-name lookup.
fine_group_by_code <- setNames(
  fine_chromophore_dictionary$fine_chromophore_group,
  fine_chromophore_dictionary$fine_chromophore_code
)
fine_chromophore_group_by_compound <- setNames(
  unname(fine_group_by_code[compound_fine_chromophore_code]),
  names(compound_fine_chromophore_code)
)
metric_ablation_results$full_fine_chromophore_group <- unname(
  fine_chromophore_group_by_compound[
    metric_ablation_results$full_nearest_compound
  ]
)
metric_ablation_results$ablated_fine_chromophore_group <- unname(
  fine_chromophore_group_by_compound[
    metric_ablation_results$ablated_nearest_compound
  ]
)
metric_ablation_results$fine_chromophore_group_changed <-
  metric_ablation_results$full_fine_chromophore_group !=
  metric_ablation_results$ablated_fine_chromophore_group
nested_validation$held_out_fine_chromophore_group <- unname(
  fine_chromophore_group_by_compound[nested_validation$held_out_compound]
)
nested_validation$nearest_fine_chromophore_group <- unname(
  fine_chromophore_group_by_compound[
    nested_validation$nearest_remaining_compound
  ]
)
nested_validation$nearest_chromophore_matches_held <-
  nested_validation$held_out_fine_chromophore_group ==
  nested_validation$nearest_fine_chromophore_group
equivalence_chromophore_code <- vapply(equivalence_groups, function(group) {
  codes <- unique(compound_chromophore_code[equivalence_groups == group])
  if (length(codes) == 1) codes else "DUV"
}, character(1))
equivalence_chromophore_family <- setNames(
  chromophore_dictionary$chromophore_family,
  chromophore_dictionary$chromophore_code
)[equivalence_chromophore_code]
equivalence_annotation_type <- setNames(
  chromophore_dictionary$annotation_type,
  chromophore_dictionary$chromophore_code
)[equivalence_chromophore_code]
equivalence_annotation_confidence <- setNames(
  chromophore_dictionary$annotation_confidence,
  chromophore_dictionary$chromophore_code
)[equivalence_chromophore_code]
compound_equivalence_chromophore_code <- setNames(
  equivalence_chromophore_code, compound_levels
)
compound_equivalence_chromophore_family <- setNames(
  equivalence_chromophore_family, compound_levels
)
library_row_chromophore_code <- compound_equivalence_chromophore_code[compounds]
chromophore_profile_summary <- do.call(rbind, lapply(
  sort(unique(library_row_chromophore_code)), function(code) {
    members <- library_spectra[library_row_chromophore_code == code, , drop = FALSE]
    data.frame(
      chromophore_code = code,
      chromophore_family = unname(setNames(
        chromophore_dictionary$chromophore_family,
        chromophore_dictionary$chromophore_code
      )[code]),
      n_compounds = length(unique(compounds[library_row_chromophore_code == code])),
      n_spectra = nrow(members), wavelength = library_spec$wavelengths,
      q10 = apply(members, 2, quantile, probs = 0.10),
      q25 = apply(members, 2, quantile, probs = 0.25),
      median = apply(members, 2, median),
      q75 = apply(members, 2, quantile, probs = 0.75),
      q90 = apply(members, 2, quantile, probs = 0.90)
    )
  }
))
library_row_fine_chromophore_code <- setNames(
  compound_fine_chromophore_code, compound_levels
)[compounds]
fine_chromophore_profile_summary <- do.call(rbind, lapply(
  sort(unique(library_row_fine_chromophore_code)), function(code) {
    members <- library_spectra[library_row_fine_chromophore_code == code, ,
                               drop = FALSE]
    data.frame(
      fine_chromophore_code = code,
      fine_chromophore_group = unname(setNames(
        fine_chromophore_dictionary$fine_chromophore_group,
        fine_chromophore_dictionary$fine_chromophore_code
      )[code]),
      n_compounds = length(unique(
        compounds[library_row_fine_chromophore_code == code]
      )),
      n_spectra = nrow(members), wavelength = library_spec$wavelengths,
      q10 = apply(members, 2, quantile, probs = 0.10),
      q25 = apply(members, 2, quantile, probs = 0.25),
      median = apply(members, 2, median),
      q75 = apply(members, 2, quantile, probs = 0.75),
      q90 = apply(members, 2, quantile, probs = 0.90)
    )
  }
))

# A sample analogue set contains every compound close enough to the best match
# to be indistinguishable under the 75th-percentile held-out-compound tie ratio.
analogue_set_limit <- minimum_compound_distance / moderate_ratio_threshold
analogue_set_indices <- lapply(seq_len(nrow(ensemble_result$ensemble)), function(i) {
  if (!library_domain_accepted[i]) integer(0) else
    which(ensemble_result$ensemble[i, ] <= analogue_set_limit[i])
})
analogue_set_names <- vapply(analogue_set_indices, function(index) {
  if (length(index)) paste(compound_levels[index], collapse = "; ") else
    "None: outside library coverage"
}, character(1))
analogue_set_classes <- vapply(analogue_set_indices, function(index) {
  if (length(index)) paste(compound_classes[index], collapse = "; ") else
    "Unresolved"
}, character(1))
analogue_set_size <- lengths(analogue_set_indices)
analogue_set_class_consistency <- vapply(analogue_set_indices, function(index) {
  if (length(index)) max(table(compound_classes[index])) / length(index) else NA_real_
}, numeric(1))
analogue_set_modal_class <- vapply(analogue_set_indices, function(index) {
  if (length(index)) names(sort(table(compound_classes[index]), decreasing = TRUE))[1]
  else "Unresolved"
}, character(1))
analogue_set_chromophore_codes <- vapply(analogue_set_indices, function(index) {
  if (length(index)) {
    paste(compound_equivalence_chromophore_code[compound_levels[index]],
          collapse = "; ")
  } else "Unresolved"
}, character(1))
analogue_set_chromophore_consistency <- vapply(
  analogue_set_indices, function(index) {
    if (!length(index)) return(NA_real_)
    codes <- compound_equivalence_chromophore_code[compound_levels[index]]
    max(table(codes)) / length(codes)
  }, numeric(1)
)
analogue_set_modal_chromophore_code <- vapply(
  analogue_set_indices, function(index) {
    if (!length(index)) return("Unresolved")
    codes <- compound_equivalence_chromophore_code[compound_levels[index]]
    names(sort(table(codes), decreasing = TRUE))[1]
  }, character(1)
)
analogue_set_modal_chromophore_family <- ifelse(
  analogue_set_modal_chromophore_code == "Unresolved", "Unresolved",
  setNames(chromophore_dictionary$chromophore_family,
           chromophore_dictionary$chromophore_code)[
             analogue_set_modal_chromophore_code
           ]
)
analogue_set_fine_chromophore_codes <- vapply(
  analogue_set_indices, function(index) {
    if (length(index)) {
      paste(compound_fine_chromophore_code[compound_levels[index]],
            collapse = "; ")
    } else "Unresolved"
  }, character(1)
)
analogue_set_fine_chromophore_consistency <- vapply(
  analogue_set_indices, function(index) {
    if (!length(index)) return(NA_real_)
    codes <- compound_fine_chromophore_code[compound_levels[index]]
    max(table(codes)) / length(codes)
  }, numeric(1)
)
analogue_set_modal_fine_chromophore_code <- vapply(
  analogue_set_indices, function(index) {
    if (!length(index)) return("Unresolved")
    codes <- compound_fine_chromophore_code[compound_levels[index]]
    names(sort(table(codes), decreasing = TRUE))[1]
  }, character(1)
)
analogue_set_modal_fine_chromophore_group <- ifelse(
  analogue_set_modal_fine_chromophore_code == "Unresolved", "Unresolved",
  setNames(fine_chromophore_dictionary$fine_chromophore_group,
           fine_chromophore_dictionary$fine_chromophore_code)[
             analogue_set_modal_fine_chromophore_code
           ]
)

# Split the former generic Ambiguous result into actionable failure modes.
ratio_failure <- minimum_to_second_ratio > moderate_ratio_threshold
agreement_failure <- method_compound_agreement < (1 / 3)
wavelength_failure <- full_nearest_compound_support < 0.60
bootstrap_disagreement <- !full_nearest_matches_bootstrap_modal
ambiguity_subtype <- rep("Not ambiguous", nrow(ensemble_result$ensemble))
ambiguous_rows <- practical_spectral_decision == "Ambiguous"
ambiguity_subtype[ambiguous_rows & mixture_results$possible_coelution] <-
  "Possible mixture/co-elution"
unassigned_ambiguous <- ambiguous_rows & ambiguity_subtype == "Not ambiguous"
ambiguity_subtype[unassigned_ambiguous & bootstrap_disagreement] <-
  "Full-data/bootstrap compound disagreement"
unassigned_ambiguous <- ambiguous_rows & ambiguity_subtype == "Not ambiguous"
ambiguity_subtype[unassigned_ambiguous & ratio_failure &
                    analogue_set_fine_chromophore_consistency >= (2 / 3)] <-
  "Compound tie, fine chromophore-consistent"
ambiguity_subtype[unassigned_ambiguous & ratio_failure &
                    analogue_set_fine_chromophore_consistency < (2 / 3) &
                    analogue_set_chromophore_consistency >= (2 / 3)] <-
  "Compound tie, broad chromophore-consistent"
ambiguity_subtype[unassigned_ambiguous & ratio_failure &
                    analogue_set_chromophore_consistency < (2 / 3)] <-
  "Compound tie, mixed chromophore families"
unassigned_ambiguous <- ambiguous_rows & ambiguity_subtype == "Not ambiguous"
ambiguity_subtype[unassigned_ambiguous & wavelength_failure] <-
  "Wavelength-sensitive"
unassigned_ambiguous <- ambiguous_rows & ambiguity_subtype == "Not ambiguous"
ambiguity_subtype[unassigned_ambiguous & agreement_failure] <-
  "Method-sensitive"
ambiguity_subtype[ambiguous_rows & ambiguity_subtype == "Not ambiguous"] <-
  "Residual spectral ambiguity"

# For tied top-two analogues, identify wavelength positions where their library
# replicate interquartile bands do not overlap. Preference on only those regions
# is supplementary evidence and never overrides the calibrated primary result.
tie_wavelength_evidence <- list()
tie_discriminating_wavelengths <- integer(nrow(ensemble_result$ensemble))
tie_preferred_analogue <- rep(NA_character_, nrow(ensemble_result$ensemble))
tie_wavelength_rmse_margin <- rep(NA_real_, nrow(ensemble_result$ensemble))
tie_rows <- which(ambiguous_rows & ratio_failure)
for (i in tie_rows) {
  first_name <- compound_levels[top_analogue_indices[i, 1]]
  second_name <- compound_levels[top_analogue_indices[i, 2]]
  first_reps <- library_spectra[compounds == first_name, , drop = FALSE]
  second_reps <- library_spectra[compounds == second_name, , drop = FALSE]
  first_median <- apply(first_reps, 2, median)
  second_median <- apply(second_reps, 2, median)
  first_q25 <- apply(first_reps, 2, quantile, probs = 0.25)
  first_q75 <- apply(first_reps, 2, quantile, probs = 0.75)
  second_q25 <- apply(second_reps, 2, quantile, probs = 0.25)
  second_q75 <- apply(second_reps, 2, quantile, probs = 0.75)
  discriminating <- (first_q25 > second_q75 | second_q25 > first_q75) &
    abs(first_median - second_median) >= 0.05
  tie_discriminating_wavelengths[i] <- sum(discriminating)
  if (any(discriminating)) {
    first_rmse <- sqrt(mean((betula_spectra[i, discriminating] -
                              first_median[discriminating])^2))
    second_rmse <- sqrt(mean((betula_spectra[i, discriminating] -
                               second_median[discriminating])^2))
    tie_preferred_analogue[i] <- if (first_rmse <= second_rmse) first_name else second_name
    tie_wavelength_rmse_margin[i] <- abs(first_rmse - second_rmse)
  }
  tie_wavelength_evidence[[length(tie_wavelength_evidence) + 1L]] <- data.frame(
    peak_id = peak_id[i], wavelength = library_spec$wavelengths,
    betula = betula_spectra[i, ], first_analogue = first_name,
    first_median = first_median, first_q25 = first_q25, first_q75 = first_q75,
    second_analogue = second_name, second_median = second_median,
    second_q25 = second_q25, second_q75 = second_q75,
    discriminating = discriminating
  )
}
tie_wavelength_evidence <- if (length(tie_wavelength_evidence)) {
  do.call(rbind, tie_wavelength_evidence)
} else data.frame()
nearest_chromophore_code <- compound_equivalence_chromophore_code[
  compound_levels[best_compound_index]
]
nearest_chromophore_family <- compound_equivalence_chromophore_family[
  compound_levels[best_compound_index]
]
nearest_fine_chromophore_code <- compound_fine_chromophore_code[
  compound_levels[best_compound_index]
]
nearest_fine_chromophore_group <- fine_chromophore_group_by_compound[
  compound_levels[best_compound_index]
]
bootstrap_modal_class <- compound_classes[bootstrap_modal_compound_index]
bootstrap_modal_fine_chromophore_code <- compound_fine_chromophore_code[
  bootstrap_modal_compound
]
bootstrap_modal_fine_chromophore_group <- fine_chromophore_group_by_compound[
  bootstrap_modal_compound
]
bootstrap_modal_chromophore_family <- compound_equivalence_chromophore_family[
  bootstrap_modal_compound
]
stable_bootstrap_alternative <- bootstrap_disagreement &
  bootstrap_modal_compound_support >= 0.60 &
  seed_sensitivity$summary$seed_modal_compound_fraction >= 0.75
chromophore_interpretation <- ifelse(
  practical_spectral_decision == "Outside library coverage",
  "No chromophore-family inference: spectrum is outside calibrated library coverage",
  ifelse(
    practical_spectral_decision == "Ambiguous",
    ifelse(analogue_set_chromophore_consistency >= (2 / 3),
           paste0("Ambiguous compound analogue; chromophore tendency: ",
                  analogue_set_modal_chromophore_family),
           "Ambiguous compound analogue; chromophore family unresolved"),
    paste0("Chromophore family of nearest supported analogue: ",
           nearest_chromophore_family)
  )
)

# Three parallel reusable labeling schemes.
class_only_label <- ifelse(
  practical_spectral_decision == "Outside library coverage", "Outside library coverage",
  ifelse(practical_spectral_decision == "Ambiguous",
         ifelse(bootstrap_disagreement, "Ambiguous",
                ifelse(analogue_set_class_consistency >= (2 / 3),
                       analogue_set_modal_class, "Ambiguous")),
         nearest_analogue_class)
)
chromophore_only_label <- ifelse(
  practical_spectral_decision == "Outside library coverage", "Outside library coverage",
  ifelse(
    practical_spectral_decision == "Ambiguous",
    ifelse(stable_bootstrap_alternative,
           bootstrap_modal_fine_chromophore_group,
           ifelse(bootstrap_disagreement, "Ambiguous",
                  ifelse(analogue_set_fine_chromophore_consistency >= (2 / 3),
                         analogue_set_modal_fine_chromophore_group,
                         ifelse(analogue_set_chromophore_consistency >= (2 / 3),
                                analogue_set_modal_chromophore_family,
                                "Ambiguous")))),
    nearest_fine_chromophore_group
  )
)
hybrid_label_type <- rep("Ambiguous", length(practical_spectral_decision))
hybrid_label_type[practical_spectral_decision == "Outside library coverage"] <-
  "Outside library coverage"
hybrid_label_type[practical_spectral_decision %in%
                    c("Strong analogue support", "Moderate analogue support")] <-
  "Class from nearest analogue"
hybrid_label_type[practical_spectral_decision == "Ambiguous" &
                    stable_bootstrap_alternative] <-
  "Bootstrap chromophore tendency"
eligible_ambiguous <- practical_spectral_decision == "Ambiguous" &
  !bootstrap_disagreement
hybrid_label_type[eligible_ambiguous &
                    analogue_set_class_consistency >= (2 / 3)] <-
  "Class tendency"
hybrid_label_type[eligible_ambiguous &
                    analogue_set_class_consistency < (2 / 3) &
                    analogue_set_fine_chromophore_consistency >= (2 / 3)] <-
  "Fine chromophore fallback"
hybrid_label_type[eligible_ambiguous &
                    analogue_set_class_consistency < (2 / 3) &
                    analogue_set_fine_chromophore_consistency < (2 / 3) &
                    analogue_set_chromophore_consistency >= (2 / 3)] <-
  "Broad chromophore fallback"

hybrid_class_then_chromophore_label <- rep(
  "Ambiguous", length(hybrid_label_type)
)
hybrid_class_then_chromophore_label[
  hybrid_label_type == "Outside library coverage"
] <- "Outside library coverage"
class_rows <- hybrid_label_type %in%
  c("Class from nearest analogue", "Class tendency")
hybrid_class_then_chromophore_label[class_rows] <- class_only_label[class_rows]
bootstrap_chromophore_rows <-
  hybrid_label_type == "Bootstrap chromophore tendency"
hybrid_class_then_chromophore_label[bootstrap_chromophore_rows] <-
  bootstrap_modal_fine_chromophore_group[bootstrap_chromophore_rows]
fine_rows <- hybrid_label_type == "Fine chromophore fallback"
hybrid_class_then_chromophore_label[fine_rows] <-
  analogue_set_modal_fine_chromophore_group[fine_rows]
broad_rows <- hybrid_label_type == "Broad chromophore fallback"
hybrid_class_then_chromophore_label[broad_rows] <-
  analogue_set_modal_chromophore_family[broad_rows]

extract_numeric_retention_time <- function(x) {
  # Accept ordinary IDs such as X4.77 and collision-safe IDs such as X4.770_A.
  # Only the first numeric token is metadata; it never enters spectral matching.
  token <- regmatches(x, regexpr("[-+]?[0-9]*\\.?[0-9]+", x))
  suppressWarnings(as.numeric(token))
}

results <- data.frame(
  peak_id = peak_id,
  analysis_input_mode = ANALYSIS_INPUT_MODE,
  source_alignment_policy = SOURCE_ALIGNMENT_POLICY,
  source_alignment_policy_display_name =
    SOURCE_ALIGNMENT_POLICY_DISPLAY_NAME,
  reference_spectrum_source = basename(DATASET_FILE),
  retention_time_metadata_only = extract_numeric_retention_time(peak_id),
  ensemble_min_distance = minimum_compound_distance,
  ensemble_mean_distance_to_other_compounds = mean_other_compound_distance,
  min_to_mean_other_distance_ratio = minimum_to_mean_other_ratio,
  min_to_second_nearest_distance_ratio = minimum_to_second_ratio,
  nearest_compound_relative_separation = 1 - minimum_to_mean_other_ratio,
  ensemble_nearest_compound = compound_levels[best_compound_index],
  nearest_library_spectral_analogue = compound_levels[best_compound_index],
  bootstrap_modal_compound = bootstrap_modal_compound,
  bootstrap_modal_compound_class = bootstrap_modal_class,
  bootstrap_modal_fine_chromophore_code =
    bootstrap_modal_fine_chromophore_code,
  bootstrap_modal_fine_chromophore_group =
    bootstrap_modal_fine_chromophore_group,
  bootstrap_modal_chromophore_family = bootstrap_modal_chromophore_family,
  full_nearest_matches_bootstrap_modal =
    full_nearest_matches_bootstrap_modal,
  full_nearest_compound_bootstrap_support = full_nearest_compound_support,
  full_nearest_compound_support_ci_lower =
    full_nearest_support_interval[, "lower"],
  full_nearest_compound_support_ci_upper =
    full_nearest_support_interval[, "upper"],
  bootstrap_modal_compound_support = bootstrap_modal_compound_support,
  bootstrap_modal_support_ci_lower =
    bootstrap_modal_support_interval[, "lower"],
  bootstrap_modal_support_ci_upper =
    bootstrap_modal_support_interval[, "upper"],
  bootstrap_support_threshold_uncertain =
    bootstrap_support_threshold_uncertain,
  seed_modal_compound = seed_sensitivity$summary$seed_modal_compound,
  seed_modal_compound_fraction =
    seed_sensitivity$summary$seed_modal_compound_fraction,
  seed_modal_support_mean = seed_sensitivity$summary$seed_modal_support_mean,
  seed_modal_support_sd = seed_sensitivity$summary$seed_modal_support_sd,
  seed_modal_compound_changed =
    seed_sensitivity$summary$seed_modal_compound_changed,
  nearest_analogue_library_class = nearest_analogue_class,
  top_three_library_analogues = apply(top_analogue_names, 1, paste,
                                      collapse = "; "),
  top_three_analogue_classes = apply(top_analogue_classes, 1, paste,
                                     collapse = "; "),
  top_three_analogue_class_consistency = top_analogue_class_consistency,
  analogue_class_interpretation = analogue_class_interpretation,
  spectral_equivalence_group = equivalence_group_by_compound[
    compound_levels[best_compound_index]
  ],
  spectral_equivalence_group_members = equivalence_members[
    compound_levels[best_compound_index]
  ],
  spectral_equivalence_group_class_composition = equivalence_class_composition[
    compound_levels[best_compound_index]
  ],
  chromophore_family_code = nearest_chromophore_code,
  chromophore_family = nearest_chromophore_family,
  chromophore_interpretation = chromophore_interpretation,
  fine_chromophore_group_code = nearest_fine_chromophore_code,
  fine_chromophore_group = nearest_fine_chromophore_group,
  class_only_label = class_only_label,
  chromophore_only_label = chromophore_only_label,
  hybrid_label_type = hybrid_label_type,
  hybrid_class_then_chromophore_label = hybrid_class_then_chromophore_label,
  accepted_analogue_set = analogue_set_names,
  analogue_set_size = analogue_set_size,
  analogue_set_classes = analogue_set_classes,
  analogue_set_modal_class = analogue_set_modal_class,
  analogue_set_class_consistency = analogue_set_class_consistency,
  analogue_set_chromophore_codes = analogue_set_chromophore_codes,
  analogue_set_modal_chromophore_code = analogue_set_modal_chromophore_code,
  analogue_set_modal_chromophore_family = analogue_set_modal_chromophore_family,
  analogue_set_chromophore_consistency = analogue_set_chromophore_consistency,
  analogue_set_fine_chromophore_codes = analogue_set_fine_chromophore_codes,
  analogue_set_modal_fine_chromophore_code =
    analogue_set_modal_fine_chromophore_code,
  analogue_set_modal_fine_chromophore_group =
    analogue_set_modal_fine_chromophore_group,
  analogue_set_fine_chromophore_consistency =
    analogue_set_fine_chromophore_consistency,
  ambiguity_subtype = ambiguity_subtype,
  tie_discriminating_wavelength_count = tie_discriminating_wavelengths,
  tie_wavelength_preferred_analogue = tie_preferred_analogue,
  tie_wavelength_rmse_margin = tie_wavelength_rmse_margin,
  compound_identity_expected_in_sample = FALSE,
  compound_result_interpretation = "Relative spectral analogue only",
  ensemble_compound_probability = best_compound_probability,
  ensemble_nearest_class = colnames(class_probability)[best_class_index],
  ensemble_class_probability = best_class_probability,
  descriptive_neighborhood_class = nearest_analogue_class,
  class_balanced_neighborhood_support = best_class_probability,
  predefined_class_used_for_primary_decision = FALSE,
  class_probability_margin = class_probability_margin,
  best_class_conformal_p_value = best_class_conformal_p,
  best_class_conformally_accepted = best_class_conformally_accepted,
  library_domain_p_value = library_domain_p,
  library_domain_accepted = library_domain_accepted,
  practical_spectral_decision = practical_spectral_decision,
  practical_class_decision = practical_spectral_decision,
  practical_decision_rule_version = "v3_exact_bootstrap_aligned",
  method_modal_compound = modal_compound,
  method_compound_agreement = method_compound_agreement,
  method_modal_class = modal_class,
  method_class_agreement = method_class_agreement,
  conformal_compound_set_75 = compound_sets,
  conformal_compound_set_size = rowSums(compound_p > COMPOUND_CONFORMAL_ALPHA),
  pooled_scaled_nearest_compound_p_value = pooled_compound_p[
    cbind(seq_len(nrow(pooled_compound_p)), best_compound_index)
  ],
  pooled_scaled_compound_set_75 = apply(
    pooled_compound_p, 1, function(p) {
      selected <- names(p)[p > COMPOUND_CONFORMAL_ALPHA]
      if (length(selected)) paste(selected, collapse = "; ") else "Unsupported"
    }
  ),
  pooled_scaled_compound_set_size =
    rowSums(pooled_compound_p > COMPOUND_CONFORMAL_ALPHA),
  exact_compound_conformal_diagnostic_used_for_decision = FALSE,
  conformal_class_set_90 = class_sets,
  conformal_class_set_size = rowSums(class_p > CONFORMAL_ALPHA),
  class_conformal_diagnostic_unsupported =
    rowSums(class_p > CONFORMAL_ALPHA) == 0,
  conformal_unsupported = !library_domain_accepted,
  wavelength_block_modal_compound = colnames(block_support$compound_support)[
    max.col(block_support$compound_support)
  ],
  wavelength_block_compound_support = bootstrap_modal_compound_support,
  wavelength_block_modal_class = colnames(block_support$class_support)[
    max.col(block_support$class_support)
  ],
  wavelength_block_class_support = apply(block_support$class_support, 1, max),
  betula_consensus_cluster = consensus$groups,
  betula_consensus_cluster_stability = consensus$cluster_stability[consensus$groups],
  betula_fine_consensus_cluster = fine_groups,
  betula_fine_consensus_cluster_stability = fine_stability[fine_groups],
  mixture_results,
  quality,
  check.names = FALSE
)

# ---------------------- AMBIGUOUS-PEAK DISCOVERY ---------------------------
# These peaks could not be assigned confidently to one library class. We
# therefore cluster them without library labels to discover recurring spectral
# shapes that can be compared with published UV spectra later. Quantile ribbons
# are preferred to standard-error bars because spectra vary across wavelength
# and the within-family distributions need not be Gaussian.
ambiguous_indices <- which(results$practical_class_decision == "Ambiguous")
ambiguous_consensus <- NULL
ambiguous_profile_summary <- data.frame()
ambiguous_cluster_annotation <- data.frame()
ambiguous_representative_maxima <- data.frame()
results$ambiguous_spectral_family <- NA_integer_
results$ambiguous_family_stability <- NA_real_

if (length(ambiguous_indices) >= 3) {
  ambiguous_consensus <- betula_consensus_clustering(
    betula_wavelet[ambiguous_indices, , drop = FALSE],
    repetitions = CONSENSUS_REPETITIONS,
    minimum_cluster_size = 3,
    seed = 20260807
  )
  ambiguous_groups <- ambiguous_consensus$groups
  results$ambiguous_spectral_family[ambiguous_indices] <- ambiguous_groups
  results$ambiguous_family_stability[ambiguous_indices] <-
    ambiguous_consensus$cluster_stability[ambiguous_groups]

  ambiguous_profile_summary <- do.call(rbind, lapply(
    sort(unique(ambiguous_groups)), function(group) {
      family_rows <- ambiguous_indices[ambiguous_groups == group]
      family_spectra <- betula_spectra[family_rows, , drop = FALSE]
      data.frame(
        spectral_family = group,
        n_peaks = length(family_rows),
        wavelength = library_spec$wavelengths,
        q10 = apply(family_spectra, 2, quantile, probs = 0.10),
        q25 = apply(family_spectra, 2, quantile, probs = 0.25),
        median = apply(family_spectra, 2, median),
        q75 = apply(family_spectra, 2, quantile, probs = 0.75),
        q90 = apply(family_spectra, 2, quantile, probs = 0.90)
      )
    }
  ))

  ambiguous_cluster_annotation <- do.call(rbind, lapply(
    sort(unique(ambiguous_groups)), function(group) {
      family_rows <- ambiguous_indices[ambiguous_groups == group]
      class_votes <- sort(table(results$ensemble_nearest_class[family_rows]),
                          decreasing = TRUE)
      compound_votes <- sort(table(results$ensemble_nearest_compound[family_rows]),
                             decreasing = TRUE)
      data.frame(
        spectral_family = group,
        n_peaks = length(family_rows),
        consensus_stability = ambiguous_consensus$cluster_stability[group],
        modal_provisional_class = names(class_votes)[1],
        modal_class_fraction = max(class_votes) / sum(class_votes),
        top_library_analogues = paste(head(names(compound_votes), 3),
                                      collapse = "; "),
        median_minimum_distance = median(results$ensemble_min_distance[family_rows]),
        median_class_probability = median(
          results$ensemble_class_probability[family_rows]
        )
      )
    }
  ))

  family_medians <- do.call(rbind, lapply(
    sort(unique(ambiguous_groups)), function(group) {
      apply(betula_spectra[ambiguous_indices[ambiguous_groups == group], ,
                           drop = FALSE], 2, median)
    }
  ))
  # Ignore the solvent/frontier-dominated extreme UV edge when preparing
  # literature-search maxima. Peaks below 210 nm are rarely shape-specific in
  # this dataset and otherwise dominate the representative spectrum.
  literature_columns <- library_spec$wavelengths >= 210
  maxima <- top_peak_positions(
    family_medians[, literature_columns, drop = FALSE],
    library_spec$wavelengths[literature_columns], number = 3
  )
  ambiguous_representative_maxima <- data.frame(
    spectral_family = sort(unique(ambiguous_groups)),
    n_peaks = as.integer(table(ambiguous_groups)),
    absorption_maximum_1_nm = maxima[, 1],
    absorption_maximum_2_nm = maxima[, 2],
    absorption_maximum_3_nm = maxima[, 3]
  )
}

results$recommended_followup <- ifelse(
  results$practical_spectral_decision == "Strong analogue support",
  "Report nearest compound as a strong spectral analogue and inherit its chromophore family; retain chemical class as metadata",
  ifelse(
    results$practical_spectral_decision == "Moderate analogue support",
    "Report nearest compound and chromophore family provisionally; inspect overlay and bootstrap support",
    ifelse(
      results$practical_spectral_decision == "Ambiguous",
      paste0(results$ambiguity_subtype, ": use chromophore composition of the accepted analogue set, spectral family ",
             ifelse(is.na(results$ambiguous_spectral_family), "unresolved",
                    results$ambiguous_spectral_family),
             ", and representative-spectrum literature search"),
      "Do not inherit a library class; compare spectra outside library coverage with literature"
    )
  )
)
results$recommended_followup[results$possible_coelution] <- paste0(
  results$recommended_followup[results$possible_coelution],
  "; also inspect possible two-component/co-elution fit"
)

# Label consensus clusters after clustering; labels did not influence clustering.
cluster_annotation <- do.call(rbind, lapply(sort(unique(consensus$groups)), function(group) {
  rows <- consensus$groups == group
  class_votes <- sort(table(results$ensemble_nearest_class[rows]), decreasing = TRUE)
  data.frame(
    consensus_cluster = group,
    n_peaks = sum(rows),
    stability = consensus$cluster_stability[group],
    modal_library_class = names(class_votes)[1],
    modal_class_fraction = max(class_votes) / sum(class_votes),
    unsupported_fraction = mean(results$conformal_unsupported[rows]),
    possible_coelution_fraction = mean(results$possible_coelution[rows])
  )
}))
fine_cluster_annotation <- do.call(rbind, lapply(sort(unique(fine_groups)), function(group) {
  rows <- fine_groups == group
  class_votes <- sort(table(results$ensemble_nearest_class[rows]), decreasing = TRUE)
  data.frame(
    fine_consensus_cluster = group,
    n_peaks = sum(rows),
    stability = fine_stability[group],
    modal_library_class = names(class_votes)[1],
    modal_class_fraction = max(class_votes) / sum(class_votes),
    unsupported_fraction = mean(results$conformal_unsupported[rows]),
    possible_coelution_fraction = mean(results$possible_coelution[rows])
  )
}))

# ------------------------------ TABLE OUTPUTS -------------------------------
# Write a complete master table first, followed by narrower interpretation,
# calibration, stability, clustering, and dictionary tables. Shared library
# files are repeated in each genus folder so every analysis is self-contained.
write.csv(results, file.path(OUTPUT_DIR_3, dataset_filename("comprehensive_results.csv")),
          row.names = FALSE)

# This compact key is intended for joining the spectral interpretation back to
# ecological, taxonomic, treatment, or other contextual data. It deliberately
# keeps the primary label, the evidence source behind that label, open-set
# support, ambiguity diagnostics, and spectral-quality fields together. The
# complete master table above remains the audit source for specialist details.
peak_annotation_metadata <- data.frame(
  peak_id = results$peak_id,
  analysis_input_mode = results$analysis_input_mode,
  reference_spectrum_source = results$reference_spectrum_source,
  retention_time_metadata_only = results$retention_time_metadata_only,
  annotation_label = results$hybrid_class_then_chromophore_label,
  annotation_basis = results$hybrid_label_type,
  annotation_evidence_tier = results$practical_spectral_decision,
  class_only_label = results$class_only_label,
  chromophore_only_label = results$chromophore_only_label,
  nearest_library_spectral_analogue =
    results$nearest_library_spectral_analogue,
  nearest_analogue_library_class = results$nearest_analogue_library_class,
  nearest_compound_bootstrap_support =
    results$full_nearest_compound_bootstrap_support,
  nearest_compound_support_ci_lower =
    results$full_nearest_compound_support_ci_lower,
  nearest_compound_support_ci_upper =
    results$full_nearest_compound_support_ci_upper,
  library_domain_p_value = results$library_domain_p_value,
  within_library_coverage = results$library_domain_accepted,
  ensemble_min_distance = results$ensemble_min_distance,
  min_to_second_nearest_distance_ratio =
    results$min_to_second_nearest_distance_ratio,
  min_to_mean_other_distance_ratio =
    results$min_to_mean_other_distance_ratio,
  method_compound_agreement = results$method_compound_agreement,
  wavelength_block_compound_support =
    results$wavelength_block_compound_support,
  accepted_analogue_set = results$accepted_analogue_set,
  analogue_set_size = results$analogue_set_size,
  analogue_set_classes = results$analogue_set_classes,
  analogue_set_class_consistency = results$analogue_set_class_consistency,
  analogue_set_modal_fine_chromophore_group =
    results$analogue_set_modal_fine_chromophore_group,
  analogue_set_fine_chromophore_consistency =
    results$analogue_set_fine_chromophore_consistency,
  ambiguity_subtype = results$ambiguity_subtype,
  ambiguous_spectral_family = results$ambiguous_spectral_family,
  ambiguous_family_stability = results$ambiguous_family_stability,
  possible_coelution = results$possible_coelution,
  maximum_raw_absorbance = results$maximum_raw_absorbance,
  raw_dynamic_range = results$raw_dynamic_range,
  raw_roughness = results$raw_roughness,
  smoothing_rmse = results$smoothing_rmse,
  number_local_maxima = results$number_local_maxima,
  fraction_signal_above_400nm = results$fraction_signal_above_400nm,
  recommended_followup = results$recommended_followup,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(
  peak_annotation_metadata,
  file.path(OUTPUT_DIR_3, dataset_filename("peak_annotation_metadata.csv")),
  row.names = FALSE, na = "NA"
)

# Preserve every value and row from the submitted or landscape-derived
# reference-spectrum table, normalize only the identifier header to `peak_id`, and append
# the compact annotation key. A strict one-to-one match prevents silent row
# recycling or accidental joins by retention time.
original_peak_data <- betula_data
names(original_peak_data)[1] <- "peak_id"
annotation_match <- match(original_peak_data$peak_id,
                          peak_annotation_metadata$peak_id)
if (anyNA(annotation_match) || anyDuplicated(original_peak_data$peak_id) ||
    anyDuplicated(peak_annotation_metadata$peak_id)) {
  stop("Could not make a unique one-to-one original-data annotation join.")
}
original_data_with_peak_annotations <- cbind(
  original_peak_data,
  peak_annotation_metadata[
    annotation_match,
    setdiff(names(peak_annotation_metadata), "peak_id"),
    drop = FALSE
  ]
)
write.csv(
  original_data_with_peak_annotations,
  file.path(
    OUTPUT_DIR_3,
    dataset_filename("original_data_with_peak_annotations.csv")
  ),
  row.names = FALSE, na = "NA"
)

compound_first_output <- results[, c(
  "peak_id", "analysis_input_mode", "reference_spectrum_source",
  "retention_time_metadata_only", "practical_spectral_decision",
  "nearest_library_spectral_analogue", "nearest_analogue_library_class",
  "bootstrap_modal_compound", "full_nearest_matches_bootstrap_modal",
  "full_nearest_compound_bootstrap_support",
  "full_nearest_compound_support_ci_lower",
  "full_nearest_compound_support_ci_upper",
  "bootstrap_modal_compound_support", "seed_modal_compound_fraction",
  "seed_modal_compound_changed",
  "analogue_class_interpretation", "library_domain_p_value",
  "ensemble_min_distance", "min_to_second_nearest_distance_ratio",
  "min_to_mean_other_distance_ratio", "method_compound_agreement",
  "wavelength_block_compound_support", "top_three_library_analogues",
  "top_three_analogue_classes", "top_three_analogue_class_consistency",
  "spectral_equivalence_group", "spectral_equivalence_group_members",
  "chromophore_family_code", "chromophore_family",
  "chromophore_interpretation", "fine_chromophore_group_code",
  "fine_chromophore_group", "class_only_label", "chromophore_only_label",
  "hybrid_label_type", "hybrid_class_then_chromophore_label",
  "accepted_analogue_set", "analogue_set_size", "analogue_set_classes",
  "analogue_set_class_consistency", "analogue_set_chromophore_codes",
  "analogue_set_modal_chromophore_family",
  "analogue_set_chromophore_consistency",
  "analogue_set_fine_chromophore_codes",
  "analogue_set_modal_fine_chromophore_group",
  "analogue_set_fine_chromophore_consistency", "ambiguity_subtype",
  "tie_discriminating_wavelength_count", "tie_wavelength_preferred_analogue",
  "tie_wavelength_rmse_margin",
  "ambiguous_spectral_family", "ambiguous_family_stability",
  "possible_coelution", "recommended_followup"
)]
write.csv(compound_first_output,
          file.path(OUTPUT_DIR_3, dataset_filename("compound_first_interpretation.csv")),
          row.names = FALSE)
parallel_label_output <- results[, c(
  "peak_id", "practical_spectral_decision",
  "nearest_library_spectral_analogue", "nearest_analogue_library_class",
  "fine_chromophore_group_code", "fine_chromophore_group",
  "chromophore_family_code", "chromophore_family",
  "class_only_label", "chromophore_only_label", "hybrid_label_type",
  "hybrid_class_then_chromophore_label", "ambiguity_subtype",
  "analogue_set_class_consistency",
  "analogue_set_fine_chromophore_consistency",
  "analogue_set_chromophore_consistency"
)]
write.csv(parallel_label_output,
          file.path(OUTPUT_DIR_3, dataset_filename("parallel_label_interpretations.csv")),
          row.names = FALSE)
write.csv(as.data.frame(table(results$hybrid_label_type),
                        responseName = "n_peaks"),
          file.path(OUTPUT_DIR_3, "hybrid_label_type_summary.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(results$practical_class_decision),
                        responseName = "n_peaks"),
          file.path(OUTPUT_DIR_3, "practical_class_decision_summary.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(results$practical_spectral_decision),
                        responseName = "n_peaks"),
          file.path(OUTPUT_DIR_3, "practical_spectral_decision_summary.csv"),
          row.names = FALSE)
write.csv(library_domain_calibration,
          file.path(OUTPUT_DIR_3,
                    "library_leave_one_compound_out_domain_calibration.csv"),
          row.names = FALSE)
write.csv(data.frame(metric = names(ensemble_model$weights),
                     reliability_weight = unname(ensemble_model$weights)),
          file.path(OUTPUT_DIR_3, "library_metric_reliability_weights.csv"),
          row.names = FALSE)
write.csv(metric_ablation_results,
          file.path(OUTPUT_DIR_3,
                    dataset_filename("metric_ablation_sensitivity.csv")),
          row.names = FALSE)
write.csv(seed_sensitivity$long[, c(
  "peak_id", "seed", "modal_compound", "modal_support"
)], file.path(OUTPUT_DIR_3,
             dataset_filename("bootstrap_seed_sensitivity.csv")),
row.names = FALSE)
write.csv(nested_validation,
          file.path(OUTPUT_DIR_3,
                    "library_nested_leave_one_compound_validation.csv"),
          row.names = FALSE)
write.csv(class_structure,
          file.path(OUTPUT_DIR_3, "library_class_structure_diagnostics.csv"),
          row.names = FALSE)
write.csv(
  library_bundle$input_manifest,
  file.path(OUTPUT_DIR_3, "library_input_validation.csv"),
  row.names = FALSE
)
write.csv(
  library_bundle$class_dictionary,
  file.path(OUTPUT_DIR_3, "compound_class_dictionary.csv"),
  row.names = FALSE
)
write.csv(
  library_bundle$chromophore_mapping,
  file.path(OUTPUT_DIR_3, "compound_chromophore_dictionary.csv"),
  row.names = FALSE
)
equivalence_output <- data.frame(
  compound = compound_levels,
  library_class = compound_classes,
  spectral_equivalence_group = equivalence_groups,
  group_members = equivalence_members,
  group_class_composition = equivalence_class_composition,
  chromophore_family_code = equivalence_chromophore_code,
  chromophore_family = equivalence_chromophore_family,
  annotation_type = equivalence_annotation_type,
  annotation_confidence = equivalence_annotation_confidence,
  structural_rationale = setNames(
    chromophore_dictionary$structural_rationale,
    chromophore_dictionary$chromophore_code
  )[equivalence_chromophore_code],
  fine_chromophore_group_code = compound_fine_chromophore_code,
  fine_chromophore_group = fine_chromophore_group_by_compound
)
write.csv(equivalence_output,
          file.path(OUTPUT_DIR_3, "library_spectral_equivalence_groups.csv"),
          row.names = FALSE)
write.csv(chromophore_dictionary,
          file.path(OUTPUT_DIR_3, "chromophore_family_dictionary.csv"),
          row.names = FALSE)
write.csv(fine_chromophore_dictionary,
          file.path(OUTPUT_DIR_3, "fine_chromophore_group_dictionary.csv"),
          row.names = FALSE)
write.csv(expanded_chromophore_dictionary,
          file.path(OUTPUT_DIR_3,
                    "expanded_chromophore_reference_dictionary.csv"),
          row.names = FALSE)
write.csv(chromophore_profile_summary,
          file.path(OUTPUT_DIR_3,
                    "library_chromophore_representative_spectra_with_uncertainty.csv"),
          row.names = FALSE)
write.csv(fine_chromophore_profile_summary,
          file.path(OUTPUT_DIR_3,
                    "library_fine_chromophore_representative_spectra_with_uncertainty.csv"),
          row.names = FALSE)
write.csv(as.data.frame(table(results$ambiguity_subtype[
  results$practical_spectral_decision == "Ambiguous"
]), responseName = "n_peaks"),
file.path(OUTPUT_DIR_3, "ambiguous_peak_subtype_summary.csv"),
row.names = FALSE)
if (nrow(tie_wavelength_evidence)) {
  write.csv(tie_wavelength_evidence,
            file.path(OUTPUT_DIR_3, "ambiguous_tie_wavelength_evidence.csv"),
            row.names = FALSE)
}
if (!is.null(ambiguous_consensus)) {
  write.csv(ambiguous_profile_summary,
            file.path(OUTPUT_DIR_3,
                      "ambiguous_family_representative_spectra.csv"),
            row.names = FALSE)
  write.csv(ambiguous_cluster_annotation,
            file.path(OUTPUT_DIR_3, "ambiguous_family_annotations.csv"),
            row.names = FALSE)
  write.csv(ambiguous_representative_maxima,
            file.path(OUTPUT_DIR_3,
                      "ambiguous_family_literature_search_maxima.csv"),
            row.names = FALSE)
  write.csv(data.frame(
    peak_id = peak_id[ambiguous_indices],
    ambiguous_spectral_family = ambiguous_consensus$groups,
    family_stability = ambiguous_consensus$cluster_stability[
      ambiguous_consensus$groups
    ]
  ), file.path(OUTPUT_DIR_3, "ambiguous_peak_family_membership.csv"),
  row.names = FALSE)
}
write.csv(data.frame(peak_id = peak_id, ensemble_probability, check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("ensemble_compound_probabilities.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, class_probability, check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("ensemble_class_probabilities.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, compound_p, check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("conformal_compound_pvalues.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, pooled_compound_p,
                     check.names = FALSE),
          file.path(OUTPUT_DIR_3,
                    dataset_filename("pooled_scaled_compound_pvalues.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, class_p, check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("conformal_class_pvalues.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, method_neighbor, check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("distance_method_neighbors.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, block_support$compound_support,
                     check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("wavelength_block_compound_support.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, block_support$class_support,
                     check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("wavelength_block_class_support.csv")),
          row.names = FALSE)
omission_support_long <- do.call(rbind, lapply(
  dimnames(block_support$omission_class_support)[[3]], function(block_name) {
    block_matrix <- block_support$omission_class_support[, , block_name]
    do.call(rbind, lapply(seq_len(nrow(block_matrix)), function(i) {
      data.frame(peak_id = peak_id[i], omitted_wavelength_block = block_name,
                 class = colnames(block_matrix), support = block_matrix[i, ])
    }))
  }
))
write.csv(omission_support_long,
          file.path(OUTPUT_DIR_3, dataset_filename("wavelength_omission_class_support.csv")),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, consensus$co_cluster,
                     check.names = FALSE),
          file.path(OUTPUT_DIR_3, dataset_filename("consensus_coclustering.csv")),
          row.names = FALSE)
write.csv(cluster_annotation,
          file.path(OUTPUT_DIR_3, dataset_filename("consensus_cluster_annotations.csv")),
          row.names = FALSE)
write.csv(fine_cluster_annotation,
          file.path(OUTPUT_DIR_3, dataset_filename("fine_consensus_cluster_annotations.csv")),
          row.names = FALSE)
write.csv(consensus$silhouette,
          file.path(OUTPUT_DIR_3, dataset_filename("consensus_cluster_number.csv")),
          row.names = FALSE)
write.csv(conformal_model$validation,
          file.path(OUTPUT_DIR_3, "library_conformal_calibration.csv"),
          row.names = FALSE)
write.csv(library_mixture_reference,
          file.path(OUTPUT_DIR_3, "library_mixture_reference.csv"),
          row.names = FALSE)
write.csv(
  betula_grid_alignment$audit,
  file.path(OUTPUT_DIR_3, "wavelength_grid_alignment.csv"),
  row.names = FALSE
)
write.csv(data.frame(
  # Record only the portable filename; the MD5 below establishes exact content
  # without embedding a collaborator-specific absolute filesystem path.
  sample_input_file = basename(DATASET_FILE),
  sample_input_md5 = unname(tools::md5sum(DATASET_FILE)),
  analysis_input_mode = ANALYSIS_INPUT_MODE,
  source_alignment_policy = SOURCE_ALIGNMENT_POLICY,
  source_alignment_policy_display_name =
    SOURCE_ALIGNMENT_POLICY_DISPLAY_NAME,
  native_wavelength_min_nm =
    betula_grid_alignment$audit$native_wavelength_min_nm,
  native_wavelength_max_nm =
    betula_grid_alignment$audit$native_wavelength_max_nm,
  native_wavelength_count =
    betula_grid_alignment$audit$native_wavelength_count,
  analysis_wavelength_min_nm =
    betula_grid_alignment$audit$analysis_wavelength_min_nm,
  analysis_wavelength_max_nm =
    betula_grid_alignment$audit$analysis_wavelength_max_nm,
  analysis_wavelength_count =
    betula_grid_alignment$audit$analysis_wavelength_count,
  wavelength_interpolation_required =
    betula_grid_alignment$audit$interpolation_required,
  wavelength_alignment_method =
    betula_grid_alignment$audit$alignment_method,
  ensemble_bandwidth = ensemble_bandwidth,
  compound_conformal_alpha = COMPOUND_CONFORMAL_ALPHA,
  class_conformal_alpha = CONFORMAL_ALPHA,
  library_domain_alpha = CONFORMAL_ALPHA,
  wavelength_block_bootstrap_repetitions = BOOTSTRAP_REPETITIONS,
  seed_sensitivity_repetitions_per_seed = SEED_SENSITIVITY_REPETITIONS,
  seed_sensitivity_number_of_seeds = length(SEED_SENSITIVITY_SEEDS),
  true_clade_bootstrap_repetitions = SEED_SENSITIVITY_REPETITIONS,
  ensemble_metric_weighting = "library_replicate_reliability",
  ensemble_metric_weights = paste(
    paste(names(ensemble_model$weights),
          round(ensemble_model$weights, 6), sep = "="), collapse = "; "
  ),
  strong_min_to_second_ratio_threshold = strong_ratio_threshold,
  moderate_min_to_second_ratio_threshold = moderate_ratio_threshold,
  strong_method_agreement_threshold = agreement_threshold,
  spectral_equivalence_distance_threshold = equivalence_distance_threshold,
  analogue_set_tie_ratio_threshold = moderate_ratio_threshold,
  mixture_improvement_threshold = mixture_improvement_threshold,
  consensus_k = consensus$k,
  fine_consensus_k = fine_k,
  retention_time_used = FALSE
), file.path(OUTPUT_DIR_3, "analysis_parameters.csv"), row.names = FALSE)

# -------------------------------- FIGURES -----------------------------------
class_levels <- sort(unique(classes))
class_colors <- library_category_colors(class_levels)
chromophore_codes <- chromophore_dictionary$chromophore_code
chromophore_colors <- library_category_colors(
  chromophore_codes,
  # Saturated qualitative colors; none uses gray, which is reserved for
  # ambiguous and outside-library-coverage decisions in chromatograms.
  base_colors = c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#FFD92F", "#A65628", "#F781BF"
  )
)
fine_chromophore_codes <- fine_chromophore_dictionary$fine_chromophore_code
fine_chromophore_colors <- setNames(
  grDevices::hcl.colors(length(fine_chromophore_codes), "Dark 3"),
  fine_chromophore_codes
)

# Exact-compound alignment is now a primary decision requirement. Horizontal
# intervals show Wilson uncertainty around support for the full-data nearest
# compound; shapes reveal whether the bootstrap modal compound agrees.
alignment_plot <- results
alignment_plot$bootstrap_alignment <- ifelse(
  alignment_plot$full_nearest_matches_bootstrap_modal,
  "Full-data and bootstrap agree", "Different bootstrap-modal compound"
)
figure <- ggplot(
  alignment_plot,
  aes(full_nearest_compound_bootstrap_support,
      bootstrap_modal_compound_support,
      color = practical_spectral_decision,
      shape = bootstrap_alignment)
) +
  geom_errorbar(aes(
    xmin = full_nearest_compound_support_ci_lower,
    xmax = full_nearest_compound_support_ci_upper
  ), orientation = "y", width = 0, alpha = 0.35, linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "grey55") +
  geom_vline(xintercept = c(0.60, 0.80), linetype = c(2, 3),
             color = "grey50") +
  geom_point(size = 2.4, alpha = 0.82) +
  scale_color_manual(values = c(
    "Strong analogue support" = "#009E73",
    "Moderate analogue support" = "#0072B2",
    "Ambiguous" = "#999999", "Outside library coverage" = "#4D4D4D"
  ), name = "Decision") +
  scale_shape_manual(values = c(
    "Full-data and bootstrap agree" = 16,
    "Different bootstrap-modal compound" = 4
  ), name = "Compound alignment") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Exact-compound support and bootstrap alignment",
    subtitle = paste0(
      "Intervals are 95% Wilson intervals; vertical lines mark moderate ",
      "and strong exact-support thresholds"
    ),
    x = "Bootstrap support for full-data nearest compound",
    y = "Support for bootstrap-modal compound"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.box = "vertical",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 dataset_filename("compound_bootstrap_alignment.pdf")),
       figure, width = 10, height = 7.5)

metric_weight_data <- data.frame(
  metric = factor(names(ensemble_model$weights),
                  levels = rev(names(sort(ensemble_model$weights)))),
  reliability_weight = unname(ensemble_model$weights)
)
figure <- ggplot(metric_weight_data,
                 aes(reliability_weight, metric)) +
  geom_col(fill = "#0072B2", width = 0.72) +
  geom_text(aes(label = sprintf("%.2f", reliability_weight)),
            hjust = -0.18, size = 3.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(title = "Unsupervised reliability weights for spectral distances",
       subtitle = "Weights are learned only from library replicate reproducibility",
       x = "Normalized ensemble weight", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "library_metric_reliability_weights.pdf"),
       figure, width = 8.5, height = 5.5)

ablation_summary <- do.call(rbind, lapply(
  split(metric_ablation_results, metric_ablation_results$omitted_metric),
  function(group) data.frame(
    omitted_metric = group$omitted_metric[1],
    outcome = c("Nearest compound changed", "Fine chromophore changed",
                "Library-domain decision changed"),
    fraction = c(mean(group$nearest_compound_changed),
                 mean(group$fine_chromophore_group_changed),
                 mean(group$library_domain_decision_changed))
  )
))
ablation_summary$omitted_metric <- factor(
  ablation_summary$omitted_metric,
  levels = c("cosine", "correlation", "shifted_correlation", "derivative",
             "wavelet", "peak_position"),
  labels = c("Cosine", "Correlation", "Shifted correlation", "Derivative",
             "Wavelet", "Peak position")
)
figure <- ggplot(ablation_summary,
                 aes(omitted_metric, fraction, fill = outcome)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.70) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  scale_fill_manual(values = c(
    "Nearest compound changed" = "#0072B2",
    "Fine chromophore changed" = "#E69F00",
    "Library-domain decision changed" = "#CC79A7"
  ), name = NULL) +
  labs(title = "Leave-one-metric-out sensitivity",
       subtitle = "Each ensemble is renormalized and library-domain calibrated without one metric",
       x = "Omitted spectral distance", y = "Fraction of sample peaks") +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        axis.text.x = element_text(angle = 25, hjust = 1),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 dataset_filename("metric_ablation_sensitivity.pdf")),
       figure, width = 10, height = 7)

nested_summary <- data.frame(
  measure = c("Accepted as library-covered", "Nearest class matches",
              "Nearest fine chromophore matches"),
  fraction = c(mean(nested_validation$nested_decision != "Outside library coverage"),
               mean(nested_validation$nearest_class_matches_held_class),
               mean(nested_validation$nearest_chromophore_matches_held))
)
figure <- ggplot(nested_summary,
                 aes(fraction, reorder(measure, fraction))) +
  geom_col(fill = "#009E73", width = 0.68) +
  geom_text(aes(label = scales::percent(fraction, accuracy = 1)),
            hjust = -0.18, size = 3.5) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.08), expand = c(0, 0)) +
  labs(title = "Nested leave-one-compound-out validation",
       subtitle = "All thresholds are re-estimated inside each outer library fold",
       x = "Fraction of held-out library replicates", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 "library_nested_leave_one_compound_validation.pdf"),
       figure, width = 9, height = 5.5)

seed_plot <- seed_sensitivity$summary
figure <- ggplot(seed_plot,
                 aes(seed_modal_support_mean, seed_modal_support_sd,
                     color = seed_modal_compound_changed)) +
  geom_point(size = 2.3, alpha = 0.82) +
  scale_color_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
                     labels = c(`FALSE` = "Stable modal compound",
                                `TRUE` = "Modal compound changed"),
                     name = "Across four seeds") +
  labs(title = "Bootstrap sensitivity to random seed",
       subtitle = paste0("Each seed used ", SEED_SENSITIVITY_REPETITIONS,
                         " wavelength-block/replicate iterations"),
       x = "Mean modal-compound support", y = "Support standard deviation") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 dataset_filename("bootstrap_seed_sensitivity.pdf")),
       figure, width = 9, height = 6.5)

# Compound-level library atlas. Class remains the main visual encoding, while
# the facet strip reports the compound's fine chromophore group. This avoids a
# second 16-color legend and keeps the structural annotation readable. Each
# panel summarizes the library replicates for one compound, so uncertainty is
# displayed consistently with all other representative-spectrum figures.
library_compound_profile_summary <- do.call(rbind, lapply(compound_levels,
  function(compound_name) {
    members <- library_spectra[compounds == compound_name, , drop = FALSE]
    data.frame(
      compound = compound_name,
      library_class = classes[match(compound_name, compounds)],
      fine_chromophore_code = unname(
        compound_fine_chromophore_code[compound_name]
      ),
      fine_chromophore_group = unname(
        fine_chromophore_group_by_compound[compound_name]
      ),
      wavelength = library_spec$wavelengths,
      median = apply(members, 2, median),
      q10 = apply(members, 2, quantile, probs = 0.10),
      q25 = apply(members, 2, quantile, probs = 0.25),
      q75 = apply(members, 2, quantile, probs = 0.75),
      q90 = apply(members, 2, quantile, probs = 0.90),
      stringsAsFactors = FALSE
    )
  }
))
library_compound_profile_summary$facet_label <- paste0(
  library_compound_profile_summary$compound, "\n",
  library_compound_profile_summary$fine_chromophore_code, " - ",
  library_compound_profile_summary$fine_chromophore_group
)
figure <- ggplot(
  library_compound_profile_summary,
  aes(wavelength, median, color = library_class, fill = library_class)
) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.10,
              color = NA, show.legend = FALSE) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25,
              color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.80) +
  facet_wrap(~ facet_label, ncol = 3, scales = "free_y") +
  scale_color_manual(values = class_colors, name = "Library class") +
  scale_fill_manual(values = class_colors) +
  labs(
    title = "Library compound representative spectra",
    subtitle = paste0(
      "Color denotes library class; each panel names its fine chromophore ",
      "group. Median line, IQR dark ribbon, and 10th-90th percentile light ribbon"
    ),
    x = "Wavelength (nm)", y = "Normalized absorbance"
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 6.8),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 1, byrow = TRUE))
ggsave(
  file.path(OUTPUT_DIR_3,
            "library_compound_representative_spectra_by_class.pdf"),
  figure, width = 14,
  height = 3.0 * ceiling(length(compound_levels) / 3)
)

chromophore_profile_summary$facet_label <- paste0(
  chromophore_profile_summary$chromophore_code, " - ",
  chromophore_profile_summary$chromophore_family
)
figure <- ggplot(chromophore_profile_summary,
                 aes(wavelength, median, color = chromophore_code,
                     fill = chromophore_code)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.10,
              color = NA, show.legend = FALSE) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25,
              color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.85, show.legend = FALSE) +
  facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +
  scale_color_manual(values = chromophore_colors) +
  scale_fill_manual(values = chromophore_colors) +
  labs(title = "Library chromophore-family representative spectra",
       subtitle = paste0("Median line; dark ribbon = interquartile range; ",
                         "light ribbon = 10th-90th percentiles"),
       x = "Wavelength (nm)", y = "Normalized absorbance") +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 7.5),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 "library_chromophore_family_representative_spectra.pdf"),
       figure, width = 12,
       height = 3.2 * ceiling(length(unique(chromophore_codes)) / 2))

fine_chromophore_profile_summary$facet_label <- paste0(
  fine_chromophore_profile_summary$fine_chromophore_code, " - ",
  fine_chromophore_profile_summary$fine_chromophore_group
)
figure <- ggplot(fine_chromophore_profile_summary,
                 aes(wavelength, median, color = fine_chromophore_code,
                     fill = fine_chromophore_code)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.10,
              color = NA, show.legend = FALSE) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.25,
              color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.80, show.legend = FALSE) +
  facet_wrap(~ facet_label, ncol = 3, scales = "free_y") +
  scale_color_manual(values = fine_chromophore_colors) +
  scale_fill_manual(values = fine_chromophore_colors) +
  labs(title = "Fine chromophore-group representative spectra",
       subtitle = paste0("Fixed library vocabulary; median line with ",
                         "interquartile and 10th-90th percentile ribbons"),
       x = "Wavelength (nm)", y = "Normalized absorbance") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 6.7),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 "library_fine_chromophore_group_representative_spectra.pdf"),
       figure, width = 14,
       height = 3.0 * ceiling(length(unique(fine_chromophore_codes)) / 3))

# Quantitative class-structure figures. These are diagnostics, not training
# targets: they show directly whether the supplied chemical classes form
# coherent neighborhoods in the compound-representative distance space.
compound_rep_names <- rownames(ensemble_model$representatives$shape)
compound_rep_classes <- classes[match(compound_rep_names, compounds)]
pair_index <- which(upper.tri(representative_ensemble_distance), arr.ind = TRUE)
class_pair_data <- data.frame(
  compound_1 = compound_rep_names[pair_index[, 1]],
  compound_2 = compound_rep_names[pair_index[, 2]],
  relationship = ifelse(
    compound_rep_classes[pair_index[, 1]] ==
      compound_rep_classes[pair_index[, 2]],
    "Within predefined class", "Between predefined classes"
  ),
  ensemble_distance = representative_ensemble_distance[pair_index]
)
write.csv(class_pair_data,
          file.path(OUTPUT_DIR_3, "library_within_between_class_distances.csv"),
          row.names = FALSE)
figure <- ggplot(class_pair_data, aes(relationship, ensemble_distance,
                                      fill = relationship)) +
  geom_violin(alpha = 0.35, trim = FALSE, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.10, alpha = 0.35, size = 1) +
  scale_fill_manual(values = c("Within predefined class" = "#4477AA",
                               "Between predefined classes" = "#999999"),
                    guide = "none") +
  labs(title = "Predefined classes overlap in spectral-distance space",
       subtitle = paste0("Permutation p = ",
                         format(class_structure$value[
                           class_structure$statistic == "permutation_p_value"
                         ], digits = 3)),
       x = NULL, y = "Ensemble distance between compound representatives") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "library_class_distance_overlap.pdf"),
       figure, width = 9, height = 7)

silhouette_distance <- representative_ensemble_distance
diag(silhouette_distance) <- 0
class_silhouette <- cluster::silhouette(
  as.integer(factor(compound_rep_classes)), as.dist(silhouette_distance)
)
silhouette_data <- data.frame(
  compound = compound_rep_names,
  class = compound_rep_classes,
  class_silhouette_width = class_silhouette[, "sil_width"]
)
write.csv(silhouette_data,
          file.path(OUTPUT_DIR_3, "library_compound_class_silhouette.csv"),
          row.names = FALSE)
silhouette_data$compound <- factor(
  silhouette_data$compound,
  levels = silhouette_data$compound[order(silhouette_data$class_silhouette_width)]
)
figure <- ggplot(silhouette_data,
                 aes(class_silhouette_width, compound, fill = class)) +
  geom_col() +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.4) +
  scale_fill_manual(values = class_colors, name = "Predefined class") +
  labs(title = "Compound agreement with predefined spectral classes",
       subtitle = "Negative silhouette widths indicate a closer other-class neighborhood",
       x = "Silhouette width using predefined class labels", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "library_class_silhouette.pdf"),
       figure, width = 10, height = 8)

library_tree <- hclust(as.dist(silhouette_distance), method = "average")
library_label_classes <- setNames(compound_rep_classes, compound_rep_names)
colored_library_tree <- color_dendrogram_labels(
  as.dendrogram(library_tree), library_label_classes, class_colors
)
pdf(file.path(OUTPUT_DIR_3, "library_compound_dendrogram_class_mixing.pdf"),
    width = 13, height = 8)
par(mar = c(9, 4, 4, 2) + 0.1)
plot(colored_library_tree,
     main = "Library compound dendrogram with descriptive class colors",
     ylab = "Ensemble spectral distance", cex = 0.70)
legend("topright", legend = names(class_colors), col = class_colors,
       lty = 1, lwd = 3, bty = "n", ncol = 2, cex = 0.8)
dev.off()

equivalence_colors <- setNames(
  grDevices::hcl.colors(length(unique(equivalence_groups)), "Dark 3"),
  sort(unique(equivalence_groups))
)
equivalence_label_groups <- setNames(as.character(equivalence_groups),
                                     compound_levels)
equivalence_tree_colored <- color_dendrogram_labels(
  as.dendrogram(equivalence_tree), equivalence_label_groups,
  equivalence_colors
)
pdf(file.path(OUTPUT_DIR_3, "library_spectral_equivalence_dendrogram.pdf"),
    width = 14, height = 8)
par(mar = c(9, 4, 4, 2) + 0.1)
plot(equivalence_tree_colored,
     main = "Library spectral-equivalence groups",
     ylab = "Complete-linkage ensemble distance", cex = 0.68)
abline(h = equivalence_distance_threshold, lty = 2, col = "grey45")
legend("topright",
       legend = c(paste("Group", sort(unique(equivalence_groups))),
                  "Grouping threshold"),
       col = c(equivalence_colors, "grey45"),
       lty = c(rep(1, length(equivalence_colors)), 2),
       lwd = c(rep(3, length(equivalence_colors)), 1),
       bty = "n", ncol = 3, cex = 0.65)
dev.off()

compound_chromophore_labels <- setNames(
  compound_equivalence_chromophore_code[compound_levels], compound_levels
)
chromophore_tree_colored <- color_dendrogram_labels(
  as.dendrogram(equivalence_tree), compound_chromophore_labels,
  chromophore_colors
)
pdf(file.path(OUTPUT_DIR_3, "library_chromophore_family_dendrogram.pdf"),
    width = 14, height = 10)
# Keep the relatively long chromophore-family key in its own panel. This
# prevents the legend from obscuring branches while preserving readable labels.
layout(matrix(c(1, 2), nrow = 2), heights = c(4.6, 1.25))
par(mar = c(9, 4, 4, 2) + 0.1)
plot(chromophore_tree_colored,
     main = "Library dendrogram with curated chromophore-family annotations",
     ylab = "Complete-linkage ensemble distance", cex = 0.68)
abline(h = equivalence_distance_threshold, lty = 2, col = "grey45")

# Dedicated legend panel below the tree, including the grouping threshold.
par(mar = c(0.5, 1, 0.5, 1))
plot.new()
legend("center",
       legend = c(
         paste0(chromophore_dictionary$chromophore_code, " - ",
                chromophore_dictionary$chromophore_family),
         "Spectral-equivalence threshold"
       ),
       col = c(
         chromophore_colors[chromophore_dictionary$chromophore_code],
         "grey45"
       ),
       lty = c(rep(1, nrow(chromophore_dictionary)), 2),
       lwd = c(rep(3, nrow(chromophore_dictionary)), 1),
       bty = "n", ncol = 3, cex = 0.70,
       x.intersp = 0.7, y.intersp = 1.05)
layout(1)
dev.off()

# One-at-a-time insertion dendrograms with support for the branch directly
# containing each sample peak. The support is neighborhood-based: it sums the
# block-bootstrap nearest-compound frequencies over every library compound in
# the peak's immediate sibling branch. Consequently, harmless movement among
# spectrally equivalent compounds in that same branch does not lower support.
insertion_branch_leaves <- function(merge_matrix, reference, labels) {
  if (reference < 0) return(labels[-reference])
  c(insertion_branch_leaves(merge_matrix, merge_matrix[reference, 1], labels),
    insertion_branch_leaves(merge_matrix, merge_matrix[reference, 2], labels))
}

insertion_branch_summary <- vector("list", nrow(betula_spectra))
insertion_trees <- vector("list", nrow(betula_spectra))
observed_branch_members <- vector("list", nrow(betula_spectra))
linkage_sensitivity <- list()
linkage_row <- 1L
for (index in seq_len(nrow(betula_spectra))) {
  query_label <- paste0(toupper(DATASET_DISPLAY_NAME), " ", peak_id[index])
  combined_labels <- c(compound_levels, query_label)
  combined_distance <- matrix(
    0, length(combined_labels), length(combined_labels),
    dimnames = list(combined_labels, combined_labels)
  )
  reference_distance <- representative_ensemble_distance
  diag(reference_distance) <- 0
  combined_distance[compound_levels, compound_levels] <- reference_distance
  combined_distance[query_label, compound_levels] <-
    ensemble_result$ensemble[index, compound_levels]
  combined_distance[compound_levels, query_label] <-
    ensemble_result$ensemble[index, compound_levels]

  insertion_tree <- hclust(as.dist(combined_distance), method = "average")
  insertion_trees[[index]] <- insertion_tree
  query_leaf_reference <- -match(query_label, insertion_tree$labels)
  query_merge_row <- which(
    insertion_tree$merge[, 1] == query_leaf_reference |
      insertion_tree$merge[, 2] == query_leaf_reference
  )[1]
  sibling_reference <- if (
    insertion_tree$merge[query_merge_row, 1] == query_leaf_reference
  ) insertion_tree$merge[query_merge_row, 2] else
    insertion_tree$merge[query_merge_row, 1]
  sibling_members <- insertion_branch_leaves(
    insertion_tree$merge, sibling_reference, insertion_tree$labels
  )
  sibling_members <- intersect(sibling_members, compound_levels)
  observed_branch_members[[index]] <- sort(sibling_members)
  branch_support <- sum(
    block_support$compound_support[index, sibling_members], na.rm = TRUE
  )
  nearest_compound <- compound_levels[best_compound_index[index]]
  insertion_branch_summary[[index]] <- data.frame(
    peak_id = peak_id[index],
    immediate_sibling_branch_members = paste(sibling_members, collapse = "; "),
    branch_neighborhood_bootstrap_support = branch_support,
    nearest_compound = nearest_compound,
    nearest_compound_bootstrap_support =
      block_support$compound_support[index, nearest_compound],
    insertion_height = insertion_tree$height[query_merge_row],
    stringsAsFactors = FALSE
  )

  # Repeat the identical insertion distance matrix under alternative linkage
  # rules. Ward linkage is intentionally excluded because the composite
  # ensemble distance is not guaranteed to be squared Euclidean.
  for (linkage in c("average", "complete", "single")) {
    linkage_tree <- hclust(as.dist(combined_distance), method = linkage)
    query_leaf <- -match(query_label, linkage_tree$labels)
    merge_row <- which(linkage_tree$merge[, 1] == query_leaf |
                         linkage_tree$merge[, 2] == query_leaf)[1]
    sibling_reference <- if (
      linkage_tree$merge[merge_row, 1] == query_leaf
    ) linkage_tree$merge[merge_row, 2] else linkage_tree$merge[merge_row, 1]
    members <- sort(intersect(insertion_branch_leaves(
      linkage_tree$merge, sibling_reference, linkage_tree$labels
    ), compound_levels))
    union_members <- union(observed_branch_members[[index]], members)
    linkage_sensitivity[[linkage_row]] <- data.frame(
      peak_id = peak_id[index], linkage = linkage,
      sibling_branch_members = paste(members, collapse = "; "),
      branch_size = length(members),
      jaccard_vs_average = if (length(union_members)) {
        length(intersect(observed_branch_members[[index]], members)) /
          length(union_members)
      } else 1,
      full_nearest_in_branch =
        compound_levels[best_compound_index[index]] %in% members,
      stringsAsFactors = FALSE
    )
    linkage_row <- linkage_row + 1L
  }
}
insertion_branch_summary <- do.call(rbind, insertion_branch_summary)
linkage_sensitivity <- do.call(rbind, linkage_sensitivity)

# True tree-resampling support for the observed insertion branch. Unlike the
# neighborhood statistic above, this reconstructs the full insertion tree in
# every replicate/wavelength-block resample and asks whether the immediate
# sibling clade is reproduced exactly or at least overlaps the observed clade.
set.seed(20260819)
exact_clade_count <- integer(nrow(betula_spectra))
overlap_clade_count <- integer(nrow(betula_spectra))
block_breaks <- c(190, 230, 270, 320, 400, 501)
bootstrap_block <- cut(library_spec$wavelengths, breaks = block_breaks,
                       include.lowest = TRUE, right = FALSE)
for (iteration in seq_len(SEED_SENSITIVITY_REPETITIONS)) {
  selected_blocks <- sample(levels(bootstrap_block),
                            length(levels(bootstrap_block)) - 1)
  selected_columns <- which(bootstrap_block %in% selected_blocks)
  selected_rows <- vapply(compound_levels, function(compound) {
    sample(which(compounds == compound), 1)
  }, integer(1))
  bootstrap_reps <- library_spectra[selected_rows, selected_columns,
                                    drop = FALSE]
  rownames(bootstrap_reps) <- compound_levels
  for (index in seq_len(nrow(betula_spectra))) {
    query_label <- paste0(toupper(DATASET_DISPLAY_NAME), " ", peak_id[index])
    query_row <- matrix(betula_spectra[index, selected_columns], nrow = 1)
    colnames(query_row) <- colnames(bootstrap_reps)
    combined <- rbind(bootstrap_reps, query_row)
    rownames(combined)[nrow(combined)] <- query_label
    tree <- hclust(as.dist(correlation_distances(combined, combined)),
                   method = "average")
    query_leaf <- -match(query_label, tree$labels)
    merge_row <- which(tree$merge[, 1] == query_leaf |
                         tree$merge[, 2] == query_leaf)[1]
    sibling_reference <- if (tree$merge[merge_row, 1] == query_leaf) {
      tree$merge[merge_row, 2]
    } else tree$merge[merge_row, 1]
    members <- sort(intersect(insertion_branch_leaves(
      tree$merge, sibling_reference, tree$labels
    ), compound_levels))
    exact_clade_count[index] <- exact_clade_count[index] +
      as.integer(identical(members, observed_branch_members[[index]]))
    overlap_clade_count[index] <- overlap_clade_count[index] +
      as.integer(length(intersect(members,
                                  observed_branch_members[[index]])) > 0)
  }
}
exact_clade_interval <- wilson_interval(
  exact_clade_count, SEED_SENSITIVITY_REPETITIONS
)
overlap_clade_interval <- wilson_interval(
  overlap_clade_count, SEED_SENSITIVITY_REPETITIONS
)
insertion_branch_summary$exact_sibling_clade_bootstrap_support <-
  exact_clade_count / SEED_SENSITIVITY_REPETITIONS
insertion_branch_summary$exact_sibling_clade_ci_lower <-
  exact_clade_interval[, "lower"]
insertion_branch_summary$exact_sibling_clade_ci_upper <-
  exact_clade_interval[, "upper"]
insertion_branch_summary$overlapping_branch_bootstrap_support <-
  overlap_clade_count / SEED_SENSITIVITY_REPETITIONS
insertion_branch_summary$overlapping_branch_ci_lower <-
  overlap_clade_interval[, "lower"]
insertion_branch_summary$overlapping_branch_ci_upper <-
  overlap_clade_interval[, "upper"]
write.csv(linkage_sensitivity,
          file.path(OUTPUT_DIR_3,
                    dataset_filename("insertion_linkage_sensitivity.csv")),
          row.names = FALSE)
write.csv(
  insertion_branch_summary,
  file.path(OUTPUT_DIR_3,
            dataset_filename("insertion_branch_bootstrap_support.csv")),
  row.names = FALSE
)

linkage_plot <- linkage_sensitivity[
  linkage_sensitivity$linkage != "average", , drop = FALSE
]
figure <- ggplot(linkage_plot,
                 aes(linkage, jaccard_vs_average, fill = linkage)) +
  geom_boxplot(width = 0.30, outlier.shape = NA, alpha = 0.72) +
  geom_jitter(width = 0.10, alpha = 0.35, size = 1.4,
              shape = 21, color = "#333333") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_fill_manual(values = c(complete = "#0072B2", single = "#E69F00"),
                    guide = "none") +
  labs(title = "Insertion-branch sensitivity to linkage rule",
       subtitle = "Jaccard overlap of sibling-branch members relative to average linkage",
       x = "Alternative linkage", y = "Branch-member overlap") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 dataset_filename("insertion_linkage_sensitivity.pdf")),
       figure, width = 8.5, height = 6.5)

pdf(file.path(
  OUTPUT_DIR_3,
  dataset_filename("individual_insertion_dendrograms_with_branch_support.pdf")
), width = 16, height = 11)
par(mfrow = c(2, 2), mar = c(8, 3.5, 4.5, 1.5) + 0.1)
for (index in seq_len(nrow(betula_spectra))) {
  insertion_tree <- insertion_trees[[index]]
  query_label <- paste0(toupper(DATASET_DISPLAY_NAME), " ", peak_id[index])
  sample_peak_group <- paste(DATASET_DISPLAY_NAME, "peak")
  label_classes <- c(
    setNames(compound_classes, compound_levels),
    setNames(sample_peak_group, query_label)
  )
  insertion_colors <- c(class_colors,
                        setNames("#000000", sample_peak_group))
  colored_insertion_tree <- color_dendrogram_labels(
    as.dendrogram(insertion_tree), label_classes, insertion_colors
  )
  support <- insertion_branch_summary$branch_neighborhood_bootstrap_support[index]
  clade_support <-
    insertion_branch_summary$exact_sibling_clade_bootstrap_support[index]
  nearest_support <-
    insertion_branch_summary$nearest_compound_bootstrap_support[index]
  plot(
    colored_insertion_tree,
    main = paste0(
      peak_id[index], " | exact clade = ", round(100 * clade_support),
      "% | neighborhood = ", round(100 * support), "%\n",
      "nearest: ", insertion_branch_summary$nearest_compound[index],
      " (", round(100 * nearest_support), "%)"
    ),
    ylab = "Ensemble spectral distance", xlab = "", sub = "", cex = 0.58
  )
  mtext(
    paste0("Sibling branch: ",
           insertion_branch_summary$immediate_sibling_branch_members[index]),
    side = 1, line = 6.2, cex = 0.62, col = "#333333"
  )
}
dev.off()

ambiguous_subtype_data <- as.data.frame(table(
  results$ambiguity_subtype[results$practical_spectral_decision == "Ambiguous"]
), responseName = "n_peaks")
names(ambiguous_subtype_data)[1] <- "ambiguity_subtype"
ambiguous_subtype_data$ambiguity_subtype <- factor(
  ambiguous_subtype_data$ambiguity_subtype,
  levels = ambiguous_subtype_data$ambiguity_subtype[
    order(ambiguous_subtype_data$n_peaks)
  ]
)
figure <- ggplot(ambiguous_subtype_data,
                 aes(n_peaks, ambiguity_subtype)) +
  geom_col(fill = "#777777") +
  geom_text(aes(label = n_peaks), hjust = -0.25, size = 3.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Why compound-level results remain ambiguous",
       subtitle = "Each ambiguous peak is assigned its primary actionable failure mode",
       x = paste("Number of", DATASET_DISPLAY_NAME, "peaks"), y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "ambiguous_peak_failure_modes.pdf"),
       figure, width = 10, height = 6)

analogue_set_plot <- results[results$practical_spectral_decision !=
                               "Outside library coverage", , drop = FALSE]
figure <- ggplot(analogue_set_plot,
                 aes(analogue_set_size, analogue_set_chromophore_consistency,
                     color = practical_spectral_decision)) +
  geom_jitter(width = 0.08, height = 0.015, alpha = 0.75, size = 2.2) +
  scale_color_manual(values = c("Strong analogue support" = "#228833",
                                "Moderate analogue support" = "#4477AA",
                                "Ambiguous" = "#888888"),
                     name = "Spectral decision") +
  scale_x_continuous(breaks = seq(1, max(analogue_set_plot$analogue_set_size))) +
  scale_y_continuous(limits = c(0, 1.02)) +
  labs(title = "Accepted analogue-set size and chromophore consistency",
       subtitle = paste0("A multi-compound set can still provide a chromophore ",
                         "tendency when its members agree"),
       x = "Number of spectrally equivalent library analogues",
       y = "Fraction of analogue set in its modal chromophore family") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3,
                 "analogue_set_size_and_chromophore_consistency.pdf"),
       figure, width = 9, height = 7)

if (nrow(tie_wavelength_evidence)) {
  tie_counts <- aggregate(discriminating ~ peak_id, tie_wavelength_evidence, sum)
  # Include every peak for which a top-pair comparison was calculated. The
  # ordering still places the most discriminating examples first, but no
  # retention-time peak is silently dropped from this individual-spectrum view.
  selected_ties <- tie_counts$peak_id[
    order(tie_counts$discriminating, decreasing = TRUE)
  ]
  selected_evidence <- tie_wavelength_evidence[
    tie_wavelength_evidence$peak_id %in% selected_ties, , drop = FALSE
  ]
  selected_evidence$comparison <- paste0(
    selected_evidence$peak_id, ": ", selected_evidence$first_analogue,
    " vs ", selected_evidence$second_analogue
  )

  # Give every panel an interpretation, rather than asking the reader to infer
  # the tie from three overlapping curves. At wavelengths where the two library
  # IQRs do not overlap, count which candidate median is closer to the sample.
  # This is descriptive evidence only; it does not replace the full-distance
  # and conformal criteria used for the final classification.
  tie_panel_labels <- do.call(rbind, lapply(selected_ties, function(current_peak) {
    current <- selected_evidence[selected_evidence$peak_id == current_peak, ]
    separating <- current[current$discriminating, , drop = FALSE]
    first_name <- as.character(current$first_analogue[1])
    second_name <- as.character(current$second_analogue[1])
    first_fraction <- if (nrow(separating)) {
      mean(abs(separating$betula - separating$first_median) <=
             abs(separating$betula - separating$second_median))
    } else {
      0.5
    }
    data.frame(
      peak_id = current_peak,
      panel_label = paste0(
        current_peak, " | ", first_name, " vs ", second_name, "\n",
        nrow(separating), " separating wavelengths; sample closer to ",
        first_name, " at ", round(100 * first_fraction), "% and ",
        second_name, " at ", round(100 * (1 - first_fraction)), "%"
      ),
      stringsAsFactors = FALSE
    )
  }))
  selected_evidence$panel_label <- tie_panel_labels$panel_label[
    match(selected_evidence$peak_id, tie_panel_labels$peak_id)
  ]
  first_long <- data.frame(
    selected_evidence[, c("peak_id", "wavelength", "comparison", "panel_label")],
    role = "Leading library candidate", median = selected_evidence$first_median,
    q25 = selected_evidence$first_q25, q75 = selected_evidence$first_q75
  )
  second_long <- data.frame(
    selected_evidence[, c("peak_id", "wavelength", "comparison", "panel_label")],
    role = "Alternative library candidate", median = selected_evidence$second_median,
    q25 = selected_evidence$second_q25, q75 = selected_evidence$second_q75
  )
  analogue_long <- rbind(first_long, second_long)

  # Six comparisons per page keeps spectra and annotations readable. Gold
  # vertical bands mark wavelengths where the candidates' replicate IQRs do
  # not overlap; unshaded regions provide little pair-specific discrimination.
  ties_per_page <- 6L
  tie_pages <- split(selected_ties,
                     ceiling(seq_along(selected_ties) / ties_per_page))
  pdf(file.path(OUTPUT_DIR_3, "ambiguous_top_pair_wavelength_evidence.pdf"),
      width = 12, height = 10.5, onefile = TRUE)
  for (page_peaks in tie_pages) {
    page_evidence <- selected_evidence[
      selected_evidence$peak_id %in% page_peaks, , drop = FALSE
    ]
    page_analogues <- analogue_long[
      analogue_long$peak_id %in% page_peaks, , drop = FALSE
    ]
    page_evidence$panel_label <- factor(
      page_evidence$panel_label,
      levels = tie_panel_labels$panel_label[
        match(page_peaks, tie_panel_labels$peak_id)
      ]
    )
    page_analogues$panel_label <- factor(
      page_analogues$panel_label, levels = levels(page_evidence$panel_label)
    )
    separating_bands <- page_evidence[page_evidence$discriminating, ]

    figure <- ggplot() +
      geom_rect(data = separating_bands,
                aes(xmin = wavelength - 2.5, xmax = wavelength + 2.5,
                    ymin = -Inf, ymax = Inf),
                inherit.aes = FALSE, fill = "#E6B800", alpha = 0.14) +
      geom_ribbon(data = page_analogues,
                  aes(wavelength, ymin = q25, ymax = q75, fill = role),
                  alpha = 0.20, color = NA) +
      geom_line(data = page_analogues,
                aes(wavelength, median, color = role), linewidth = 0.75) +
      geom_line(data = page_evidence,
                aes(wavelength, betula, color = "Sample spectrum"),
                linewidth = 1.0) +
      facet_wrap(~ panel_label, ncol = 2) +
      scale_color_manual(
        values = c("Sample spectrum" = "#111111",
                   "Leading library candidate" = "#4477AA",
                   "Alternative library candidate" = "#CC6677"),
        breaks = c("Sample spectrum", "Leading library candidate",
                   "Alternative library candidate"), name = NULL
      ) +
      scale_fill_manual(
        values = c("Leading library candidate" = "#4477AA",
                   "Alternative library candidate" = "#CC6677"),
        guide = "none"
      ) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                         expand = expansion(mult = c(0, 0.02))) +
      coord_cartesian(ylim = c(0, 1.03)) +
      labs(
        title = "Why the leading compound assignments remain ambiguous",
        subtitle = paste0(
          DATASET_DISPLAY_NAME,
          " spectrum compared with its two strongest library candidates; ",
          "ribbons show replicate IQR"
        ),
        caption = paste0(
          "Gold bands: candidate IQRs do not overlap. Percentages summarize ",
          "which median is closer only within those bands; final decisions ",
          "also use full-spectrum distance and conformal support."
        ),
        x = "Wavelength (nm)", y = "Normalized absorbance"
      ) +
      theme_bw(base_size = 9.5) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        legend.position = "bottom",
        legend.box = "horizontal",
        strip.background = element_rect(fill = "grey94", color = "grey75"),
        strip.text = element_text(face = "bold", size = 7.5, lineheight = 1.05),
        plot.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, color = "grey30", size = 8),
        plot.margin = margin(8, 10, 8, 8)
      )
    print(figure)
  }
  dev.off()
}

heatmap_long <- do.call(rbind, lapply(seq_len(nrow(ensemble_probability)), function(i) {
  data.frame(peak_id = peak_id[i], compound = compound_levels,
             probability = ensemble_probability[i, ],
             class = compound_classes)
}))
peak_order <- peak_id[order(results$conformal_unsupported,
                            -results$ensemble_compound_probability)]
heatmap_long$peak_id <- factor(heatmap_long$peak_id, levels = rev(peak_order))
heatmap_long$compound <- factor(heatmap_long$compound,
                                levels = compound_levels[order(compound_classes)])
figure <- ggplot(heatmap_long, aes(compound, peak_id, fill = probability)) +
  geom_tile() +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      name = "Relative analogue support") +
  labs(title = paste0(DATASET_DISPLAY_NAME, "-to-library compound support"),
       subtitle = "Consensus across six spectral distance representations",
       x = "Library compound", y = paste(DATASET_DISPLAY_NAME, "peak")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_text(size = 5), plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, dataset_filename("library_support_heatmap.pdf")),
       figure, width = 13, height = 16)

conformal_long <- do.call(rbind, lapply(seq_len(nrow(class_p)), function(i) {
  data.frame(peak_id = peak_id[i], class = colnames(class_p), p_value = class_p[i, ])
}))
conformal_long$peak_id <- factor(conformal_long$peak_id, levels = rev(peak_order))
figure <- ggplot(conformal_long, aes(class, peak_id, fill = p_value)) +
  geom_tile() +
  geom_hline(yintercept = integer(0)) +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      limits = c(0, 1), name = "Conformal p-value") +
  labs(title = "Diagnostic predefined-class compatibility",
       subtitle = paste0("Shown for interpretation only; predefined classes ",
                         "are not used for the primary decision"),
       x = "Library class", y = paste(DATASET_DISPLAY_NAME, "peak")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 5), plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, dataset_filename("conformal_class_heatmap.pdf")),
       figure, width = 9, height = 16)

co_long <- as.data.frame(as.table(consensus$co_cluster))
names(co_long) <- c("row", "column", "support")
order_indices <- consensus$tree$order
co_long$row <- factor(peak_id[co_long$row], levels = rev(peak_id[order_indices]))
co_long$column <- factor(peak_id[co_long$column], levels = peak_id[order_indices])
figure <- ggplot(co_long, aes(column, row, fill = support)) +
  geom_raster() +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      limits = c(0, 1), name = "Co-clustering support") +
  labs(title = paste0(DATASET_DISPLAY_NAME, "-only consensus clustering"),
       subtitle = paste0(CONSENSUS_REPETITIONS,
                         " feature-resampling trees; k = ", consensus$k),
       x = paste(DATASET_DISPLAY_NAME, "peak"),
       y = paste(DATASET_DISPLAY_NAME, "peak")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, dataset_filename("consensus_matrix.pdf")),
       figure, width = 10, height = 9)

agreement_long <- do.call(rbind, lapply(seq_len(nrow(method_neighbor)), function(i) {
  data.frame(peak_id = peak_id[i], method = colnames(method_neighbor),
             neighbor = method_neighbor[i, ],
             class = method_class[i, ])
}))
agreement_long$peak_id <- factor(agreement_long$peak_id, levels = rev(peak_order))
figure <- ggplot(agreement_long, aes(method, peak_id, fill = class)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = class_colors, name = "Neighbor class") +
  labs(title = "Agreement among spectral distance methods",
       subtitle = "Each cell shows the class of the nearest compound",
       x = "Distance representation", y = paste(DATASET_DISPLAY_NAME, "peak")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 5), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "distance_method_agreement.pdf"),
       figure, width = 10, height = 16)

figure <- ggplot(results, aes(mixture_relative_improvement,
                              mixture_minor_fraction,
                              color = possible_coelution,
                              size = best_single_rmse)) +
  geom_vline(xintercept = mixture_improvement_threshold,
             color = "grey60", linetype = 2) +
  geom_hline(yintercept = 0.15, color = "grey60", linetype = 2) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(`FALSE` = "#4477AA", `TRUE` = "#EE6677"),
                     name = "Possible co-elution") +
  labs(title = "Single-compound versus two-compound reconstruction",
       subtitle = "Upper-right peaks gain materially from two library spectra",
       x = "Relative reduction in reconstruction error",
       y = "Minor component fraction", size = "Single-model RMSE") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "mixture_coelution_diagnostics.pdf"),
       figure, width = 9, height = 7)

figure <- ggplot(results, aes(smoothing_rmse, ensemble_compound_probability,
                              color = conformal_unsupported,
                              size = raw_dynamic_range)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(`FALSE` = "#4477AA", `TRUE` = "#EE6677"),
                     name = "Conformally unsupported") +
  labs(title = "Spectral quality versus library support",
       subtitle = "Distinguishes noisy spectra from structured unsupported peaks",
       x = "Smoothing residual (RMSE)",
       y = "Best ensemble compound probability",
       size = "Raw dynamic range") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "spectral_quality_vs_support.pdf"),
       figure, width = 9, height = 7)

block_long <- do.call(rbind, lapply(seq_len(nrow(block_support$class_support)),
                                   function(i) {
  data.frame(peak_id = peak_id[i], class = colnames(block_support$class_support),
             support = block_support$class_support[i, ])
}))
block_long$peak_id <- factor(block_long$peak_id, levels = rev(peak_order))
figure <- ggplot(block_long, aes(class, peak_id, fill = support)) +
  geom_tile() +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      limits = c(0, 1), name = "Block-bootstrap support") +
  labs(title = "Wavelength-region class support",
       subtitle = "Contiguous wavelength blocks are omitted across bootstrap runs",
       x = "Nearest-compound class", y = paste(DATASET_DISPLAY_NAME, "peak")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 5), plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "wavelength_block_class_support.pdf"),
       figure, width = 9, height = 16)

omission_modal <- do.call(rbind, lapply(split(
  omission_support_long,
  interaction(omission_support_long$peak_id,
              omission_support_long$omitted_wavelength_block, drop = TRUE)
), function(group) group[which.max(group$support), ]))
omission_modal$peak_id <- factor(omission_modal$peak_id,
                                 levels = rev(peak_order))
figure <- ggplot(omission_modal,
                 aes(omitted_wavelength_block, peak_id,
                     fill = class, alpha = support)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_manual(values = class_colors, name = "Modal class") +
  scale_alpha_continuous(range = c(0.25, 1), limits = c(0, 1),
                         name = "Conditional support") +
  labs(title = "Sensitivity to omitted wavelength regions",
       subtitle = "Color is the modal class when each region is withheld",
       x = "Omitted wavelength block", y = paste(DATASET_DISPLAY_NAME, "peak")) +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = 5), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "wavelength_omission_sensitivity.pdf"),
       figure, width = 10, height = 16)

figure <- ggplot(results,
                 aes(ensemble_min_distance,
                     library_domain_p_value,
                     color = conformal_unsupported,
                     size = method_compound_agreement)) +
  geom_point(alpha = 0.82) +
  scale_color_manual(values = c(`FALSE` = "#4477AA", `TRUE` = "#EE6677"),
                     name = "Conformally unsupported") +
  labs(title = "Library distance and conformal support",
       subtitle = "Calibrated by withholding complete library compounds",
       x = "Minimum ensemble distance to a library compound",
       y = "Spectral-library coverage p-value",
       size = "Analogue-method agreement") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, "library_distance_vs_conformal_support.pdf"),
       figure, width = 9, height = 7)

# Visualization-only pseudochromatogram suite. The class, chromophore, and
# class-first hybrid views use the same peaks and heights and differ only in
# reporting labels and colors. Retention-time metadata determines horizontal
# placement but is never supplied to a classifier or distance model.
# For the chromophore chromatogram, display deployed broad families rather than
# the more numerous fine groups. Fine assignments remain in the exported table.
fine_to_broad_code <- setNames(
  fine_chromophore_dictionary$broad_chromophore_code,
  fine_chromophore_dictionary$fine_chromophore_group
)
broad_family_by_code <- setNames(
  chromophore_dictionary$chromophore_family,
  chromophore_dictionary$chromophore_code
)
chromophore_display_label <- results$chromophore_only_label
fine_label_rows <- chromophore_display_label %in% names(fine_to_broad_code)
chromophore_display_label[fine_label_rows] <- unname(broad_family_by_code[
  fine_to_broad_code[chromophore_display_label[fine_label_rows]]
])

# In the hybrid view, classes retain their class colors. All chromophore-based
# interpretations are collapsed to broad families for a concise legend. The
# displayed label records whether the color represents a class or chromophore;
# this avoids making a chromophore-based interpretation look like a class call.
hybrid_display_label <- results$hybrid_class_then_chromophore_label
fallback_rows <- results$hybrid_label_type %in%
  c("Fine chromophore fallback", "Broad chromophore fallback")
hybrid_display_label[fallback_rows] <-
  results$analogue_set_modal_chromophore_family[fallback_rows]
bootstrap_chromophore_rows <-
  results$hybrid_label_type == "Bootstrap chromophore tendency"
hybrid_display_label[bootstrap_chromophore_rows] <-
  results$bootstrap_modal_chromophore_family[bootstrap_chromophore_rows]
hybrid_plot_label <- hybrid_display_label
hybrid_class_rows <- results$hybrid_label_type %in%
  c("Class from nearest analogue", "Class tendency")
hybrid_chromophore_rows <- results$hybrid_label_type %in% c(
  "Bootstrap chromophore tendency", "Fine chromophore fallback",
  "Broad chromophore fallback"
)
hybrid_plot_label[hybrid_class_rows] <- paste0(
  "Class: ", hybrid_display_label[hybrid_class_rows]
)
hybrid_plot_label[hybrid_chromophore_rows] <- paste0(
  "Chromophore: ", hybrid_display_label[hybrid_chromophore_rows]
)
pseudochromatogram <- data.frame(
  peak_id = results$peak_id,
  retention_time = results$retention_time_metadata_only,
  relative_height = 100 * results$maximum_raw_absorbance /
    max(results$maximum_raw_absorbance, na.rm = TRUE),
  spectral_decision = results$practical_spectral_decision,
  class_only_label = results$class_only_label,
  chromophore_only_label = results$chromophore_only_label,
  chromophore_display_label = chromophore_display_label,
  hybrid_label_type = results$hybrid_label_type,
  hybrid_class_then_chromophore_label =
    results$hybrid_class_then_chromophore_label,
  hybrid_display_label = hybrid_display_label,
  hybrid_plot_label = hybrid_plot_label,
  analogue_support = results$ensemble_compound_probability
)
pseudochromatogram <- pseudochromatogram[
  is.finite(pseudochromatogram$retention_time) &
    is.finite(pseudochromatogram$relative_height), , drop = FALSE
]
write.csv(pseudochromatogram,
          file.path(OUTPUT_DIR_3, dataset_filename("pseudochromatogram_data.csv")),
          row.names = FALSE)
peak_width <- 0.035
trace_time <- sort(unique(c(
  seq(min(pseudochromatogram$retention_time) - 4 * peak_width,
      max(pseudochromatogram$retention_time) + 4 * peak_width,
      length.out = 5000),
  pseudochromatogram$retention_time
)))
trace_intensity <- vapply(trace_time, function(current_time) {
  sum(pseudochromatogram$relative_height *
        exp(-0.5 * ((current_time - pseudochromatogram$retention_time) /
                      peak_width)^2))
}, numeric(1))
trace_scale <- 100 / max(trace_intensity)
trace_data <- data.frame(retention_time = trace_time,
                         intensity = trace_intensity * trace_scale)

make_pseudochromatogram <- function(display_group, colors, title, subtitle,
                                    legend_title, filename,
                                    mirror_filename,
                                    legend_rows = 3, height = 9) {
  peak_curves <- do.call(rbind, lapply(seq_len(nrow(pseudochromatogram)), function(i) {
    center <- pseudochromatogram$retention_time[i]
    curve_time <- seq(center - 4 * peak_width, center + 4 * peak_width,
                      length.out = 81)
    data.frame(
      peak_id = pseudochromatogram$peak_id[i], retention_time = curve_time,
      intensity = pseudochromatogram$relative_height[i] *
        exp(-0.5 * ((curve_time - center) / peak_width)^2) * trace_scale,
      display_group = display_group[i]
    )
  }))
  # Numerical invariant behind both figure designs: after the shared scaling,
  # the black envelope must remain the pointwise sum of all component curves.
  interpolated_components <- lapply(
    split(peak_curves, peak_curves$peak_id), function(component_curve) {
      stats::approx(
        component_curve$retention_time, component_curve$intensity,
        xout = trace_data$retention_time, yleft = 0, yright = 0
      )$y
    }
  )
  reconstructed_trace <- Reduce(`+`, interpolated_components)
  # The analytic Gaussian used to choose the initial scale has infinitesimal
  # tails beyond the displayed +/-4-width component support. Rebuild the black
  # curve from the exact displayed components, then rescale both together so
  # its maximum is exactly 100%.
  display_correction <- 100 / max(reconstructed_trace, na.rm = TRUE)
  peak_curves$intensity <- peak_curves$intensity * display_correction
  interpolated_components <- lapply(
    interpolated_components, `*`, display_correction
  )
  trace_data$intensity <- reconstructed_trace * display_correction
  trace_scale <- trace_scale * display_correction
  reconstructed_trace <- Reduce(`+`, interpolated_components)
  figure <- ggplot() +
    geom_area(data = peak_curves,
              aes(retention_time, intensity, group = peak_id,
                  fill = display_group),
              position = "identity", color = NA, alpha = 0.72) +
    geom_line(data = trace_data, aes(retention_time, intensity),
              color = "#222222", linewidth = 0.55) +
    geom_hline(yintercept = 0, color = "#222222", linewidth = 0.4) +
    scale_fill_manual(
      values = colors, name = legend_title,
      labels = function(labels) vapply(
        labels, function(label) paste(strwrap(label, width = 34),
                                      collapse = "\n"),
        character(1), USE.NAMES = FALSE
      )
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = title, subtitle = subtitle,
         x = "Retention time (display metadata only)",
         y = "Relative summed apex signal (%)") +
    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    theme_classic(base_size = 11) +
    theme(legend.position = "bottom", legend.box = "vertical",
          legend.text = element_text(size = 9),
          legend.key.height = grid::unit(0.45, "cm"),
          plot.title = element_text(face = "bold", margin = margin(b = 3)),
          plot.subtitle = element_text(margin = margin(b = 8)),
          plot.title.position = "plot",
          plot.margin = margin(12, 14, 12, 12))
  ggsave(file.path(OUTPUT_DIR_3, filename), figure,
         width = 14, height = height)

  # Companion mirror view. The upper trace is the exact pointwise sum of the
  # Gaussian pseudopeaks. The lower sticks show each peak's apex contribution
  # and interpretation on the identical numerical scale.
  stick_rows <- vapply(
    split(seq_len(nrow(peak_curves)), peak_curves$peak_id),
    function(rows) rows[which.max(peak_curves$intensity[rows])], integer(1)
  )
  sticks <- peak_curves[stick_rows, c(
    "peak_id", "retention_time", "intensity", "display_group"
  )]
  mirror_limit <- max(c(trace_data$intensity, sticks$intensity), na.rm = TRUE)
  mirror <- ggplot() +
    geom_area(data = trace_data, aes(retention_time, intensity),
              fill = "#D9D9D9", alpha = 0.42, color = NA) +
    geom_line(data = trace_data, aes(retention_time, intensity),
              color = "#000000", linewidth = 0.72) +
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
      limits = c(-1.05 * mirror_limit, 1.05 * mirror_limit),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = paste0(title, " - signal/composition mirror"),
      subtitle = paste(
        "Top: summed black pseudochromatogram; bottom: individual peak",
        "apex contributions colored by interpretation"
      ),
      caption = paste(
        "Both halves use the same scale. Lower sticks are discrete peak",
        "contributions; retention time is display metadata only."
      ),
      x = "Retention time (display metadata only)",
      y = "Relative summed apex signal (%)"
    ) +
    guides(color = guide_legend(nrow = legend_rows, byrow = TRUE)) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom", legend.box = "vertical",
      legend.text = element_text(size = 9),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 8),
      plot.title.position = "plot"
    )
  ggsave(file.path(OUTPUT_DIR_3, mirror_filename), mirror,
         width = 14, height = height)

  data.frame(
    figure_file = filename,
    mirror_file = mirror_filename,
    shared_scale_factor = trace_scale,
    black_trace_maximum = max(trace_data$intensity, na.rm = TRUE),
    largest_component_apex = max(sticks$intensity, na.rm = TRUE),
    maximum_reconstruction_error = max(
      abs(reconstructed_trace - trace_data$intensity), na.rm = TRUE
    ),
    black_sum_not_below_any_component = all(vapply(
      interpolated_components,
      function(values) all(values <= trace_data$intensity + 1e-8),
      logical(1)
    )),
    stringsAsFactors = FALSE
  )
}

neutral_colors <- c(Ambiguous = "#BDBDBD",
                    `Outside library coverage` = "#4D4D4D")
pseudochromatogram_scaling_audit <- list()
pseudochromatogram_scaling_audit[["class"]] <- make_pseudochromatogram(
  pseudochromatogram$class_only_label,
  c(class_colors, neutral_colors),
  paste(DATASET_DISPLAY_NAME, "class-annotated pseudochromatogram"),
  paste0("Class is inherited from supported analogue evidence; unresolved peaks ",
         "are gray and peaks outside library coverage are dark gray"),
  "Class interpretation", dataset_filename("class_pseudochromatogram.pdf"),
  dataset_filename("class_signal_composition_mirror.pdf"),
  legend_rows = 2, height = 8
)

fine_color_by_name <- setNames(
  fine_chromophore_colors,
  fine_chromophore_dictionary$fine_chromophore_group
)
broad_color_by_name <- setNames(
  chromophore_colors,
  chromophore_dictionary$chromophore_family
)
pseudochromatogram_scaling_audit[["chromophore"]] <- make_pseudochromatogram(
  pseudochromatogram$chromophore_display_label,
  c(broad_color_by_name, neutral_colors),
  paste(DATASET_DISPLAY_NAME, "chromophore-family pseudochromatogram"),
  paste0(length(setdiff(chromophore_codes, "DUV")),
         " deployed broad families provide concise colors; fine-group ",
         "assignments remain in the results table; ambiguous peaks are gray"),
  "Broad chromophore family", dataset_filename("chromophore_pseudochromatogram.pdf"),
  dataset_filename("chromophore_signal_composition_mirror.pdf"),
  legend_rows = 3, height = 9
)

hybrid_class_colors <- setNames(
  unname(class_colors), paste0("Class: ", names(class_colors))
)
hybrid_chromophore_colors <- setNames(
  unname(broad_color_by_name),
  paste0("Chromophore: ", names(broad_color_by_name))
)
hybrid_colors <- c(
  hybrid_class_colors, hybrid_chromophore_colors, neutral_colors
)
pseudochromatogram_scaling_audit[["hybrid"]] <- make_pseudochromatogram(
  pseudochromatogram$hybrid_plot_label,
  hybrid_colors,
  paste(DATASET_DISPLAY_NAME, "combined class/chromophore pseudochromatogram"),
  NULL,
  "Interpretation shown",
  dataset_filename("hybrid_class_chromophore_pseudochromatogram.pdf"),
  dataset_filename("hybrid_signal_composition_mirror.pdf"),
  legend_rows = 4, height = 10
)
write.csv(
  do.call(rbind, pseudochromatogram_scaling_audit),
  file.path(
    OUTPUT_DIR_3, dataset_filename("pseudochromatogram_scaling_audit.csv")
  ),
  row.names = FALSE
)

fine_profile <- do.call(rbind, lapply(sort(unique(fine_groups)), function(group) {
  members <- betula_spectra[fine_groups == group, , drop = FALSE]
  data.frame(wavelength = library_spec$wavelengths,
             absorbance = apply(members, 2, median),
             q10 = apply(members, 2, quantile, probs = 0.10),
             q25 = apply(members, 2, quantile, probs = 0.25),
             q75 = apply(members, 2, quantile, probs = 0.75),
             q90 = apply(members, 2, quantile, probs = 0.90),
             fine_cluster = paste("Fine cluster", group))
}))
write.csv(fine_profile,
          file.path(OUTPUT_DIR_3,
                    dataset_filename("consensus_representative_spectra_with_uncertainty.csv")),
          row.names = FALSE)
figure <- ggplot(fine_profile, aes(wavelength, absorbance)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), fill = "#4477AA",
              alpha = 0.12, color = NA) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#4477AA",
              alpha = 0.28, color = NA) +
  geom_line(color = "#4477AA", linewidth = 0.8) +
  facet_wrap(~ fine_cluster, ncol = 2, scales = "free_y") +
  labs(title = paste("Representative spectra of", DATASET_DISPLAY_NAME,
                     "consensus groups"),
       subtitle = paste0("Median line; dark ribbon = interquartile range; ",
                         "light ribbon = 10th-90th percentiles"),
       x = "Wavelength (nm)", y = "Normalized absorbance") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR_3, dataset_filename("consensus_cluster_profiles.pdf")),
       figure, width = 10, height = 3 * ceiling(fine_k / 2))

pdf(file.path(OUTPUT_DIR_3, dataset_filename("consensus_dendrogram.pdf")),
    width = 15, height = 8)
group_palette <- grDevices::hcl.colors(fine_k, "Dark 3")
label_groups <- setNames(as.character(fine_groups), peak_id)
group_colors <- setNames(group_palette, as.character(seq_len(fine_k)))
colored_tree <- color_dendrogram_labels(as.dendrogram(consensus$tree),
                                        label_groups, group_colors)
par(mar = c(10, 4, 4, 2) + 0.1)
plot(colored_tree, main = paste0(DATASET_DISPLAY_NAME,
                                 "-only consensus dendrogram (fine k = ",
                                 fine_k, ")"),
     ylab = "Consensus dissimilarity", cex = 0.45)
legend("topright", legend = paste("Fine cluster", seq_len(fine_k)),
       col = group_palette, lty = 1, lwd = 3, cex = 0.7, bty = "n", ncol = 2)
dev.off()

unsupported_indices <- which(results$conformal_unsupported)
if (length(unsupported_indices) >= 2) {
  unknown_tree <- hclust(dist(betula_wavelet[unsupported_indices, , drop = FALSE]),
                         method = "average")
  pdf(file.path(OUTPUT_DIR_3,
                dataset_filename("unsupported_family_dendrogram.pdf")),
      width = 10, height = 7)
  par(mar = c(8, 4, 4, 2) + 0.1)
  plot(unknown_tree, labels = peak_id[unsupported_indices],
       main = paste("Relationships among conformally unsupported",
                    DATASET_DISPLAY_NAME, "peaks"),
       xlab = "", sub = "", ylab = "Wavelet-feature distance", cex = 0.75)
  dev.off()
}

if (!is.null(ambiguous_consensus)) {
  family_levels <- sort(unique(ambiguous_consensus$groups))
  family_colors <- setNames(grDevices::hcl.colors(length(family_levels), "Dark 3"),
                            family_levels)
  ambiguous_profile_summary$spectral_family <- factor(
    ambiguous_profile_summary$spectral_family, levels = family_levels,
    labels = paste("Family", family_levels)
  )
  figure <- ggplot(ambiguous_profile_summary,
                   aes(wavelength, median, group = spectral_family)) +
    geom_ribbon(aes(ymin = q10, ymax = q90), fill = "#4477AA",
                alpha = 0.12, color = NA) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#4477AA",
                alpha = 0.28, color = NA) +
    geom_line(color = "#1F4E79", linewidth = 0.9) +
    facet_wrap(~ spectral_family, ncol = 2, scales = "free_y") +
    labs(title = paste("Representative spectra of ambiguous", DATASET_DISPLAY_NAME,
                       "families"),
         subtitle = paste0("Line = median; dark ribbon = interquartile range; ",
                           "light ribbon = 10th-90th percentiles"),
         x = "Wavelength (nm)", y = "Normalized absorbance") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTPUT_DIR_3,
                   "ambiguous_family_representative_spectra.pdf"),
         figure, width = 11,
         height = 3.1 * ceiling(length(family_levels) / 2))
  ambiguous_co_long <- as.data.frame(as.table(ambiguous_consensus$co_cluster))
  names(ambiguous_co_long) <- c("row", "column", "support")
  ambiguous_peak_ids <- peak_id[ambiguous_indices]
  ambiguous_order <- ambiguous_consensus$tree$order
  ambiguous_co_long$row <- factor(
    ambiguous_peak_ids[ambiguous_co_long$row],
    levels = rev(ambiguous_peak_ids[ambiguous_order])
  )
  ambiguous_co_long$column <- factor(
    ambiguous_peak_ids[ambiguous_co_long$column],
    levels = ambiguous_peak_ids[ambiguous_order]
  )
  figure <- ggplot(ambiguous_co_long, aes(column, row, fill = support)) +
    geom_raster() +
    scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                        limits = c(0, 1), name = "Co-clustering support") +
    labs(title = paste("Consensus structure among ambiguous", DATASET_DISPLAY_NAME,
                       "peaks"),
         subtitle = paste0(CONSENSUS_REPETITIONS,
                           " feature-resampling trees; k = ",
                           ambiguous_consensus$k),
         x = "Ambiguous peak", y = "Ambiguous peak") +
    theme_minimal(base_size = 8) +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTPUT_DIR_3, "ambiguous_family_consensus_matrix.pdf"),
         figure, width = 9, height = 8)

  pdf(file.path(OUTPUT_DIR_3, "ambiguous_family_dendrogram.pdf"),
      width = 13, height = 7)
  ambiguous_label_groups <- setNames(
    as.character(ambiguous_consensus$groups), ambiguous_peak_ids
  )
  ambiguous_tree <- color_dendrogram_labels(
    as.dendrogram(ambiguous_consensus$tree), ambiguous_label_groups,
    family_colors
  )
  par(mar = c(9, 4, 4, 2) + 0.1)
  plot(ambiguous_tree,
       main = paste0("Ambiguous-peak spectral families (k = ",
                     ambiguous_consensus$k, ")"),
       ylab = "Consensus dissimilarity", cex = 0.65)
  legend("topright", legend = paste("Family", family_levels),
         col = family_colors, lty = 1, lwd = 3, bty = "n", ncol = 2)
  dev.off()
}

# Complete spectral-overlay atlas. Earlier exploratory versions selected nine
# examples; the reusable workflow now plots every retention-time peak exactly
# once so no individual spectrum is omitted by a subjective example filter.
overlay_indices <- seq_len(nrow(betula_spectra))
pdf(file.path(OUTPUT_DIR_3, "representative_spectral_overlays.pdf"),
    width = 12, height = 10)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1) + 0.1)
for (index in overlay_indices) {
  nearest_order <- order(ensemble_result$ensemble[index, ])[1:3]
  analogue_colors <- c("#4477AA", "#EE6677", "#228833")
  plot(
    library_spec$wavelengths, betula_spectra[index, ],
    type = "l", lty = 1, lwd = 2.5, col = "black",
    xlab = "Wavelength (nm)", ylab = "Normalized absorbance",
    main = paste0(peak_id[index], " | ",
                  results$practical_spectral_decision[index]),
    ylim = range(c(0, betula_spectra[index, ],
                                          shape_representatives[nearest_order, ]))
  )
  for (analogue_index in seq_along(nearest_order)) {
    compound_name <- rownames(shape_representatives)[nearest_order[analogue_index]]
    replicate_spectra <- library_spectra[compounds == compound_name, , drop = FALSE]
    lower <- apply(replicate_spectra, 2, quantile, probs = 0.25)
    upper <- apply(replicate_spectra, 2, quantile, probs = 0.75)
    polygon(c(library_spec$wavelengths, rev(library_spec$wavelengths)),
            c(lower, rev(upper)), border = NA,
            col = grDevices::adjustcolor(analogue_colors[analogue_index],
                                         alpha.f = 0.18))
    lines(library_spec$wavelengths,
          shape_representatives[nearest_order[analogue_index], ],
          col = analogue_colors[analogue_index], lwd = 1.2)
  }
  lines(library_spec$wavelengths, betula_spectra[index, ],
        col = "black", lwd = 2.5)
  legend("topright",
         legend = c(DATASET_DISPLAY_NAME,
                    rownames(shape_representatives)[nearest_order]),
         col = c("black", analogue_colors),
         lty = 1, lwd = c(2.5, 1.2, 1.2, 1.2), cex = 0.55, bty = "n")
}
dev.off()

cat("Comprehensive analysis complete.\n")
cat("Consensus clusters:", consensus$k, "\n")
cat("Outside library coverage", DATASET_NAME, "peaks:",
    sum(results$conformal_unsupported), "\n")
cat("Possible co-elution peaks:", sum(results$possible_coelution), "\n")
cat("Practical spectral decisions:\n")
print(table(results$practical_spectral_decision))
cat("Results written to:", normalizePath(OUTPUT_DIR_3), "\n")
