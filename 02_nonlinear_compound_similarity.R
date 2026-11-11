## Nonlinear, compound-first spectral matching
##
## The spectral model is learned without chemical-class labels. Library classes
## are attached only after a Betula spectrum has been compared with individual
## compound representatives.
##
## Role in the project
##   This file supplies reusable preprocessing, wavelet, representative,
##   bootstrap-insertion, and plotting functions. The current comprehensive
##   raw-DAD annotation engine sources these functions. It is a helper module,
##   not a production entry point; use RUN_AGILENT_DAD.R instead.
##
## Interpretation
##   Similarity identifies a closest library analogue; it does not establish
##   chemical identity. Bootstrap neighbor support and outgroup frequency show
##   whether that relative match is stable.
##
## Data contract
##   The active engines supply validated library and sample data. When sourced
##   by stages 03/04, the reusable functions in this module are used. The
##   production entry point is RUN_AGILENT_DAD.R.

library(ggplot2)

source("00_library_inputs.R")

BETULA_FILE <- file.path("data", "Betula_MaxAbsorbances.csv")
OUTPUT_DIR <- "nonlinear_similarity_results"
WAVELENGTH_MIN <- 190
WAVELENGTH_MAX <- 500
SMOOTHING_WINDOW <- 7
SMOOTHING_DEGREE <- 3


# ------------------------- INPUT AND PREPROCESSING --------------------------
# Shared input API. Spectra are rows and wavelengths are columns. Normalization
# and smoothing emphasize wavelength shape rather than concentration-dependent
# raw absorbance; neither operation uses sample class labels.

# Preserve literal wavelength headers (for example, `190nm`) and metadata names.
read_spectral_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# Standardize display labels in the retained standalone workflow. These edits
# affect legend text only; all spectra and analytical calculations are unchanged.
format_display_compound <- function(x) {
  vapply(x, function(value) {
    value <- trimws(value)
    if (tolower(value) == "p-cresol") return("p-Cresol")
    first_letter <- regexpr("[[:alpha:]]", value)[1]
    if (first_letter > 0) substr(value, first_letter, first_letter) <-
      toupper(substr(value, first_letter, first_letter))
    value
  }, character(1), USE.NAMES = FALSE)
}
format_display_class <- function(x) {
  ifelse(toupper(trimws(x)) == "NFP", "NFP",
         tools::toTitleCase(tolower(trimws(x))))
}


# Return integer indices and numeric wavelengths within the requested range;
# non-wavelength metadata columns are ignored.
spectral_columns <- function(data, minimum, maximum) {
  wavelengths <- suppressWarnings(as.numeric(sub("nm$", "", names(data))))
  keep <- !is.na(wavelengths) & wavelengths >= minimum & wavelengths <= maximum
  if (!any(keep)) stop("No wavelength columns were found.")
  list(columns = which(keep), wavelengths = wavelengths[keep])
}


# Derive centered Savitzky-Golay coefficients without an extra package.
sgolay_coefficients <- function(window, degree) {
  if (window %% 2 != 1 || degree >= window) {
    stop("The smoothing window must be odd and larger than the polynomial degree.")
  }
  offsets <- seq(-(window %/% 2), window %/% 2)
  design <- outer(offsets, 0:degree, `^`)
  as.numeric(solve(crossprod(design), t(design))[1, ])
}


# Smooth one spectrum and retain the original normalized edge values where a
# centered convolution cannot be evaluated.
smooth_spectrum <- function(x, window, degree) {
  coefficients <- sgolay_coefficients(window, degree)
  smoothed <- as.numeric(stats::filter(x, coefficients, sides = 2))
  edge <- is.na(smoothed)
  smoothed[edge] <- x[edge]
  smoothed
}


# Validate, row-range-normalize to [0, 1], and smooth a spectral matrix. Missing,
# non-finite, or nearly constant spectra stop the analysis explicitly. The
# matrix method is also used after a denser or offset sample grid has been
# interpolated onto the fixed library wavelength coordinates.
preprocess_spectral_matrix <- function(spectra, window, degree) {
  spectra <- as.matrix(spectra)
  storage.mode(spectra) <- "double"
  if (any(!is.finite(spectra))) stop("Spectra contain missing/non-finite values.")

  normalized <- t(apply(spectra, 1, function(x) {
    span <- max(x) - min(x)
    if (span <= sqrt(.Machine$double.eps)) stop("A spectrum is nearly constant.")
    (x - min(x)) / span
  }))
  t(apply(normalized, 1, smooth_spectrum, window = window, degree = degree))
}

preprocess_spectra <- function(data, columns, window, degree) {
  preprocess_spectral_matrix(
    data[, columns, drop = FALSE], window = window, degree = degree
  )
}


# ------------------------- NONLINEAR FEATURE SPACE --------------------------
# Haar features retain broad-band energy and localized changes without assuming
# that spectra form linearly separable groups. They supplement direct curve
# distances rather than replacing wavelength-resolved spectra.

# Pad one spectrum to a power of two and summarize detail coefficients at every
# scale by energy, mean absolute magnitude, and standard deviation.
haar_features_one <- function(spectrum) {
  target_length <- 2^ceiling(log2(length(spectrum)))
  approximation <- c(spectrum, rep(tail(spectrum, 1), target_length - length(spectrum)))
  features <- numeric()
  level <- 1L

  while (length(approximation) >= 2L) {
    odd <- approximation[seq(1, length(approximation), by = 2)]
    even <- approximation[seq(2, length(approximation), by = 2)]
    detail <- (odd - even) / sqrt(2)
    approximation <- (odd + even) / sqrt(2)
    features <- c(
      features,
      setNames(sum(detail^2), paste0("detail_energy_L", level)),
      setNames(mean(abs(detail)), paste0("detail_mean_abs_L", level)),
      setNames(if (length(detail) > 1) stats::sd(detail) else 0,
               paste0("detail_sd_L", level))
    )
    level <- level + 1L
  }
  c(features, approximation_energy = sum(approximation^2))
}


# Apply the Haar summary row-wise and retain stable feature names for scaling.
extract_haar_features <- function(spectra) {
  features <- t(apply(spectra, 1, haar_features_one))
  colnames(features) <- names(haar_features_one(spectra[1, ]))
  features
}


# ---------------------- REPRESENTATIVES AND MATCHING ------------------------
# Replicates are collapsed to robust compound medians. Class labels are attached
# only after spectra have been compared with individual compounds.

# Fit median/MAD scaling on library features only; invariant features are dropped.
fit_robust_scaler <- function(features) {
  centers <- apply(features, 2, median)
  scales <- apply(features, 2, mad, constant = 1.4826)
  keep <- is.finite(scales) & scales > sqrt(.Machine$double.eps)
  list(center = centers[keep], scale = scales[keep], keep = keep)
}


# Apply a frozen library scaler without learning from the sample dataset.
apply_robust_scaler <- function(features, scaler) {
  scaled <- sweep(features[, scaler$keep, drop = FALSE], 2, scaler$center, "-")
  sweep(scaled, 2, scaler$scale, "/")
}


# Return one pointwise median feature vector per compound.
make_compound_representatives <- function(features, compounds) {
  compound_names <- sort(unique(compounds))
  representatives <- t(vapply(compound_names, function(compound) {
    apply(features[compounds == compound, , drop = FALSE], 2, median)
  }, numeric(ncol(features))))
  rownames(representatives) <- compound_names
  representatives
}


# Compute all query-to-reference squared Euclidean distances efficiently.
squared_distances <- function(x, representatives) {
  x2 <- rowSums(x^2)
  r2 <- rowSums(representatives^2)
  distances <- outer(x2, r2, "+") - 2 * tcrossprod(x, representatives)
  pmax(distances, 0)
}


# Use median library nearest-neighbor spacing as the Gaussian-kernel bandwidth.
choose_kernel_bandwidth <- function(representatives) {
  distances <- sqrt(squared_distances(representatives, representatives))
  diag(distances) <- Inf
  bandwidth <- median(apply(distances, 1, min))
  if (!is.finite(bandwidth) || bandwidth <= 0) stop("Could not estimate bandwidth.")
  bandwidth
}


# Convert distances to relative support over the supplied standards. Rows sum
# to one, but values are not probabilities of molecular identity.
compound_probabilities <- function(features, representatives, bandwidth) {
  distances_squared <- squared_distances(features, representatives)
  log_weights <- -distances_squared / (2 * bandwidth^2)
  log_weights <- log_weights - apply(log_weights, 1, max)
  weights <- exp(log_weights)
  probabilities <- weights / rowSums(weights)
  colnames(probabilities) <- rownames(representatives)
  probabilities
}


# Aggregate compound support after matching. Averaging prevents classes with
# more standards from receiving support solely because they are larger.
aggregate_class_probabilities <- function(compound_probabilities, compound_classes) {
  class_levels <- sort(unique(compound_classes))
  output <- sapply(class_levels, function(class_name) {
    # Average rather than sum so a class with more library compounds does not
    # receive greater support solely because it has more spectral anchors.
    rowMeans(compound_probabilities[, compound_classes == class_name, drop = FALSE])
  })
  if (is.null(dim(output))) output <- matrix(output, ncol = 1)
  colnames(output) <- class_levels
  output / rowSums(output)
}


# Withhold every replicate of one compound at a time. This tests transfer to an
# unseen compound without leaving sibling replicates in the reference set.
grouped_leave_one_compound_out <- function(features, compounds, classes) {
  results <- list()
  support_distances <- numeric()
  result_row <- 1L

  for (held_compound in sort(unique(compounds))) {
    training <- compounds != held_compound
    representatives <- make_compound_representatives(features[training, , drop = FALSE],
                                                     compounds[training])
    bandwidth <- choose_kernel_bandwidth(representatives)
    held_rows <- which(compounds == held_compound)
    probabilities <- compound_probabilities(features[held_rows, , drop = FALSE],
                                            representatives, bandwidth)
    representative_classes <- classes[match(rownames(representatives), compounds)]
    class_probabilities <- aggregate_class_probabilities(probabilities,
                                                         representative_classes)
    predicted <- colnames(class_probabilities)[max.col(class_probabilities)]
    distances <- sqrt(squared_distances(features[held_rows, , drop = FALSE],
                                       representatives))
    support_distances <- c(support_distances, apply(distances, 1, min))

    for (i in seq_along(held_rows)) {
      observed <- classes[held_rows[i]]
      ordered <- sort(class_probabilities[i, ], decreasing = TRUE)
      results[[result_row]] <- data.frame(
        library_row = held_rows[i],
        held_out_compound = held_compound,
        observed_class = observed,
        predicted_class = predicted[i],
        predicted_probability = ordered[1],
        observed_class_probability = class_probabilities[i, observed],
        probability_margin = ordered[1] - ordered[2],
        nearest_training_distance = min(distances[i, ]),
        correct = predicted[i] == observed
      )
      result_row <- result_row + 1L
    }
  }

  list(results = do.call(rbind, results), support_distances = support_distances)
}


# ---------------------- DENDROGRAM INSERTION SUPPORT ------------------------
# Each query is inserted into a frozen library tree. Bootstrap runs resample
# compound replicates and wavelet features to measure how often its nearest
# compound, nearest class, or outgroup status remains unchanged.

# Insert one query independently and report its nearest standard, distance
# relative to typical library spacing, first merge height, and sibling size.
insertion_metrics_one <- function(query, representatives) {
  reference_distances <- sqrt(squared_distances(representatives, representatives))
  diag(reference_distances) <- Inf
  reference_scale <- median(apply(reference_distances, 1, min))
  query_distances <- sqrt(squared_distances(matrix(query, nrow = 1),
                                            representatives))[1, ]
  nearest_index <- which.min(query_distances)

  combined <- rbind(representatives, query)
  tree <- hclust(as.dist(sqrt(squared_distances(combined, combined))),
                 method = "average")
  query_leaf <- -nrow(combined)
  merge_row <- which(tree$merge[, 1] == query_leaf |
                     tree$merge[, 2] == query_leaf)[1]
  partner <- if (tree$merge[merge_row, 1] == query_leaf) {
    tree$merge[merge_row, 2]
  } else {
    tree$merge[merge_row, 1]
  }
  cluster_sizes <- integer(nrow(tree$merge))
  for (i in seq_len(nrow(tree$merge))) {
    cluster_sizes[i] <- sum(vapply(tree$merge[i, ], function(code) {
      if (code < 0) 1L else cluster_sizes[code]
    }, integer(1)))
  }
  partner_size <- if (partner < 0) 1L else cluster_sizes[partner]

  data.frame(
    nearest_compound = rownames(representatives)[nearest_index],
    nearest_distance = query_distances[nearest_index],
    nearest_distance_ratio = query_distances[nearest_index] / reference_scale,
    insertion_height = tree$height[merge_row],
    insertion_height_ratio = tree$height[merge_row] / reference_scale,
    joined_cluster_size = partner_size
  )
}


# Row-wise wrapper returning one insertion summary for every query.
insertion_metrics <- function(features, representatives) {
  do.call(rbind, lapply(seq_len(nrow(features)), function(i) {
    insertion_metrics_one(features[i, ], representatives)
  }))
}


# Estimate expected insertion behavior for a known-but-withheld compound by
# removing all of its replicates before inserting them into the reduced library.
grouped_known_insertion_metrics <- function(features, compounds) {
  output <- list()
  output_row <- 1L
  for (held_compound in sort(unique(compounds))) {
    training <- compounds != held_compound
    fold_representatives <- make_compound_representatives(
      features[training, , drop = FALSE], compounds[training]
    )
    held_rows <- which(compounds == held_compound)
    metrics <- insertion_metrics(features[held_rows, , drop = FALSE],
                                 fold_representatives)
    metrics$library_row <- held_rows
    metrics$held_out_compound <- held_compound
    output[[output_row]] <- metrics
    output_row <- output_row + 1L
  }
  do.call(rbind, output)
}


# Repeatedly resample feature columns and select one replicate per compound.
# Frequencies summarize exact-neighbor stability, class-level stability,
# within-class switching, singleton attachment, and outgroup behavior.
bootstrap_neighbor_support <- function(query_features, library_features,
                                       compounds, classes,
                                       full_representatives,
                                       outgroup_threshold,
                                       repetitions = 250,
                                       fraction = 0.8,
                                       seed = 20260806) {
  set.seed(seed)
  compound_levels <- rownames(full_representatives)
  compound_classes <- classes[match(compound_levels, compounds)]
  class_levels <- sort(unique(compound_classes))
  full_nearest <- insertion_metrics(query_features,
                                    full_representatives)$nearest_compound

  compound_counts <- matrix(
    0L, nrow(query_features), length(compound_levels),
    dimnames = list(NULL, compound_levels)
  )
  class_counts <- matrix(
    0L, nrow(query_features), length(class_levels),
    dimnames = list(NULL, class_levels)
  )
  outgroup_counts <- integer(nrow(query_features))
  singleton_counts <- integer(nrow(query_features))

  for (iteration in seq_len(repetitions)) {
    selected_features <- sort(sample(
      seq_len(ncol(query_features)),
      ceiling(fraction * ncol(query_features))
    ))

    # Choosing one replicate per compound propagates injection/technical
    # variability into the neighbor support rather than treating the median
    # representative as perfectly known.
    selected_rows <- vapply(compound_levels, function(compound) {
      sample(which(compounds == compound), 1)
    }, integer(1))
    bootstrap_representatives <- library_features[
      selected_rows, selected_features, drop = FALSE
    ]
    rownames(bootstrap_representatives) <- compound_levels

    for (i in seq_len(nrow(query_features))) {
      metrics <- insertion_metrics_one(
        query_features[i, selected_features], bootstrap_representatives
      )
      nearest_compound <- metrics$nearest_compound
      nearest_class <- compound_classes[match(nearest_compound, compound_levels)]
      compound_counts[i, nearest_compound] <-
        compound_counts[i, nearest_compound] + 1L
      class_counts[i, nearest_class] <- class_counts[i, nearest_class] + 1L
      outgroup_counts[i] <- outgroup_counts[i] +
        as.integer(metrics$insertion_height_ratio > outgroup_threshold)
      singleton_counts[i] <- singleton_counts[i] +
        as.integer(metrics$joined_cluster_size == 1L)
    }
  }

  compound_support <- compound_counts / repetitions
  class_support <- class_counts / repetitions
  modal_compound_index <- max.col(compound_support, ties.method = "first")
  modal_class_index <- max.col(class_support, ties.method = "first")
  full_nearest_class <- compound_classes[match(full_nearest, compound_levels)]

  entropy <- function(probabilities) {
    probabilities <- probabilities[probabilities > 0]
    -sum(probabilities * log(probabilities))
  }

  summary <- data.frame(
    bootstrap_full_nearest_compound_support = compound_support[
      cbind(seq_len(nrow(query_features)), match(full_nearest, compound_levels))
    ],
    bootstrap_modal_compound = compound_levels[modal_compound_index],
    bootstrap_modal_compound_support = compound_support[
      cbind(seq_len(nrow(query_features)), modal_compound_index)
    ],
    bootstrap_full_nearest_class_support = class_support[
      cbind(seq_len(nrow(query_features)), match(full_nearest_class, class_levels))
    ],
    bootstrap_modal_class = class_levels[modal_class_index],
    bootstrap_modal_class_support = class_support[
      cbind(seq_len(nrow(query_features)), modal_class_index)
    ],
    bootstrap_within_class_switching = pmax(
      0,
      class_support[cbind(seq_len(nrow(query_features)), modal_class_index)] -
        compound_support[cbind(seq_len(nrow(query_features)), modal_compound_index)]
    ),
    bootstrap_compound_entropy = apply(compound_support, 1, entropy),
    bootstrap_class_entropy = apply(class_support, 1, entropy),
    bootstrap_outgroup_frequency = outgroup_counts / repetitions,
    bootstrap_single_compound_attachment_frequency = singleton_counts / repetitions
  )

  list(summary = summary,
       compound_support = compound_support,
       class_support = class_support)
}


# --------------------------- FIGURE UTILITIES -------------------------------
# Dendrogram label coloring is kept separate from analysis so presentation
# choices cannot affect distances or cluster assignments.

# Recursively color leaf text from a named label-to-group mapping. Internal
# branches remain neutral because annotation color is not branch support.
color_dendrogram_labels <- function(dendrogram, label_classes, class_colors) {
  dendrapply(dendrogram, function(node) {
    if (is.leaf(node)) {
      label <- attr(node, "label")
      node_parameters <- attr(node, "nodePar")
      node_parameters$lab.col <- unname(class_colors[label_classes[label]])
      node_parameters$lab.cex <- 0.8
      node_parameters$pch <- NA
      attr(node, "nodePar") <- node_parameters
    }
    node
  })
}


# Standalone legacy comparison workflow. Sourcing this file loads only the
# reusable functions above; the block runs only when this file itself is called
# with Rscript. The comprehensive pipeline is the recommended entry point.
run_nonlinear_analysis <- function() {
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------- FIT USING LIBRARY ONLY ----------------------------
library_bundle <- load_library_inputs()
library_data <- library_bundle$spectra
library_data$Compound <- format_display_compound(library_data$Compound)
library_data$Short_Class <- format_display_class(library_data$Short_Class)
library_spec <- spectral_columns(library_data, WAVELENGTH_MIN, WAVELENGTH_MAX)
library_spectra <- preprocess_spectra(library_data, library_spec$columns,
                                     SMOOTHING_WINDOW, SMOOTHING_DEGREE)
library_features_raw <- extract_haar_features(library_spectra)
scaler <- fit_robust_scaler(library_features_raw)
library_features <- apply_robust_scaler(library_features_raw, scaler)

library_compounds <- library_data$Compound
library_classes <- library_data$Short_Class
representatives <- make_compound_representatives(library_features, library_compounds)
representative_classes <- library_classes[match(rownames(representatives), library_compounds)]
bandwidth <- choose_kernel_bandwidth(representatives)
validation <- grouped_leave_one_compound_out(
  library_features, library_compounds, library_classes
)

known_insertion <- grouped_known_insertion_metrics(library_features,
                                                   library_compounds)
outgroup_threshold <- unname(quantile(
  known_insertion$insertion_height_ratio, 0.95
))

# The spectral model is frozen before Betula is opened.
# -------------------------- APPLY TO BETULA --------------------------------
betula_data <- read_spectral_csv(BETULA_FILE)
betula_native_spec <- spectral_columns(betula_data, -Inf, Inf)
betula_grid_alignment <- harmonize_wavelength_matrix(
  betula_data[, betula_native_spec$columns, drop = FALSE],
  betula_native_spec$wavelengths,
  library_spec$wavelengths,
  input_id = basename(BETULA_FILE)
)
betula_spec <- list(wavelengths = library_spec$wavelengths)
betula_spectra <- preprocess_spectral_matrix(
  betula_grid_alignment$values, SMOOTHING_WINDOW, SMOOTHING_DEGREE
)
betula_features <- apply_robust_scaler(extract_haar_features(betula_spectra), scaler)

compound_probs <- compound_probabilities(betula_features, representatives, bandwidth)
conditional_class_probs <- aggregate_class_probabilities(compound_probs,
                                                         representative_classes)
betula_insertion <- insertion_metrics(betula_features, representatives)
bootstrap_support <- bootstrap_neighbor_support(
  betula_features,
  library_features,
  library_compounds,
  library_classes,
  representatives,
  outgroup_threshold,
  repetitions = 250
)
betula_stability <- bootstrap_support$summary
nearest_distances <- betula_insertion$nearest_distance
unknown_prob <- betula_stability$bootstrap_outgroup_frequency
known_prob <- 1 - unknown_prob
open_class_probs <- conditional_class_probs * known_prob
open_probs <- cbind(open_class_probs, Unknown = unknown_prob)

best_compound <- max.col(compound_probs, ties.method = "first")
best_known_class <- max.col(conditional_class_probs, ties.method = "first")
best_joint_outcome <- max.col(open_probs, ties.method = "first")
inherited_class <- representative_classes[best_compound]
compound_sorted <- t(apply(compound_probs, 1, sort, decreasing = TRUE))
class_sorted <- t(apply(open_probs, 1, sort, decreasing = TRUE))
hierarchical_label <- ifelse(
  known_prob >= unknown_prob,
  colnames(conditional_class_probs)[best_known_class],
  "Unknown"
)
hierarchical_joint_probability <- ifelse(
  hierarchical_label == "Unknown",
  unknown_prob,
  open_class_probs[cbind(seq_len(nrow(open_class_probs)), best_known_class)]
)

peak_id <- betula_data[[1]]
assignments <- data.frame(
  peak_id = peak_id,
  retention_time = suppressWarnings(as.numeric(sub("^X", "", peak_id))),
  provisional_compound = colnames(compound_probs)[best_compound],
  compound_inherited_class = inherited_class,
  compound_probability = compound_sorted[, 1],
  second_compound_probability = compound_sorted[, 2],
  compound_probability_margin = compound_sorted[, 1] - compound_sorted[, 2],
  mixture_most_likely_class = colnames(conditional_class_probs)[best_known_class],
  conditional_class_probability = conditional_class_probs[
    cbind(seq_len(nrow(conditional_class_probs)), best_known_class)
  ],
  inherited_class_probability = conditional_class_probs[
    cbind(seq_len(nrow(conditional_class_probs)),
          match(inherited_class, colnames(conditional_class_probs)))
  ],
  library_support_probability = known_prob,
  unknown_probability = unknown_prob,
  provisional_open_set_label = hierarchical_label,
  provisional_open_set_probability = hierarchical_joint_probability,
  library_vs_unknown_decision_probability = pmax(known_prob, unknown_prob),
  seven_way_joint_argmax_label = colnames(open_probs)[best_joint_outcome],
  seven_way_joint_argmax_probability = class_sorted[, 1],
  second_open_set_probability = class_sorted[, 2],
  open_set_probability_margin = class_sorted[, 1] - class_sorted[, 2],
  nearest_library_distance = nearest_distances,
  nearest_distance_ratio = betula_insertion$nearest_distance_ratio,
  dendrogram_insertion_height = betula_insertion$insertion_height,
  insertion_height_ratio = betula_insertion$insertion_height_ratio,
  joined_cluster_size = betula_insertion$joined_cluster_size,
  betula_stability,
  check.names = FALSE
)

write.csv(assignments, file.path(OUTPUT_DIR, "betula_provisional_labels.csv"),
          row.names = FALSE)
write.csv(cbind(betula_data, assignments[, -c(1:2), drop = FALSE]),
          file.path(OUTPUT_DIR, "betula_spectra_with_provisional_labels.csv"),
          row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, compound_probs, check.names = FALSE),
          file.path(OUTPUT_DIR, "betula_compound_probabilities.csv"), row.names = FALSE)
write.csv(data.frame(peak_id = peak_id, open_probs, check.names = FALSE),
          file.path(OUTPUT_DIR, "betula_class_and_unknown_probabilities.csv"),
          row.names = FALSE)
write.csv(validation$results,
          file.path(OUTPUT_DIR, "grouped_leave_one_compound_out_validation.csv"),
          row.names = FALSE)
write.csv(data.frame(
  peak_id = peak_id,
  bootstrap_support$compound_support,
  check.names = FALSE
), file.path(OUTPUT_DIR, "betula_bootstrap_compound_neighbor_support.csv"),
row.names = FALSE)
write.csv(data.frame(
  peak_id = peak_id,
  bootstrap_support$class_support,
  check.names = FALSE
), file.path(OUTPUT_DIR, "betula_bootstrap_class_neighbor_support.csv"),
row.names = FALSE)
write.csv(known_insertion,
          file.path(OUTPUT_DIR, "held_out_compound_insertion_reference.csv"),
          row.names = FALSE)
write.csv(data.frame(
  metric = c("grouped_class_accuracy", "outgroup_height_ratio_threshold"),
  value = c(mean(validation$results$correct), outgroup_threshold)
), file.path(OUTPUT_DIR, "model_diagnostics.csv"), row.names = FALSE)
write.csv(data.frame(
  compound = rownames(representatives),
  class = representative_classes,
  representatives,
  check.names = FALSE
), file.path(OUTPUT_DIR, "compound_wavelet_representatives.csv"), row.names = FALSE)
write.csv(data.frame(
  bandwidth = bandwidth,
  bootstrap_repetitions = 250,
  bootstrap_feature_fraction = 0.8,
  outgroup_height_ratio_threshold = outgroup_threshold,
  wavelength_min = WAVELENGTH_MIN,
  wavelength_max = WAVELENGTH_MAX
), file.path(OUTPUT_DIR, "model_parameters.csv"), row.names = FALSE)

# ------------------------------ FIGURES ------------------------------------
class_levels <- sort(unique(representative_classes))
class_colors <- setNames(
  c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377")
    [seq_along(class_levels)], class_levels
)

compound_tree <- hclust(as.dist(sqrt(squared_distances(representatives,
                                                       representatives))),
                        method = "average")
compound_label_classes <- setNames(representative_classes, rownames(representatives))
pretty_dendrogram <- color_dendrogram_labels(
  as.dendrogram(compound_tree), compound_label_classes, class_colors
)
dendrogram_pdf <- file.path(OUTPUT_DIR, "nonlinear_compound_dendrogram.pdf")
pdf(dendrogram_pdf, width = 12, height = 8)
par(mar = c(10, 4, 4, 2) + 0.1)
plot(pretty_dendrogram,
     main = "Compound similarity in multiscale wavelet-feature space",
     ylab = "Robust-scaled Euclidean distance", xlab = "", sub = "")
mtext("Library compound", side = 1, line = 7.5)
mtext("Class colors are annotations only and were not used to build the tree",
      side = 1, line = 9, cex = 0.85)
legend("topright", legend = class_levels, col = class_colors[class_levels],
       pch = 15, title = "Class annotation", bty = "n")
dev.off()

probability_long <- do.call(rbind, lapply(seq_len(nrow(open_probs)), function(i) {
  data.frame(peak_id = peak_id[i], outcome = colnames(open_probs),
             probability = open_probs[i, ])
}))
top_peaks <- assignments$peak_id[order(assignments$unknown_probability)[1:25]]
probability_long <- probability_long[probability_long$peak_id %in% top_peaks, ]
probability_long$peak_id <- factor(probability_long$peak_id, levels = rev(top_peaks))
probability_figure <- ggplot(probability_long,
                            aes(probability, peak_id, fill = outcome)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(class_colors, Unknown = "#777777"),
                    name = "Outcome") +
  labs(
    title = "Provisional Betula class probabilities",
    subtitle = "Twenty-five peaks with the strongest support from the library",
    x = "Probability", y = "Betula peak"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "betula_open_set_probabilities.pdf"),
       probability_figure, width = 10, height = 9)

stability_figure <- ggplot(
  assignments,
  aes(bootstrap_outgroup_frequency,
      bootstrap_modal_class_support,
      color = bootstrap_modal_class,
      size = bootstrap_within_class_switching)
) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = class_colors, name = "Modal neighbor class") +
  scale_size_continuous(range = c(1.5, 5), name = "Within-class switching") +
  labs(
    title = "Class-neighbor support versus outgroup behavior",
    subtitle = "Each Betula peak was inserted separately in 250 bootstrap trees",
    x = "Frequency classified as a dendrogram outgroup",
    y = "Bootstrap support for the modal neighboring class"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "betula_insertion_stability.pdf"),
       stability_figure, width = 10, height = 7)

switching_figure <- ggplot(
  assignments,
  aes(bootstrap_modal_compound_support,
      bootstrap_modal_class_support,
      color = bootstrap_modal_class,
      size = bootstrap_within_class_switching)
) +
  geom_abline(slope = 1, intercept = 0, color = "grey70", linewidth = 0.7) +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = class_colors, name = "Modal neighbor class") +
  scale_size_continuous(range = c(1.5, 5), name = "Within-class switching") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "Does an unstable compound neighbor remain in the same class?",
    subtitle = "Points above the diagonal switch compounds while retaining class support",
    x = "Bootstrap support for modal compound neighbor",
    y = "Bootstrap support for modal class neighbor"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "compound_vs_class_neighbor_support.pdf"),
       switching_figure, width = 9, height = 8)

class_support_long <- do.call(rbind, lapply(seq_len(nrow(bootstrap_support$class_support)),
                                           function(i) {
  data.frame(peak_id = peak_id[i],
             class = colnames(bootstrap_support$class_support),
             support = bootstrap_support$class_support[i, ])
}))
peak_order <- assignments$peak_id[order(
  assignments$bootstrap_outgroup_frequency,
  -assignments$bootstrap_modal_class_support
)]
class_support_long$peak_id <- factor(class_support_long$peak_id,
                                     levels = rev(peak_order))
class_support_figure <- ggplot(
  class_support_long, aes(class, peak_id, fill = support)
) +
  geom_tile() +
  scale_fill_gradient(low = "#F7FBFF", high = "#2166AC",
                      limits = c(0, 1), name = "Bootstrap support") +
  labs(
    title = "Betula class-neighbor bootstrap support",
    subtitle = "Rows are ordered from library-like to persistent outgroup behavior",
    x = "Class of nearest library compound",
    y = "Betula peak"
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 5.5),
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "betula_class_neighbor_support_heatmap.pdf"),
       class_support_figure, width = 9, height = 16)

# PCA is used only to display the learned nonlinear feature space.
pca_display_features <- rbind(
  representatives,
  betula_features
)
pca_display <- prcomp(pca_display_features, center = TRUE, scale. = FALSE)
pca_variance <- 100 * pca_display$sdev^2 / sum(pca_display$sdev^2)
pca_data <- data.frame(
  PC1 = pca_display$x[, 1],
  PC2 = pca_display$x[, 2],
  source = c(
    rep("Library compounds", nrow(representatives)),
    rep("Betula", nrow(betula_features))
  )
)
pca_colors <- c(
  "Library compounds" = "#4477AA",
  "Betula" = "#228833"
)
pca_figure <- ggplot(pca_data, aes(PC1, PC2, color = source, shape = source)) +
  geom_point(alpha = 0.65, size = 2) +
  scale_color_manual(values = pca_colors, name = "Spectrum source") +
  labs(
    title = "Library and Betula spectra in wavelet-feature space",
    subtitle = "PCA is used for visualization only; matching uses the full feature space",
    x = paste0("PC1 (", round(pca_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_variance[2], 1), "%)")
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold")) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = "none")
ggsave(file.path(OUTPUT_DIR, "wavelet_feature_space.pdf"),
       pca_figure, width = 10, height = 7)

# Six representative one-at-a-time insertions: high support, boundary, outgroup.
example_indices <- unique(c(
  order(assignments$library_support_probability, decreasing = TRUE)[1:2],
  order(abs(assignments$library_support_probability - 0.5))[1:2],
  order(assignments$library_support_probability)[1:2]
))
pdf(file.path(OUTPUT_DIR, "example_betula_insertions.pdf"), width = 16, height = 9)
par(mfrow = c(2, 4), mar = c(7, 3, 4, 1) + 0.1)
for (index in example_indices) {
  combined <- rbind(representatives, betula_features[index, ])
  rownames(combined)[nrow(combined)] <- paste0("BETULA ", peak_id[index])
  tree <- hclust(as.dist(sqrt(squared_distances(combined, combined))),
                 method = "average")
  label_classes <- c(
    setNames(representative_classes, rownames(representatives)),
    setNames("Betula", paste0("BETULA ", peak_id[index]))
  )
  example_colors <- c(class_colors, Betula = "#000000")
  dendrogram <- color_dendrogram_labels(as.dendrogram(tree),
                                        label_classes, example_colors)
  plot(dendrogram,
       main = paste0(peak_id[index], "  |  library support = ",
                     round(assignments$library_support_probability[index], 2)),
       ylab = "Distance", xlab = "", sub = "")
}
plot.new()
legend("center",
       legend = c(class_levels, "Betula peak"),
       col = c(class_colors[class_levels], "#000000"),
       pch = 15, ncol = 2, bty = "n", title = "Label type")
dev.off()

cat("Nonlinear compound-first model complete.\n")
cat("Grouped leave-one-compound-out class accuracy:",
    mean(validation$results$correct), "\n")
cat("Betula open-set label counts:\n")
print(table(assignments$provisional_open_set_label))
cat("Peaks with stronger class than compound bootstrap support:",
    sum(assignments$bootstrap_within_class_switching >= 0.20), "\n")
cat("Results written to:", normalizePath(OUTPUT_DIR), "\n")
}

if (sys.nframe() == 0L) {
  message(
    "02_nonlinear_compound_similarity.R is a shared helper module. ",
    "Run RUN_AGILENT_DAD.R instead."
  )
}
