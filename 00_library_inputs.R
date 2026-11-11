## ========================================================================== ##
## VALIDATED ACTIVE LIBRARY INPUT CONTRACT
## ========================================================================== ##
## Production analyses use exactly three authoritative CSV inputs:
##
##   1. Absorbance_StandardLibrary_RawDAD.csv - three active replicate spectra
##   2. compound_class_dictionary.csv        - one descriptive class per compound
##   3. compound_chromophore_dictionary.csv  - one fine and broad chromophore
##                                              definition per compound
##
## The spectral model is still compound-first and unsupervised with respect to
## class/chromophore labels. Dictionaries are joined only for interpretation,
## validation summaries, and figure colors. Embedded class columns in an older
## spectral library are ignored and overwritten by the class dictionary.
##
## The repository ships the validated active spectral table so deployment does
## not require raw authenticated-standard `.D` files. Maintainers may rebuild
## it explicitly with `Rscript VALIDATE_LIBRARY.R --rebuild-library`; raw peak-
## consensus spectra are admitted only as auxiliary uncertainty/QC evidence
## after extraction, peak-purity, anchor-agreement, and replicate-concordance
## checks. Raw standard retention times and peak widths are metadata only.
##
## Deployment
##   Ordinary users supply the active spectral CSV and two dictionaries in
##   library/. Raw standards and identity anchors are needed only to rebuild or
##   extend that library.
## Individual file overrides are also supported through the environment
## variables documented below. All joins are case-insensitive but preserve the
## compound spelling supplied by the spectral library.
## ========================================================================== ##

LIBRARY_DIR <- Sys.getenv("UV_LIBRARY_DIR", unset = "library")

library_input_path <- function(environment_variable, default_filename) {
  override <- Sys.getenv(environment_variable, unset = "")
  if (nzchar(override)) override else file.path(LIBRARY_DIR, default_filename)
}

LIBRARY_FILE <- library_input_path(
  "UV_SPECTRAL_LIBRARY_FILE", "Absorbance_StandardLibrary_RawDAD.csv"
)
CLASS_DICTIONARY_FILE <- library_input_path(
  "UV_CLASS_DICTIONARY_FILE", "compound_class_dictionary.csv"
)
CHROMOPHORE_DICTIONARY_FILE <- library_input_path(
  "UV_CHROMOPHORE_DICTIONARY_FILE", "compound_chromophore_dictionary.csv"
)

library_required_files <- function() {
  c(
    spectral_library = LIBRARY_FILE,
    class_dictionary = CLASS_DICTIONARY_FILE,
    chromophore_dictionary = CHROMOPHORE_DICTIONARY_FILE
  )
}

library_compound_key <- function(x) {
  tolower(gsub("[[:space:]]+", " ", trimws(as.character(x))))
}

format_library_compound_name <- function(x) {
  vapply(x, function(value) {
    value <- trimws(as.character(value))
    if (tolower(value) == "p-cresol") return("p-Cresol")
    first_letter <- regexpr("[[:alpha:]]", value)[1]
    if (first_letter > 0) {
      substr(value, first_letter, first_letter) <-
        toupper(substr(value, first_letter, first_letter))
    }
    value
  }, character(1), USE.NAMES = FALSE)
}

# ---------------------- WAVELENGTH-GRID HARMONIZATION ---------------------
# The authenticated library defines the coordinates used by every spectral
# model. Sample instruments may record a denser or offset grid (for example,
# every 3 nm rather than every 5 nm), so exact column equality is unnecessarily
# restrictive. These helpers interpolate raw sample absorbance onto the fixed
# library grid before normalization, smoothing, distances, or clustering.
#
# Interpolation is deliberately conservative:
#   * coordinates must be unique and strictly ordered after numeric sorting;
#   * the native grid cannot contain a gap wider than twice the library step;
#   * interior values use ordinary linear interpolation;
#   * a target endpoint may lie just beyond the native endpoint by no more than
#     one native/library interval, in which case the nearest measured endpoint
#     is carried over; and
#   * larger coverage gaps stop the analysis rather than extrapolating a shape.
#
# `values` must have spectra in rows and native wavelengths in columns. The
# returned audit row is written by each production route so grid conversion is
# visible and reproducible rather than a hidden preprocessing step.
harmonize_wavelength_matrix <- function(values, native_wavelengths,
                                        target_wavelengths, input_id) {
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  native_wavelengths <- as.numeric(native_wavelengths)
  target_wavelengths <- as.numeric(target_wavelengths)
  tolerance <- 1e-8

  if (!length(native_wavelengths) || !length(target_wavelengths) ||
      ncol(values) != length(native_wavelengths)) {
    stop("Wavelength coordinates do not match the spectral matrix for: ",
         input_id)
  }
  if (any(!is.finite(native_wavelengths)) ||
      any(!is.finite(target_wavelengths)) || any(!is.finite(values))) {
    stop("Non-finite wavelength or absorbance value in: ", input_id)
  }

  native_order <- order(native_wavelengths)
  native_wavelengths <- native_wavelengths[native_order]
  values <- values[, native_order, drop = FALSE]
  target_order <- order(target_wavelengths)
  if (!identical(target_order, seq_along(target_wavelengths))) {
    stop("The fixed library wavelength grid is not increasing.")
  }
  if (length(native_wavelengths) < 2L ||
      any(diff(native_wavelengths) <= tolerance)) {
    stop("Native wavelengths must be unique for: ", input_id)
  }
  if (length(target_wavelengths) < 2L ||
      any(diff(target_wavelengths) <= tolerance)) {
    stop("The fixed library wavelength grid must contain unique coordinates.")
  }

  native_steps <- diff(native_wavelengths)
  target_steps <- diff(target_wavelengths)
  native_step <- stats::median(native_steps)
  target_step <- stats::median(target_steps)
  maximum_allowed_native_gap <- 2 * target_step
  if (max(native_steps) > maximum_allowed_native_gap + tolerance) {
    stop(
      "Native wavelength grid is too sparse for reliable interpolation in: ",
      input_id, ". Largest gap = ", round(max(native_steps), 4),
      " nm; allowed maximum = ", round(maximum_allowed_native_gap, 4), " nm."
    )
  }

  lower_edge_carry <- max(native_wavelengths[1] - target_wavelengths[1], 0)
  upper_edge_carry <- max(tail(target_wavelengths, 1) -
                            tail(native_wavelengths, 1), 0)
  edge_tolerance <- max(native_step, target_step)
  if (lower_edge_carry > edge_tolerance + tolerance ||
      upper_edge_carry > edge_tolerance + tolerance) {
    stop(
      "Native wavelength coverage is insufficient for the fixed library grid in: ",
      input_id, ". Native range = ", native_wavelengths[1], "-",
      tail(native_wavelengths, 1), " nm; required range = ",
      target_wavelengths[1], "-", tail(target_wavelengths, 1), " nm."
    )
  }

  target_match <- vapply(target_wavelengths, function(target) {
    candidates <- which(abs(native_wavelengths - target) <= tolerance)
    if (length(candidates) == 1L) candidates else NA_integer_
  }, integer(1))
  exact_target_coordinates <- !anyNA(target_match)
  exact_grid <- exact_target_coordinates &&
    length(native_wavelengths) == length(target_wavelengths) &&
    identical(target_match, seq_along(target_wavelengths))
  aligned <- if (exact_target_coordinates) {
    values[, target_match, drop = FALSE]
  } else {
    t(vapply(seq_len(nrow(values)), function(row_index) {
      stats::approx(
        x = native_wavelengths,
        y = values[row_index, ],
        xout = target_wavelengths,
        method = "linear",
        rule = 2,
        ties = "ordered"
      )$y
    }, numeric(length(target_wavelengths))))
  }
  rownames(aligned) <- rownames(values)
  colnames(aligned) <- as.character(target_wavelengths)

  audit <- data.frame(
    input_id = as.character(input_id),
    native_wavelength_min_nm = native_wavelengths[1],
    native_wavelength_max_nm = tail(native_wavelengths, 1),
    native_wavelength_count = length(native_wavelengths),
    native_median_step_nm = native_step,
    native_maximum_gap_nm = max(native_steps),
    analysis_wavelength_min_nm = target_wavelengths[1],
    analysis_wavelength_max_nm = tail(target_wavelengths, 1),
    analysis_wavelength_count = length(target_wavelengths),
    analysis_median_step_nm = target_step,
    interpolation_required = !exact_target_coordinates,
    lower_endpoint_carry_nm = lower_edge_carry,
    upper_endpoint_carry_nm = upper_edge_carry,
    alignment_method = if (exact_grid) {
      "Exact fixed-library grid"
    } else if (exact_target_coordinates) {
      "Exact fixed-library coordinates selected from native grid"
    } else if (lower_edge_carry > 0 || upper_edge_carry > 0) {
      "Linear interpolation with bounded nearest-endpoint carry"
    } else {
      "Linear interpolation to fixed-library grid"
    },
    stringsAsFactors = FALSE
  )
  list(values = aligned, audit = audit)
}

require_library_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
}

require_nonempty_library_values <- function(data, columns, label) {
  for (column in columns) {
    values <- trimws(as.character(data[[column]]))
    bad <- is.na(data[[column]]) | !nzchar(values)
    if (any(bad)) {
      stop(
        label, " contains blank values in `", column,
        "` at row(s): ", paste(which(bad), collapse = ", ")
      )
    }
  }
}

# A code must always mean the same thing. This catches common deployment errors
# such as reusing one fine code for two names or mapping one broad code to two
# incompatible families.
validate_code_definition <- function(data, code_column, definition_columns,
                                     label) {
  codes <- unique(data[[code_column]])
  for (code in codes) {
    local <- data[data[[code_column]] == code, , drop = FALSE]
    for (column in definition_columns) {
      values <- unique(local[[column]])
      if (length(values) != 1L) {
        stop(
          label, " maps `", code, "` to multiple values in `", column,
          "`: ", paste(values, collapse = "; ")
        )
      }
    }
  }
}

collapse_library_compounds <- function(mapping, code_column, code) {
  members <- sort(unique(mapping$Compound[mapping[[code_column]] == code]))
  paste(members, collapse = "; ")
}

derive_broad_chromophore_dictionary <- function(mapping) {
  codes <- unique(mapping$broad_chromophore_code)
  output <- do.call(rbind, lapply(codes, function(code) {
    local <- mapping[mapping$broad_chromophore_code == code, , drop = FALSE]
    data.frame(
      chromophore_code = code,
      chromophore_family = local$broad_chromophore_family[1],
      annotation_type = local$annotation_type[1],
      annotation_confidence = local$broad_annotation_confidence[1],
      structural_rationale = local$broad_structural_rationale[1],
      library_compounds = collapse_library_compounds(
        mapping, "broad_chromophore_code", code
      ),
      n_library_compounds = length(unique(local$Compound)),
      stringsAsFactors = FALSE
    )
  }))

  # Mixed spectral-equivalence groups need a conservative nonstructural label.
  # It is generated automatically and is never assigned to an authenticated
  # standard in the input mapping.
  if (!"DUV" %in% output$chromophore_code) {
    output <- rbind(
      output,
      data.frame(
        chromophore_code = "DUV",
        chromophore_family = "Deep-UV-dominant mixed-structure",
        annotation_type = "Spectral-only",
        annotation_confidence = "Low",
        structural_rationale = paste0(
          "Spectral similarity is dominated by nonspecific deep-UV shape ",
          "despite mixed dictionary chromophores"
        ),
        library_compounds =
          "None (automatic fallback for mixed equivalence groups)",
        n_library_compounds = 0L,
        stringsAsFactors = FALSE
      )
    )
  }
  rownames(output) <- NULL
  output
}

derive_fine_chromophore_dictionary <- function(mapping) {
  codes <- unique(mapping$fine_chromophore_code)
  output <- do.call(rbind, lapply(codes, function(code) {
    local <- mapping[mapping$fine_chromophore_code == code, , drop = FALSE]
    data.frame(
      fine_chromophore_code = code,
      fine_chromophore_group = local$fine_chromophore_group[1],
      broad_chromophore_code = local$broad_chromophore_code[1],
      annotation_confidence = local$fine_annotation_confidence[1],
      structural_rationale = local$fine_structural_rationale[1],
      library_compounds = collapse_library_compounds(
        mapping, "fine_chromophore_code", code
      ),
      n_library_compounds = length(unique(local$Compound)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(output) <- NULL
  output
}

# Return deterministic qualitative colors for any number of user-defined
# classes or chromophore groups. The first six preserve the established project
# palette; larger dictionaries automatically receive a complete HCL palette.
library_category_colors <- function(levels, base_colors = c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"
), palette = "Dark 3") {
  levels <- as.character(levels)
  number <- length(levels)
  if (!number) return(setNames(character(0), character(0)))
  colors <- if (number <= length(base_colors)) {
    base_colors[seq_len(number)]
  } else {
    grDevices::hcl.colors(number, palette)
  }
  setNames(colors, levels)
}

load_library_inputs <- function() {
  paths <- library_required_files()
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "The three-file library contract is incomplete. Missing: ",
      paste(unname(missing), collapse = ", ")
    )
  }

  spectra <- read.csv(
    paths[["spectral_library"]], check.names = FALSE,
    stringsAsFactors = FALSE
  )
  class_dictionary <- read.csv(
    paths[["class_dictionary"]], check.names = FALSE,
    stringsAsFactors = FALSE
  )
  chromophore_mapping <- read.csv(
    paths[["chromophore_dictionary"]], check.names = FALSE,
    stringsAsFactors = FALSE
  )

  require_library_columns(spectra, c("Compound", "Rep"), "Spectral library")
  require_nonempty_library_values(spectra, c("Compound", "Rep"),
                                  "Spectral library")
  wavelength_values <- suppressWarnings(as.numeric(sub(
    "nm$", "", names(spectra), ignore.case = TRUE
  )))
  wavelength_columns <- which(is.finite(wavelength_values))
  if (length(wavelength_columns) < 5L) {
    stop("Spectral library must contain at least five numeric wavelength columns.")
  }
  spectral_matrix <- as.matrix(spectra[, wavelength_columns, drop = FALSE])
  storage.mode(spectral_matrix) <- "double"
  if (any(!is.finite(spectral_matrix))) {
    stop("Spectral library wavelength cells must all be finite numeric values.")
  }

  spectra$Compound <- format_library_compound_name(spectra$Compound)
  spectral_keys <- library_compound_key(spectra$Compound)
  compound_levels <- unique(spectra$Compound)
  compound_keys <- library_compound_key(compound_levels)
  replicate_counts <- table(spectral_keys)
  if (any(replicate_counts < 2L)) {
    bad <- names(replicate_counts)[replicate_counts < 2L]
    stop(
      "Every compound needs at least two replicate spectra; missing replicate ",
      "support for: ", paste(compound_levels[match(bad, compound_keys)],
                              collapse = ", ")
    )
  }
  if (length(compound_levels) < 3L) {
    stop("At least three distinct library compounds are required.")
  }

  require_library_columns(
    class_dictionary, c("Compound", "class_name"), "Class dictionary"
  )
  if (!"class_detail" %in% names(class_dictionary)) {
    class_dictionary$class_detail <- class_dictionary$class_name
  }
  require_nonempty_library_values(
    class_dictionary, c("Compound", "class_name", "class_detail"),
    "Class dictionary"
  )
  class_dictionary$.compound_key <-
    library_compound_key(class_dictionary$Compound)
  if (anyDuplicated(class_dictionary$.compound_key)) {
    stop("Class dictionary must contain exactly one row per compound.")
  }
  extra_class_rows <- setdiff(class_dictionary$.compound_key, compound_keys)
  if (length(extra_class_rows)) {
    stop(
      "Class dictionary contains compound(s) absent from the spectral library: ",
      paste(class_dictionary$Compound[
        match(extra_class_rows, class_dictionary$.compound_key)
      ], collapse = ", ")
    )
  }
  class_match <- match(compound_keys, class_dictionary$.compound_key)
  if (anyNA(class_match)) {
    stop(
      "Class dictionary is missing library compound(s): ",
      paste(compound_levels[is.na(class_match)], collapse = ", ")
    )
  }
  class_dictionary <- class_dictionary[class_match, , drop = FALSE]
  class_dictionary$Compound <- compound_levels
  class_dictionary$.compound_key <- NULL
  if (length(unique(class_dictionary$class_name)) < 2L) {
    stop("At least two descriptive class categories are required.")
  }

  required_chromophore_columns <- c(
    "Compound", "fine_chromophore_code", "fine_chromophore_group",
    "broad_chromophore_code", "broad_chromophore_family",
    "annotation_type", "fine_annotation_confidence",
    "broad_annotation_confidence", "fine_structural_rationale",
    "broad_structural_rationale"
  )
  require_library_columns(
    chromophore_mapping, required_chromophore_columns,
    "Chromophore dictionary"
  )
  require_nonempty_library_values(
    chromophore_mapping, required_chromophore_columns,
    "Chromophore dictionary"
  )
  chromophore_mapping$.compound_key <-
    library_compound_key(chromophore_mapping$Compound)
  if (anyDuplicated(chromophore_mapping$.compound_key)) {
    stop("Chromophore dictionary must contain exactly one row per compound.")
  }
  extra_chromophore_rows <- setdiff(
    chromophore_mapping$.compound_key, compound_keys
  )
  if (length(extra_chromophore_rows)) {
    stop(
      "Chromophore dictionary contains compound(s) absent from the spectral ",
      "library: ",
      paste(chromophore_mapping$Compound[
        match(extra_chromophore_rows, chromophore_mapping$.compound_key)
      ], collapse = ", ")
    )
  }
  chromophore_match <- match(compound_keys, chromophore_mapping$.compound_key)
  if (anyNA(chromophore_match)) {
    stop(
      "Chromophore dictionary is missing library compound(s): ",
      paste(compound_levels[is.na(chromophore_match)], collapse = ", ")
    )
  }
  chromophore_mapping <- chromophore_mapping[
    chromophore_match, required_chromophore_columns, drop = FALSE
  ]
  chromophore_mapping$Compound <- compound_levels

  validate_code_definition(
    chromophore_mapping, "fine_chromophore_code",
    c(
      "fine_chromophore_group", "broad_chromophore_code",
      "fine_annotation_confidence", "fine_structural_rationale"
    ),
    "Chromophore dictionary fine-group definition"
  )
  validate_code_definition(
    chromophore_mapping, "broad_chromophore_code",
    c(
      "broad_chromophore_family", "annotation_type",
      "broad_annotation_confidence", "broad_structural_rationale"
    ),
    "Chromophore dictionary broad-family definition"
  )

  # Dictionary files are authoritative. These assignments intentionally
  # replace any legacy class columns embedded in the spectral measurements.
  class_by_key <- setNames(class_dictionary$class_name, compound_keys)
  detail_by_key <- setNames(class_dictionary$class_detail, compound_keys)
  spectra$Short_Class <- unname(class_by_key[spectral_keys])
  spectra$Known_Class <- unname(detail_by_key[spectral_keys])

  broad_dictionary <- derive_broad_chromophore_dictionary(
    chromophore_mapping
  )
  fine_dictionary <- derive_fine_chromophore_dictionary(chromophore_mapping)
  manifest <- data.frame(
    input_role = names(paths),
    file = unname(paths),
    md5 = unname(tools::md5sum(paths)),
    n_rows = c(nrow(spectra), nrow(class_dictionary),
               nrow(chromophore_mapping)),
    n_library_compounds = length(compound_levels),
    n_categories = c(
      NA_integer_, length(unique(class_dictionary$class_name)),
      length(unique(chromophore_mapping$fine_chromophore_code))
    ),
    validation_status = "Passed",
    stringsAsFactors = FALSE
  )

  list(
    spectra = spectra,
    class_dictionary = class_dictionary,
    chromophore_mapping = chromophore_mapping,
    broad_chromophore_dictionary = broad_dictionary,
    fine_chromophore_dictionary = fine_dictionary,
    input_manifest = manifest,
    paths = paths
  )
}
