## Validate the compact active library without sample analyses.
##
## Default:
##   Rscript VALIDATE_LIBRARY.R
##
## Maintainer rebuild from raw authenticated-standard `.D` files:
##   Rscript VALIDATE_LIBRARY.R --rebuild-library
##
## Another library folder:
##   UV_LIBRARY_DIR=/path/to/library Rscript VALIDATE_LIBRARY.R

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 1L ||
    (length(arguments) == 1L && arguments != "--rebuild-library")) {
  stop("Usage: Rscript VALIDATE_LIBRARY.R [--rebuild-library]")
}
rebuild_library <- identical(arguments, "--rebuild-library")
if (rebuild_library) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"), "00_build_raw_dad_library.R"
  )
  if (!identical(status, 0L)) stop("Raw DAD library build failed.")
} else {
  message("Validating the packaged active library without rebuilding raw standards.")
}

# Identity conflicts are warnings during the build so collaborators can review
# the audit table, but a deployment validation must not silently certify a raw
# library that still contains an unresolved strong label conflict.
identity_audit_file <- file.path(
  "library", "raw_dad_derived", "raw_standard_identity_conflict_audit.csv"
)
if (file.exists(identity_audit_file)) {
  identity_audit <- read.csv(identity_audit_file, stringsAsFactors = FALSE,
                             check.names = FALSE)
  if (any(identity_audit$identity_conflict_flag)) {
    stop(
      "Unresolved raw-standard identity conflict(s): ",
      paste(identity_audit$Compound[identity_audit$identity_conflict_flag],
            collapse = ", "),
      ". Review the raw manifest and authenticated anchors."
    )
  }
} else if (rebuild_library) {
  stop("Raw-standard identity conflict audit was not created.")
}

standard_curve_status_file <- file.path(
  "library", "raw_dad_derived", "standard_curve_readiness.csv"
)
source("00_library_inputs.R")

bundle <- load_library_inputs()
cat("Library validation passed.\n")
cat("Spectral standards:", nrow(bundle$spectra), "replicate spectra\n")
cat("Compounds:", length(unique(bundle$spectra$Compound)), "\n")
cat("Classes:", length(unique(bundle$class_dictionary$class_name)), "\n")
cat(
  "Fine chromophore groups:",
  length(unique(bundle$chromophore_mapping$fine_chromophore_code)), "\n"
)
cat(
  "Broad chromophore families:",
  length(unique(bundle$chromophore_mapping$broad_chromophore_code)), "\n\n"
)
if (file.exists(identity_audit_file)) {
  cat("Unresolved raw-standard identity conflicts: 0\n")
}
if (file.exists(standard_curve_status_file)) {
  standard_curve_status <- read.csv(
    standard_curve_status_file, stringsAsFactors = FALSE, check.names = FALSE
  )
  cat(
    "Standard-curve rows available:",
    standard_curve_status$calibration_rows_available[1], "\n"
  )
  cat(
    "Quantitative calibration enabled:",
    standard_curve_status$quantitative_calibration_enabled[1], "\n\n"
  )
}
print(bundle$input_manifest, row.names = FALSE)
