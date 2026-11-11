## Lightweight release check that requires no private raw `.D` files.
## Run from the repository root with: Rscript tests/smoke_test.R

required_files <- c(
  "RUN_AGILENT_DAD.R",
  "00_library_inputs.R",
  "01_raw_agilent_dad_alignment.R",
  "02_nonlinear_compound_similarity.R",
  "03_comprehensive_unsupervised_analysis.R",
  "08_raw_dad_genus_pseudochromatograms.R",
  "09_organize_raw_dad_outputs.R",
  "VALIDATE_LIBRARY.R",
  "VALIDATE_OUTPUTS.R",
  "library/Absorbance_StandardLibrary_RawDAD.csv",
  "library/Absorbance_StandardLibrary_IdentityAnchors.csv",
  "library/compound_class_dictionary.csv",
  "library/compound_chromophore_dictionary.csv"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing release files: ", paste(missing_files, collapse = ", "))
}

r_scripts <- list.files(
  ".", pattern = "[.]R$", recursive = TRUE, full.names = TRUE
)
r_scripts <- r_scripts[!grepl("(^|/)library/Library_Dfiles/", r_scripts)]
invisible(lapply(r_scripts, parse))
cat("Parsed", length(r_scripts), "R scripts successfully.\n")

validation_status <- system2(
  file.path(R.home("bin"), "Rscript"), "VALIDATE_LIBRARY.R"
)
if (!identical(validation_status, 0L)) {
  stop("Packaged library validation failed.")
}

cat("Release smoke test passed.\n")
