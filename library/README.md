# Authenticated raw-DAD spectral library

Ordinary users need the three compact validated inputs below. Maintainers need
the additional raw/anchor components only when rebuilding or extending the
library.

| Source | Role |
|---|---|
| `Absorbance_StandardLibrary_RawDAD.csv` | Packaged active spectral library: three validated replicate spectra per compound. |
| `compound_class_dictionary.csv` | One descriptive chemical class per compound. |
| `compound_chromophore_dictionary.csv` | One fine and broad chromophore definition per compound. |
| `Library_Dfiles/<Compound>/**/*.D/DAD1.uv` | Maintainer-only: three authenticated raw standard injections per compound; excluded from Git. |
| `Absorbance_StandardLibrary_Test.csv` | Maintainer input: curated identity anchors and point-estimate spectra. |
| `standard_curve_manifest.csv` | Optional concentration-series injection manifest; currently header-only and calibration is disabled. |

`raw_dad_library_manifest.csv` is the authoritative raw-file inventory. Only
rows marked `standard_replicate` may train the library; method blanks are never
treated as standards.

Chromophore assignments are curated a priori from the known structure of each
authenticated compound, not inferred from sample spectra or dendrogram
membership. Each dictionary row identifies the dominant UV-active structural
motif, its broader family, separate confidence levels, and written structural
rationales. The label describes the motif expected to dominate spectral shape
rather than every functional group present in the molecule. These mappings are
attached only after compound-level spectral matching and therefore do not
supervise the similarity model.

`00_build_raw_dad_library.R` creates
`Absorbance_StandardLibrary_RawDAD.csv` when a maintainer explicitly requests a
rebuild. That active table always contains
exactly three point-estimate rows per compound. Raw peak consensus is retained
as auxiliary within-standard evidence only when all three raw replicates pass:

- adequate chromatographic signal;
- stable wavelength shape across the selected peak;
- agreement with the authenticated identity anchor;
- agreement among the three raw replicate spectra; and
- a consistent selected peak location across the three raw standard runs.

Absolute agreement with the curated CSV retention time is deliberately not an
admission criterion because the curated and raw standards may have been
acquired with different HPLC methods.

The curated CSV spectra remain active for every compound because they provided
clearer leave-one-replicate-out separation than the raw consensus spectra.
Passing raw spectra remain in the QC and uncertainty outputs; failed sets are
quarantined. Raw and CSV versions are never appended because they likely
represent the same injections and would otherwise be counted twice.

The builder also compares each compound-level raw consensus with every
authenticated CSV anchor. A large assigned-anchor disagreement combined with a
close, at least threefold better alternative is an identity conflict and blocks
raw-evidence admission. Reciprocal A-to-B and B-to-A conflicts are explicitly
reported as a possible naming switch. This check detected the historical
luteolin/naringenin raw-folder switch, which is documented in
`raw_dad_identity_corrections.csv` and has been corrected in the active
manifest.

Standard HPLC methods may differ. As a transparent extraction aid, the builder
first evaluates candidates within
±0.45 min of each curated authenticated RT. If that tight region lacks a
credible peak, it permits a wider ±1.30-min recovery search for a reproducible
method shift; the complete chromatogram remains a final diagnostic fallback.
Every selected candidate must be an observed chromatographic local maximum with
adequate profile SNR. The curated RT is shown separately in the atlas and can
never be plotted or analyzed as though it were an observed peak. Its absolute
offset from a raw candidate is diagnostic provenance only and cannot reject
otherwise reproducible raw auxiliary evidence.
Candidate wavelength shape is compared with both the corresponding explicit
replicate spectrum and the compound-median curated spectrum. Retention time and
peak width are used only to extract and audit raw standards and never enter
biological-sample classification.

## Validate or rebuild

```bash
Rscript VALIDATE_LIBRARY.R
```

This validates the packaged active library and does not require raw standards.
To rebuild first, use:

```bash
Rscript VALIDATE_LIBRARY.R --rebuild-library
```

Review these files under `raw_dad_derived/`:

| Output | Use |
|---|---|
| `raw_standard_compound_admission_decisions.csv` | See which compounds have admitted raw auxiliary evidence and which raw sets were quarantined. |
| `raw_standard_replicate_extraction_qc.csv` | Audit selected RT, width, SNR, purity, anchor distance, and nearest anchor for every raw replicate. |
| `raw_standard_peak_candidate_audit.csv` | Inspect every RT/spectral candidate considered, its search tier, explicit-replicate distance, and whether it was selected. |
| `raw_standard_identity_conflict_audit.csv` | Compare every raw consensus with its assigned anchor and nearest alternative; review any identity conflict or possible reciprocal naming switch. |
| `raw_standard_identity_conflict_audit.pdf` | Visualize assigned-anchor distance against the nearest alternative. Points above the diagonal favor the assigned identity; strong flagged points below it require review. |
| `standard_curve_readiness.csv` | Confirm whether quantitative calibration data are available and enabled. Calibration is currently disabled. |
| `all_raw_standard_peak_extractions.pdf` | Inspect every chromatogram. Red dashed lines show detected candidate apices; purple dotted lines show curated reference RTs. A separation between them is an explicit QC warning, not a second observed peak. |
| `raw_vs_authenticated_standard_spectra_atlas.pdf` | Compare raw and curated median spectra with replicate IQR for every compound. |
| `leave_one_replicate_out_library_comparison.csv` | Replicate-level compound/class recovery and separation statistics. |
| `library_version_validation_summary.csv` | Compact baseline-versus-active validation summary. |
| `library_version_validation_comparison.pdf` | Compare held-out compound and class accuracy. |
| `library_version_separation_comparison.pdf` | Compare nearest-versus-second-nearest separation; lower is clearer. |

## Add or replace a standard

1. Add three authenticated `.D` folders beneath a folder whose name exactly
   matches the curated compound name.
2. Add one row per raw file to `raw_dad_library_manifest.csv` and mark the role
   `standard_replicate`.
3. Add three curated identity-anchor/point-estimate rows to
   `Absorbance_StandardLibrary_Test.csv` on the existing wavelength grid.
4. Add exactly one matching row to each dictionary. For the chromophore
   dictionary, document one dominant fine motif, one broad family, confidence
   at both levels, and a concise structural rationale. Reuse an existing group
   only when the dominant UV-active π-system is genuinely comparable; create
   a new fine group when the structural motif is materially different.
5. Run `Rscript VALIDATE_LIBRARY.R` and inspect the admission and atlas outputs.
6. Rerun every genus, because changing the active library changes its empirical
   distances, thresholds, and fingerprint.

## Future concentration calibration

`standard_curve_manifest.csv` is a header-only template until concentration-
series data are supplied. Add one row per injection with the compound, raw data
path, concentration value and unit, injection volume, replicate, role, and any
notes. Standard curves must be compound-specific and acquired under the same
HPLC-DAD method used for the samples being quantified. A future calibration
module must use the same response definition as the sample matrix and validate
working range, residuals, and back-calculated error before enabling absolute
concentrations. The present pipeline does not silently treat relative DAD area
as concentration.

Class and chromophore dictionaries annotate compound-level results after
spectral matching; they do not supervise the distance model. Retention-time
metadata columns are allowed but are ignored by every library similarity step.
