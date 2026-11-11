# Raw HPLC-DAD spectral analogue pipeline

[![Validate compact pipeline release](https://github.com/jordandowell/UVLibraryAnalysis/actions/workflows/validate.yml/badge.svg)](https://github.com/jordandowell/UVLibraryAnalysis/actions/workflows/validate.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`UVLibraryAnalysis` converts raw Agilent HPLC-DAD acquisitions into aligned,
quantitative spectral features and compares their wavelength shapes with a
fixed authenticated compound library. It is designed for chemically diverse
multi-sample studies where real presence–absence variation is expected.

Key features include:

- joint retention-time and full-spectrum alignment across biological samples;
- conservative handling of rare peaks, split alignments, resolved doublets,
  and sample-level assignment uncertainty;
- compound-first, open-set spectral analogue assessment without using sample
  retention time for library identification;
- class, chromophore, and class-first hybrid interpretations;
- NA-free primary abundance matrices with explicit uncertainty flags and
  feature-wise upper bounds; and
- publication-oriented QC tables, spectra, dendrograms, and
  pseudochromatograms.

This project now has one production workflow: original Agilent `.D` files from
biological samples are peak-picked and aligned across samples, then their
representative spectra are compared with an authenticated standard library.
The repository includes a compact validated spectral library. Raw standard `.D`
files are not required for ordinary use and are intentionally excluded from the
GitHub release; maintainers can rebuild the library explicitly when standards
are added or replaced.

The model is compound-first and open-set. A result names the closest spectral
analogue represented by the library; it does **not** confirm that the biological
feature is that exact compound. Chemical class and chromophore terms are
interpretations inherited from supported analogues.

## Documentation map

| Document | Purpose |
|---|---|
| [README.md](README.md) | Installation, quick start, output order, and practical interpretation. |
| [Technical methods reference](METHODS.md) | Detailed computational methods, terminology, thresholds, assumptions, and reproducibility information. |
| [Library guide](library/README.md) | Active-library contract, chromophore curation, raw-standard rebuilding, and adding standards. |
| [Data guide](data/README.md) | Required Agilent `.D` folder layout and wavelength-grid behavior. |
| [Contributing guide](CONTRIBUTING.md) | Validation expectations and scientific safeguards for changes. |

## Run the analysis

### Software installation

The production workflow requires R, `ggplot2`, and `chromConverter`:

```r
install.packages(c("ggplot2", "remotes"))
remotes::install_github("ethanbass/chromConverter")
```

Run commands from the repository root. Paths are relative to that folder.

Place biological sample acquisitions here:

```text
data/Genus_Agilent/
  Sample001.D/DAD1.uv
  Sample002.D/DAD1.uv
```

The packaged library under `library/` is ready for use. Run one of the two
explicitly named alignment policies:

```bash
Rscript RUN_AGILENT_DAD.R Genus --alignment-policy=retained_uniqueness
Rscript RUN_AGILENT_DAD.R Genus --alignment-policy=expanded_rescue
```

For a genus-level study, `retained_uniqueness` is the primary analysis because
it better protects genuinely rare chemistry. `expanded_rescue` is a sensitivity
analysis that attempts more low-prevalence reconciliations. Neither policy
changes the evidence required for a merge.

| Policy | Recommended role | Main tradeoff |
|---|---|---|
| `retained_uniqueness` | Primary genus-level analysis | Better protects genuinely rare chemistry; may retain more possible alignment fragments. |
| `expanded_rescue` | Alignment-sensitivity analysis | Reconciles more low-prevalence detections; may be less protective of true rare features. |

For the Rhododendron example used throughout the Methods guide:

```bash
Rscript RUN_AGILENT_DAD.R Rhododendron --alignment-policy=retained_uniqueness
Rscript RUN_AGILENT_DAD.R Rhododendron --alignment-policy=expanded_rescue
```

Raw biological `.D` files are intentionally not committed. The commands above
work after the corresponding `data/Rhododendron_Agilent/` folder is supplied.

Maintainers with authenticated raw standard acquisitions can rebuild and audit
the active library before a run:

```bash
Rscript VALIDATE_LIBRARY.R --rebuild-library
Rscript RUN_AGILENT_DAD.R Genus \
  --alignment-policy=retained_uniqueness --rebuild-library
```

Ordinary collaborators should omit `--rebuild-library`.

The output is written to `Genus_RawDAD_<Policy>/`. Open
`00_START_HERE.csv` first; the folders are numbered in interpretation order.

## What happens in one run

1. `00_build_raw_dad_library.R` searches chromatographic local maxima first
   within ±0.45 min of the curated authenticated retention time, then within
   ±1.30 min only when the tight region lacks a credible peak. Both the matching
   explicit replicate spectrum and the compound-median spectrum contribute to
   peak selection. A reference RT can guide the search but can never substitute
   for a detected chromatographic apex.
2. A normalized consensus is extracted across the high-signal peak region.
   This characterizes broad peaks without making peak width a predictor.
3. Raw replicates must pass signal, within-peak purity, anchor-agreement,
   replicate-concordance, and within-method peak-location-concordance checks
   before they are retained as auxiliary within-standard evidence. Absolute
   agreement with the curated reference RT is not an admission requirement.
   A separate all-against-all anchor audit blocks large identity conflicts and
   reports reciprocal conflicts as possible two-compound naming switches.
4. The curated spectra remain active point estimates because they showed
   clearer leave-one-replicate-out separation. Raw spectra are never appended
   as duplicate observations of the same injections.
5. Biological peaks are detected independently, aligned using retention time
   and spectral shape, and summarized into sample-by-aligned-feature matrices.
   Every binary DAD cube is decoded in a preflight check. Unreadable files are
   listed in `raw_agilent_input_exclusions.csv` and omitted, not zero-filled;
   the run stops if more than 5% of input folders are unreadable.
6. Retention time is then excluded. Representative wavelength shapes alone are
   compared with the active library using the nonlinear distance ensemble,
   leave-one-replicate-out conformity, leave-one-compound-out library-coverage
   calibration, wavelength bootstraps, and sensitivity analyses.
7. Final class, broad-chromophore, and class-first hybrid pseudochromatograms
   are drawn from the aligned quantitative data.

The standard runs may use different HPLC methods. Their retention times and
peak widths are therefore recorded only as extraction-QC metadata. The three
raw injections must locate a consistent peak relative to one another because
they share a method, but their absolute RT need not match the curated CSV
reference. Standard RTs are never compared with sample RTs and never enter
analogue, class, or chromophore support.

## Interpretation order

| Order | Output | Use |
|---:|---|---|
| 1 | `03_Final_Results/00_READ_FIRST/READ_FIRST_output_order.csv` | File-by-file interpretation guide. |
| 2 | `03_Final_Results/01_Primary_Tables/` | NA-free strict confident-abundance matrix, matching uncertainty flags, and annotated feature table. |
| 3 | `03_Final_Results/02_Primary_Figures/` | Sample-by-peak map and genus class, chromophore, and hybrid pseudochromatograms. |
| 4 | `03_Final_Results/03_Sensitivity_and_Review/` | Feature-wise abundance upper bounds, complete matrices, and family-summed sensitivity tables. |
| 5 | `01_Alignment_QC/` | Peak picking, RT correction, alignment, reconciliation, and spectral uncertainty. |
| 6 | `02_Library_Annotation/` | Compound analogue, open-set, bootstrap, clustering, and ambiguity evidence. |
| 7 | `library/raw_dad_derived/` | Raw-standard peak extractions, admissions, spectra, and library-version validation. |

The primary annotated table is deliberately concise: it contains the essential
alignment and library-interpretation fields plus every sample abundance. The
full conditional diagnostics remain in the complete review table and numbered
QC folders, so “not applicable” diagnostics do not fill the primary handoff
with `NA` values.

The primary abundance matrix contains confidently assigned integrated DAD area
and has no missing cells. A positive value is assigned signal; zero means no
signal was confidently assigned, not proven biological absence. When a local
peak could plausibly belong to multiple features, the same-shaped uncertainty-
flag matrix marks each candidate cell with `1`, and the feature-wise upper-bound
matrix adds that unresolved area to each candidate separately. Upper-bound
columns must never be summed because the same unresolved peak can appear as a
possible contribution to more than one feature. The complete matrix retains all
discovery features, while the strict matrix excludes features whose alignment
evidence is inadequate for a primary model.

## Practical naming rules

- Use the named compound only as a **closest spectral analogue**, never as a
  confirmed identity, unless independent chemical evidence is available.
- Use a **class tendency** when compound support is ambiguous but supported
  analogues converge on a chemical class.
- Use a **chromophore tendency** when the wavelength shape supports a shared
  UV-active motif more clearly than a broad class.
- Use the class-first hybrid label for concise figures: class when supported,
  otherwise chromophore, otherwise ambiguous or outside library coverage.
- Report **outside library coverage** when the peak fails the grouped open-set
  criterion. This means that the supplied standards do not adequately cover
  its shape; it does not mean the peak is noise.

## Core glossary

| Term | Meaning |
|---|---|
| Peak detection | One chromatographic peak in one biological sample. |
| Aligned feature | Detections judged to be one homologous RT-spectral unit across samples; one output matrix column. |
| Confident abundance | Integrated DAD area assigned safely to one aligned feature; the NA-free primary matrix contains these values. |
| Assignment-uncertainty flag | `1` when a detected local peak could belong to that sample-feature cell but could not be assigned uniquely; otherwise `0`. |
| Feature-wise upper bound | Confident abundance plus all unresolved candidate signal that could belong to one feature; use one feature at a time and never sum columns. |
| Active library replicate | One curated point-estimate spectrum, accompanied by admitted raw-DAD auxiliary evidence when available. Exactly three training rows are retained per compound. |
| Identity anchor | Curated spectrum used to locate the authenticated compound within its raw standard chromatogram. |
| Input exclusion | A sample whose binary DAD cube cannot be decoded. It is absent from all matrices and documented in the run manifest; it is never encoded as a zero. |
| Closest spectral analogue | Library compound with the most similar UV spectrum; not a confirmed molecular identity. |
| Relative similarity weight | Support relative to the supplied standards; not a posterior probability of identity. |
| Within/outside library coverage | Whether the feature resembles the spectral domain represented by held-out standards. |
| Reconciliation | Evidence-gated reassignment of a low-prevalence fragment to a recurrent aligned feature. |
| Unresolved neighboring-feature family | Nearby, mutually exclusive features retained separately in the complete matrix and summed only for sensitivity analysis. |

The full controlled vocabulary and computational details are documented in the
[technical methods reference](METHODS.md).

## Validation

```bash
Rscript VALIDATE_LIBRARY.R
Rscript RUN_AGILENT_DAD.R Rhododendron --alignment-policy=retained_uniqueness
Rscript VALIDATE_OUTPUTS.R Rhododendron_RawDAD_RetainedUniqueness
```

The runner validates the active library and the complete output before and
after organizing it. Every annotation folder stores the active library's
MD5 fingerprint so outputs cannot silently be interpreted with another library
version.

## If final organization fails

On Windows or inside OneDrive, an open CSV/PDF can prevent replacement during
the final organization step. The runner now checks existing organized folders
before sample processing and reports the exact locked path. Close Excel, PDF
and RStudio viewers (and pause cloud sync if necessary), then resume without
reprocessing samples:

```bash
Rscript 09_organize_raw_dad_outputs.R Genus_RawDAD_RetainedUniqueness
Rscript VALIDATE_OUTPUTS.R Genus_RawDAD_RetainedUniqueness
```

The run is complete only when validation passes and `00_START_HERE.csv` exists.

## Repository contents and data policy

The GitHub release contains the analysis scripts, compact validated library,
class/chromophore dictionaries, library QC summaries, documentation, tests, and
continuous-integration configuration. Biological raw acquisitions,
authenticated-standard raw `.D` files, and generated genus result folders are
excluded through `.gitignore`. Users retain responsibility for permissions,
de-identification, and appropriate sharing of their own data.

## Citation and license

Citation metadata are provided in [CITATION.cff](CITATION.cff). Please cite the
repository version or commit used in an analysis and separately cite the
software used to decode Agilent files. The code and documentation are released
under the [MIT License](LICENSE); library spectra and compound annotations
should also be attributed when reused.
