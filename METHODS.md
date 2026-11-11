# Raw HPLC-DAD spectral alignment and library-analogue workflow

## Technical methods reference

This document describes the computational methods implemented by the production
workflow when original Agilent ChemStation HPLC-DAD `.D` folders are available
for every biological sample. It is a technical reference for understanding,
auditing, and reproducing the implemented analysis.

The workflow aligns homologous sample peaks using retention time and spectral
shape, quantifies each aligned feature using integrated DAD response, and only
then compares its representative UV spectrum with authenticated library
standards. Retention time helps align peaks measured under the common sample
HPLC method, but it is never used to select a library compound, class, or
chromophore. Because sample compounds are not expected to be identical to the
limited set of standards, chemical results are reported as spectral analogues
or structural tendencies rather than confirmed identities.

The active library contains 23 authenticated compounds with three independent
standard runs per compound. Raw standard `.D` files were processed during
library construction and validation, while the curated CSV spectra serve as
identity anchors and safe active point estimates. The compact validated active
library is packaged for ordinary deployment; raw standards are needed only for
an explicit maintainer rebuild. Raw spectra are retained as auxiliary within-
standard evidence after QC; raw and CSV spectra from the same injections are
never appended as independent observations.

## Controlled vocabulary and glossary

### The three analytical levels

Most apparent terminology conflicts disappear when every term is assigned to
one of three levels:

```text
one-sample observation          cross-sample quantitative unit          chemical interpretation
peak detection          ->      aligned feature                 ->      spectral analogue
```

- A **peak detection** exists in one sample.
- An **aligned feature** links compatible detections across samples and becomes
  one column in a quantitative matrix.
- A **spectral analogue** is assigned only after alignment and describes which
  authenticated standard has the most similar UV spectrum.

These terms must not be substituted for one another. In particular, an aligned
feature is not automatically one known compound, and the nearest analogue does
not prove compound identity.

### Peak detection and alignment terms

| Preferred term | Definition | Important distinction |
|---|---|---|
| **Chromatographic profile** | The baseline-corrected root-mean-square DAD signal across wavelengths as a function of retention time. | It is the trace used for peak picking, not an individual wavelength chromatogram. |
| **Peak detection** | One local chromatographic maximum passing the permissive profile-height and local-prominence gates in one sample. | Use “detection” when referring to a sample-level row. Avoid calling it an aligned feature or compound. |
| **Core detection** | A detection also passing the stricter chromatographic and spectral-quality criteria required to seed alignment. | “Core” describes suitability for alignment, not chemical confidence. Preferred documentation term: **landmark-eligible detection**. Existing tables retain `core_detection`. |
| **Retention-time landmark** | A high-confidence cross-sample peak correspondence used to estimate sample RT correction. | A landmark supports correction; it is not necessarily a final feature or library match. |
| **Aligned feature** | A cross-sample group containing at most one detection per sample and satisfying RT-spectral alignment rules. | This is the quantitative column and preferred term for the final homologous unit. “Homologous peak” may be used descriptively, but **aligned feature** is less likely to imply proven chemical identity. |
| **Feature representative spectrum** | Median normalized spectrum of the detections assigned to one aligned feature, with uncertainty intervals retained separately. | This is sample-derived and is distinct from a library compound representative. |
| **Recurrent feature** | An aligned feature observed in at least the minimum number of samples required to act as a reconciliation target. | Recurrence measures sample occurrence, not annotation certainty. |
| **Rare feature** | An aligned feature observed in few samples. | Rare does not mean erroneous; it may represent real presence-absence variation or limited detection sensitivity. |
| **Unresolved local detection** | A peak detected in a sample that could not be assigned safely to one aligned feature. | It is excluded from confident abundance, marked in the matching assignment-uncertainty matrix, and retained quantitatively in the detection audit and feature-wise upper bounds. |
| **No confidently assigned signal** | No signal was assigned uniquely to that sample-feature cell. | It is represented by zero in the primary matrix; consult the uncertainty flag to distinguish an unresolved candidate from no qualifying local detection. It is not proof of biological absence. |
| **Prevalence** | Fraction of all input samples containing an assigned detection for an aligned feature. | Prevalence combines biological occurrence and analytical detectability. |

### Alignment reconciliation terms

| Preferred term | Definition | Important distinction |
|---|---|---|
| **Retained-uniqueness policy** | Primary genus-level policy that reviews only groups found in one to three samples as possible displaced duplicates. | It more strongly preserves rare features. It is not a claim that every retained feature is unique chemistry. |
| **Expanded-rescue policy** | Sensitivity policy that reviews a broader range of low-prevalence groups using the same evidence gates. | It expands review scope; it does not weaken the assignment thresholds. |
| **Reconciled detection** | Preferred prose for a detection moved from a low-prevalence candidate group into a larger recurrent feature after all gates pass. | Existing code and files use “rescued detection.” “Reconciled” better communicates evidence-based reassignment rather than recovery of a known truth. |
| **Initial alignment subgroup** | A provisional group produced by strict complete-linkage alignment that may later be reunited with an equivalent group. | Prefer this over the vague phrase “complete-linkage fragment.” |
| **Probable split-alignment pair** | Two mutually exclusive, RT-close, spectrally equivalent features judged to be an algorithmic split of one homologous unit and merged automatically. | This is an alignment correction, not evidence for two chemical variants. Existing files retain “split feature.” |
| **Resolved chromatographic doublet** | Two close features observed together in at least one sample and therefore protected from automatic merging. | Co-occurrence provides direct evidence for two resolved peaks, although their molecular identities may remain unknown. |
| **Unresolved neighboring-feature family** | Close features whose assignment cannot be resolved confidently. They remain separate in the complete matrix and are summed only for sensitivity analysis. | It is not a chemical family, class, chromophore family, or confirmed set of variants. |
| **Provisional alignment fragment** | A low-prevalence feature that may duplicate a recurrent feature but did not pass every reassignment gate. | Retained for review and excluded from the strict-primary matrix; not synonymous with a rare feature. |

### Library interpretation terms

| Preferred term | Definition | Important distinction |
|---|---|---|
| **Library standard** | An authenticated compound measurement supplied as reference data. | A standard is known; a sample feature is not. |
| **Library replicate spectrum** | One spectrum or raw-file-derived spectrum from one standard replicate. | Replicates characterize within-standard variation. |
| **Identity-anchor spectrum** | Curated spectrum used to find the authenticated compound peak anywhere in its raw standard chromatogram. | It guides extraction but is not added as a second observation when the raw spectrum is accepted. |
| **Active library replicate** | Curated point-estimate spectrum accompanied by admitted raw-DAD auxiliary evidence when available. | Exactly three training rows per compound prevent duplicate weighting. |
| **Compound representative spectrum** | Median processed spectrum of the accepted replicates for one authenticated compound. | It is library-derived and distinct from the feature representative spectrum. |
| **Closest spectral analogue** | Library compound with the smallest ensemble spectral distance to a sample feature. | Always mathematically available, even for a poor match; therefore it is not by itself an accepted annotation. |
| **Accepted analogue set** | Standards close enough to the best match to remain plausible under the calibrated tie rule. | It represents unresolved library neighbors, not multiple compounds proven to occur in the sample. |
| **Spectral-equivalence group** | Library standards that cannot be separated reliably by the current UV spectra and distance resolution. | It is a property of the library evidence, not a structural class or sample cluster. |
| **Chemical class** | Broad descriptive chemical category assigned to authenticated standards. | It is attached after compound comparison and is not used to train spectral similarity. |
| **Fine chromophore group** | Specific structure-informed UV-active motif associated with library compounds. | More specific than a broad chromophore family, but still an interpretation rather than a full molecular identity. |
| **Broad chromophore family** | Higher-level grouping of related UV-active motifs. | Used when fine-group or class evidence is not sufficiently consistent. |
| **Class-first hybrid label** | Reporting hierarchy that uses supported class evidence first and supported chromophore evidence as fallback. | It is not a fourth predictive model. |
| **Sample-only spectral family** | An unsupervised cluster of ambiguous or outside-coverage sample features with similar spectral shapes. | It is a discovery aid and must not be confused with a chemical class, chromophore family, or spectral-equivalence group. |

### Evidence and uncertainty terms

| Preferred term | Definition | Important distinction |
|---|---|---|
| **Ensemble spectral distance** | Weighted combination of six complementary shape distances; smaller values indicate greater similarity. | It is a dissimilarity measure, not a probability. |
| **Relative similarity weight** | Kernel-normalized support for one standard relative to the other supplied standards. | Existing CSV fields use “probability,” but these values are not posterior probabilities of identity. Use **relative similarity weight** in prose. |
| **Library coverage** | Open-set decision asking whether a feature resembles the spectral domain represented by the authenticated standards. | Separate from which standard is nearest. A feature can have a clear nearest standard and still be outside coverage. |
| **Calibrated compatibility p-value** | Conformal or empirical p-value measuring compatibility with a reference distribution. | It is not the probability that the compound or class is correct. Existing fields retain `conformal_p_value` for reproducibility. |
| **Bootstrap support** | Fraction of wavelength/replicate perturbations returning the same compound, class, chromophore, or alignment result. | It measures stability to the specified perturbation, not truth or biological prevalence. |
| **Branch-neighborhood support** | Frequency with which an inserted sample spectrum returns to the same library tree neighborhood. | A feature may move among spectrally equivalent standards yet retain high neighborhood support. |
| **Assignment margin** | Difference between the best and next-best candidate assignment scores. | It measures uniqueness of the preferred assignment, not absolute fit. |
| **Strong analogue support** | Highest practical tier after coverage, separation, method agreement, and stability gates pass. | Still an analogue-level result, not confirmed identity. |
| **Moderate analogue support** | Library-covered result passing the calibrated moderate evidence gates. | More cautious than strong support and not equivalent to a class prediction. |
| **Ambiguous** | Within library coverage, but no unique analogue or structural interpretation passes the required stability and separation rules. | Ambiguous is different from outside library coverage. |
| **Outside library coverage** | Insufficient compatibility with the current standard domain. | It means “not represented adequately by this library,” not “noise” or “not a real compound.” |

### Quantitative output terms

| Preferred term | Definition | Recommended use |
|---|---|---|
| **Complete matrix** | All aligned features, including provisional and unresolved-review columns. | Audit, rare-feature review, and sensitivity analysis. |
| **Primary-filtered matrix** | Complete matrix excluding provisional alignment fragments. | Intermediate view retained for transparency; not the default biological matrix. |
| **Strict-primary matrix** | Primary-filtered matrix also excluding members of unresolved neighboring-feature families. | Recommended input for primary biological models. |
| **Family-summed sensitivity matrix** | Stable feature columns plus one summed column for each unresolved neighboring-feature family. | Sensitivity analysis when uncertain neighboring assignments may represent one underlying signal. |
| **Assignment-uncertainty flag matrix** | Matrix matching the abundance table exactly, with `1` when an unresolved detected peak could belong to a cell and `0` otherwise. | Pair with the confident matrix to distinguish unresolved candidates from cells with no qualifying detection. |
| **Feature-wise upper-bound matrix** | Confident abundance plus every unresolved signal that could individually belong to a feature. | Analyze one feature at a time as a sensitivity limit; never sum columns because one unresolved peak can contribute to multiple candidate upper bounds. |
| **Integrated DAD-area response** | Background-corrected absorbance integrated across the chromatographic peak and wavelength range. | Relative quantitative response; do not call it concentration without compound-specific calibration. |

### Machine-readable field names

Some filenames and CSV columns use shorter internal labels. The table below
maps those machine-readable names to the technical terms used in this guide.

| Machine-readable wording | Preferred technical term |
|---|---|
| `peak` when referring to a matrix column | aligned feature |
| `core_detection` / core peak | landmark-eligible detection |
| rescued peak or rescued detection | reconciled detection |
| split feature | probable split-alignment pair |
| competition family | unresolved neighboring-feature family |
| compound probability / class probability | relative compound/class similarity weight |
| conformal support | calibrated compatibility p-value or conformal prediction-set membership, whichever the field actually reports |
| homologous peak | aligned feature, unless explaining the biological intention of alignment |

This vocabulary changes interpretation and writing, not the underlying
algorithms or machine-readable schema.

## Reproducing the analysis

### Required inputs

Store one Agilent DAD folder per biological sample using this layout:

```text
data/
  Genus_Agilent/
    Sample_001.D/
      DAD1.uv
    Sample_002.D/
      DAD1.uv
```

The public `library/` folder contains the validated active spectral table,
curated identity-anchor spectra, a compound class dictionary, a compound-to-
chromophore mapping, and raw-standard QC summaries. Authenticated raw standard
`.D` files are maintained separately and are required only to rebuild or extend
the library. The sample folder name before `_Agilent` is supplied as the dataset
name. No maximum-absorbance table or exported landscape CSV is used in this
route.

### Commands

For a genus-level study with substantial true presence-absence variation, run
retained uniqueness as the primary analysis:

```bash
Rscript VALIDATE_LIBRARY.R
Rscript RUN_AGILENT_DAD.R Rhododendron \
  --alignment-policy=retained_uniqueness
```

Then run expanded rescue as an alignment-sensitivity analysis:

```bash
Rscript RUN_AGILENT_DAD.R Rhododendron \
  --alignment-policy=expanded_rescue
```

The runner executes raw-data alignment, library-only spectral annotation,
pseudochromatogram construction, output validation, and organization. A run is
complete only when both validators report `Collaborator-handoff validation
passed`. Begin with `Genus_RawDAD_<Policy>/00_START_HERE.csv`, followed by
`03_Final_Results/00_READ_FIRST/READ_FIRST_output_order.csv`.

# Computational methods

The following sections document the implemented calculations. Instrument model,
column, mobile phases, gradient, flow rate, injection volume, DAD acquisition
interval, and biological sampling design are external experimental metadata and
cannot be inferred from these scripts.

## Analysis design and separation of alignment from annotation

Raw HPLC-DAD data were analyzed in two distinct stages. First, chromatographic
peaks were detected and aligned among biological samples using retention time
and full-spectrum similarity. Second, one representative spectrum from each
aligned feature was compared with a fixed library of authenticated standards.
Library compounds and their class or chromophore labels did not influence peak
detection, retention-time correction, or cross-sample feature construction.
Conversely, retention time was excluded from all library-similarity,
probability, clustering, and chemical-label decisions. This separation avoided
using retention times collected under a sample HPLC method as evidence for
standards measured under other chromatographic conditions.

## Raw-data import and wavelength harmonization

Agilent ChemStation `DAD1.uv` files were decoded with `chromConverter`. Each
file was represented as a retention-time-by-wavelength absorbance matrix. Input
coordinates and absorbance values were required to be finite, and retention
times were required to increase strictly. All binary cubes were decoded in a
preflight step before alignment. Files that could not be decoded were recorded
in `raw_agilent_input_exclusions.csv` and omitted from every matrix rather than
encoded as biological zeros; analyses stopped when more than 5% of discovered
sample folders were unreadable. Spectra were harmonized to the fixed
library wavelength coordinates from 190 to 900 nm using bounded linear
interpolation. Endpoint carrying was allowed only when a boundary mismatch did
not exceed one native or target wavelength interval; inputs with wider gaps or
inadequate wavelength resolution were rejected and recorded in the grid audit.

The 190-900 nm range was retained for sample peak detection, extraction,
alignment, and integrated DAD quantification. Chemical analogue comparisons
used the information-bearing 190-500 nm portion shared with the current
authenticated library. Signal above 500 nm therefore contributed to detecting
and quantifying a peak but did not create evidence for a library identity.

## Baseline correction and chromatographic peak detection

A running-median baseline with a 0.80-min window was estimated independently at
each wavelength and subtracted; negative corrected values were truncated to
zero. A chromatographic profile was formed as the root-mean-square corrected
absorbance across wavelengths and smoothed over 0.035 min. Profile noise was
estimated robustly from the median absolute deviation of successive profile
differences.

Candidate local maxima required chromatographic profile signal-to-noise ratio
(SNR) of at least 4 and local-prominence SNR of at least 2.5. Nonmaximum
suppression retained the stronger candidate when maxima were less than 0.020
min apart. These permissive thresholds were selected from a prespecified local
sensitivity grid. In the 91-sample Rhododendron audit, they retained every peak
found by the earlier 5/3 thresholds while adding 88 candidates, a 3.45%
increase. They were not selected by maximizing peak count.

For each apex, chromatographic boundaries were estimated from full width at
half maximum and capped at 0.12 min on either side. A local spectral background
was estimated outside the peak within a 0.18-min neighborhood. The non-negative
background-corrected spectrum was integrated across the observed peak width.
The quantitative response was the further integral across 190-900 nm. It is an
integrated DAD-area response rather than apex height and remains relative rather
than molar because extinction coefficients can differ among compounds.

## Landmark-eligible detections and retention-time correction

Peak detection was more permissive than landmark selection. A detection could
seed alignment only when profile SNR was at least 8, local-prominence SNR was at
least 4, spectral-shape SNR was at least 4, and at least 2.5% of wavelengths
were informative above local background. High-quality peak pairs within 0.14
min and composite spectral distance 0.10 formed a landmark network.

Retention-time corrections were estimated relative to the network barycenter,
not a selected reference sample. A constrained affine correction was eligible
only when at least five landmarks spanned at least 1.5 min and it reduced error
relative to a shift-only model. Affine slopes were restricted to 0.97-1.03 and
intercept magnitude to at most 0.12 min. Other samples retained the pairwise
shift; no unconstrained global warping was applied.

## Presence-absence-aware feature alignment

Landmark-eligible detections were grouped by complete linkage using the maximum of normalized
retention-time and spectral dissimilarities. Within a core feature, corrected
retention-time separation could not exceed 0.060 min, composite spectral
distance could not exceed 0.15, and no feature could contain two detections
from the same sample. Disjoint fragments were reunited only when their sample
sets were nonoverlapping, consensus retention times differed by at most 0.035
min, consensus spectral distance was at most 0.08, and combined retention-time
span was at most 0.090 min.

Lower-quality detections could attach to a fixed core consensus only within
0.065 min, at spectral distance no greater than 0.21, and with a best-versus-
second assignment margin of at least 0.12. An unsafe local detection remained
unresolved rather than being converted to absence or forced into its nearest
feature.

## Low-prevalence reconciliation and alignment policies

Low-prevalence groups were reviewed for possible displacement from recurrent
features using a policy chosen before analysis. Under `retained_uniqueness`,
only source groups observed in one to three samples were eligible. Under
`expanded_rescue`, groups observed in up to 20 samples were eligible when at
least 20 samples were present; for smaller datasets the limit was
`ceiling(0.20 × N)`. Both policies used identical evidence thresholds.

A detection could enter a larger feature only when the target had been
recurrent in at least four samples before iteration, corrected RT difference
was at most 0.12 min, spectral distance was at most 0.07, normalized joint score
was at most 1, assignment margin was at least 0.12, wavelength-block bootstrap
support was at least 85%, and the target did not already contain that sample.
Consensus spectra were rebuilt after accepted moves until a complete iteration
accepted no reconciliation. A detection could move only once and newly
reconciled detections could not
create a target, preventing snowball merging.

Retained uniqueness is recommended as the primary policy across a chemically
diverse genus because rare metabolites and true species-level absences are
plausible. Expanded rescue asks how conclusions change when a broader set of
groups is tested as possible displaced fragments. It is a sensitivity analysis,
not a more accurate result by definition.

## Neighboring-feature reconciliation

Very close, mutually exclusive features were evaluated as possible fragments of
one aligned feature. Automatic merging required consensus RT difference no
greater than 0.035 min, spectral distance no greater than 0.12, combined RT
span no greater than 0.10 min, and at least 85% wavelength-block support. If any
sample contained both peaks, the pair was protected as a resolved doublet.
Inconclusive pairs stayed separate in the complete matrix, were excluded from
the strict matrix, and were supplied as a summed family for sensitivity work.

## Quantitative matrices and explicit assignment uncertainty

The complete matrix contains samples in rows, aligned features in columns, and
confidently assigned integrated DAD-area responses in cells. It contains no
missing abundance values. A positive value is signal assigned safely to that
feature. Zero means that no signal was confidently assigned to that cell; it is
not, by itself, proof of biological absence.

Assignment uncertainty is represented explicitly rather than encoded as `NA`.
A matching flag matrix contains `1` when a detected local peak could plausibly
belong to the sample-feature cell but could not be assigned uniquely, and `0`
otherwise. A feature-wise upper-bound matrix equals confident abundance plus
all unresolved candidate area that could belong to that feature. These upper
bounds support one-feature-at-a-time sensitivity analyses and must not be
summed across columns because one unresolved peak is represented in every
plausible destination column. The detection-level ambiguity audit retains each
unresolved peak's area and complete candidate-feature list.

The strict-primary matrix excludes provisional alignment fragments and members
of unresolved neighboring-feature families and is recommended for primary
biological models. The complete matrix supports review of rare chemistry, while
the family-summed matrix provides a sensitivity analysis for unresolved
neighbors. Feature filtering changes which columns enter a model; it does not
retroactively change the alignment.

The strict annotated handoff contains only universally populated interpretation
fields and the NA-free sample abundances. Conditional diagnostics that are not
applicable to every feature—such as tie-specific wavelength tests or a nearest
recurrent-feature distance—remain in the complete audit and numbered QC tables
rather than appearing as unexplained missing cells in the primary table.

## Representative feature spectra

For every detected peak, absorbance above a local wavelength-specific
background was integrated across the chromatographic peak boundaries and then
range-normalized across wavelength. Consequently, representative construction
used spectral shape rather than absolute abundance. A detected peak had already
passed the permissive profile-SNR and local-prominence-SNR gates of 4 and 2.5,
respectively, before it could contribute to an aligned feature. Stricter
profile-SNR, prominence-SNR, spectral-shape-SNR, and informative-wavelength
criteria were used to determine whether a detection could seed an alignment
core, although lower-quality detections could subsequently attach to a fixed
core only when both retention-time and spectral-distance requirements were
satisfied.

Each aligned feature was represented by the pointwise median of the normalized
spectra from all sample detections assigned to that feature. Thus, the method
did not select one sample as the chemical representative and did not weight the
median toward the most concentrated sample. With one detected sample, the
representative was necessarily that single normalized spectrum; with two
samples, the pointwise median equaled their pointwise mean; and with three or
more samples, the median reduced the influence of an atypical spectrum at an
individual wavelength. The detection with the largest composite quality score
was recorded as the representative source and supplied only a magnitude for
compatibility plots. Library comparison range-normalized the feature again, so
that magnitude and the identity of the strongest sample did not influence
spectral matching.

Pointwise interquartile and 10th-90th percentile intervals across all assigned
spectra were retained to display within-feature variation. Prevalence was the
number of samples containing an aligned detection divided by all input samples.
Only the median sample-derived spectrum entered the current library-distance
calculation; individual member spectra and the uncertainty interval did not
contribute independent votes to the reported compound support. This avoids
pseudoreplication and prevents frequently observed features from dominating the
library calibration, but it means that the current compound-support statistic
does not increase merely because many high-quality, concordant observations are
available.

## Authenticated spectral library and preprocessing

The current reference library contains 23 authenticated compounds with three
standard runs per compound. Within each standard chromatogram, peak selection
first evaluated a primary region within ±0.45 min of the explicit curated
retention time. When that region contained no credible candidate, or only a
very weak candidate relative to an overwhelmingly stronger shifted peak, a
±1.30-min recovery region was evaluated. The complete chromatogram remained a
diagnostic fallback. Candidates were required to be chromatographic local
maxima with adequate profile SNR; an explicit RT could prioritize candidates
but could not create or replace an observed peak. Candidate wavelength shapes
were scored equally against the corresponding explicit replicate spectrum and
the compound-median curated spectrum. Several running-median profile windows
were evaluated so narrow and broad standard peaks could be represented. A
normalized consensus spectrum was calculated across the high-signal portion of
the selected peak; a baseline-corrected apex spectrum was retained when
integration degraded anchor agreement.

All three raw replicates for a compound were retained as auxiliary evidence
only when they passed
signal, within-peak spectral-stability, identity-anchor agreement, between-
replicate concordance, and within-method peak-location-concordance criteria.
As an additional identity safeguard, the median raw consensus assigned to each
manifest label was compared with every authenticated CSV compound anchor. A
raw set was flagged when its assigned-anchor distance was at least 0.15, another
anchor was within 0.08, and the alternative improved the distance by at least
threefold. If standard A instead matched anchor B and standard B reciprocally
matched anchor A, the audit reported a possible two-compound naming switch.
Flagged raw evidence was withheld from admission pending review; close
off-diagonal matches without a large assigned-anchor disagreement were retained
as spectral-equivalence ties rather than treated as labeling errors.
Absolute agreement with the curated CSV retention time was not required because
the two standard sources could have been acquired using different HPLC methods.
The curated CSV
spectra remained active point estimates because their leave-one-replicate-out
nearest-versus-second separation was clearer. Raw spectra were not appended,
preventing the same injections from being counted twice. Retention time, peak width, and the
standard HPLC method were retained as extraction-QC metadata only and were not
used in biological-sample annotation.

Library standards alone estimated compound
representatives, normalization parameters, metric scaling and weights,
acceptance thresholds, and open-set coverage. Sample spectra did not update or
retrain the model.

Absorbance from 190 to 500 nm was range-normalized per spectrum to emphasize
shape rather than concentration and smoothed with a second-degree, seven-point
Savitzky-Golay filter. Haar-wavelet features summarized nonlinear multiscale
structure and were centered and scaled using library medians and median
absolute deviations. Replicates were collapsed within compounds by pointwise
medians. Chemical class was not used to fit compound similarity; class and
chromophore dictionaries were attached only after compound evidence was found.

## A priori chromophore annotation of library compounds

Chromophore labels were assigned to the authenticated library before analysis
and were not learned from the biological samples or from the library
dendrogram. For each standard, the known molecular structure was inspected to
identify the dominant UV-active π-system expected to determine its absorbance
shape between 190 and 500 nm. Assignment considered aromatic ring architecture,
the extent of conjugation, coupling of aromatic rings to alkene or carbonyl
groups, nitrogen-containing heteroaromatic systems, oxygenated flavonoid
scaffolds, furan-containing diterpenoids, and weak nonconjugated aliphatic
unsaturation. The assignment represents the dominant spectrally informative
motif rather than an inventory of every functional group in the molecule.

Each of the 23 compounds received exactly one fine chromophore group and one
broad chromophore family in `compound_chromophore_dictionary.csv`. Sixteen fine
groups described relatively specific motifs, such as hydroxycinnamate
conjugation, a flavanone aromatic carbonyl, or a xanthine carbonyl heterocycle.
Structurally related fine groups were then collapsed into seven broad families:
nitrogen-containing heteroaromatics, simple single-ring aromatics,
hydroxybenzoates, hydroxycinnamates, polyoxygenated flavonoid aromatics,
diterpenoid-associated chromophores, and weak unsaturated
aliphatic/lipid-like chromophores. Fine- and broad-level confidence were
recorded as high or moderate, and a written structural rationale was retained
for every assignment. Lower confidence was used where absorption was expected
to be weak, localized, or influenced by more than one structural element.

These labels were interpretive metadata rather than supervised training
classes. Sample spectra were first compared with individual compound standards.
Only after compound-level spectral support was calculated were the fixed
compound-to-chromophore mappings used to summarize fine- or broad-chromophore
tendencies. Thus, a reported chromophore label means that the sample spectrum
resembled standards carrying that structural motif; it does not establish the
presence of the exact library compound or independently confirm molecular
structure.

## Future quantitative standard-curve extension

No concentration-series standards were available for the present analysis;
therefore, DAD responses were interpreted as relative integrated areas and no
absolute concentration, recovery, limit of detection, or limit of quantitation
was estimated. The empty `standard_curve_manifest.csv` establishes the planned
input contract without enabling calibration. Future rows will record compound,
raw data path, replicate identifier, known concentration and unit, injection
volume, analytical role, and notes.

Quantitative calibration should be fitted separately for each authenticated
compound because molar absorptivity and detector response differ among
chromophores. Concentration levels, blanks, and replicate injections should be
acquired under the same HPLC-DAD method as the samples to be quantified. The
calibration response must use the same peak-integration and wavelength-response
definition as the sample matrix. Curve form and working range should be chosen
from residual behavior rather than assumed; linearity, heteroscedasticity,
back-calculated concentration error, and lack of fit should be reviewed before
reporting concentrations. Until that module is implemented and validated, the
pipeline writes `standard_curve_readiness.csv` with calibration disabled.

## Nonlinear spectral-distance ensemble

Every feature was compared with every compound representative using cosine
distance, mean-centered correlation distance, shifted correlation allowing two
wavelength intervals, first-derivative correlation distance, distance in robust
Haar-wavelet space, and displacement of the three principal absorbance maxima.
Metrics were scaled by library-derived nearest-neighbor distance and truncated
to prevent domination. Unsupervised weights were derived from between-compound
separation relative to within-compound replicate variation and normalized to
sum to one.

Kernel-normalized compound and class values describe relative similarity among
the supplied standards; they are not posterior probabilities of molecular
identity. A high relative value can coexist with rejection from the library
domain when all supplied standards are poor matches.

## Calibration, open-set coverage, and robustness

Compound conformity was calibrated by leaving out individual standard
replicates and rebuilding the compound representative. Class conformity
excluded the held replicate's entire compound, preventing same-compound
recognition from inflating class performance. Open-set coverage left out all
replicates of one compound at a time and treated it as unseen. These folds
calibrated minimum distance, best-to-second and best-to-mean-other ratios, and
agreement among distance representations. Features with library-domain upper-
tail empirical p-value no greater than 0.10 were outside library coverage.

Robustness used 1,000 iterations that resampled library replicates and omitted
one contiguous wavelength block: 190-230, 230-270, 270-320, 320-400, or 400-500
nm. Wilson intervals summarized compound support. Four additional seeds with
250 iterations each evaluated Monte Carlo stability. Sample spectra were also
inserted into library dendrograms, with wavelength-block resampling used to
quantify branch-neighborhood stability.

## Practical analogue and structural labels

Strong analogue support required library-domain acceptance, calibrated
best-versus-second separation, method agreement, agreement between full-data
and bootstrap-modal compounds, and at least 80% exact bootstrap support.
Moderate support used less stringent calibrated separation and agreement with
at least 60% exact support. Library-covered features failing these rules were
ambiguous. Outside-coverage features inherited no class or chromophore label.

Compound names were interpreted as closest authenticated spectral analogues.
When a unique analogue was unstable, accepted-set or bootstrap agreement could
support a class, fine-chromophore, or broad-chromophore tendency. The class-first
hybrid label displays supported class evidence and otherwise falls back to
supported chromophore information; it is not another classifier.

## Exploratory analyses and pseudochromatograms

Feature spectra were clustered without library labels using scaled wavelet
features and feature-resampling consensus clustering. Ambiguous and outside-
coverage features were additionally grouped into sample-only spectral families
with uncertainty envelopes for literature comparison. These groups are
hypothesis-generating. Two-standard nonnegative mixture fits only flagged
possible co-elution and were not evidence that the named components occurred.

Genus pseudochromatograms placed features at consensus retention times and
scaled their shapes using observed widths and mean integrated DAD area. Class,
chromophore, and class-first hybrid views used the same features. Ambiguous and
outside-coverage features retained neutral colors. Pseudochromatograms are
descriptive summaries; exact abundance values remain in the tables.

## Reproducibility

Randomized library procedures used fixed seeds beginning with 20260806. Raw
alignment used 100 wavelength-block bootstrap repetitions; library annotation
used 1,000 wavelength-block/replicate iterations, four additional 250-iteration
seed checks, and 250 tree-resampling iterations. Every run recorded input and
grid audits, alignment policy, parameters, and downstream handoff. Software
versions should be added from the environment used for the final frozen run.
