# Contributing

Thank you for improving the raw HPLC-DAD spectral analogue pipeline.

## Before opening a change

1. Do not commit biological raw `.D` acquisitions, authenticated-standard raw
   files, generated genus results, or identifiable sample metadata.
2. Keep retention time restricted to biological peak alignment. Library
   analogue, class, and chromophore support must continue to use wavelength
   shape only.
3. Preserve the distinction between aligned-feature confidence and chemical
   annotation confidence.
4. Add comments explaining scientific thresholds and any change to an output
   contract.

## Validation

From the repository root, run:

```bash
Rscript tests/smoke_test.R
```

For changes that affect alignment or annotation, also run both policies on a
local raw dataset and validate each output folder:

```bash
Rscript RUN_AGILENT_DAD.R Genus --alignment-policy=retained_uniqueness
Rscript RUN_AGILENT_DAD.R Genus --alignment-policy=expanded_rescue
Rscript VALIDATE_OUTPUTS.R Genus_RawDAD_RetainedUniqueness
Rscript VALIDATE_OUTPUTS.R Genus_RawDAD_ExpandedRescue
```

Raw-standard changes require the maintainer-only rebuild and visual review of
all raw-standard QC outputs:

```bash
Rscript VALIDATE_LIBRARY.R --rebuild-library
```

## Pull requests

Describe the scientific reason for the change, affected thresholds or output
files, datasets used for stress testing, and validation results. Do not describe
a nearest library analogue as a confirmed compound identity.
