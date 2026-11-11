# Biological sample inputs

The production workflow reads original Agilent ChemStation HPLC-DAD folders.
Create one folder per genus:

```text
data/Rhododendron_Agilent/
  JAD010603.D/
    DAD1.uv
  JAD010604.D/
    DAD1.uv
```

The `.D` folder name without `.D` becomes the sample identifier. Auxiliary
ChemStation files may remain in each folder; `DAD1.uv` is required. Run from the
project root:

```bash
Rscript RUN_AGILENT_DAD.R Rhododendron --alignment-policy=retained_uniqueness
```

Before peak picking, every `DAD1.uv` is decoded and checked for finite time,
wavelength, and absorbance values. An unreadable or truncated cube is written
to `01_Alignment_QC/00_Run_Manifest/raw_agilent_input_exclusions.csv` and is
absent from all quantitative matrices rather than being treated as a sample
with zero peaks. The analysis stops if more than 5% of input folders fail this
check, because that indicates an upstream export or transfer problem.

Sample wavelength grids need not exactly match the active library. Dense or
offset grids are conservatively interpolated to the library coordinates, and
the conversion is written to the output audit. Coverage gaps that would require
unsafe extrapolation stop the analysis.

Retention time is used to align peaks across biological samples. It is never
used to decide which library compound, class, or chromophore a feature resembles.

Only raw Agilent `.D` folders belong in `data/` for the current production
pipeline. The retired maximum-absorbance tables and exported wavelength-by-
retention-time landscapes are deliberately excluded so that users cannot
accidentally select an unsupported input route.
