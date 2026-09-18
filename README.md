# GCCFP food plant distributions

Geographic distribution of food plant taxa from the World Checklist of Food
Plants (WCFP) 2026, per country and per TDWG level 3 / level 4 area.

A taxon is counted as present in a place when **2 or more of 3 independent
sources** agree: raw occurrence points, GRIN-Global distributions, and
WCVP/TDWG distributions.

## Layout

```
GCCFP_food_plant_distributions/
  Spatial_distribution_of_wcfp_taxa_methods.R   the workflow - run this
  inputs/wcfp_taxa_distribution/                all input data (528 MB)
  outputs/                                      everything the run produces
```

## Running it

> **TODO — write this section.** Should cover: how to launch, expected runtime
> and output size, and the two operational gotchas (redirect R's temp directory
> off the C: drive; run from a copy of the script, since Rscript reads the file
> incrementally and editing it mid-run corrupts the run).
>
> Until then the same material is in `WORKFLOW_NOTES_open_issues.txt` under
> "ENVIRONMENT FOR A CLEAN RUN", and the script's own header block lists the
> packages, inputs and outputs.

## What the outputs contain

Figures below are from the run of 2026-09-18 in this project.

**Inputs reconciled**

| | |
|---|---|
| WCFP taxa in the plantlist | 26,633 |
| matched to a WCVP `plant_name_id` | 25,445 |
| with any WCVP TDWG level 3 distribution | 25,445 |
| GRIN records | 222,446 (199,390 native / 23,056 non-native) |
| GRIN region codes expanded to countries | 18 of 21 → 22,788 records from 3,893 rows |
| GRIN codes with no country mapping, dropped | 007, 014, 035 (282 rows) |
| countries with occurrence data | 235 |

*(coverage table — countries / L4 / L3 areas and row counts — added when the
run completes)*

**The 17 files**

*Workbooks* — one sheet per area, listing every taxon confirmed there by ≥2
sources, with its authority, IUCN Red List code and category, and which
sources confirmed it. A `Summary` sheet gives the per-area taxon count.

*Interactive maps* — choropleth shaded by taxon richness. Hover an area for its
name and count; click for the full taxa table, with a button to download that
table as `.xlsx`. Variants:

- `country_map`, `L4area_map`, `L3area_map` — one layer, taxa confirmed by ≥2
  sources. The two TDWG maps carry the other resolution as a toggleable
  boundary overlay.
- `by_source_*` — four mutually exclusive layers: each source alone, plus the
  confirmed-by-≥2 result. Selecting one at a time is deliberate; the layers are
  not additive.
- `by_redlist_*` — three schemes (grouped / all categories / threatened-only)
  over the ≥2-confirmed taxa. **Each layer has its own colour scale**, so
  shading is not comparable between layers — Least Concern has thousands of
  taxa and Extinct has two, and a shared scale would render every threatened
  layer as a uniform pale wash.

*Static figures* — publication-style PNG/PDF of the country and level 3
choropleths, magma palette.

## Running your own taxa list

Set `USER_TAXA_FILE` near the top of the script to an `.xlsx`/`.csv` holding a
column of taxon names. Names may carry authors, odd spacing or a plain `x` for
the hybrid sign — they are normalised before matching. Outputs go to their own
subfolder under `outputs/`, and a match report lists anything that could not be
matched.

**Your taxa must be in the WCFP list.** Two of the three input sources are
WCFP-scoped exports, so a taxon outside WCFP can only ever be seen by one
source and can never meet the 2-of-3 rule. The matcher reports unmatched names
rather than letting them disappear.

## Known issues

- **`WORKFLOW_NOTES_open_issues.txt`**
- **`WORKFLOW_NOTES_ARCHIVE_fixed-issues.txt`**
