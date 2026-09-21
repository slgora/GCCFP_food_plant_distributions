# GCCFP food plant distributions

Where the food plant taxa of the World Checklist of Food Plants (WCFP) 2026
occur — per country, per TDWG level 3 area and per TDWG level 4 area.

A taxon counts as present in a place when **2 or more of 3 independent sources**
agree.

## Inputs

All local, under `inputs/`. Nothing is downloaded at run time.

### 1. Occurrence points — `data/occurrences_data/occurrences_land.rds`

5,064,165 georeferenced records of 26,126 taxa, already filtered to land.

| source | records |
|---|---|
| GBIF observations | 3,839,495 |
| Genesys | 1,012,558 |
| WIEWS (FAO World Information and Early Warning System) | 195,235 |
| GBIF living collections | 16,877 |

Of these, 1,210,790 come from genebanks and 13,880 from botanic gardens.

### 2. GRIN distributions — `data/GRIN_distributions_data/WebQueryGRIN-Global_distributions_combined.xlsx`

222,446 records from a GRIN-Global web query (USDA-ARS Germplasm Resources
Information Network). Each carries a status code: native (199,390), introduced
(18,190), cultivated (3,910), and smaller counts for absent, uncertain and
other.

### 3. WCFP 2026 release — `data/…figshare_repo…/`

Khoury, Colin; Gora, Sarah L.; et al. (2026). *World Checklist of Food Plants
2026 (WCFP): Data.* figshare. <https://doi.org/10.6084/m9.figshare.31441615.v2>

| file | supplies |
|---|---|
| `WCFP.xlsx` | the taxon list — 26,622 taxa with `WCFP_ID`, accepted names, authorities and IUCN Red List categories |
| `WCFP_distribution_native_introduced.csv` | 376,373 taxon × TDWG level 3 records with native / introduced / extinct / doubtful status |

These two are located by **filename** anywhere under `inputs/data/`, so renaming
the folder they sit in will not break the run.

### 4. Boundaries — `TDWG_L3_L4/`

`level3/level3.shp` and `level4/level4.geojson`, the TDWG World Geographical
Scheme for Recording Plant Distributions. Country outlines come from the
`rnaturalearth` package, not a file.

**How the three sources are kept comparable.** Occurrence and GRIN taxa are
matched to `WCFP.xlsx` before use, so all three are scoped to the same taxon
list. A handful of synonyms are absent from the release and are excluded; the
run names every taxon it drops.

## Outputs

**Runtime about 65 minutes**, writing 26 files, roughly 2.0 GB. Existing
outputs are overwritten. Roughly the last 20 minutes is the Equal Earth set.

```
outputs/
  taxa_lists_and_summaries/        4 workbooks
  interactive_html_maps/
    natural_earth_maps/            9 maps, Web Mercator
    equal_earth_maps/              the same 9, equal-area
  figures_static_maps/             4 PNG / PDF
```

Filenames follow one pattern:

```
<kind>_food-plants_[IUCN-<scheme>_]distributions_by-<unit>_<qualifier>.<ext>

kind       map | taxa-lists | counts | figure
unit       country | L3 | L4
qualifier  2-or-more-sources | by-data-source
```

**Workbooks** — one sheet per area, listing every taxon confirmed there with
`WCFP_ID`, authority, IUCN category and which sources confirmed it. The level 3
workbook also gives each taxon's level 4 sub-areas. A `Summary` sheet holds
per-area counts. `counts_…` holds per-source totals.

**Interactive maps** — choropleth shaded by taxon richness. Hover for name and
count; click for the full taxa table and an `.xlsx` download. Each has a **Find
a taxon** box, top right: type three or more letters, click a result, and every
area holding that taxon is outlined.

- `…_2-or-more-sources` — one shaded layer, with the other TDWG resolution as a
  toggleable overlay. On the **L4 maps** that overlay is selectable: click
  inside a level 3 area for its level 3 taxa table. While it is ticked on it
  takes every click in its area — untick it to interact with level 4.
- `…_by-data-source` — four mutually exclusive layers: each source alone, plus
  the confirmed-by-≥2 result. They are not additive, hence one at a time.
- `…IUCN-…` — grouped, all categories, or threatened only. **Each layer has its
  own colour scale**, so shading is not comparable between layers.

**Equal Earth maps** are the same nine files under identical names. Equal Earth
is equal-area; Web Mercator is not, and inflates high latitudes so much that
Greenland reads as comparable to Africa when it is about a fourteenth of it.
Those maps carry a graticule and have **no scale bar** — Leaflet cannot project
to Equal Earth, so the geometry is projected in R and the coordinates are
neither degrees nor metres.

**Antarctica** is omitted from every map — the continent has no confirmed taxa.
The subantarctic islands are not: the Falkland Is. (22 taxa), Macquarie,
Kerguelen, South Georgia and four more hold 47 records between them and stay in
both the maps and the workbooks.

## Running it

Open `GCCFP_food_plant_distributions.Rproj`, then:

```r
source("Spatial_distribution_of_wcfp_taxa_methods.R")
```

`inputs_dir` and `outputs_root` are absolute paths near the top of the script;
edit them only if you move the project.

Two things worth knowing before a long run:

- **Point R's temp directory off the C: drive.** Windows Automatic Maintenance
  can empty `%TEMP%` mid-run. Set `TMP`/`TEMP` to something like `D:/Rtemp`.
- **Don't edit the script while it runs.** `Rscript` reads the file
  incrementally, so an edit part-way through changes what still executes.

## Running your own taxa list

Set `USER_TAXA_FILE` near the top of the script to an `.xlsx` or `.csv` holding
a column of taxon names, and `USER_TAXA_COLUMN` if it is not the first one:

```r
USER_TAXA_FILE   <- "D:/my_taxa.xlsx"
USER_TAXA_COLUMN <- NULL     # NULL = first column
```

Names may carry authors, odd spacing, curly quotes or a plain `x` for the
hybrid sign — they are normalised before matching, which runs in tiers from
exact to authors-stripped. Every input name is reported with the tier that
matched it, and unmatched names are written out rather than dropped silently.

Outputs go to their own subfolder under `outputs/`, so a full run is never
overwritten and the two can be compared.

**Your taxa must be in the WCFP list.** Two of the three sources are
WCFP-scoped, so a taxon outside WCFP can only ever be seen by one source and can
never meet the 2-of-3 rule.

## Results

Run of 2026-09-21.

| | areas | (area, taxon) rows |
|---|---|---|
| countries | 235 | 190,731 |
| TDWG level 4 | 591 (within 360 level 3 areas) | 460,437 |
| TDWG level 3 | 360 of 368 | 274,739 |

**25,192** of 26,633 taxa are confirmed somewhere at level 3. The level 3 result
applies the 2-of-3 rule at level 3 rather than rolling up level 4, which
recovers 208 (area, taxon) pairs that no single level 4 area confirms on its
own — the rows the level 3 workbook marks `confirmed at level 3 only`.

**IUCN Red List** over those 25,192 — NE 14,521 · LC 8,677 · VU 582 · NT 475 ·
EN 443 · DD 360 · CR 129 · EW 3 · EX 2.

## Known issues

- `WORKFLOW_NOTES_open_issues.txt`
- `WORKFLOW_NOTES_ARCHIVE_fixed-issues.txt`
