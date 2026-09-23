# GCCFP food plant distributions

Where the food plant taxa of the World Checklist of Food Plants (WCFP) 2026
occur — per country, per TDWG level 3 area and per TDWG level 4 area.

A taxon counts as present in a place when **2 or more of 3 independent sources**
agree.

## Inputs

All local, under `inputs/`. Nothing is downloaded at run time.

### 1. Occurrence points (accessed 2025) — `data/occurrences_data/occurrences_land.rds`

5,064,165 georeferenced records of 26,126 taxa, already filtered to land.

| source | records | citation |
|---|---|---|
| GBIF observations | 3,839,495 | GBIF.org (2025). <https://www.gbif.org> |
| Genesys | 1,012,558 | Genesys PGR (2025). <https://www.genesys-pgr.org> |
| WIEWS | 195,235 | FAO (2025). *WIEWS — World Information and Early Warning System on Plant Genetic Resources for Food and Agriculture.* [Cited 4 August 2025]. <https://www.fao.org/wiews/> |
| GBIF living collections | 16,877 | GBIF.org (2025). <https://www.gbif.org> |

Of these, 1,210,790 come from genebanks and 13,880 from botanic gardens.

### 2. GRIN distributions (accessed 2026) — `data/GRIN_distributions_data/WebQueryGRIN-Global_distributions_combined.xlsx`

222,446 records from a GRIN-Global web query. Each carries a status code:
native (199,390), introduced (18,190), cultivated (3,910), and smaller counts
for absent, uncertain and other.

> Global Crop Diversity Trust, Bioversity International & USDA Agricultural
> Research Service (2026). *GRIN-Global Server* (Germplasm Resources
> Information Network). <https://www.grin-global.org/>

### 3. WCFP 2026 release — `data/…figshare_repo…/`

Khoury, Colin; Gora, Sarah L.; et al. (2026). *World Checklist of Food Plants
2026 (WCFP): Data.* figshare. <https://doi.org/10.6084/m9.figshare.31441615.v2>

| file | supplies |
|---|---|
| `WCFP.xlsx` | the taxon list — 26,622 taxa with `WCFP_ID`, accepted names, authorities and IUCN Red List categories |
| `WCFP_distribution_native_introduced.csv` | 376,373 taxon × TDWG level 3 records with native / introduced / extinct / doubtful status |

These two are located by **filename** anywhere under `inputs/data/`, so renaming
the folder they sit in will not break the run.

### 4. Boundaries (WGSRPD Edition 2, 2001) — `TDWG_L3_L4/`

`level3/level3.shp` (369 areas) and `level4/level4.geojson` (610 areas) — 368
and 609 once Antarctica is dropped, which is what the maps and the Results
table below count.

> **Scheme:** Brummitt, R.K. (2001). *World Geographical Scheme for Recording
> Plant Distributions*, Edition 2. Plant Taxonomic Database Standards No. 2.
> Hunt Institute for Botanical Documentation, for TDWG (Biodiversity
> Information Standards).
>
> **GIS dataset:** Royal Botanic Gardens, Kew (2001). *WGSRPD level 3 and
> level 4 boundaries.* Distributed by TDWG via
> <https://github.com/tdwg/wgsrpd> (last updated 2021). File metadata records
> a 2000 release date and a 2011 revision. The underlying administrative
> boundaries are © ESRI and GMI 1992–1997, used by permission.
>
> **Country outlines:** Natural Earth 1:50m Admin 0 (public domain), via the
> `rnaturalearth` 1.2.0 and `rnaturalearthdata` 1.0.0 packages, not a file.

Edition 2 is still the current version of the standard. TDWG's Geographical
Schemes Interest Group has task groups working on a revision, with no release
date announced.

**How TDWG areas become countries.** Each level 4 area in the Kew dataset
(Royal Botanic Gardens, Kew, 2001) carries an ISO 3166-1 code; those codes are
joined to Natural Earth country names to resolve a level 3 area to one country,
or to several where its level 4 children span more than one. Areas with no
level 4 children fall back to a direct name match against Natural Earth, then to
a manual override table for naming variants and multi-country composites.

**Area names shown on the maps come from the WCFP release, not the shapefile.**
The Kew build still carries Edition 2's 2001 political vocabulary, so 15 areas
differ — Zaïre, Swaziland, Czechoslovakia, Yugoslavia and Surinam among them.
The code is what every join uses, so this affects labels only. Set
`L3_NAME_OVERRIDES` near the level 3 read to correct an individual name.

**How the three sources are kept comparable.** Occurrence and GRIN taxa are
matched to `WCFP.xlsx` before use, so all three are scoped to the same taxon
list. A handful of synonyms are absent from the release and are excluded; the
run names every taxon it drops.

## Outputs

Writes **32 files** into a dated folder on the shared GCCFP Drive:

```
…/Food plant distributions/outputs/outputs_<date>/
  taxa_lists_and_summaries/        4 workbooks
  interactive_html_maps/
    equal_earth_maps/              10 maps, equal-area
  figures_static_maps/             3 figures, PNG + PDF
    colour_scheme_options/         the same 3 in 4 other palettes, PNG
```

Each run gets its own dated folder, so **no run overwrites an earlier one**.
`outputs_root` is a Google Drive shortcut-target path: Drive must be mounted
and the folder synced. The script proves the destination exists and is writable
before it reads a single input, so a missing mount costs seconds rather than an
hour.

**Runtime under an hour.** The exact figure will be filled in after the next
full run — the map set has just changed, so any number quoted now would be
stale.

Filenames follow one pattern:

```
<kind>_food-plants_[IUCN-<scheme>_]distributions_by-<unit>_<qualifier>.<ext>

kind       map | taxa-lists | counts | figure
unit       country | L3 | L4
qualifier  2-or-more-sources | by-data-source
```

Alternative-palette figures carry the palette name as a final suffix, e.g.
`figure_food-plants_distributions_by-L3_2-or-more-sources_magma.png`.

**Workbooks** — one sheet per area, listing every taxon confirmed there with
`WCFP_ID`, authority, IUCN category and which sources confirmed it. Every sheet
also names its own area on every row, so a sheet copied out of the workbook
still says where it is from: the country workbook carries `country` and
`country_iso3`, the level 3 workbook `l3_area_code` and `l3_area_name` plus each
taxon's level 4 sub-areas. A `Summary` sheet holds per-area counts. `counts_…`
holds per-source totals.

**Interactive maps** — choropleth shaded by taxon richness. Hover for name and
count; click for the full taxa table and an `.xlsx` download. Each has a **Find
a taxon** box, top right: type three or more letters, click a result, and every
area holding that taxon is outlined.

Popup titles carry the area's code in the scheme that map uses — *India
(ISO3: IND)* on the country maps, *India (L3: IND)* at level 3, *Austria
(L3-L4: AUT-AU)* at level 4 — and the downloaded sheet repeats it in its own
`Area code` column. The download also carries columns the on-screen table
leaves out, because they are long and repeat on every row: a level 3 area's
level 4 sub-areas, and a level 4 area's parent level 3 area.

- `…_2-or-more-sources` — one shaded layer, with the other TDWG resolution as a
  toggleable overlay. On the **L4 maps** that overlay is selectable: click
  inside a level 3 area for its own level 3 taxa table and its own download,
  separate from the level 4 one underneath. Ticked on, it lays a light neutral
  veil over the map so you can see which level 3 area you are in while the
  level 4 shading still reads through it, and it takes every click in its area
  — untick it to interact with level 4 again.
- `…_by-data-source` — four mutually exclusive layers: each source alone, plus
  the confirmed-by-≥2 result. They are not additive, hence one at a time.
- `…_by-L3-L4_…` — **two merged maps, kept alongside the separate by-L3 and
  by-L4 ones so the two arrangements can be compared.** Each carries both TDWG
  resolutions as selectable layers rather than one shaded level with the other
  as an overlay: `2-or-more-sources` has two layers, `by-data-source` has eight
  (four sources × two levels). The levels are on **separate colour scales**, so
  the legend follows the selected layer and shading is not comparable across
  them. Decide which arrangement you prefer and the other set can be dropped.
- `…IUCN-…` — two maps: categories grouped, or all nine separately. Layers run
  in IUCN severity order, most severe first: EX, EW, CR, EN, VU, NT, LC, DD, NE.
  On the grouped map **Threatened is CR + EN + VU**, the IUCN definition;
  Extinct and Extinct in the Wild are outcomes rather than threat categories, so
  they get their own layers.

  Each layer has **its own colour ramp, anchored on the official IUCN category
  colour** — black for Extinct, purple for Extinct in the Wild, red, orange and
  yellow for the threatened categories, green for Least Concern, grey for Data
  Deficient. Within a layer, darker means more taxa. The grouped map's Threatened
  ramp runs yellow → orange → red across its three categories. Because every
  layer has its own scale, **shading is not comparable between layers**.

  The lightest shade of every ramp is a visible tint, never white: white is
  reserved for "no taxa in this category", which is most of the map on the
  Extinct layer. Each map opens on its most informative layer — Threatened for
  the grouped map, Critically Endangered for the per-category one — rather than
  on Extinct, which holds 2 taxa.

**Projection.** All nine are Equal Earth (EPSG:8857), which is equal-area — the
right property for a map whose whole subject is how much richness sits where.
Leaflet's own default, Web Mercator, is not equal-area: it inflates high
latitudes so much that Greenland reads as comparable to Africa when it is about
a fourteenth of it. No Web Mercator maps are produced.

The maps carry a graticule and a **scale bar that re-measures as you zoom and
pan**, bottom left. Leaflet's own scale control cannot be used here: the map is
drawn with `L.CRS.Simple`, so Leaflet would measure raw map units and label them
metres. The conversion is exact instead — the geometry is projected to EPSG:8857,
whose units are metres, then divided by `EE_SCALE`, so one map unit is 100 km.
`EE_SCALEBAR_MAX_PX` sets how long the bar is allowed to get.

The bar is marked **approx.** because Equal Earth is equal-*area*, not
conformal: linear distance is not uniform across the map, so no single bar is
exactly right at every latitude. It is measured at the centre of the current
view.

The maps open slightly wider than the data, so the world is not pressed against
the edges of the pane; `EE_ZOOM_OUT` is that margin.

**Static figures** — three, one per spatial unit: country, level 3 and level 4.
Each shows the same thing its interactive `…_2-or-more-sources` twin shows, as
a single shaded layer with no toggles, and each is written as both PNG (300 dpi)
and PDF at 12 × 7 inches.

They are Equal Earth like the interactive maps, and use the same **YlGn**
palette, so a reader moving between the HTML and the figures sees one colour
scheme. The legend is a plain horizontal colour bar on a white background,
bottom left; a scale bar sits bottom right, marked approximate for the same
reason as the interactive one. Areas with no confirmed taxa are light grey,
distinct from white.

`colour_scheme_options/` holds the same three maps in **magma**, **viridis**,
**plasma** and **cividis** — twelve PNGs for picking a different scheme. All
four are colourblind-safe. They are PNG only, since they exist to be looked at
rather than placed in a document. To switch the main figures to one of them,
change `STATIC_PALETTE` near the bottom of the script.

**Antarctica** is omitted from every map — the continent has no confirmed taxa.
The subantarctic islands are not: the Falkland Is. (22 taxa), Macquarie,
Kerguelen, South Georgia and four more hold 47 records between them and stay in
both the maps and the workbooks.

## Running it

Open `GCCFP_food_plant_distributions.Rproj`, then:

```r
source("Spatial_distribution_of_food_taxa_methods.R")
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
