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

Open the `.Rproj`, then run the script. A full run is **~37 minutes** and
writes **17 files, ~1.1 GB** to `outputs/`. Nothing is fetched from the
internet.

The script's own header block lists the required packages (with an
`install.packages()` line to paste), every input file and what it feeds, and
every output. Read that first.

From a terminal, run it detached and redirect R's temp directory off the C:
drive — Windows disk cleanup has killed a run mid-flight before, and this
workflow writes a scratch folder on every workbook save:

```powershell
$env:TMP = $env:TEMP = $env:TMPDIR = 'D:\Rtemp'
Start-Process 'C:\Program Files\R\R-4.5.3\bin\x64\Rscript.exe' `
  -ArgumentList '"D:\GCCFP_food_plant_distributions\Spatial_distribution_of_wcfp_taxa_methods.R"' `
  -RedirectStandardOutput run.log -RedirectStandardError run.err
```

Run from a **copy** of the script, not the file you are editing: Rscript reads
the file incrementally, so editing it mid-run corrupts the run.

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

## Provenance and known issues

The plantlist is a patched copy — two IUCN Red List categories present in the
published figshare release were blank in the working copy. See
`inputs/wcfp_taxa_distribution/WCFP_plantlist_final_2026-06-02_redlist-patched_PROVENANCE.txt`.

- **`WORKFLOW_NOTES_open_issues.txt`** — read this before using the numbers.
  Open data-quality issues, and a section on claims the data does *not*
  support.
- **`WORKFLOW_NOTES_ARCHIVE_fixed-issues.txt`** — issues already fixed. Read
  before "tidying" anything in the script that looks odd: several non-obvious
  choices are deliberate and were expensive to arrive at, and undoing one
  silently loses data.

Three caveats that most often get stated wrongly:

- **France ranks 3rd** only because French Guiana, Guadeloupe, Martinique,
  Réunion and Mayotte fold into the France polygon — Natural Earth has no
  separate polygon for any of them. Tropical species from three continents
  count as French. A ranked table needs a footnote.
- **Two thirds of GRIN records are not placed at a specific level 4 area.**
  They are smeared across every L4 area in their country, so L4 GRIN evidence
  is largely inferred rather than observed.
- **57% of confirmed taxa have never been IUCN-assessed.** Any statement of the
  form "N taxa in region X are threatened" describes the assessed minority.
