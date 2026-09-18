# Project: GCCFP:Species Distribution (in countries, L3, L4)
# Main workflow
# Aims (2)                                                                           #
#     * 1. Be able to recognize how many food plant species are in any given         #
#          country, and which food plant species those are.                          #
#          Food plants = WCFP 2026                                                   #
#           (To do later (include CWR): Food plants = WCFP + CWR)                    #
# *    2. Be able to perform these counts/queries for any given list of food plants  #


## ------------------------------------------------------------------------------ ##
# Methods:
#
#  3 methods below.
#    Derive combined method (e.g. if species is in country in 2 or more datasets?)
#
# What comes from the WCFP 2026 paper's own code vs. what's custom to this
# script:
#  - Method 3 (WCVP/TDWG level 3 join) reproduces data_summaries.R section 8
#    "NATIVE / INTRODUCED FLAG" from the paper's data-and-code archive
#    (https://doi.org/10.6084/m9.figshare.31441615.v2): the join of
#    wcfp_plantlist -> WCVP plant_name_id -> WCVP distribution table, using
#    the same wcvp_names.csv/wcvp_distribution.csv snapshot files the paper
#    ships, and the same native/introduced/extinct/doubtful classification
#    logic. The "count all occurrence_status values, not just native" rule
#    used for richness counts also mirrors the paper's section 9a.
#  - Everything else is custom to this script, not from the paper: Method 1
#    (raw occurrence points), Method 2 (GRIN Taxonomy distributions), the
#    3-source "combined method" synthesis (>=2 of 3 sources), the static
#    ggplot country-level richness map, and the interactive leaflet map.
## ------------------------------------------------------------------------------ ##


## ============================================================================== ##
##  HOW TO RUN THIS WORKFLOW
##  Runtime ~37 min. Writes 17 files, ~1.1 GB. Nothing is read from the internet.
## ============================================================================== ##
##
## PACKAGES
##   Attached below : dplyr stringr tidyr readxl readr sf ggplot2 viridis scales
##                    grid cowplot rnaturalearth leaflet htmlwidgets
##   Attached later : openxlsx, at the first workbook write - left there on
##                    purpose so it cannot change which functions are masked in
##                    the several hundred lines that run before it
##   Used as pkg::  : base64enc htmltools ragg tibble  (installed, not attached)
##                    openxlsx:: is ALSO called before it is attached, by the
##                    user-taxa matcher - namespaced deliberately so it works there
##   NOT needed     : expowo, s2 - they appear in comments only.
##
##   install.packages(c("dplyr","stringr","tidyr","readxl","readr","sf","ggplot2",
##     "viridis","scales","cowplot","rnaturalearth","leaflet","htmlwidgets",
##     "openxlsx","base64enc","htmltools","ragg","tibble"))
##
## RUN YOUR OWN TAXA LIST
##   Set USER_TAXA_FILE just below inputs_dir to an .xlsx/.csv holding a column
##   of taxon names. Outputs go to their own subfolder, and a match report names
##   anything that could not be matched. See the block there for the one limit:
##   your taxa must be in the WCFP list.
##
## INPUTS - all in inputs_dir, all local
##   occurrences_land.rds                              Method 1: occurrence points
##   WebQueryGRIN-Global_distributions_combined.xlsx   Method 2: GRIN distributions
##   wcvp_names.csv  +  wcvp_distribution.csv          Method 3: WCVP/TDWG (pipe-sep)
##   WCFP_plantlist_final_2026-06-02.xlsx              the taxon list + Red List codes
##   level3/level3.shp  +  level4.geojson              TDWG WGSRPD boundaries
##   Country polygons come from the rnaturalearth PACKAGE, not a file.
##
## OUTPUTS - all written to outputs_dir. 17 files:
##   Workbooks (4)
##     WCFP_data_source_summary_counts.xlsx
##     WCFP_taxa-lists_by_country_2plus_sources.xlsx      236 sheets
##     WCFP_taxa-lists_by_L4area_2plus_sources.xlsx       592 sheets
##     WCFP_taxa-lists_by_L3area_2plus_sources.xlsx
##   Interactive maps (9)  ..._map_interactive.html
##     country_map / by_source_map                        country resolution
##     L4area_map  / by_source_L4area                     TDWG level 4
##     L3area_map  / by_source_L3area                     TDWG level 3
##     by_redlist_grouped_L3area / _categories_ / _threatened_
##   Static figures (4)    ..._v1_simple.png / .pdf
##     country_map, L3area_map
##
## ORDER OF PLAY
##   1 setup  2 Method 1 occurrences  3 Method 2 GRIN  4 Method 3 WCVP/TDWG
##   5 combine (>=2 of 3) + country outputs  6 level 4  7 level 3  8 Red List
##
## KNOWN ISSUES: outputs/WORKFLOW_NOTES_issues_and_fixes_2026-09-08.txt
## ============================================================================== ##


# Load required packages
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
if (!requireNamespace("cowplot", quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("leaflet", quietly = TRUE)) install.packages("leaflet")
if (!requireNamespace("htmlwidgets", quietly = TRUE)) install.packages("htmlwidgets")
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(readr)
library(sf)
library(ggplot2)
library(viridis)
library(scales)
library(grid)
library(cowplot)
library(rnaturalearth)
library(leaflet)
library(htmlwidgets)

# All input data used by this workflow is consolidated here (see
# copy_workflow_inputs.R, which copies everything from its original
# scattered locations into this one folder).
inputs_dir <- "D:/GCCFP_food_plant_distributions/inputs/wcfp_taxa_distribution"

## ============================================================================== ##
##  RUN YOUR OWN TAXA LIST - set USER_TAXA_FILE and run the script unchanged
## ============================================================================== ##
##
## Leave as NULL to run the whole WCFP list, exactly as before.
##
## Give it a path and the workflow runs on YOUR taxa instead: same three
## methods, same >=2-of-3 rule, same 17 outputs, written to their own folder
## under outputs/ so nothing of the full run is overwritten.
##
##   USER_TAXA_FILE <- "D:/my_project/my_taxa.xlsx"
##
## THE FILE: .xlsx or .csv, with a column of taxon names - that is all. Names
## may carry authors ("Musa acuminata Colla"), odd spacing, or "x" for the
## hybrid sign; they are normalised before matching. Set USER_TAXA_COLUMN if
## the column is not found automatically.
##
## THE ONE LIMIT, AND WHY. Your taxa must be IN the WCFP list. Two of the three
## input sources are WCFP-scoped exports - occurrences_land.rds is 100% WCFP
## names and the GRIN export is 99.9% - so a taxon outside WCFP can only ever
## be seen by Method 3 (WCVP) and can NEVER meet the >=2-of-3 rule. It would
## vanish from every output with no error. The matcher therefore reports
## unmatched names loudly and writes them to a file, rather than letting them
## disappear. Supplying wider occurrence and GRIN extracts is the only way to
## lift this.
## ============================================================================== ##

USER_TAXA_FILE   <- NULL   # path to your .xlsx/.csv, or NULL for the full WCFP list
USER_TAXA_COLUMN <- NULL   # name of the taxon-name column, or NULL to auto-detect

# Everything is written under here. A user-list run gets its own subfolder.
outputs_root <- "D:/GCCFP_food_plant_distributions/outputs"


# ---------------------------------------- #
#     Filter method 1. Occurrence data     #
# ---------------------------------------- #

# filter using occurrence data
# How many species have geo data?
# ---------------------------
# Load occurrences dataset
occurrences_data <- readRDS(file.path(inputs_dir, "occurrences_land.rds"))

## ------------------------------------------------------------------ ##
## Fold breakaway / unrecognised polygons into their parent state.
##
## Natural Earth maps de facto control, so it draws Somaliland as its own
## admin-0 polygon, separate from Somalia. None of the data sources used
## here do: GRIN reports "SOM", and WCVP/TDWG reports Somalia. The
## Somaliland polygon therefore never received any data and rendered as a
## white "no data" hole in the Horn of Africa - which on a richness map
## reads as "we looked and found nothing" rather than "this is Somalia".
##
## Folding unions the breakaway geometry into the parent's and drops the
## breakaway row, keeping the parent's attributes (ISO codes, continent).
## 30 of Natural Earth's 31 worldviews already code Somaliland as SOM, so
## this matches the majority treatment rather than inventing one.
##
## Add rows here to fold others - Northern Cyprus -> Cyprus and Siachen
## Glacier -> India are the obvious candidates, both currently left as-is.
## ------------------------------------------------------------------ ##
NE_FOLD_POLYGONS <- tribble(
  ~breakaway,   ~parent,
  "Somaliland", "Somalia"
)

fold_breakaway_polygons <- function(x) {
  for (i in seq_len(nrow(NE_FOLD_POLYGONS))) {
    b <- NE_FOLD_POLYGONS$breakaway[i]
    p <- NE_FOLD_POLYGONS$parent[i]
    if (!(b %in% x$admin) || !(p %in% x$admin)) next
    merged <- st_union(st_geometry(x[x$admin %in% c(b, p), ]))
    x <- x[x$admin != b, ]
    g <- st_geometry(x)
    g[which(x$admin == p)] <- merged
    st_geometry(x) <- g
  }
  x
}

# ---------------------------
# Spatial Data Preparation
# Load world land polygons (exclude Antarctica)
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") %>%
  fold_breakaway_polygons()
world_4326 <- st_transform(world, 4326)

# Convert occurrences to sf
occurrences_sf <- st_as_sf(
  occurrences_data,
  coords = c("longitude", "latitude"),
  crs = 4326, remove = FALSE)

# Spatial join: assign each occurrence to a country (land only)
occurrences_land <- st_join(occurrences_sf, world_4326, join = st_within, left = FALSE)

# ---------------------------
# Calculate country-level species richness
# The spatial join adds the `admin` column (country name) from the world polygons
country_richness <- occurrences_land %>%
  st_drop_geometry() %>%
  filter(!is.na(admin)) %>%
  group_by(admin) %>%
  summarise(richness = n_distinct(wcfp_name_match), .groups = "drop")

# ---------------------------
# Join richness data with the world map for plotting
map_data <- world_4326 %>%
  left_join(country_richness, by = "admin")

# ---------------------------
# Check for unmatched countries
unmatched <- country_richness %>%
  filter(!admin %in% world_4326$admin)

if (nrow(unmatched) > 0) {
  message("WARNING: ", nrow(unmatched), " countries in your data did not match the world map:")
  print(unmatched)
} else {
  message("All countries matched successfully.")
}

# ---------------------------
# Summary stats
message("Max species richness per country: ", max(map_data$richness, na.rm = TRUE))
message("Total distinct species: ", n_distinct(occurrences_land$wcfp_name_match))
message("Total countries with data: ", nrow(country_richness))




# --------------------------------------------------------- #
#     Filter method 2. GRIN Taxonomy distributions data     #
# --------------------------------------------------------- #

# 2. Filter using GRIN Taxonomy distributions data - use query
# How many species have this data?

# Web query output. Replaces the earlier two-part CSV export
# (WebQuery_GRINGlobal_countries_wcfp_p1/p2.csv), which turned out to be a
# NATIVE-ONLY extract: its 199,390 rows are exactly the geography_status_code
# == "n" subset of this file. That is why the old status column looked
# uninformative - every non-native record had been filtered out upstream.
#
# This export adds 23,056 non-native rows and 116 species that GRIN records
# only as introduced/cultivated - including staples such as Allium cepa,
# Allium porrum and Abelmoschus esculentus, which previously had no GRIN
# evidence at all. It also adds country_name, continent and subcontinent
# columns (not used below, but available).
GRIN_countries <- read_excel(
  file.path(inputs_dir, "WebQueryGRIN-Global_distributions_combined.xlsx"),
  guess_max = 100000
) %>%
  # Status collapsed to a two-level field: GRIN's "n" is native, every other
  # code (i introduced, c cultivated, a, u, o, and 5 blank rows) is
  # non-native. Carried through the pipeline below and reported in the
  # popups, workbooks and combined tables alongside the WCFP status.
  mutate(grin_status = if_else(geography_status_code == "n", "native", "non-native"))

length(unique(GRIN_countries$name))
length(unique(GRIN_countries$country_code))
cat("GRIN records:", nrow(GRIN_countries), "-",
    sum(GRIN_countries$grin_status == "native"), "native,",
    sum(GRIN_countries$grin_status == "non-native"), "non-native\n")



# ------------------------------------------------------------ #
#     Filter method 3. Filter using WCVP                       #
#                     World checklist of Food Vascular Plants  #
# ------------------------------------------------------------ #


# 3. Filter using World Checklist of Vascular Plants (WCVP) / Taxonomic Databases Working Group -Biodiversity Information Standards  (TDWG) area (level 3, countries or subnational units). - (as we did in WCFP paper)

# How many species have this data?

# Reproduces the exact join used in the WCFP paper's own code
# (data_summaries.R, section 8: "NATIVE / INTRODUCED FLAG") - join your
# species list's accepted-name IPNI LSID to WCVP's plant_name_id, then to
# WCVP's distribution table - rather than scraping POWO species-by-species
# via expowo::powoSpDist (which hit persistent connection errors and
# deterministic genus mismatches against Kew's live server).
#
# https://powo.science.kew.org/about-wcvp#geographicaldistribution
# https://www.nature.com/articles/s41597-021-00997-6
# https://powo.science.kew.org/
# https://matildabrown.github.io/rWCVP/
# WCFP paper data-and-code archive: https://doi.org/10.6084/m9.figshare.31441615.v2
#
# Citation for the data-and-code archive:
# Diazgranados, Mauricio; Gianella, Maraeva; Kor, Laura; Gori, Benedetta;
# Khoury, Colin; Gora, Sarah L.; et al. (2026). World Checklist of Food
# Plants 2026 (WCFP): Data. figshare. Dataset.
# https://doi.org/10.6084/m9.figshare.31441615.v2

# load wcfp
wcfp_plantlist <- read_excel(
  file.path(inputs_dir, "WCFP_plantlist_final_2026-06-02.xlsx"))


## ------------------------------------------------------------------ ##
## Normalise a user-supplied taxon name to WCFP's form.
##
## wcfp_name_match is the accepted binomial with NO authors, and is identical
## to taxon_name_accepted for all 26,632 rows. Shapes present: 26,428 two-word
## names and 204 three-word ones - infix hybrids ("Achillea x serrata"),
## leading hybrids ("x Pyraria irregularis"), leading "+" graft chimaeras
## ("+ Pyrocydonia danielii") and "sect." names. All four must survive
## normalisation intact, which is why authors are stripped by TOKEN COUNT
## rather than by pattern - a pattern that strips capitalised trailing words
## would eat "Taraxacum sect. Taraxacum".
## ------------------------------------------------------------------ ##
normalise_taxon_name <- function(x) {
  x <- str_squish(as.character(x))
  x <- gsub("×", "x", x, fixed = TRUE)   # unify the hybrid sign to plain x
  x <- gsub("’|‘|`", "'", x)         # curly quotes -> straight
  x <- sub("^([A-Za-z])", "\\U\\1", tolower(x), perl = TRUE)  # sentence case
  x
}

# How many leading tokens are the NAME, the rest being authors.
.name_token_count <- function(x) {
  tok <- strsplit(x, " ", fixed = TRUE)
  vapply(tok, function(t) {
    if (length(t) >= 2 && t[1] %in% c("x", "+")) return(3L)   # x Genus species
    if (length(t) >= 3 && t[2] %in% c("x", "+", "sect.")) return(3L)
    2L
  }, integer(1))
}

## ------------------------------------------------------------------ ##
## Subset the WCFP plantlist to a user-supplied list of taxon names.
##
## Matching runs in tiers, most exact first, and every input row is reported
## with the tier that matched it so the result is auditable rather than a bare
## count. Unmatched names are written out - silently dropping them is the
## failure mode this whole function exists to prevent.
## ------------------------------------------------------------------ ##
subset_plantlist_to_user_taxa <- function(plantlist, path, column = NULL,
                                          report_dir = NULL) {
  if (!file.exists(path)) stop("USER_TAXA_FILE does not exist: ", path)

  ext <- tolower(tools::file_ext(path))
  user <- switch(ext,
    xlsx = ,
    xls  = readxl::read_excel(path, guess_max = 50000),
    csv  = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    tsv  = ,
    txt  = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    stop("USER_TAXA_FILE must be .xlsx, .xls, .csv, .tsv or .txt - got '", ext, "'"))
  if (nrow(user) == 0) stop("USER_TAXA_FILE has no rows: ", path)

  # Pick the name column: the one asked for, else the first whose name looks
  # taxonomic, else the first column that actually contains binomial-shaped text.
  if (is.null(column)) {
    cand <- grep("^(taxon|taxa|name|species|scientific|binomial)",
                 names(user), ignore.case = TRUE, value = TRUE)
    if (!length(cand)) {
      looks <- vapply(user, function(v) {
        v <- str_squish(as.character(v)); v <- v[!is.na(v) & nzchar(v)]
        if (!length(v)) return(0)
        mean(grepl("^[A-Za-z+x×'.-]+ [A-Za-z'.-]+", v))
      }, numeric(1))
      if (max(looks) < 0.5) {
        stop("could not find a taxon-name column in ", basename(path),
             " - set USER_TAXA_COLUMN. Columns present: ",
             paste(names(user), collapse = ", "))
      }
      cand <- names(user)[which.max(looks)]
    }
    column <- cand[1]
  }
  if (!column %in% names(user)) {
    stop("USER_TAXA_COLUMN '", column, "' is not in ", basename(path),
         ". Columns present: ", paste(names(user), collapse = ", "))
  }
  message("Reading user taxa from '", basename(path), "', column '", column, "'")

  raw  <- as.character(user[[column]])
  keep <- !is.na(raw) & nzchar(str_squish(raw))
  raw  <- raw[keep]
  norm <- normalise_taxon_name(raw)

  # Reference forms, normalised the same way.
  ref_match <- normalise_taxon_name(plantlist$wcfp_name_match)
  ref_gs    <- normalise_taxon_name(plantlist$genus_species)

  # Tier 1: the name as given.  Tier 2: against genus_species.
  # Tier 3: authors stripped by token count.
  hit  <- match(norm, ref_match)
  tier <- ifelse(is.na(hit), NA_character_, "exact name")

  todo <- is.na(hit)
  if (any(todo)) {
    h2 <- match(norm[todo], ref_gs)
    hit[todo]  <- h2
    tier[todo] <- ifelse(is.na(h2), NA_character_, "genus_species")
  }
  todo <- is.na(hit)
  if (any(todo)) {
    stripped <- vapply(seq_along(norm[todo]), function(i) {
      t <- strsplit(norm[todo][i], " ", fixed = TRUE)[[1]]
      paste(t[seq_len(min(length(t), .name_token_count(norm[todo][i])))], collapse = " ")
    }, character(1))
    h3 <- match(stripped, ref_match)
    hit[todo]  <- h3
    tier[todo] <- ifelse(is.na(h3), NA_character_, "authors stripped")
  }

  audit <- tibble::tibble(
    supplied   = raw,
    normalised = norm,
    matched_to = ifelse(is.na(hit), NA_character_, plantlist$wcfp_name_match[hit]),
    matched_by = tier
  )
  n_in   <- nrow(audit)
  n_hit  <- sum(!is.na(hit))
  n_taxa <- length(unique(hit[!is.na(hit)]))

  cat("\n--- user taxa list ---------------------------------------------\n")
  cat("  names supplied        :", n_in, "\n")
  cat("  matched to WCFP       :", n_hit, sprintf("(%.1f%%)", 100 * n_hit / n_in), "\n")
  cat("  distinct WCFP taxa    :", n_taxa, "\n")
  if (n_hit) print(table(audit$matched_by[!is.na(audit$matched_by)]))

  if (!is.null(report_dir)) {
    if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE)
    f <- file.path(report_dir, "user_taxa_match_report.xlsx")
    openxlsx::write.xlsx(list(all = audit, unmatched = audit[is.na(hit), ]), f,
                         overwrite = TRUE)
    cat("  match report          :", f, "\n")
  }

  if (any(is.na(hit))) {
    miss <- audit$supplied[is.na(hit)]
    cat("  NOT MATCHED           :", length(miss),
        "- these contribute NOTHING to any output\n")
    cat("    ", paste(utils::head(miss, 8), collapse = " | "),
        if (length(miss) > 8) paste0(" ... +", length(miss) - 8, " more") else "", "\n")
  }
  cat("----------------------------------------------------------------\n\n")

  if (n_taxa == 0) {
    stop("none of the supplied names matched the WCFP list - nothing to map. ",
         "See the match report for what was read.")
  }
  plantlist[sort(unique(hit[!is.na(hit)])), , drop = FALSE]
}

# Output folder: a user-list run gets its own, so a full WCFP run is never
# overwritten and the two can be compared side by side.
if (is.null(USER_TAXA_FILE)) {
  outputs_dir <- outputs_root
} else {
  outputs_dir <- file.path(outputs_root,
                           paste0("user_", tools::file_path_sans_ext(basename(USER_TAXA_FILE))))
  if (!dir.exists(outputs_dir)) dir.create(outputs_dir, recursive = TRUE)
  message("User-list run. All outputs go to: ", outputs_dir)
  wcfp_plantlist <- subset_plantlist_to_user_taxa(
    wcfp_plantlist, USER_TAXA_FILE, USER_TAXA_COLUMN, report_dir = outputs_dir)
}

wcvp_data_dir <- inputs_dir

# ----------------------------------------------------------------- #
# 3a) Load the WCVP snapshot shipped with the WCFP paper - the same
#     files used to produce the published TDWG level 3 distribution
#     counts (Govaerts et al. 2021, Scientific Data 8, 215). Shipped
#     rather than freshly downloaded because WCVP is revised
#     continuously and a later download wouldn't reproduce the same
#     counts.
# ----------------------------------------------------------------- #
wcvp_names <- read.table(file.path(wcvp_data_dir, "wcvp_names.csv"), sep = "|",
                          header = TRUE, quote = "", fill = TRUE, encoding = "UTF-8")
wcvp_distribution <- read.table(file.path(wcvp_data_dir, "wcvp_distribution.csv"), sep = "|",
                                 header = TRUE, quote = "", fill = TRUE, encoding = "UTF-8")

# A handful of records carry a lower-case TDWG code, which fails to join
# to any TDWG L3 lookup/shapefile downstream. Normalise, same as the
# paper's own script.
wcvp_distribution$area_code_l3 <- toupper(trimws(wcvp_distribution$area_code_l3))

# Accepted WCVP names with a usable (non-blank) IPNI identifier - blank
# ids would otherwise all match each other and any of your species that
# also has no LSID.
wcvp_ids <- wcvp_names %>%
  filter(taxon_status == "Accepted") %>%
  mutate(ipni_id = na_if(trimws(as.character(ipni_id)), "")) %>%
  filter(!is.na(ipni_id))

# ----------------------------------------------------------------- #
# 3b) Join wcfp_plantlist -> WCVP plant_name_id (by IPNI LSID) -> WCVP
#     distribution table (one taxon can have many TDWG L3 areas, hence
#     the one-to-many join).
# ----------------------------------------------------------------- #
wcfp_plantlist_ids <- wcfp_plantlist %>%
  mutate(LSID_accepted = na_if(trimws(as.character(LSID_accepted)), "")) %>%
  filter(!is.na(LSID_accepted)) %>%
  left_join(wcvp_ids %>% select(ipni_id, plant_name_id),
            by = c("LSID_accepted" = "ipni_id"), na_matches = "never") %>%
  filter(!is.na(plant_name_id))

cat("wcfp_plantlist species with a usable IPNI LSID matched to a WCVP plant_name_id:",
    n_distinct(wcfp_plantlist_ids$WCFP_ID), "of", n_distinct(wcfp_plantlist$WCFP_ID), "\n")

# A handful of WCFP_ID rows share the same plant_name_id (e.g. duplicate or
# synonymous entries in wcfp_plantlist pointing to the same accepted WCVP
# taxon) - not a strict one-row-per-taxon list. Flag them for review, and
# join as many-to-many below so every WCFP_ID still gets its distribution
# data rather than the join erroring out.
duplicate_plant_name_ids <- wcfp_plantlist_ids %>%
  count(plant_name_id, name = "n_wcfp_ids") %>%
  filter(n_wcfp_ids > 1)

if (nrow(duplicate_plant_name_ids) > 0) {
  cat(nrow(duplicate_plant_name_ids), "WCVP plant_name_id(s) are shared by more than one WCFP_ID",
      "(duplicate/synonymous entries in wcfp_plantlist) - see `duplicate_plant_name_ids`",
      "and `wcfp_plantlist_dupes` for the affected rows.\n")
}

wcfp_plantlist_dupes <- wcfp_plantlist_ids %>%
  semi_join(duplicate_plant_name_ids, by = "plant_name_id") %>%
  select(WCFP_ID, LSID_accepted, taxon_name_accepted, plant_name_id) %>%
  arrange(plant_name_id)

wcvp_dist_tdwg3 <- wcfp_plantlist_ids %>%
  select(WCFP_ID, LSID_accepted, taxon_name_accepted, family, plant_name_id) %>%
  inner_join(wcvp_distribution, by = "plant_name_id", relationship = "many-to-many") %>%
  filter(area_code_l3 != "") %>%
  # Collapsed to the same two-level native / non-native scheme used for GRIN
  # above, so the two sources' status fields are directly comparable. WCVP's
  # native records stay "native"; introduced, extinct and location-doubtful
  # all become "non-native". The raw introduced / extinct / location_doubtful
  # flags are still carried in the columns selected below, so the original
  # four-way distinction is not lost - only the summary field changes.
  #
  # case_when (not if_else) to preserve the original NA handling: a record
  # with NA flags falls through to "native" rather than becoming NA.
  mutate(occurrence_status = case_when(
    location_doubtful == 1 ~ "non-native",
    extinct           == 1 ~ "non-native",
    introduced        == 1 ~ "non-native",
    TRUE                   ~ "native"
  )) %>%
  select(WCFP_ID, LSID_accepted, taxon_name_accepted, family,
         continent_code_l1, continent, region_code_l2, region,
         area_code_l3, area, occurrence_status,
         introduced, extinct, location_doubtful) %>%
  distinct()

cat("wcfp_plantlist species with at least one TDWG level 3 distribution record:",
    n_distinct(wcvp_dist_tdwg3$WCFP_ID), "of", n_distinct(wcfp_plantlist$WCFP_ID), "\n")

# ----------------------------------------------------------------- #
# 3c) TDWG level 3 presence, as a long-format table (one row per
#     species x area, any occurrence_status), plus a per-area count of
#     distinct species - mirrors the paper's own "food-plant richness
#     per TDWG3" step (data_summaries.R section 9a), which also counts
#     all statuses rather than native-only.
# ----------------------------------------------------------------- #
wcvp_dist_tdwg3_long <- wcvp_dist_tdwg3 %>%
  distinct(WCFP_ID, taxon_name_accepted, area_code_l3, area)

wcvp_distinct_taxa_by_tdwg3 <- wcvp_dist_tdwg3_long %>%
  distinct(area_code_l3, area, WCFP_ID) %>%
  count(area_code_l3, area, sort = TRUE, name = "n_distinct_taxa")

wcvp_distinct_taxa_by_tdwg3

cat("wcfp_plantlist species with any WCVP TDWG level 3 distribution:",
    n_distinct(wcvp_dist_tdwg3_long$WCFP_ID), "of", n_distinct(wcfp_plantlist$WCFP_ID), "\n")




## ============================================================================ ##
## Combined method: distinct taxa per country, by data source
## Author: Sarah Gora
##
## For each country, lists every taxon confirmed present by >=2 of the 3
## methods above (occurrences, GRIN, POWO/WCVP), with a Y/NA flag showing
## which source(s) found it there. One sheet per qualifying country in a
## single .xlsx workbook. Runs right after Method 3 (rather than at the
## end of the script) because the country-level richness map below is
## driven by this section's `combined_confirmed` result, not raw
## single-source occurrence data.
##
## Source 3 (POWO/WCVP) below is adapted from the original expowo-scraped
## version: it now reads from `wcvp_dist_tdwg3` (the local WCVP-snapshot
## join built in Method 3 above) instead of expowo's scraped
## `native_to_country` field. The two are the same underlying data - WCVP's
## TDWG level 3 area names, which POWO's own species pages label loosely as
## "country" - so the country-name alias/expansion tables below are reused
## unchanged. The join back to a taxon name now uses `WCFP_ID` directly
## (carried through `wcvp_dist_tdwg3` from wcfp_plantlist) rather than
## reconstructing a sentence-cased binomial to match against, which is more
## robust for hybrids/infraspecifics where `taxon_name_accepted` differs
## from the plain binomial.
## ============================================================================ ##

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

combined_output_dir <- outputs_dir   # set near the top; user runs get their own folder
if (!dir.exists(combined_output_dir)) dir.create(combined_output_dir, recursive = TRUE)

## Self-contained: recomputes its own country polygons and occurrence-
## country join rather than reusing `world_4326`/`occurrences_land` from
## Method 1 above, so this section doesn't depend on execution order.
country_polys <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") %>%
  fold_breakaway_polygons() %>%
  st_transform(4326)

# ---------------------------------------- #
#   Source 1: occurrences                  #
# ---------------------------------------- #
occ_country_sf <- st_as_sf(
  readRDS(file.path(inputs_dir, "occurrences_land.rds")),
  coords = c("longitude", "latitude"), crs = 4326, remove = FALSE
)
occ_country_joined <- st_join(occ_country_sf, country_polys, join = st_within, left = FALSE)

# A species can have occurrence records from more than one raw data
# source (e.g. GBIF, herbarium) in a given country - aggregated here as a
# comma-separated list of distinct values per (country, taxon), same
# pattern as status_WCFP below.
src_occurrences <- occ_country_joined %>%
  st_drop_geometry() %>%
  filter(!is.na(admin), !is.na(wcfp_name_match)) %>%
  distinct(country = admin, taxon = wcfp_name_match, data_source) %>%
  group_by(country, taxon) %>%
  summarise(occurrences_data_source = paste(sort(unique(na.omit(data_source))), collapse = ", "),
            .groups = "drop")

# ---------------------------------------- #
#   Taxon name lookup: genus_species -> wcfp_name_match
# ---------------------------------------- #
# Used below by GRIN (source 2) as a fallback when GRIN's species name
# doesn't match `wcfp_name_match` directly -- most commonly hybrids, where
# `wcfp_name_match` carries the "x" hybrid symbol (e.g. "Amelanchier x
# spicata") but GRIN gives the plain binomial ("Amelanchier spicata");
# `genus_species` is that same plain binomial (lowercase, no hybrid
# symbol), so sentence-casing it recovers the match.
taxon_name_lookup <- wcfp_plantlist %>%
  transmute(
    taxon_name      = str_to_sentence(str_squish(as.character(genus_species))),
    wcfp_name_match = str_squish(wcfp_name_match)
  ) %>%
  filter(!is.na(taxon_name), !is.na(wcfp_name_match)) %>%
  distinct(taxon_name, .keep_all = TRUE)

# ---------------------------------------- #
#   Source 2: GRIN Taxonomy distributions  #
# ---------------------------------------- #
# Only rows with a real ISO3 country code -- GRIN also has a small share of
# rows coded with non-ISO regional/former-USSR-subdivision codes (e.g.
# "002" = Africa, "155" = "European part" of Russia) with no reliable
# crosswalk to specific countries, so those are dropped rather than guessed.
grin_alpha3 <- GRIN_countries %>%
  filter(str_detect(country_code, "^[A-Z]{3}$"))

iso3_to_country <- country_polys %>%
  st_drop_geometry() %>%
  transmute(
    iso3 = if_else(is.na(iso_a3) | iso_a3 == "-99", iso_a3_eh, iso_a3),
    admin
  ) %>%
  filter(!is.na(iso3), nzchar(iso3)) %>%
  distinct(iso3, .keep_all = TRUE)

# A handful of GRIN ISO3 codes still don't match even with the iso_a3_eh
# fallback above, because the underlying territory isn't a separate
# Natural Earth admin-0 polygon at all (it's folded into a parent
# country's polygon under a different code) or uses an obsolete ISO3 code
# no longer present in either iso_a3 field. Mapped directly to the parent
# country's `admin` name here, bypassing the ISO join entirely for these:
#  - GLP/GUF/MTQ/REU/XYO (Guadeloupe, French Guiana, Martinique, Réunion,
#    Mayotte): French overseas departments. These are legally integral parts
#    of France rather than dependencies, so Natural Earth draws no separate
#    admin-0 polygon for any of them at EITHER 1:50m or 1:10m - verified.
#    They are therefore folded into "France", which is the only polygon that
#    can carry their records.
#  - SJM (Svalbard and Jan Mayen): folded into "Norway"
#  - CXR (Christmas Island): Natural Earth's separate "Indian Ocean
#    Territories" polygon, not under Australia's own code
#  - TKL (Tokelau): folded into "New Zealand"
#  - TMP: an obsolete ISO3 code for East Timor, superseded by TLS, so it
#    never matches either current iso_a3 field. NOTE the target must be
#    "East Timor", which is Natural Earth's `admin` value - "Timor-Leste"
#    is its `name` value and matches no polygon. The alias previously said
#    "Timor-Leste", which passed the !is.na(admin) filter below and so put
#    14 GRIN rows into src_grin under a country that does not exist on the
#    map: invisible in every map, and worse, filed separately from the
#    occurrence/WCFP evidence for the same place, so those taxa could fail
#    to reach the >=2-source threshold.
#  - XYO: identified as MAYOTTE by GRIN's own GEOGRAPHY_COUNTRY_CODE lookup
#    (WebQueryGRIN-Global_geo_codes_2026-09-09.xlsx). Previously dropped as
#    "not a genuine ISO 3166-1 alpha-3 code". It is a French overseas
#    department and is not a separate Natural Earth admin-0 polygon, so it
#    folds into "France" exactly like GLP/GUF/MTQ/REU above.
#  - GIB (Gibraltar): a British Overseas Territory not present as its own
#    Natural Earth admin-0 polygon at medium scale, folded into "United
#    Kingdom"
#  - ANT: the retired code for the Netherlands Antilles, superseded by CUW
#    (Curaçao). 147 rows, previously dropped entirely. Note that GRIN's
#    `state` field names a more specific island on 71 of those rows (Saba 28,
#    Curacao 17, St. Eustatius 15, Bonaire 6, South St. Martin 5); mapping the
#    code as a whole to Curaçao therefore also brings in rows that GRIN
#    attributes to Saba, St. Eustatius, Bonaire and Sint Maarten. Applied per
#    explicit user instruction.
grin_country_code_alias <- c(
  GLP = "France",
  GUF = "France",
  MTQ = "France",
  REU = "France",
  XYO = "France",
  SJM = "Norway",
  CXR = "Indian Ocean Territories",
  TKL = "New Zealand",
  TMP = "East Timor",
  GIB = "United Kingdom",
  ANT = "Curaçao"
)

# Guard: every alias target must be a real Natural Earth admin name, or the
# rows silently become a phantom country (as TMP did). Fail loudly instead.
.bad_alias <- setdiff(unname(grin_country_code_alias), country_polys$admin)
if (length(.bad_alias) > 0) {
  stop("grin_country_code_alias targets not present in country_polys$admin: ",
       paste(.bad_alias, collapse = ", "))
}

# ---------------------------------------- #
#   State promotion                        #
# ---------------------------------------- #
# GRIN files some territories as a `state` of their parent country, e.g.
# Puerto Rico is country_code = "USA", state = "Puerto Rico". The
# country-level pipeline keys on country_code alone, so those records were
# being counted as the parent - Puerto Rico's 496 species counted as United
# States - while Natural Earth's separate Puerto Rico polygon received
# nothing and rendered as a white "no data" hole.
#
# These are the only (country_code, state) pairs where GRIN names a state
# AND Natural Earth draws that place as its own admin-0 polygon. For them,
# the state wins over the country code.
#
# KEYED ON THE PAIR, DELIBERATELY. Two GRIN states share a name with a
# country: USA/"Georgia" (729 rows) and NGA/"Niger" (1 row, Niger State,
# Nigeria's largest state). Matching on the state name alone would move all
# 729 American records to the country of Georgia. Both are correctly coded
# as-is and must NOT be promoted.
grin_state_promotion <- tribble(
  ~country_code, ~state,                        ~promoted_admin,
  "USA",         "Puerto Rico",                 "Puerto Rico",
  "USA",         "Virgin Islands, U.S.",        "United States Virgin Islands",
  "USA",         "Guam",                        "Guam",
  "USA",         "Northern Mariana Islands",    "Northern Mariana Islands",
  "USA",         "American Samoa",              "American Samoa",
  "AUS",         "Ashmore and Cartier Islands", "Ashmore and Cartier Islands",
  # Redundant while ANT aliases to Curaçao wholesale, but kept explicit so
  # these rows still land correctly if that alias is ever changed.
  "ANT",         "Curacao",                     "Curaçao"
)

.bad_promo <- setdiff(grin_state_promotion$promoted_admin, country_polys$admin)
if (length(.bad_promo) > 0) {
  stop("grin_state_promotion targets not present in country_polys$admin: ",
       paste(.bad_promo, collapse = ", "))
}

# ---------------------------------------- #
#   Historical / regional numeric codes    #
# ---------------------------------------- #
# GRIN's numeric codes are its legacy geography namespace - historical or
# aggregate entities (Soviet Union, Rhodesia, Ceylon, Former Yugoslavia).
# They are not ISO3 and were previously dropped SILENTLY, because the
# `^[A-Z]{3}$` filter removes them before the alias check ever runs.
#
# Each is expanded to its modern constituent countries: a taxon recorded in
# "Former Yugoslavia" is counted present in ALL seven successor states,
# mirroring how manual_country_overrides already treats multi-country TDWG
# areas. Targets are Natural Earth `admin` names, which differ from the
# obvious country name in two places here - Serbia is "Republic of Serbia",
# and Congo-Brazzaville is "Republic of the Congo".
#
# Membership lists supplied by the user. Note SLO (IOC/FIFA) and XKX
# (user-assigned) are not ISO 3166-1 codes; they resolve to the Natural
# Earth admin names "Slovenia" and "Kosovo" respectively.
grin_region_members <- list(
  "002" = c("Saudi Arabia", "Yemen", "Oman", "United Arab Emirates", "Qatar",
            "Bahrain", "Kuwait"),                                        # Arabia
  "004" = c("Turkey"),                                                   # Asia Minor
  "012" = c("Chad", "Central African Republic", "Republic of the Congo",
            "Gabon"),                                                    # French Equatorial Africa
  "015" = c("North Korea", "South Korea"),                               # Korea
  "018" = c("Bahrain", "Cyprus", "Egypt", "Iran", "Iraq", "Israel", "Jordan",
            "Kuwait", "Lebanon", "Oman", "Qatar", "Saudi Arabia", "Syria",
            "Turkey", "United Arab Emirates", "Yemen"),                  # Middle East
  "019" = c("Papua New Guinea", "Indonesia"),                            # New Guinea
  "020" = c("Algeria", "Egypt", "Morocco", "Sudan", "Tunisia",
            "Western Sahara"),                                           # North Africa
  "023" = c("Israel", "Palestine", "Jordan"),                            # Ancient Palestine
  "024" = c("Zimbabwe"),                                                 # Rhodesia
  "026" = c("Somalia"),                                                  # Somaliland
  "031" = c("Kazakhstan", "Uzbekistan", "Turkmenistan", "Kyrgyzstan",
            "Tajikistan", "Afghanistan", "China"),                       # Turkistan
  "036" = c("Armenia", "Azerbaijan", "Bahrain", "Cyprus", "Georgia", "Iraq",
            "Israel", "Jordan", "Kuwait", "Lebanon", "Oman", "Qatar",
            "Saudi Arabia", "Palestine", "Syria", "Turkey",
            "United Arab Emirates", "Yemen"),                            # Western Asia
  "108" = c("Armenia", "Azerbaijan", "Belarus", "Estonia", "Georgia",
            "Kazakhstan", "Kyrgyzstan", "Latvia", "Lithuania", "Moldova",
            "Russia", "Tajikistan", "Turkmenistan", "Ukraine",
            "Uzbekistan"),                                               # Soviet Union
  "155" = c("Russia"),                                                   # RF European part
  "156" = c("Russia"),                                                   # RF Western Siberia
  "157" = c("Russia"),                                                   # RF Eastern Siberia
  "158" = c("Russia"),                                                   # RF Far East
  "159" = c("Russia"),                                                   # RF Ciscaucasia
  "161" = c("Slovenia", "Croatia", "Bosnia and Herzegovina",
            "Republic of Serbia", "Montenegro", "North Macedonia",
            "Kosovo"),                                                   # Former Yugoslavia
  # CZE / SVK. Natural Earth's admin for CZE is "Czechia", not
  # "Czech Republic" - the latter matches no polygon.
  "163" = c("Czechia", "Slovakia"),                                      # Czechoslovakia
  "013" = c("Haiti", "Dominican Republic"),                              # Hispaniola
  "005" = c("Indonesia", "Malaysia", "Brunei")                           # Borneo
)

# ---------------------------------------- #
#   Continent-level codes                  #
# ---------------------------------------- #
# 001 Africa, 003 Asia, 009 Europe, 021 North America, 027 South America.
# Expanded per explicit user instruction: a continent record cannot confirm
# anything on its own, because the combined method still requires >=2 of the
# 3 sources - so it acts only as corroboration where occurrences or WCFP
# independently place the taxon in that country.
#
# Membership comes from NATURAL EARTH's continent/subregion fields, NOT from
# GRIN's own continent column. GRIN assigns continent per RECORD, based on
# the territory involved, so a country appears in whichever continents its
# territories fall in: Madeira makes Portugal "African", the Canaries make
# Spain "African", Puerto Rico makes the USA "South American", Sinai makes
# Egypt "Asian". Using that directly would count an "Africa" taxon as
# present in Spain and Portugal. Natural Earth's continent is a property of
# the country, so it has no such ambiguity.
#
# Derived from country_polys (i.e. AFTER fold_breakaway_polygons), so
# Somaliland is not referenced - it no longer exists as its own admin.
#
# Two deliberate departures from Natural Earth, to match GRIN's TDWG-style
# split of the Americas:
#   - Mexico is placed in Northern America (TDWG region 79, and GRIN's own
#     data files it there); Natural Earth puts it in Central America.
#   - Southern America absorbs Central America and the Caribbean, as TDWG
#     region 8 does.
.ne_meta <- country_polys %>% st_drop_geometry() %>% select(admin, continent, subregion)

grin_continent_members <- list(
  "001" = .ne_meta$admin[.ne_meta$continent == "Africa"],
  "009" = .ne_meta$admin[.ne_meta$continent == "Europe"],
  "003" = .ne_meta$admin[.ne_meta$continent == "Asia"],
  "021" = union(.ne_meta$admin[.ne_meta$subregion == "Northern America"], "Mexico"),
  "027" = setdiff(.ne_meta$admin[.ne_meta$continent == "South America" |
                                 .ne_meta$subregion %in% c("Central America", "Caribbean")],
                  "Mexico")
)

grin_region_members <- c(grin_region_members, grin_continent_members)

grin_region_expansion <- tibble(
  country_code = rep(names(grin_region_members), lengths(grin_region_members)),
  admin        = unlist(grin_region_members, use.names = FALSE)
)

# Same guard as the alias table: a target that is not a real Natural Earth
# admin would silently create a phantom country.
.bad_region <- setdiff(unique(grin_region_expansion$admin), country_polys$admin)
if (length(.bad_region) > 0) {
  stop("grin_region_members targets not present in country_polys$admin: ",
       paste(.bad_region, collapse = ", "))
}

grin_matched <- grin_alpha3 %>%
  left_join(iso3_to_country, by = c("country_code" = "iso3")) %>%
  mutate(admin = if_else(is.na(admin) & country_code %in% names(grin_country_code_alias),
                          grin_country_code_alias[country_code], admin)) %>%
  # str_squish on both sides so stray whitespace in `state` cannot cause a
  # silent miss.
  mutate(.state_norm = str_squish(as.character(state))) %>%
  left_join(grin_state_promotion %>% mutate(.state_norm = str_squish(state)) %>%
              select(country_code, .state_norm, promoted_admin),
            by = c("country_code", ".state_norm")) %>%
  mutate(admin = coalesce(promoted_admin, admin)) %>%
  select(-.state_norm, -promoted_admin)

grin_dropped_codes <- setdiff(unique(grin_alpha3$country_code),
                               c(iso3_to_country$iso3, names(grin_country_code_alias)))
if (length(grin_dropped_codes) > 0) {
  message(length(grin_dropped_codes), " GRIN ISO3 country codes have no matching Natural Earth ",
          "country and no manual alias, so are dropped: ",
          paste(sort(grin_dropped_codes), collapse = ", "))
}

# Expand the numeric/region rows: one row per constituent country, so a
# taxon recorded in e.g. "Former Yugoslavia" is counted present in all seven
# successor states. inner_join is intentionally many-to-many - that IS the
# expansion. Region codes with no mapping (the continent-level ones, plus
# Hispaniola, West Indies etc.) simply do not join and stay dropped.
grin_regions_expanded <- GRIN_countries %>%
  inner_join(grin_region_expansion, by = "country_code",
             relationship = "many-to-many")

.numeric_codes  <- unique(GRIN_countries$country_code[!grepl("^[A-Z]{3}$", GRIN_countries$country_code)])
.mapped_codes   <- intersect(.numeric_codes, names(grin_region_members))
.unmapped_codes <- setdiff(.numeric_codes, names(grin_region_members))
message(length(.mapped_codes), " of ", length(.numeric_codes),
        " GRIN numeric/region codes expanded to constituent countries (",
        nrow(grin_regions_expanded), " rows from ",
        sum(GRIN_countries$country_code %in% .mapped_codes), " source rows).")
if (length(.unmapped_codes) > 0) {
  message(length(.unmapped_codes), " GRIN numeric/region codes have no country mapping and are ",
          "dropped (", sum(GRIN_countries$country_code %in% .unmapped_codes), " rows): ",
          paste(sort(.unmapped_codes), collapse = ", "))
}

grin_with_country <- bind_rows(
  grin_matched %>% filter(!is.na(admin)),
  grin_regions_expanded
)

# GRIN's geography_status_code is now carried through as `grin_status`
# (native / non-native, mapped at load time). Under the previous native-only
# CSV export the field was constant "n" and therefore useless; the current
# xlsx export contains the introduced/cultivated records too, so it is
# genuinely informative and is reported alongside the WCFP status.

# 1) Direct match against wcfp_name_match
grin_direct <- grin_with_country %>%
  filter(name %in% wcfp_plantlist$wcfp_name_match) %>%
  transmute(country = admin, taxon = str_squish(name), grin_status)

# 2) Fallback: match the remainder against genus_species (catches hybrids,
#    where wcfp_name_match carries the "x" symbol but GRIN/genus_species
#    don't)
grin_via_genus_species <- grin_with_country %>%
  filter(!(name %in% wcfp_plantlist$wcfp_name_match)) %>%
  inner_join(taxon_name_lookup, by = c("name" = "taxon_name")) %>%
  transmute(country = admin, taxon = wcfp_name_match, grin_status)

# 3) Fallback: manual alias for GRIN names missing a hyphen that WCFP's
#    specific epithet carries (e.g. GRIN "belladonna" vs WCFP
#    "bella-donna"), verified individually against wcfp_plantlist.
grin_hyphen_alias <- c(
  "Atropa belladonna"   = "Atropa bella-donna",
  "Lanxangia tsaoko"    = "Lanxangia tsao-ko",
  "Mezilaurus itauba"   = "Mezilaurus ita-uba",
  "Oxalis sanmiguelii"  = "Oxalis san-miguelii"
)

grin_via_hyphen_alias <- grin_with_country %>%
  filter(!(name %in% wcfp_plantlist$wcfp_name_match), name %in% names(grin_hyphen_alias)) %>%
  transmute(country = admin, taxon = recode(name, !!!grin_hyphen_alias), grin_status) %>%
  filter(taxon %in% wcfp_plantlist$wcfp_name_match)

grin_still_unmatched <- setdiff(
  unique(grin_with_country$name[!(grin_with_country$name %in% wcfp_plantlist$wcfp_name_match)]),
  c(taxon_name_lookup$taxon_name, names(grin_hyphen_alias))
)
if (length(grin_still_unmatched) > 0) {
  message(length(grin_still_unmatched), " GRIN species names match neither wcfp_name_match, ",
          "genus_species, nor a known hyphenation alias, and are dropped from the GRIN source: ",
          paste(sort(grin_still_unmatched), collapse = ", "))
}

# A taxon can be native in one GRIN record for a country and introduced in
# another (different state/subdivision) - aggregated as a comma-separated
# list of distinct statuses per (country, taxon), e.g. "native, non-native",
# the same pattern src_wcfp uses below.
src_grin <- bind_rows(grin_direct, grin_via_genus_species, grin_via_hyphen_alias) %>%
  distinct(country, taxon, grin_status) %>%
  group_by(country, taxon) %>%
  summarise(grin_status = paste(sort(unique(grin_status)), collapse = ", "), .groups = "drop")

# ---------------------------------------- #
#   Source 3: WCFP plantlist                #
#   (via WCVP-snapshot distribution join)  #
# ---------------------------------------- #
# Uses `wcvp_dist_tdwg3$area` (Method 3's local WCVP-snapshot join, TDWG
# level 3 area names) resolved to Natural Earth `admin` country names via a
# three-tier fallback:
#  1) Derive from the WGSRPD level-3/level-4 hierarchy: an area's level-4
#     children (level4.geojson) each carry an ISO-3166-1 country code;
#     level3.shp supplies the LEVEL3_COD <-> LEVEL3_NAM crosswalk used to
#     join the two. An area resolves to a single country when all its
#     level-4 children share one ISO code, or to multiple countries when
#     they span more than one (e.g. "Windward Is.", "Czechia-Slovakia") -
#     a taxon recorded under a multi-country area is counted present in
#     ALL of its resolved countries, since the source data doesn't say
#     which specific one. This covers the bulk of areas, including every
#     US state, Canadian province, Russian oblast, and
#     Brazilian/Mexican/Chinese/Argentine/Chilean/Australian region.
#  2) Areas that map onto exactly one whole country are often NOT
#     subdivided any further, so they have zero level-4 children at all
#     (tier 1 finds nothing for them) - e.g. "Suriname", "Taiwan". For
#     these, fall back to a direct match against `country_polys$admin`.
#  3) A residual set are still genuine naming variants (e.g. "Turkiye" vs
#     Natural Earth's "Turkey", "DR Congo" vs "Democratic Republic of the
#     Congo") or true multi-country composites with no level-4 breakdown
#     in level4.geojson at all (e.g. "Leeward Is.", "Netherlands
#     Antilles") - these go through `manual_country_overrides` below,
#     which also covers "South China Sea" (disputed reefs, no assigned
#     ISO code), "NW. Balkan Pen." (no corresponding TDWG level-3 area at
#     all), and "Sudan-South Sudan" (a pre-2011, pre-split legacy label
#     with no current level-3 area), assigned per explicit user
#     instruction rather than derived.
# ENCODING=LATIN1, not UTF-8. The shapefile's .dbf is Latin-1 encoded and has
# no .cpg sidecar declaring otherwise. Passing ENCODING=UTF-8 asserts the
# bytes are ALREADY UTF-8, so GDAL passes them through untranscoded and they
# get tagged UTF-8 while actually being Latin-1 - i.e. invalid UTF-8. The
# effect was that every accented area name silently failed to compare equal
# to the (correctly UTF-8) name in wcvp_distribution, so Føroyar, Galápagos,
# Juan Fernández Is. and Québec were dropped from the WCFP source entirely.
# Verified byte-for-byte: Québec read as 51 75 e9 62 65 63 under the old
# option vs 51 75 c3 a9 62 65 63 in wcvp_distribution; with LATIN1 the two
# match exactly and all four areas resolve (to Canada, Ecuador, Chile and
# the Faroe Islands respectively).
level3_shp <- st_read(file.path(inputs_dir, "level3/level3.shp"), quiet = TRUE,
                      options = "ENCODING=LATIN1") %>%
  st_drop_geometry() %>%
  distinct(LEVEL3_COD, LEVEL3_NAM)

level4 <- st_read(file.path(inputs_dir, "level4.geojson"), quiet = TRUE) %>%
  st_drop_geometry() %>%
  distinct(Level3_cod, ISO_Code)

# ISO 3166-1 alpha-2 -> Natural Earth admin country name. iso_a2_eh (the
# "everyone happy" variant) is used as a fallback for disputed territories
# that carry the placeholder "-99" in iso_a2, matching how iso_a3_eh is
# used as a fallback for GRIN's country codes in Source 2 above.
stopifnot("iso_a2" %in% names(country_polys))
iso2_to_country <- country_polys %>%
  st_drop_geometry() %>%
  transmute(
    iso2 = if_else(is.na(iso_a2) | iso_a2 == "-99", iso_a2_eh, iso_a2),
    admin
  ) %>%
  filter(!is.na(iso2), nzchar(iso2)) %>%
  distinct(iso2, .keep_all = TRUE)

# One row per level-3 area name, with a list-column of every country its
# level-4 children resolve to (length 1 for single-country areas, >1 for
# genuinely multi-country ones). Areas with no level-4 children at all
# simply don't appear here, and fall through to tiers 2/3 below.
level3_country_map <- level3_shp %>%
  inner_join(level4, by = c("LEVEL3_COD" = "Level3_cod")) %>%
  left_join(iso2_to_country, by = c("ISO_Code" = "iso2")) %>%
  filter(!is.na(admin)) %>%
  distinct(LEVEL3_NAM, admin) %>%
  group_by(LEVEL3_NAM) %>%
  summarise(countries = list(unique(admin)), .groups = "drop")

manual_country_overrides <- list(
  # Multi-country composites with no level-4 breakdown in level4.geojson:
  "Central American Pacific Is." = c("French Polynesia", "Costa Rica", "Colombia"),
  "Christmas I."                 = c("Indian Ocean Territories"),
  "Cocos (Keeling) Is."          = c("Indian Ocean Territories"),
  "Czechia-Slovakia"             = c("Czechia", "Slovakia"),
  "Gulf of Guinea Is."           = c("Equatorial Guinea", "São Tomé and Principe"),
  "Howland-Baker Is."            = c("United States of America"),
  "Leeward Is."                  = c("Antigua and Barbuda", "Anguilla", "Venezuela",
                                      "British Virgin Islands", "Montserrat",
                                      "Saint Kitts and Nevis", "United States Virgin Islands", "France"),
  "Marianas"                     = c("Guam", "Northern Mariana Islands"),
  "Mozambique Channel Is."       = c("France"),
  "Netherlands Antilles"         = c("Netherlands", "Curaçao"),
  "Southwest Caribbean"          = c("Colombia", "Honduras", "Nicaragua"),
  "Tokelau-Manihiki"             = c("Cook Islands", "American Samoa", "New Zealand"),
  "Wake I."                      = c("United States of America"),
  "Windward Is."                 = c("Barbados", "Dominica", "Grenada", "Saint Lucia",
                                      "Saint Vincent and the Grenadines", "France"),
  # Naming variants with no level-4 breakdown to derive from instead:
  "DR Congo"          = "Democratic Republic of the Congo",
  "Republic of Congo" = "Republic of the Congo",
  "Eswatini"          = "eSwatini",
  "French Guiana"     = "France",
  "Great Britain"     = "United Kingdom",
  "Kirgizstan"        = "Kyrgyzstan",
  "Panamá"            = "Panama",
  "Réunion"           = "France",
  "Svalbard"          = "Norway",
  "Türkey"            = "Turkey",
  "Türkey-in-Europe"  = "Turkey",
  "Turkiye"           = "Turkey",
  "Turkiye-in-Europe" = "Turkey",
  "Yakutiya"          = "Russia",
  # Manually assigned per user instruction (no TDWG level-3 area/level-4
  # breakdown exists for these):
  "NW. Balkan Pen."   = c("Croatia", "Slovenia", "Bosnia and Herzegovina"),
  "South China Sea"   = c("Philippines", "China", "Vietnam", "Malaysia", "Brunei"),
  "Sudan-South Sudan" = "Sudan"
)

wcvp_country_raw <- wcvp_dist_tdwg3 %>%
  filter(!is.na(area), nzchar(area)) %>%
  transmute(WCFP_ID, country_raw = str_squish(area), occurrence_status)

wcvp_country_lookup <- wcvp_country_raw %>%
  distinct(country_raw) %>%
  rowwise() %>%
  mutate(countries = list({
    matched <- level3_country_map$countries[level3_country_map$LEVEL3_NAM == country_raw]
    if (length(matched) == 1) {
      matched[[1]]
    } else if (country_raw %in% country_polys$admin) {
      country_raw
    } else if (country_raw %in% names(manual_country_overrides)) {
      manual_country_overrides[[country_raw]]
    } else {
      character(0)
    }
  })) %>%
  ungroup()

unresolved_areas <- wcvp_country_lookup$country_raw[lengths(wcvp_country_lookup$countries) == 0]
if (length(unresolved_areas) > 0) {
  message(length(unresolved_areas), " TDWG level-3 area names have no level-4 country ",
          "breakdown, no direct country name match, and no manual override, so are ",
          "dropped from the WCFP source entirely:\n",
          paste(sort(unresolved_areas), collapse = ", "))
}

# Joined back to a taxon name via WCFP_ID (a stable identifier already
# carried through from wcfp_plantlist), rather than reconstructing a
# sentence-cased binomial - more robust for hybrids/infraspecifics.
wcvp_country_long <- wcvp_country_raw %>%
  left_join(wcvp_country_lookup, by = "country_raw") %>%
  tidyr::unnest(countries) %>%
  rename(country = countries) %>%
  left_join(wcfp_plantlist %>% distinct(WCFP_ID, wcfp_name_match), by = "WCFP_ID") %>%
  filter(!is.na(wcfp_name_match))

# A species can be native in one TDWG L3 area of a country and introduced
# in another (or both) - aggregated here as a comma-separated list of
# distinct statuses per (country, taxon), e.g. "introduced, native",
# rather than collapsing to a single value.
src_wcfp <- wcvp_country_long %>%
  filter(country %in% country_polys$admin) %>%
  distinct(country, taxon = wcfp_name_match, occurrence_status) %>%
  group_by(country, taxon) %>%
  summarise(occurrence_status = paste(sort(unique(occurrence_status)), collapse = ", "), .groups = "drop")

# ---------------------------------------- #
#   Summary: distinct taxa & countries     #
#   per data source                        #
# ---------------------------------------- #
source_summary <- bind_rows(
  tibble(data_source = "occurrences",
         n_distinct_taxa = n_distinct(src_occurrences$taxon),
         n_distinct_countries = n_distinct(src_occurrences$country)),
  tibble(data_source = "GRIN",
         n_distinct_taxa = n_distinct(src_grin$taxon),
         n_distinct_countries = n_distinct(src_grin$country)),
  tibble(data_source = "WCFP",
         n_distinct_taxa = n_distinct(src_wcfp$taxon),
         n_distinct_countries = n_distinct(src_wcfp$country))
)

summary_wb <- createWorkbook()
addWorksheet(summary_wb, "Summary")
writeData(summary_wb, "Summary", source_summary, headerStyle = createStyle(textDecoration = "bold"))
setColWidths(summary_wb, "Summary", cols = 1:3, widths = c(18, 18, 22))
summary_xlsx_path <- file.path(combined_output_dir, "WCFP_data_source_summary_counts.xlsx")
saveWorkbook(summary_wb, summary_xlsx_path, overwrite = TRUE)
message("Data source summary counts saved to: ", summary_xlsx_path)

# ---------------------------------------- #
#   Combine the 3 sources                  #
# ---------------------------------------- #
combined_long <- bind_rows(
  src_occurrences %>% mutate(data_source_occurrences = "Y"),
  src_grin        %>% mutate(data_source_GRIN = "Y") %>%
    rename(status_GRIN = grin_status),
  src_wcfp        %>% mutate(data_source_WCFP = "Y") %>%
    rename(status_WCFP = occurrence_status)
)

# Native / non-native status is now tracked for BOTH GRIN and WCFP (the
# occurrence records still carry no such field - see Source 1). Each
# collapses to "" (then NA) when that source has no record for the pair, or
# to a comma-separated list of distinct statuses ("native, non-native") when
# the source reports both within the same country.
combined <- combined_long %>%
  group_by(country, taxon) %>%
  summarise(
    data_source_occurrences  = if (any(!is.na(data_source_occurrences))) "Y" else NA_character_,
    data_source_GRIN         = if (any(!is.na(data_source_GRIN)))        "Y" else NA_character_,
    data_source_WCFP         = if (any(!is.na(data_source_WCFP)))        "Y" else NA_character_,
    occurrences_data_source  = paste(sort(unique(na.omit(occurrences_data_source))), collapse = ", "),
    status_GRIN              = paste(sort(unique(na.omit(status_GRIN))), collapse = ", "),
    status_WCFP              = paste(sort(unique(na.omit(status_WCFP))), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    occurrences_data_source = na_if(occurrences_data_source, ""),
    status_GRIN = na_if(status_GRIN, ""),
    status_WCFP = na_if(status_WCFP, ""),
    n_sources = (!is.na(data_source_occurrences)) + (!is.na(data_source_GRIN)) + (!is.na(data_source_WCFP))
  )

# Keep only taxa confirmed by >=2 of the 3 methods (the "combined method"
# from the header aim), and only countries with at least one such taxon.
combined_confirmed <- combined %>%
  filter(n_sources >= 2) %>%
  select(-n_sources) %>%
  arrange(country, taxon)

# Look up each confirmed species' accepted name + author string from
# wcfp_plantlist (matched via wcfp_name_match), for display in the output
# workbook and leaflet map below - kept as separate columns rather than
# replacing `taxon`, which is still used for grouping/joins throughout
# this section.
wcfp_taxon_lookup <- wcfp_plantlist %>%
  transmute(
    wcfp_name_match = str_squish(wcfp_name_match),
    taxon_name_accepted,
    taxon_authors_accepted
  ) %>%
  filter(!is.na(wcfp_name_match)) %>%
  distinct(wcfp_name_match, .keep_all = TRUE)

combined_confirmed <- combined_confirmed %>%
  left_join(wcfp_taxon_lookup, by = c("taxon" = "wcfp_name_match"))

## ------------------------------------------------------------------ ##
## IUCN Red List category, from wcfp_plantlist$red_list_category_code.
##
## Defined here rather than in the Red List map section further down,
## because all three taxa-list workbooks carry a Red List column and are
## written long before that section runs. The maps then reuse this lookup
## instead of rebuilding it, so the workbooks and the maps cannot disagree
## about what category a taxon is in.
## ------------------------------------------------------------------ ##

# Pre-2001 "Lower Risk" subcategories, mapped to their modern equivalents:
# LR/lc -> LC and LR/nt -> NT are direct, and LR/cd (conservation dependent) is
# conventionally treated as NT since the category was retired. IUCN_LC is a
# single malformed cell, plainly meant to be LC.
RL_RECODE <- c("LR/LC" = "LC", "LR/NT" = "NT", "LR/CD" = "NT", "IUCN_LC" = "LC")

# Most severe first - used to break ties if a future plantlist ever gives one
# taxon two categories, and to validate that every code is one we recognise.
RL_SEVERITY <- c("EX", "EW", "CR", "EN", "VU", "NT", "LC", "DD", "NE")

# Category NAMES only. The code is reported in its own adjacent field in both
# the workbooks and the popup tables, so repeating it here as "Vulnerable (VU)"
# would duplicate the neighbouring cell. The map LAYER names further down keep
# the combined form, since a layer has no second field to carry the code.
RL_LABEL <- c(
  EX = "Extinct", EW = "Extinct in the Wild",
  CR = "Critically Endangered", EN = "Endangered",
  VU = "Vulnerable", NT = "Near Threatened",
  LC = "Least Concern", DD = "Data Deficient",
  NE = "Not Evaluated"
)

wcfp_redlist_lookup <- wcfp_plantlist %>%
  transmute(
    wcfp_name_match = str_squish(wcfp_name_match),
    # as.character() is deliberate: wcfp_plantlist is read with read_excel's
    # default guess_max = 1000, and this column is sparse (10,909 non-NA of
    # 26,632), so a run of empty leading cells could type it as logical.
    .raw    = toupper(str_squish(as.character(red_list_category_code))),
    rl_code = coalesce(unname(RL_RECODE[.raw]), .raw),
    # Blank or missing becomes NE (Not Evaluated), IUCN's own code for it.
    rl_code = if_else(is.na(rl_code) | !nzchar(rl_code), "NE", rl_code)
  ) %>%
  filter(!is.na(wcfp_name_match), nzchar(wcfp_name_match)) %>%
  mutate(.sev = match(rl_code, RL_SEVERITY)) %>%
  arrange(wcfp_name_match, .sev) %>%
  distinct(wcfp_name_match, .keep_all = TRUE) %>%
  select(wcfp_name_match, rl_code)

# Fails loudly rather than silently bucketing an unknown code into "Not
# assessed", which would understate threat.
.rl_unknown <- setdiff(unique(wcfp_redlist_lookup$rl_code), RL_SEVERITY)
if (length(.rl_unknown) > 0) {
  stop("unrecognised Red List category code(s) in the plantlist: ",
       paste(sort(.rl_unknown), collapse = ", "),
       " - add to RL_RECODE or RL_SEVERITY before continuing.")
}

# A taxon with no row in the plantlist lookup is NE, same as a blank
# cell. attach_redlist() is applied to each of the three confirmed tables so
# the country, L4 and L3 workbooks all resolve the category identically.
attach_redlist <- function(df) {
  df %>%
    left_join(wcfp_redlist_lookup, by = c("taxon" = "wcfp_name_match")) %>%
    mutate(rl_code = coalesce(rl_code, "NE"))
}

combined_confirmed <- attach_redlist(combined_confirmed)

qualifying_countries <- sort(unique(combined_confirmed$country))
cat(length(qualifying_countries), "countries have >=1 taxon confirmed by >=2 of the 3 data sources.\n")
cat(nrow(combined_confirmed), "total (country, taxon) rows across all qualifying countries.\n")

# ---------------------------------------- #
#   Write one workbook, one sheet/country  #
# ---------------------------------------- #
make_sheet_name <- function(country, used) {
  nm <- gsub("[\\\\/\\?\\*\\[\\]:]", "", country)
  nm <- substr(nm, 1, 31)
  base <- nm
  i <- 1
  while (nm %in% used) {
    suffix <- paste0("_", i)
    nm <- paste0(substr(base, 1, 31 - nchar(suffix)), suffix)
    i <- i + 1
  }
  nm
}

wb <- createWorkbook()
used_sheet_names <- character(0)

# ---- Summary sheet (first sheet): country x number of distinct taxa ----
country_summary <- combined_confirmed %>%
  count(country, name = "n_distinct_taxa") %>%
  arrange(country)

addWorksheet(wb, "Summary")
writeData(wb, "Summary", country_summary, headerStyle = createStyle(textDecoration = "bold"))
freezePane(wb, "Summary", firstRow = TRUE)
setColWidths(wb, "Summary", cols = 1:2, widths = c(38, 20))
used_sheet_names <- c(used_sheet_names, "Summary")

for (ctry in qualifying_countries) {
  sheet_name <- make_sheet_name(ctry, used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)

  country_data <- combined_confirmed %>%
    filter(country == ctry) %>%
    transmute(
      taxa = coalesce(taxon_name_accepted, taxon),
      authority = taxon_authors_accepted,
      iucn_red_list_category_code = rl_code,
      iucn_red_list_category      = unname(RL_LABEL[rl_code]),
      data_source_occurrences, occurrences_data_source,
      data_source_GRIN, status_GRIN, data_source_WCFP, status_WCFP
    )

  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, country_data, headerStyle = createStyle(textDecoration = "bold"))
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:10,
               widths = c(38, 22, 16, 24, 22, 22, 14, 20, 14, 20))
}

combined_xlsx_path <- file.path(combined_output_dir, "WCFP_taxa-lists_by_country_2plus_sources.xlsx")
saveWorkbook(wb, combined_xlsx_path, overwrite = TRUE)
message("Combined per-country taxa workbook saved to: ", combined_xlsx_path)




# ---------------------------------------------------------------------#


## Map: Species richness map - species confirmed by >=2 of the 3 data
## sources (occurrences, GRIN, WCFP), by country
## Author: Sarah Gora
## Three versions:
##   V1 — Simple style (original)
##   V2 — Green gradient with graticules (styled after grid maps)
##   V3 — Magma palette with graticules



# ==============================================================================
# OUTPUT DIRECTORY
# ==============================================================================
output_dir <- outputs_dir            # set near the top; user runs get their own folder
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 1) LOAD AND PREPARE DATA
# ==============================================================================
# Richness here comes from `combined_confirmed` (the "Combined method"
# section above) - taxa confirmed present in a country by >=2 of the 3
# data sources - not raw single-source occurrence records.
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") %>%
  fold_breakaway_polygons()
world_4326 <- st_transform(world, 4326)

# Calculate country-level species richness from the combined method
country_richness <- combined_confirmed %>%
  count(country, name = "richness")

# Join richness to map polygons
map_data <- world_4326 %>%
  left_join(country_richness, by = c("admin" = "country"))

# Check for unmatched countries
unmatched <- country_richness %>%
  filter(!country %in% world_4326$admin)

if (nrow(unmatched) > 0) {
  message("WARNING: ", nrow(unmatched), " countries in your data did not match the world map:")
  print(unmatched)
} else {
  message("All countries matched successfully.")
}

# Summary stats
max_richness <- max(country_richness$richness, na.rm = TRUE)
top_countries <- country_richness %>%
  filter(richness == max_richness)

message("Max taxa richness per country: ", max_richness)
message("Country/countries with max richness:")
for (i in seq_len(nrow(top_countries))) {
  message("  ", top_countries$country[i], ": ", top_countries$richness[i], " taxa")
}
message("Total distinct taxa (confirmed by >=2 sources): ", n_distinct(combined_confirmed$taxon))
message("Total countries with data: ", nrow(country_richness))





# ==============================================================================
# VERSION 1 — Simple style (magma, no graticules)
# ==============================================================================
p1 <- ggplot(map_data) +
  geom_sf(aes(fill = richness), color = "gray40", size = 0.15) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    na.value = "gray95",
    name = "Number of food plant taxa"
  ) +
  theme(
    panel.grid = element_line(color = "transparent"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.03, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.box.margin = margin(0, 0, 0, 0)
  )

print(p1)

# Save PNG
output_file_v1_png <- file.path(output_dir, "WCFP_taxa_distribution_country_map_v1_simple.png")
ragg::agg_png(output_file_v1_png, width = 12, height = 7, units = "in", res = 300)
print(p1)
dev.off()
message("V1 PNG saved to: ", output_file_v1_png)

# Save PDF
output_file_v1_pdf <- file.path(output_dir, "WCFP_taxa_distribution_country_map_v1_simple.pdf")
ggsave(
  filename = output_file_v1_pdf,
  plot     = p1,
  width    = 12,
  height   = 7,
  device   = cairo_pdf,
  bg       = "white"
)
message("V1 PDF saved to: ", output_file_v1_pdf)

message("V1 map saved to: ", output_dir)




## ============================================================================ ##
## Interactive leaflet map: species richness by country (>=2 data sources),
## click a country to see its confirmed species list and which source(s)
## (occurrences / GRIN / WCFP) and native/introduced status confirmed each one.
## Author: Sarah Gora
## ============================================================================ ##

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

# Escapes a string for safe embedding inside a single-quoted JS string
# literal (e.g. within an inline onclick="..." HTML attribute) - distinct
# from escape_html(), which is for HTML text content, not JS string syntax.
# Needed because some country names contain an apostrophe (e.g. "Cote
# d'Ivoire"), which would otherwise break the JS string.
escape_js <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("'",  "\\'",  x, fixed = TRUE)
  x
}

## ---------------------------------------------------------------- ##
## saveWidget() replacement that self-contains without pandoc.
##
## htmlwidgets::saveWidget(selfcontained = TRUE) shells out to pandoc,
## which parses the whole document into memory. The L4 maps further down
## embed ~420k (area, taxon) rows, so their raw HTML is ~135 MB, and
## pandoc dies on it with "osCommitMemory ... the paging file is too
## small" (error 251) on a 16 GB machine. Using saveWidget also makes the
## whole script depend on a pandoc install being discoverable, which it
## is not by default outside RStudio (plain Rscript fails with "Saving a
## widget with selfcontained = TRUE requires pandoc").
##
## The dependencies themselves are only ~1 MB - it is purely the embedded
## data that is large - so they are inlined here directly instead: write
## the widget with selfcontained = FALSE, then splice each <script src>/
## <link href> pointing into libdir back in as an inline <script>/<style>.
## Any url() asset a stylesheet refers to (e.g. Leaflet's layers.png, used
## by the layer-toggle control on the by-source maps) is base64-encoded
## into a data: URI, so the finished file has nothing left pointing at the
## lib folder once it is deleted.
##
## Only the handful of short <head> tag lines are rewritten; the single
## ~141 MB line holding the widget's JSON payload is passed straight
## through, which is what keeps memory use flat.
## ---------------------------------------------------------------- ##
if (!requireNamespace("base64enc", quietly = TRUE)) install.packages("base64enc")

save_widget_selfcontained <- function(widget, file, libdir) {
  htmlwidgets::saveWidget(widget, file = file, selfcontained = FALSE, libdir = libdir)

  base_dir <- dirname(normalizePath(file))
  lib_name <- basename(libdir)
  lines    <- readr::read_lines(file)

  as_data_uri <- function(path) {
    mime <- switch(tolower(tools::file_ext(path)),
                   png  = "image/png",   gif   = "image/gif",
                   jpg  = "image/jpeg",  jpeg  = "image/jpeg",
                   svg  = "image/svg+xml",
                   woff = "font/woff",   woff2 = "font/woff2",
                   ttf  = "font/ttf",    eot   = "application/vnd.ms-fontobject",
                   "application/octet-stream")
    paste0("data:", mime, ";base64,", base64enc::base64encode(path))
  }

  inline_css <- function(rel) {
    p <- file.path(base_dir, rel)
    if (!file.exists(p)) return(NULL)
    css <- paste(readr::read_lines(p), collapse = "\n")
    for (u in unique(unlist(regmatches(css, gregexpr("url\\([^)]+\\)", css))))) {
      ref <- trimws(sub("\\)$", "", sub("^url\\(", "", u)))
      ref <- gsub('^["\']|["\']$', "", ref)
      # Skip anything that isn't a local file reference. The "#" case is not
      # hypothetical: leaflet.css carries the legacy IE hack url(#default#VML),
      # which strips to an empty path and would otherwise resolve to the lib
      # directory itself.
      if (!nzchar(ref) || grepl("^(data:|https?:|//|#)", ref)) next
      asset <- file.path(dirname(p), sub("[?#].*$", "", ref))
      # file.exists() is TRUE for directories too, hence the dir.exists() guard.
      if (nzchar(basename(asset)) && file.exists(asset) && !dir.exists(asset)) {
        css <- gsub(u, paste0("url(", as_data_uri(asset), ")"), css, fixed = TRUE)
      }
    }
    c("<style>", css, "</style>")
  }

  inline_js <- function(rel) {
    p <- file.path(base_dir, rel)
    if (!file.exists(p)) return(NULL)
    c("<script>", readr::read_lines(p), "</script>")
  }

  css_re <- sprintf('<link[^>]+href="(%s/[^"]+\\.css)"', lib_name)
  js_re  <- sprintf('<script[^>]+src="(%s/[^"]+\\.js)"',  lib_name)

  out <- vector("list", length(lines))
  for (i in seq_along(lines)) {
    ln <- lines[[i]]
    out[[i]] <- ln
    # Only the short <head> tag lines can match. nchar(type = "bytes") is
    # what makes skipping the giant payload line cheap.
    if (nchar(ln, type = "bytes") > 1000L) next
    m_css <- regmatches(ln, regexec(css_re, ln))[[1]]
    m_js  <- regmatches(ln, regexec(js_re,  ln))[[1]]
    if (length(m_css) == 2L) {
      rep <- inline_css(m_css[2]); if (!is.null(rep)) out[[i]] <- rep
    } else if (length(m_js) == 2L) {
      rep <- inline_js(m_js[2]);   if (!is.null(rep)) out[[i]] <- rep
    }
  }

  readr::write_lines(unlist(out), file)
  if (dir.exists(libdir)) unlink(libdir, recursive = TRUE)
  invisible(file)
}

## ---------------------------------------------------------------- ##
## Shared onRender JS body, used by every interactive map below.
##
## Two things, both of which can only happen once the widget has rendered:
##  1. set the blue ocean map background;
##  2. move Leaflet's +/- zoom control into the bottom-left corner, sitting
##     above the scale bar.
##
## The zoom control is relocated by moving its DOM node into the
## .leaflet-bottom.leaflet-left corner container, rather than by calling
## zoomControl.setPosition(): setPosition() APPENDS to the target corner,
## which would place the zoom control BELOW the scale bar (already added
## there during render). insertBefore(..., corner.firstChild) puts it above
## instead. Doing it on the DOM also avoids depending on Leaflet's private
## _controlCorners API, and on `this` being bound to the map object.
##
## The one-shot setTimeout retry covers the case where the controls have
## not been attached yet at the moment onRender fires.
## ---------------------------------------------------------------- ##
map_onrender_common_js <- "
  el.style.background = '#9dc3da';
  // Zero Leaflet's default 10px control margin for the map-source strip
  // only, so it sits flush against the bottom edge of the map, and strip
  // its border/shadow. Injected as a stylesheet rather than inline styles
  // because .leaflet-control's margin is set by Leaflet's own CSS and an
  // inline style on our inner div cannot reach the wrapper element.
  if (!document.getElementById('wcfp-map-css')) {
    var st = document.createElement('style');
    st.id = 'wcfp-map-css';
    st.textContent =
      '.leaflet-control.wcfp-mapsource{margin:0 !important;' +
      'border:none !important;box-shadow:none !important;' +
      'border-radius:0 !important;}' +
      // Leaflet's scale bar is 11px text inside a 2px rule by default, which
      // looks undersized next to the enlarged title and legend. Widening it
      // via scaleBarOptions(maxWidth) alone would not touch the type, so the
      // text and rule are scaled here.
      '.leaflet-control-scale-line{font-size:15px !important;' +
      'line-height:1.5 !important;padding:4px 9px 3px !important;' +
      'border-width:3px !important;}' +
      // Restores Leaflet's own rule for the lifted popup pane - see
      // liftPopupPane() below. Leaflet fades popups in with
      //   .leaflet-fade-anim .leaflet-popup           {opacity:0}
      //   .leaflet-fade-anim .leaflet-map-pane .leaflet-popup {opacity:1}
      // and the second selector stops matching once the pane is moved out
      // of .leaflet-map-pane, which would leave EVERY popup at opacity 0.
      '.leaflet-fade-anim .wcfp-popup-pane .leaflet-popup{opacity:1;}';
    document.head.appendChild(st);
  }
  function moveZoomAboveScale() {
    var zc     = el.querySelector('.leaflet-control-zoom');
    var corner = el.querySelector('.leaflet-bottom.leaflet-left');
    if (!zc || !corner) { return false; }
    if (corner.firstChild !== zc) { corner.insertBefore(zc, corner.firstChild); }
    return true;
  }
  if (!moveZoomAboveScale()) { setTimeout(moveZoomAboveScale, 200); }

  // Sit the map-source strip on the bottom border immediately LEFT of
  // Leaflet's attribution link. Both live in the .leaflet-bottom.leaflet-right
  // corner, but every Leaflet control is a block element, so by default they
  // stack vertically. Wrapping just these two in a flex row puts them side by
  // side without disturbing the legend, which stays a block above the row.
  // The row is appended last, so it sits at the very bottom of the corner;
  // Leaflet already gives .leaflet-control-attribution margin:0, and the
  // strip's own margin is zeroed by the CSS above, so the row is flush with
  // the bottom edge.
  function placeSourceStrip() {
    var strip  = el.querySelector('.leaflet-control.wcfp-mapsource');
    var attrib = el.querySelector('.leaflet-control-attribution');
    if (!strip || !attrib || !attrib.parentNode) { return false; }
    if (strip.parentNode && strip.parentNode.classList &&
        strip.parentNode.classList.contains('wcfp-attrib-row')) { return true; }
    var corner = attrib.parentNode;
    var row = document.createElement('div');
    row.className = 'wcfp-attrib-row';
    row.style.display = 'flex';
    row.style.alignItems = 'stretch';
    row.style.justifyContent = 'flex-end';
    row.style.clear = 'both';
    corner.appendChild(row);
    row.appendChild(strip);
    row.appendChild(attrib);
    return true;
  }
  if (!placeSourceStrip()) { setTimeout(placeSourceStrip, 200); }

  // Popups render ON TOP of the legend, hint box, layers control and
  // map-source strip rather than sliding underneath them.
  //
  // Leaflet puts popups in .leaflet-popup-pane (z-index 700), which lives
  // inside .leaflet-map-pane. That pane is itself z-index 400 AND carries a
  // translate3d transform, either of which makes it a stacking context - so
  // nothing inside it can ever paint above the control corners
  // (.leaflet-top / .leaflet-bottom, z-index 1000). Raising the popup pane's
  // own z-index does nothing, because a descendant cannot escape an
  // ancestor's stacking context. The pane has to leave the map pane.
  //
  // So it is re-parented to the map container, as a sibling of the map pane,
  // and given a z-index above the corners. Leaflet still positions popups in
  // LAYER coordinates - i.e. relative to the map pane's transform - so the
  // lifted pane mirrors that transform on every change and the popups land
  // exactly where they did before.
  //
  // Two of Leaflet's own rules key off the pane's ancestry and break when it
  // moves. The fade rule is restored in the stylesheet above; the zoom
  // transition is driven by .leaflet-zoom-anim, which Leaflet adds to the MAP
  // PANE for the 250ms of a zoom, so that class is mirrored too - without it
  // an open popup jumps to its new position instead of gliding with the map.
  // Captured here rather than inside liftPopupPane(), which may run from a
  // setTimeout where `this` is no longer the map.
  var wcfpMap = (this && typeof this.on === 'function') ? this : null;

  function liftPopupPane() {
    var container = (el.classList && el.classList.contains('leaflet-container'))
                      ? el : el.querySelector('.leaflet-container');
    var mapPane   = el.querySelector('.leaflet-map-pane');
    var popupPane = el.querySelector('.leaflet-popup-pane');
    if (!container || !mapPane || !popupPane) { return false; }
    if (popupPane.parentNode === container) { return true; }
    popupPane.className += ' wcfp-popup-pane';
    popupPane.style.zIndex = '1100';
    container.appendChild(popupPane);
    function syncPopupPane() {
      popupPane.style.transform       = mapPane.style.transform;
      popupPane.style.webkitTransform = mapPane.style.webkitTransform;
      var zooming = mapPane.className.indexOf('leaflet-zoom-anim') > -1;
      var marked  = popupPane.className.indexOf('leaflet-zoom-anim') > -1;
      if (zooming && !marked) {
        popupPane.className += ' leaflet-zoom-anim';
      } else if (!zooming && marked) {
        popupPane.className =
          popupPane.className.replace(/ ?leaflet-zoom-anim/g, '');
      }
    }
    syncPopupPane();
    // Leaflet's own events fire synchronously, in the same frame that moves
    // the map pane, so the popup tracks a drag exactly. A MutationObserver
    // alone would deliver a microtask later and leave the popup a frame
    // behind the map; it is kept as a backstop for any transform change that
    // does not come with an event.
    if (wcfpMap) {
      wcfpMap.on('move zoom zoomanim viewreset moveend zoomend', syncPopupPane);
    }
    if (window.MutationObserver) {
      new MutationObserver(syncPopupPane).observe(
        mapPane, { attributes: true, attributeFilter: ['style', 'class'] });
    }
    return true;
  }
  if (!liftPopupPane()) { setTimeout(liftPopupPane, 200); }
"

# Plain version for the two combined maps (no layers control to label).
map_base_onrender_js <- paste0("function(el, x) {", map_onrender_common_js, "}")

## ---------------------------------------------------------------- ##
## Slightly tighter initial view: contract a data bounding box toward its
## centre, so each map opens zoomed in a little rather than fitting the
## full data extent edge to edge.
##
## Must be paired with a fractional `zoomSnap` on the leaflet() call.
## Leaflet's fitBounds() snaps to WHOLE zoom levels by default, and zoom
## levels are powers of two, so a modest contraction like this would
## otherwise round straight back to the same zoom and have no visible
## effect at all.
##
## Returns a plain named numeric (not an sf bbox), which is all the
## as.numeric(bounds["xmin"]) style calls below need.
## ---------------------------------------------------------------- ##
# frac is the single knob for the opening zoom: the fraction of the data
# extent trimmed off each side. Higher = tighter. 0.12 opened too close,
# so it is dialled back to 0.06 (retaining 94% of the extent rather than
# 88%). Paired with zoomSnap = 0.05 on the leaflet() calls, so a change
# this small still lands on a different zoom level instead of rounding
# back to the untrimmed fit.
zoom_in_bbox <- function(bb, frac = 0.06) {
  v <- as.numeric(bb[c("xmin", "ymin", "xmax", "ymax")])
  names(v) <- c("xmin", "ymin", "xmax", "ymax")
  cx <- (v[["xmin"]] + v[["xmax"]]) / 2
  cy <- (v[["ymin"]] + v[["ymax"]]) / 2
  hw <- (v[["xmax"]] - v[["xmin"]]) / 2 * (1 - frac)
  hh <- (v[["ymax"]] - v[["ymin"]]) / 2 * (1 - frac)
  c(xmin = cx - hw, ymin = cy - hh, xmax = cx + hw, ymax = cy + hh)
}

# Shared <h4> heading for all four interactive maps. Defined once so the
# wording cannot drift between them; each map's own subtitle underneath is
# what distinguishes it.
MAP_HEADING <- paste0(
  "<h4 style='margin:0; font-size:20px; line-height:1.25; color:#000;'>",
  "Geographic distribution of food plant taxa in the World Checklist ",
  "of Food Plants (2026)</h4>"
)

## ------------------------------------------------------------------ ##
## "Number of food plant taxa" legend, drawn as a continuous gradient bar.
##
## Replaces the four near-identical blocks of discrete swatch rows that
## each map used to build for itself. Still built by hand rather than with
## addLegend(), for the original reason: addLegend() always orders a
## continuous palette ascending top-to-bottom, and these read greatest at
## the top.
##
## The bar is a CSS linear-gradient sampled from the SAME colorNumeric
## palette the polygons are filled with, so the legend cannot drift away
## from the map.
##
## Labels are placed by value (absolute top:%) rather than spaced evenly
## down the bar. pretty() breaks do not necessarily sit at even fractions
## of the bar's range, so spacing them evenly would put the numbers next
## to the wrong colours.
## ------------------------------------------------------------------ ##
## `subtitle` is for the Red List legends. Those swap a different legend in
## per layer, so they need the layer name on the legend - but the bold title
## stays "Number of food plant taxa" as on every other map, and the layer name goes
## underneath in smaller grey text.
# Sizing is deliberately compact - the legend is a reference scale, not a
# feature of the map, and at the original 150x22 bar with 17px type it was
# covering a visible slice of ocean on every map. Every dimension below is
# driven from bar_h/bar_w and the two font sizes, so the whole block can be
# rescaled from here without hunting through the markup.
build_gradient_legend <- function(pal, breaks, title = "Number of food plant taxa",
                                  subtitle = NULL, bar_h = 105, bar_w = 15) {
  breaks <- sort(unique(breaks[is.finite(breaks)]), decreasing = TRUE)
  if (length(breaks) < 2) return("")
  hi <- max(breaks)
  lo <- min(breaks)

  # 12 stops is a visually smooth ramp at this bar height.
  stops <- pal(seq(hi, lo, length.out = 12))
  bar <- paste0(
    "<div style='width:", bar_w, "px; height:", bar_h, "px; border:1px solid #999; ",
    "background:linear-gradient(to bottom,", paste(stops, collapse = ","), ");'></div>"
  )

  pct  <- (hi - breaks) / (hi - lo) * 100
  labs <- paste0(
    "<span style='position:absolute; left:0; top:", round(pct, 2), "%; ",
    "transform:translateY(-50%); white-space:nowrap;'>",
    scales::comma(breaks), "</span>",
    collapse = ""
  )

  paste0(
    "<div style='background:white; padding:7px 9px; border-radius:4px; ",
    "box-shadow:0 0 4px rgba(0,0,0,0.3); font-size:11px; color:#333;'>",
    "<div style='font-weight:bold; font-size:13px; margin-bottom:",
    if (is.null(subtitle)) "5px" else "1px", ";'>", title, "</div>",
    if (is.null(subtitle)) "" else
      paste0("<div style='font-size:10px; color:#555; margin-bottom:5px;'>",
             subtitle, "</div>"),
    "<div style='display:flex; gap:7px; align-items:flex-start;'>",
    bar,
    "<div style='position:relative; height:", bar_h, "px; min-width:42px;'>", labs, "</div>",
    "</div></div>"
  )
}

# Concatenates the sources confirming a single (country, taxon) row into
# one readable "Distribution data source" string, e.g. "Occurrence records (GBIF);
# GRIN-Global distribution data; WCFP distribution data (status: native,
# introduced)" - only the sources that actually confirmed it are
# included. Shared between the popup table and the downloadable per-
# country .xlsx so the two stay in sync.
build_data_source_label <- function(has_occ, occ_src, has_grin, has_wcfp,
                                    wcfp_status, grin_status = NA_character_) {
  parts <- character(0)
  if (!is.na(has_occ)) {
    parts <- c(parts, if (is.na(occ_src) || !nzchar(occ_src)) {
      "Occurrence records"
    } else {
      paste0("Occurrence records (", occ_src, ")")
    })
  }
  if (!is.na(has_grin)) {
    parts <- c(parts, if (is.na(grin_status) || !nzchar(grin_status)) {
      "GRIN-Global distribution data"
    } else {
      paste0("GRIN-Global distribution data (status: ", grin_status, ")")
    })
  }
  if (!is.na(has_wcfp)) {
    parts <- c(parts, if (is.na(wcfp_status) || !nzchar(wcfp_status)) {
      "WCFP distribution data"
    } else {
      paste0("WCFP distribution data (status: ", wcfp_status, ")")
    })
  }
  paste(parts, collapse = "; ")
}

## ================================================================== ##
## Client-side popup + downloadable taxa table (shared by all 4 maps)
##
## Every map gives each country / TDWG L4 area a clickable popup listing
## its taxa, with a button to download that list as .xlsx. On the two
## by-source maps this exists PER LAYER, so clicking a country on the
## GRIN layer lists what GRIN reports there, and downloads only that.
##
## Why the table is built in the browser rather than pre-rendered in R:
## the original design embedded each taxa list TWICE - once as popup HTML
## and once as the JSON behind the Download button. The combined L4 map
## alone was 135 MB that way; doing it for four source layers (WCFP is
## smeared across every L4 child of its L3 area, and ~66% of GRIN records
## across every L4 area in their country) would have run to several
## hundred MB. Here each list is embedded ONCE as JSON, and the popup
## table is generated on first open, roughly halving the payload.
##
## Payload shape passed to onRender(data = ...):
##   list(
##     areaNames = list(<areaKey> = <display name>),
##     headers   = list(<layer>   = <detail header> OR list(<h1>, <h2>, ...)),
##     layers    = list(<layer>   = list(<areaKey> =
##                        list(c(taxa, authority, <detail1>, <detail2>, ...), ...)))
##   )
## A layer's `headers` entry sets how many detail columns its rows carry: a
## bare string means one, a list means that many. The Red List maps pass two
## (category code, category name); every other map passes one.
## Named lists (not named vectors) throughout - jsonlite serialises those
## to JSON objects without the "named vector" deprecation warning that
## previously flooded the run log.
## ================================================================== ##

# escape_html() is for text content; HTML attribute values also need the
# double quote escaped, since the popup stub below puts values inside
# data-layer="..." / data-area="...".
escape_attr <- function(x) gsub('"', "&quot;", escape_html(x), fixed = TRUE)

# `df` must supply: area_key, taxa, authority, and one column per name in
# `detail_cols`. Every map but the Red List ones has a single detail column
# ("info"); the Red List maps pass two, so the category code and the category
# name are separate fields rather than one "Vulnerable (VU)" string. A row is
# emitted as a flat array - taxa, authority, then the detail values in order -
# and the popup JS takes its column count from the layer's `headers` entry, so
# R and the browser cannot drift out of step.
build_layer_export <- function(df, detail_cols = "info") {
  val_cols <- c("taxa", "authority", detail_cols)
  d <- df %>%
    select(all_of(c("area_key", val_cols))) %>%
    filter(!is.na(area_key), nzchar(as.character(area_key)), !is.na(taxa)) %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(all_of(c("authority", detail_cols)),
                  ~ ifelse(is.na(.x), "", .x))) %>%
    distinct() %>%
    arrange(area_key, taxa)

  lapply(split(d, d$area_key), function(x) {
    # do.call(Map, ...) stays vectorised over the whole area; an lapply over
    # seq_len(nrow(x)) would go row by row and is far slower on big areas.
    unname(do.call(Map, c(list(f = function(...) list(...)),
                          unname(as.list(x[val_cols])))))
  })
}

# The per-source tables (src_occurrences, src_grin, src_wcfp and their L4
# equivalents) carry `taxon` = wcfp_name_match, not the accepted name, so
# resolve display name + authority the same way the combined tables already do.
attach_taxon_names <- function(df) {
  df %>%
    left_join(wcfp_taxon_lookup, by = c("taxon" = "wcfp_name_match")) %>%
    mutate(taxa      = coalesce(taxon_name_accepted, taxon),
           authority = coalesce(taxon_authors_accepted, ""))
}

## ------------------------------------------------------------------ ##
## "Map data source" strip, below the "Number of food plant taxa" legend.
##
## Names ONLY the source of the map polygons - the taxon data sources are
## already reported by the by-source layers, the subtitle counts and the
## popups, so repeating them here would be redundant.
##
## Single line, semi-transparent grey fill, no border or shadow, sitting on
## the bottom border of the map immediately left of Leaflet's attribution
## link. Two things are needed for that, both keyed off the
## className = "wcfp-mapsource" passed to addControl() below:
##  - CSS injected in map_onrender_common_js zeroes Leaflet's default 10px
##    .leaflet-control margin (and its border/shadow/radius);
##  - placeSourceStrip() in the same JS moves this control and the
##    attribution control into a shared flex row, since Leaflet controls
##    are block elements and would otherwise stack vertically.
##
## It is still added to the bottomright corner after the legend, so the
## legend remains a block above the row.
## ------------------------------------------------------------------ ##
# Small instruction box telling the reader the polygons are clickable, shown
# on all four maps. Placed in the top-left corner below the title (and, on
# the by-source maps, below the "Distribution data source" layer toggle).
build_click_hint <- function(txt) {
  paste0("<div style='background:white; padding:4px 10px; border-radius:4px; ",
         "font-size:12px; color:#333; box-shadow:0 0 4px rgba(0,0,0,0.3);'>",
         txt, "</div>")
}

click_hint_country    <- build_click_hint(
  "Click map to view/download food plants taxa list per country")
click_hint_l4         <- build_click_hint(
  "Click map to view/download food plants taxa list per TDWG Level-4 botanical area")
# The by-source maps get two lines in the one box: how to switch layers,
# then what clicking does. The second line spells out "for the selected data
# source" because those popups list only the active layer's taxa - without
# that, the counts look inconsistent with the combined map.
click_hint_country_bs <- build_click_hint(paste0(
  "Select one distribution data source - layers are exclusive, not stacked<br>",
  "Click map to view/download food plants taxa list per country, for the selected data source"))
click_hint_l4_bs      <- build_click_hint(paste0(
  "Select one distribution data source - layers are exclusive, not stacked<br>",
  "Click map to view/download food plants taxa list per TDWG Level-4 botanical area, for the selected data source"))

build_map_source_strip <- function(src) {
  paste0(
    "<div style='background:rgba(110,110,110,0.35); color:#1a1a1a; ",
    "font-size:12px; padding:3px 10px; white-space:nowrap;'>",
    "<b>Map data source:</b> ", src, "</div>"
  )
}

# Country maps draw country_polys/continent_polys (both Natural Earth);
# the L4 maps draw level4_sf/level3_sf (both TDWG WGSRPD).
map_source_strip_country <- build_map_source_strip("Natural Earth")
map_source_strip_l4      <- build_map_source_strip("TDWG WGSRPD Level-3 / Level-4")

# Leaflet's popup maxWidth defaults to 300px, which would clip the taxa
# table however wide the table's own CSS is - so the width has to be raised
# here as well as in buildHtml(). minWidth stops narrow tables (a country
# with two taxa) collapsing to a cramped box. Shared by all ten polygon
# layers so every popup on every map is sized the same.
WCFP_POPUP_OPTIONS <- popupOptions(maxWidth = 820, minWidth = 500,
                                   autoPan = TRUE,
                                   # Generous vertical padding: the table is
                                   # tall, and after the re-pan triggered in
                                   # wcfp_popup_js this keeps it clear of the
                                   # top edge rather than flush against it.
                                   autoPanPadding = c(24, 40))

# The popup actually bound to each polygon: a tiny stub carrying its layer
# and area in data- attributes. Filled in on open by wcfp_popup_js.
popup_stub <- function(layer_label, area_key, area_name, has_data) {
  # ifelse() returns a value shaped like `test`, NOT like `yes`/`no`. So a
  # scalar has_data - which the three combined maps pass, since every polygon
  # they stub is known to have data - collapsed the whole vector to a SINGLE
  # stub. Leaflet then recycled that one stub across every polygon, so all
  # 235/592/368 areas carried the first area's data-area and every popup
  # opened the same taxa list. Silent: the maps rendered, the popups worked,
  # they were just all Alberta.
  #
  # Recycling has_data to the data length first makes the call correct
  # whatever shape the caller passes, and the length check below turns any
  # future recurrence into a hard stop rather than 17 plausible-looking
  # outputs.
  n  <- max(length(area_key), length(area_name), length(has_data))
  hd <- rep_len(as.logical(has_data), n)
  # ifelse(NA, ...) yields NA, which would reach the map as a literal "NA"
  # popup. An unknown is treated as "no data" - the safe direction.
  hd[is.na(hd)] <- FALSE
  out <- ifelse(
    hd,
    paste0("<div class=\"wcfp-popup\" data-layer=\"", escape_attr(layer_label),
           "\" data-area=\"", escape_attr(area_key), "\"></div>"),
    paste0("<b>", escape_html(area_name), "</b><br>No taxa recorded by this source")
  )
  if (length(out) != n) {
    stop("popup_stub() produced ", length(out), " stub(s) for ", n,
         " areas - every polygon would get the same popup.")
  }
  out
}

# XLSX writing uses the same CDN library as before, so downloading (but
# not viewing) needs an internet connection when the saved file is opened.
#
# Popups are filled via Leaflet's own popupopen event where possible, so
# that e.popup.update() can resize the popup around the injected table.
# A MutationObserver on the popup pane is kept as a fallback in case
# `this` is not bound to the map object; both paths guard on data-filled
# so a popup can never be built twice.
wcfp_popup_js <- "
function(el, x, data) {
  if (!window.XLSX) {
    var s = document.createElement('script');
    s.src = 'https://cdn.jsdelivr.net/npm/xlsx@0.18.5/dist/xlsx.full.min.js';
    document.head.appendChild(s);
  }
  var D      = data || {};
  var LAYERS = D.layers || {};
  var NAMES  = D.areaNames || {};
  var HDRS   = D.headers || {};
  // Layer labels that should NOT be repeated as the grey sub-heading in the
  // popup. The combined layer is listed because 'Confirmed by >=2 data
  // sources' is already stated in the map's own title box, so printing it
  // again above every taxa table is noise.
  var HIDELBL = D.hideLabels || [];
  function showLabel(layer) {
    for (var i = 0; i < HIDELBL.length; i++) {
      if (HIDELBL[i] === layer) { return false; }
    }
    return true;
  }

  function esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }
  function rowsFor(layer, area) {
    var L = LAYERS[layer] || {};
    return L[area] || [];
  }
  function nameFor(area) { return NAMES[area] || area; }
  // A layer's headers entry is either one string (single detail column) or an
  // array of them. Normalising to an array here is what lets the table below
  // render 1 or 2 detail columns from the same code.
  function headersFor(layer) {
    var h = HDRS[layer];
    if (h == null) { return ['Detail']; }
    return Object.prototype.toString.call(h) === '[object Array]' ? h : [h];
  }

  window.wcfpDownload = function(layer, area) {
    if (!window.XLSX) {
      alert('Download library is still loading - please try again in a moment.');
      return;
    }
    var rows = rowsFor(layer, area);
    var hdrs = headersFor(layer);
    var nm   = nameFor(area);
    var aoa  = [['Area', 'Layer', 'Taxa', 'Authority'].concat(hdrs)];
    for (var i = 0; i < rows.length; i++) {
      aoa.push([nm, layer].concat(rows[i]));
    }
    var ws = XLSX.utils.aoa_to_sheet(aoa);
    var wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, ws, 'Taxa');
    var safe = (nm + '_' + layer).replace(/[^A-Za-z0-9]+/g, '_');
    XLSX.writeFile(wb, safe + '_wcfp_taxa.xlsx');
  };

  function buildHtml(layer, area) {
    var rows = rowsFor(layer, area);
    var hdrs = headersFor(layer);
    // white-space:nowrap on the Taxa column makes the table's automatic
    // layout size that column to the LONGEST taxon name, so names never
    // wrap onto a second line; Authority and the detail column wrap into
    // whatever width is left. vertical-align:top keeps a wrapped
    // Authority/detail cell aligned with the top of its taxon name instead
    // of floating to the vertical middle of the row.
    var tdTaxa = \"style='white-space:nowrap; vertical-align:top; padding:3px 12px 3px 0;'\";
    var tdRest = \"style='vertical-align:top; padding:3px 12px 3px 0;'\";
    var tdLast = \"style='vertical-align:top; padding:3px 0;'\";
    var out  = [];
    out.push(\"<div style='min-width:500px; font-family:sans-serif;'>\");
    out.push(\"<div style='display:flex; justify-content:space-between; align-items:center; gap:8px;'>\");
    out.push('<span><b>' + esc(nameFor(area)) + '</b> &mdash; ' +
             rows.length.toLocaleString() + ' food plants taxa' +
             (showLabel(layer)
                ? \"<br><span style='font-size:10px; color:#666;'>\" + esc(layer) + '</span>'
                : '') +
             '</span>');
    out.push(\"<button class='wcfp-dl' style='font-size:10px; padding:2px 6px; \" +
             \"cursor:pointer; white-space:nowrap;'>Download</button>\");
    out.push('</div>');
    // overflow:auto (not just overflow-y) so a pathologically long name
    // scrolls horizontally rather than bursting the popup.
    out.push(\"<div style='max-height:300px; overflow:auto; margin-top:6px;'>\");
    out.push(\"<table style='font-size:11px; border-collapse:collapse; width:100%;'>\");
    // Only the final detail column drops the right padding, whatever the
    // column count, so a 2-detail Red List table is spaced like a 1-detail one.
    function tdFor(i) { return (i === hdrs.length - 1) ? tdLast : tdRest; }
    var head = '<td ' + tdTaxa + '>Taxa</td>' + '<td ' + tdRest + '>Authority</td>';
    for (var h = 0; h < hdrs.length; h++) {
      head += '<td ' + tdFor(h) + '>' + esc(hdrs[h]) + '</td>';
    }
    out.push(\"<tr style='font-weight:bold; border-bottom:1px solid #999; text-align:left;'>\" +
             head + '</tr>');
    for (var j = 0; j < rows.length; j++) {
      var tr = '<td ' + tdTaxa + '><i>' + esc(rows[j][0]) + '</i></td>' +
               '<td ' + tdRest + '>' + esc(rows[j][1]) + '</td>';
      for (var h2 = 0; h2 < hdrs.length; h2++) {
        tr += '<td ' + tdFor(h2) + '>' + esc(rows[j][2 + h2]) + '</td>';
      }
      out.push(\"<tr style='border-bottom:1px solid #eee;'>\" + tr + '</tr>');
    }
    out.push('</table></div></div>');
    return out.join('');
  }

  function fill(host) {
    if (!host || host.getAttribute('data-filled') === '1') { return false; }
    var layer = host.getAttribute('data-layer');
    var area  = host.getAttribute('data-area');
    host.innerHTML = buildHtml(layer, area);
    host.setAttribute('data-filled', '1');
    var btn = host.querySelector('.wcfp-dl');
    if (btn) {
      btn.addEventListener('click', function() { window.wcfpDownload(layer, area); });
    }
    return true;
  }

  var map = this;
  if (map && typeof map.on === 'function') {
    map.on('popupopen', function(e) {
      var p = e.popup;
      var root = (p && p.getElement) ? p.getElement() : null;
      if (!root) { return; }
      if (!fill(root.querySelector('.wcfp-popup'))) { return; }
      // Re-measure and re-pan now that the table actually exists.
      //
      // Deliberately NOT p.update(): update() calls _updateContent(), which
      // resets the popup's innerHTML back to the empty stub and discards the
      // table we just built. These three do the layout / position / pan work
      // without touching the content.
      //
      // _adjustPan() is the important one. Leaflet decides whether to pan the
      // map when the popup OPENS, at which point this popup is still an empty
      // stub and looks tiny, so it concludes no pan is needed. The table is
      // injected immediately afterwards and the popup grows upward - so a
      // country near the top of the map ends up with its table clipped off
      // the top edge. Re-running _adjustPan() after filling pans the map down
      // so the full table is visible.
      if (p._updateLayout)   { p._updateLayout(); }
      if (p._updatePosition) { p._updatePosition(); }
      if (p._adjustPan)      { p._adjustPan(); }
    });
  }

  var pane = el.querySelector('.leaflet-popup-pane') || el;
  if (window.MutationObserver) {
    new MutationObserver(function() {
      var hosts = el.querySelectorAll('.wcfp-popup:not([data-filled])');
      for (var k = 0; k < hosts.length; k++) { fill(hosts[k]); }
    }).observe(pane, { childList: true, subtree: true });
  }
}"

# Label used for this map's single layer, in the popup header, the .xlsx
# filename and the payload key below.
COMBINED_LAYER <- "Confirmed by ≥2 data sources"

# Toggle labels for the two combined choropleths, which each carry the OTHER
# resolution as a boundary-only overlay. Overlays (checkboxes), not base
# groups - unlike the data-source and Red List layers these are meant to be
# stacked, since the point is seeing how the two grids relate. Boundary-only:
# no fill, no popup, no taxa payload, so the second layer costs geometry
# alone rather than doubling the file.
LYR_L3 <- "TDWG Level-3 botanical countries"
LYR_L4 <- "TDWG Level-4 botanical areas"

# One row per (country, taxon) with the concatenated "Distribution data source" string
# that the popup's third column shows. vapply over rows because
# build_data_source_label() is inherently row-wise.
combined_country_rows <- combined_confirmed %>%
  mutate(
    area_key  = country,
    taxa      = coalesce(taxon_name_accepted, taxon),
    authority = coalesce(taxon_authors_accepted, ""),
    info      = vapply(seq_len(n()), function(i) build_data_source_label(
      data_source_occurrences[i], occurrences_data_source[i],
      data_source_GRIN[i], data_source_WCFP[i], status_WCFP[i], status_GRIN[i]
    ), character(1))
  )

country_payload <- list(
  areaNames = as.list(setNames(sort(unique(combined_confirmed$country)),
                               sort(unique(combined_confirmed$country)))),
  hideLabels = list(COMBINED_LAYER),
  headers   = setNames(list("Distribution data source"), COMBINED_LAYER),
  layers    = setNames(list(build_layer_export(combined_country_rows)), COMBINED_LAYER)
)

# Built directly from `country_polys`/`combined_confirmed` (both from the
# Combined method section above) rather than reusing `map_data` from the
# static V1 map below - keeps the leaflet map runnable on its own even if
# the static ggplot map section is skipped or hasn't run yet.
leaflet_country_richness <- combined_confirmed %>%
  count(country, name = "richness")

leaflet_map_data <- country_polys %>%
  left_join(leaflet_country_richness, by = c("admin" = "country")) %>%
  mutate(
    popup_html = ifelse(
      is.na(richness),
      paste0("<b>", escape_html(admin), "</b><br>No taxa confirmed by ≥2 data sources"),
      popup_stub(COMBINED_LAYER, admin, admin, TRUE)
    )
  )

# YlGn (ColorBrewer, colorblind-safe) - pale yellow at low values, dark
# green at high values. Unlike magma, this is already dark = high in its
# natural order, so no reverse is needed here.
richness_pal <- colorNumeric(
  palette = "YlGn",
  domain  = leaflet_map_data$richness,
  # White rather than the original grey90 - reads cleanly against the blue
  # ocean background. The legend no longer carries a "no data" swatch, so
  # this only has to work against the map itself.
  na.color = "#FFFFFF",
  reverse = FALSE
)

# Dissolve country polygons into continent-level outlines (thicker black
# border) drawn on top of the thin black country borders below, matching
# the mockup's country-border / coastline-outline distinction.
#
# st_wrap_dateline() after the union is required for continents that span
# the antimeridian (e.g. Oceania, with Pacific islands on both sides of
# +/-180 longitude) - without it, the dissolved polygon isn't split at the
# dateline, and Leaflet draws a straight line connecting a point near
# -180 to one near +180, which shows up as a spurious horizontal line
# crossing the entire map.
stopifnot("continent" %in% names(country_polys))
continent_polys <- country_polys %>%
  filter(continent != "Antarctica") %>%
  group_by(continent) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE)

legend_breaks <- pretty(leaflet_map_data$richness, n = 6)
legend_breaks <- sort(
  legend_breaks[legend_breaks >= min(leaflet_map_data$richness, na.rm = TRUE) &
                  legend_breaks <= max(leaflet_map_data$richness, na.rm = TRUE)],
  decreasing = TRUE
)

legend_html <- build_gradient_legend(richness_pal, legend_breaks)

# Centroid of each country's single largest polygon part (rather than the
# whole multipolygon's centroid), used below for the country name labels
# so archipelago/multi-part countries (e.g. Indonesia, Chile) get a label
# on their main landmass instead of potentially over open water between
# islands.
largest_part_centroids <- leaflet_map_data %>%
  st_cast("POLYGON", warn = FALSE) %>%
  mutate(part_area = st_area(.)) %>%
  group_by(admin) %>%
  slice_max(part_area, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  st_centroid()

# Scale each label's font size to its country's land area, so small
# countries (which get badly overlapped/cluttered at world-view zoom when
# every label is shown at a fixed size) shrink down while large countries
# (Russia, Canada, China) stay prominent. Uses sqrt(area) rather than raw
# area - area scales with the square of visual size, so a linear map from
# area to font-size would make every country except a handful of giants
# collapse to the minimum.
area_sqrt <- sqrt(as.numeric(largest_part_centroids$part_area))
label_font_px <- scales::rescale(area_sqrt, to = c(6, 16))
largest_part_centroids$label_html <- lapply(
  sprintf("<span style='font-size:%.0fpx;'>%s</span>",
          label_font_px, htmltools::htmlEscape(largest_part_centroids$admin)),
  htmltools::HTML
)

map_title_html <- sprintf(
  paste0(MAP_HEADING,
         "<div style='font-size:13px; color:#000; margin-top:3px;'>",
         "%s taxa across %s countries<br><span style='font-style:italic;'>Distribution of taxa confirmed in &ge;2 data sources</span></div>"),
  scales::comma(n_distinct(combined_confirmed$taxon)),
  scales::comma(length(qualifying_countries))
)

# Extent of all country polygons, used below to explicitly fit the map's
# initial view to the data on open - without this, the map can otherwise
# load at Leaflet's default center/zoom (or an unhelpful auto-fit based on
# just one layer, e.g. the label point markers) rather than showing the
# whole map area.
map_bounds <- zoom_in_bbox(st_bbox(leaflet_map_data))

species_richness_leaflet <- leaflet(leaflet_map_data,
                                    options = leafletOptions(zoomSnap = 0.05, zoomDelta = 0.5)) %>%
  fitBounds(
    lng1 = as.numeric(map_bounds["xmin"]),
    lat1 = as.numeric(map_bounds["ymin"]),
    lng2 = as.numeric(map_bounds["xmax"]),
    lat2 = as.numeric(map_bounds["ymax"])
  ) %>%
  # No tile basemap - every country is already drawn as a filled polygon
  # below (richness color, or white for no confirmed species), so a tile
  # layer would only be supplying the ocean color anyway, and tile
  # providers keep hitting API-key walls (CartoDB.Positron) or missing
  # options (no Esri equivalent). Setting the container background
  # directly gives an exact, dependency-free ocean color instead.
  addPolygons(
    fillColor    = ~richness_pal(richness),
    fillOpacity  = 0.8,
    color        = "black",
    weight       = 1,
    label        = ~lapply(paste0("<b>", escape_html(admin), "</b>: ",
                                  scales::comma(ifelse(is.na(richness), 0, richness)),
                                  " taxa"), htmltools::HTML),
    popup        = ~popup_html,
    popupOptions = WCFP_POPUP_OPTIONS,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data    = continent_polys,
    fill    = FALSE,
    color   = "black",
    weight  = 1,
    opacity = 1,
    options = pathOptions(interactive = FALSE)
  ) %>%
  # Country name labels from the base tile layer get covered by the
  # choropleth fill above. No free labels-only reference tile layer is
  # available (Esri has none, CartoDB's requires the same API key
  # CartoDB.Positron does) - added directly as our own permanent labels
  # instead, which as a vector layer we control always renders on top,
  # with no external tile dependency.
  addLabelOnlyMarkers(
    data = largest_part_centroids,
    label = largest_part_centroids$label_html,
    labelOptions = labelOptions(
      noHide = TRUE,
      textOnly = TRUE,
      direction = "center",
      style = list(
        "font-weight" = "bold",
        "color" = "#000000",
        "text-shadow" = "-1px -1px 0 #fff, 1px -1px 0 #fff, -1px 1px 0 #fff, 1px 1px 0 #fff, 0 0 3px #fff"
      )
    )
  ) %>%
  addControl(
    html = legend_html,
    position = "bottomright"
  ) %>%
  addControl(
    html = map_source_strip_country,
    position = "bottomright",
    className = "wcfp-mapsource"
  ) %>%
  # Title added before the "Click map to view/download" hint below so that, with
  # both in the same corner, the title stacks above it.
  addControl(
    html = map_title_html,
    position = "topleft"
  ) %>%
  addControl(
    html = click_hint_country,
    position = "topleft"
  ) %>%
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(maxWidth = 260, metric = TRUE, imperial = FALSE, updateWhenIdle = TRUE)
  ) %>%
  htmlwidgets::onRender(map_base_onrender_js) %>%
  # Builds each popup's taxa table and wires its Download button. The
  # payload is embedded once as JSON and the table rendered on first open
  # - see the wcfp_popup_js block above for why.
  htmlwidgets::onRender(wcfp_popup_js, data = country_payload)

leaflet_html_path <- file.path(combined_output_dir, "WCFP_taxa_distribution_country_map_interactive.html")
leaflet_lib_tmp <- file.path(combined_output_dir, "lib_tmp")
save_widget_selfcontained(
  widget = species_richness_leaflet,
  file   = leaflet_html_path,
  libdir = leaflet_lib_tmp
)
message("Interactive leaflet map saved to: ", leaflet_html_path)




## ============================================================================ ##
## Interactive leaflet map (by data source): richness per country for each
## of the 3 sources individually (occurrences, GRIN, WCFP), shown as
## toggleable checkbox layers, plus a fourth layer holding the strict
## ">=2 of 3 sources" combined result so it can be compared directly
## against each individual source. Built from the same src_occurrences/
## src_grin/src_wcfp tables and the same `combined_confirmed` (via
## `leaflet_country_richness`) as the combined map above, using the same
## country_polys/continent_polys/largest_part_centroids - no new input data.
##
## Note: since all 4 layers are opaque country polygons covering the same
## area, checking more than one at once shows only whichever is drawn on
## top (no blending) - this is meant for toggling one at a time to compare,
## not overlaying. The three individual sources start hidden, so the map
## opens showing only "Confirmed by >=2 data sources".
## Author: Sarah Gora
## ============================================================================ ##

richness_occurrences <- src_occurrences %>% count(country, name = "richness")
richness_grin        <- src_grin        %>% count(country, name = "richness")
richness_wcfp        <- src_wcfp        %>% count(country, name = "richness")

map_data_occurrences <- country_polys %>% left_join(richness_occurrences, by = c("admin" = "country"))
map_data_grin        <- country_polys %>% left_join(richness_grin,        by = c("admin" = "country"))
map_data_wcfp        <- country_polys %>% left_join(richness_wcfp,        by = c("admin" = "country"))

# Fourth layer: the ">=2 of 3 sources" combined result, reusing the same
# `leaflet_country_richness` (derived from `combined_confirmed`) that drives
# the main richness map above - so the strict combined view can be toggled
# against each individual source without leaving this map.
map_data_combined <- country_polys %>%
  left_join(leaflet_country_richness, by = c("admin" = "country"))

# One shared palette/domain across all 4 layers (rather than each scaled
# independently) so toggling between them is a fair visual comparison. The
# combined layer is included in the domain so its values can never fall
# outside the palette range.
by_source_max_richness <- max(
  c(map_data_occurrences$richness, map_data_grin$richness, map_data_wcfp$richness,
    map_data_combined$richness),
  na.rm = TRUE
)
by_source_pal <- colorNumeric(
  palette  = "YlGn",
  domain   = c(0, by_source_max_richness),
  # White rather than the original grey90 - reads cleanly against the blue
  # ocean background. The legend no longer carries a "no data" swatch, so
  # this only has to work against the map itself.
  na.color = "#FFFFFF",
  reverse  = FALSE
)

by_source_legend_breaks <- pretty(c(0, by_source_max_richness), n = 6)
by_source_legend_breaks <- sort(by_source_legend_breaks[by_source_legend_breaks <= by_source_max_richness],
                                 decreasing = TRUE)

by_source_legend_html <- build_gradient_legend(by_source_pal, by_source_legend_breaks)

# Used by both by-source maps (country and L4 below). Everything the other
# maps do (grey background + zoom control above the scale bar, via
# map_onrender_common_js), plus a "Distribution data source" heading on the layers
# control - addLayersControl() has no title argument of its own, so the
# heading is inserted into the control's DOM here instead. It goes inside
# .leaflet-control-layers-list (falling back to the control itself) so it
# sits above the checkboxes rather than above the collapsed-mode toggle.
by_source_onrender_js <- paste0("function(el, x) {", map_onrender_common_js, "
  function addLayersTitle() {
    var lc = el.querySelector('.leaflet-control-layers');
    if (!lc) { return false; }
    if (lc.querySelector('.wcfp-layers-title')) { return true; }
    var h = document.createElement('div');
    h.className = 'wcfp-layers-title';
    h.textContent = 'Distribution data source';
    h.style.fontWeight = 'bold';
    h.style.fontSize = '13px';
    h.style.marginBottom = '5px';
    var list = lc.querySelector('.leaflet-control-layers-list') || lc;
    list.insertBefore(h, list.firstChild);
    return true;
  }
  if (!addLayersTitle()) { setTimeout(addLayersTitle, 200); }
}")

# ---------------------------------------- #
#   Per-layer popups + downloadable tables  #
# ---------------------------------------- #
# One taxa list per (layer, country), so clicking a country while the GRIN
# layer is active lists what GRIN reports there and downloads only that.
LBL_OCC  <- "Occurrence records"
LBL_GRIN <- "GRIN-Global distribution data"
LBL_WCFP <- "WCFP distribution data"

bs_rows_occ  <- src_occurrences %>%
  transmute(area_key = country, taxon, info = occurrences_data_source) %>% attach_taxon_names()
bs_rows_grin <- src_grin %>%
  transmute(area_key = country, taxon, info = grin_status) %>% attach_taxon_names()
bs_rows_wcfp <- src_wcfp %>%
  transmute(area_key = country, taxon, info = occurrence_status) %>% attach_taxon_names()
bs_rows_comb <- combined_country_rows %>% select(area_key, taxa, authority, info)

bs_layer_names <- c(LBL_OCC, LBL_GRIN, LBL_WCFP, COMBINED_LAYER)

# ---------------------------------------- #
#   Title + per-source subtitle             #
# ---------------------------------------- #
# Four subtitle lines, one per toggleable layer, reporting what that source
# actually contributes. Counts come from the same src_* / combined_confirmed
# tables the layers are drawn from, so the numbers can't drift from the map.
bs_taxa_counts <- c(
  n_distinct(src_occurrences$taxon),
  n_distinct(src_grin$taxon),
  n_distinct(src_wcfp$taxon),
  n_distinct(combined_confirmed$taxon)
)
bs_area_counts <- c(
  n_distinct(src_occurrences$country),
  n_distinct(src_grin$country),
  n_distinct(src_wcfp$country),
  n_distinct(combined_confirmed$country)
)

by_source_title_html <- paste0(
  MAP_HEADING,
  "<div style='font-size:13px; color:#000; margin-top:3px;'>",
  paste0(bs_layer_names, ": ", scales::comma(bs_taxa_counts), " taxa across ",
         scales::comma(bs_area_counts), " countries", collapse = "<br>"),
  "</div>"
)

by_source_payload <- list(
  areaNames = as.list(setNames(country_polys$admin, country_polys$admin)),
  hideLabels = list(COMBINED_LAYER),
  headers   = setNames(as.list(c("Occurrence source", "Native/non-native status",
                                 "Native/non-native status", "Distribution data source")),
                       bs_layer_names),
  layers    = setNames(list(build_layer_export(bs_rows_occ),
                            build_layer_export(bs_rows_grin),
                            build_layer_export(bs_rows_wcfp),
                            build_layer_export(bs_rows_comb)),
                       bs_layer_names)
)

by_source_leaflet <- leaflet(options = leafletOptions(zoomSnap = 0.05, zoomDelta = 0.5)) %>%
  fitBounds(
    lng1 = as.numeric(map_bounds["xmin"]), lat1 = as.numeric(map_bounds["ymin"]),
    lng2 = as.numeric(map_bounds["xmax"]), lat2 = as.numeric(map_bounds["ymax"])
  ) %>%
  addPolygons(
    data        = map_data_occurrences,
    fillColor   = ~by_source_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1,
    label       = ~lapply(paste0("<b>", escape_html(admin), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_OCC, admin, admin, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_OCC,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_grin,
    fillColor   = ~by_source_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1,
    label       = ~lapply(paste0("<b>", escape_html(admin), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_GRIN, admin, admin, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_GRIN,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_wcfp,
    fillColor   = ~by_source_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1,
    label       = ~lapply(paste0("<b>", escape_html(admin), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_WCFP, admin, admin, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_WCFP,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_combined,
    fillColor   = ~by_source_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1,
    label       = ~lapply(paste0("<b>", escape_html(admin), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(COMBINED_LAYER, admin, admin, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = COMBINED_LAYER,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data    = continent_polys,
    fill    = FALSE,
    color   = "black",
    weight  = 1,
    opacity = 1,
    options = pathOptions(interactive = FALSE)
  ) %>%
  addLabelOnlyMarkers(
    data = largest_part_centroids,
    label = largest_part_centroids$label_html,
    labelOptions = labelOptions(
      noHide = TRUE,
      textOnly = TRUE,
      direction = "center",
      style = list(
        "font-weight" = "bold",
        "color" = "#000000",
        "text-shadow" = "-1px -1px 0 #fff, 1px -1px 0 #fff, -1px 1px 0 #fff, 1px 1px 0 #fff, 0 0 3px #fff"
      )
    )
  ) %>%
  # Title is added BEFORE the layers control, and both sit in the same
  # corner, so the title stacks above the "Distribution data source" layer toggle.
  addControl(
    html = by_source_title_html,
    position = "topleft"
  ) %>%
  addLayersControl(
    baseGroups = c("Occurrence records", "GRIN-Global distribution data",
                   "WCFP distribution data", "Confirmed by ≥2 data sources"),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  # Base groups, not overlays: the layers are mutually exclusive, so exactly
  # one source is ever drawn and they cannot be stacked over one another. The
  # combined >=2-source layer is the headline result, so it is the one selected
  # on open; hideGroup() removes the other three, which leaves the radio button
  # for the combined layer as the checked one.
  hideGroup(c(LBL_OCC, LBL_GRIN, LBL_WCFP)) %>%
  # Added after the layers control so it sits below it, keeping the title
  # and the "Distribution data source" toggle adjacent.
  addControl(
    html = click_hint_country_bs,
    position = "topleft"
  ) %>%
  addControl(
    html = by_source_legend_html,
    position = "bottomright"
  ) %>%
  addControl(
    html = map_source_strip_country,
    position = "bottomright",
    className = "wcfp-mapsource"
  ) %>%
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(maxWidth = 260, metric = TRUE, imperial = FALSE, updateWhenIdle = TRUE)
  ) %>%
  htmlwidgets::onRender(by_source_onrender_js) %>%
  htmlwidgets::onRender(wcfp_popup_js, data = by_source_payload)

by_source_html_path <- file.path(combined_output_dir, "WCFP_taxa_distribution_by_source_map_interactive.html")
by_source_lib_tmp <- file.path(combined_output_dir, "lib_tmp_by_source")
save_widget_selfcontained(
  widget = by_source_leaflet,
  file   = by_source_html_path,
  libdir = by_source_lib_tmp
)
message("Interactive by-source leaflet map saved to: ", by_source_html_path)




## ============================================================================ ##
## TDWG Level 4 versions of both interactive maps above: same 3 data
## sources, same >=2-source combined logic, but keyed on TDWG level 4
## areas instead of countries.
##
## Data-resolution caveat: only occurrence records have genuine L4
## granularity (each point spatially joins to exactly one L4 polygon).
## GRIN and WCFP do not natively go finer than country/L3 level:
##  - GRIN: matched to one specific L4 area via its own `state` field
##    (exact name match against Level_4_Na, within the correct ISO
##    country) where possible; records with no usable state match are
##    smeared across every L4 area in that country instead, per explicit
##    user instruction.
##  - WCFP/WCVP (`wcvp_dist_tdwg3`): fundamentally a TDWG level 3
##    dataset, so every record is smeared across all L4 children of its
##    L3 area - a direct join on area_code_l3 == Level3_cod, no name
##    matching or guessing needed, since WCVP already stores the level-3
##    code itself.
## Author: Sarah Gora
## ============================================================================ ##

# Antarctica (TDWG level 3 code "ANT", level 4 "ANT-OO") is dropped from both
# layers, to match the country-level maps above - `country_polys` there is
# already filtered with admin != "Antarctica", so without this the L4 maps
# would be the only ones still showing it. It also matters for the initial
# view: keeping it pushes the data bounding box down to -90 latitude, which
# forces both L4 maps to open zoomed far out. Sub-antarctic islands (South
# Georgia, Kerguelen, Heard-McDonald etc.) are level 3 areas in their own
# right and are deliberately NOT removed, again matching the country maps.
level4_sf <- st_read(file.path(inputs_dir, "level4.geojson"), quiet = TRUE) %>%
  filter(Level3_cod != "ANT")
level4_lookup_full <- level4_sf %>%
  st_drop_geometry() %>%
  transmute(ISO_Code, Level4_cod, Level_4_Na = str_squish(Level_4_Na), Level3_cod)

# ENCODING=LATIN1 - see the level3_shp read in Method 3 above for why UTF-8
# is wrong here.
level3_sf <- st_read(file.path(inputs_dir, "level3/level3.shp"), quiet = TRUE,
                      options = "ENCODING=LATIN1") %>%
  filter(LEVEL3_COD != "ANT")

# ---------------------------------------- #
#   L4 Source 1: occurrences (true L4      #
#   resolution via direct spatial join)    #
# ---------------------------------------- #
# level4.geojson does not pass sf's spherical-geometry (s2) validity check:
# three features - FIN-OO (Finland), MAG-OO (Magadan), NUN-OO (Nunavut) -
# have crossing edges, which makes st_join() abort with
# "Loop 24 edge 0 crosses loop 31 edge 46". Neither st_make_valid() nor
# s2::s2_rebuild(check = FALSE) repairs them (both leave 2-3 features still
# s2-invalid and the join still fails), so this one join is run in planar
# mode, which GEOS handles without complaint. s2 is restored immediately
# afterwards so no other section is affected.
#
# Planar is safe for this particular operation: it is a point-in-polygon
# test against densely-vertexed L4 polygons, where the great-circle vs
# straight-line edge difference is negligible. The four features whose
# bounding box spans the antimeridian (Antarctica, Fiji, Magadan, Aleutian
# Is.) are multipolygons with separate parts either side of +/-180, not
# single polygons wrapping across it, so they do not degenerate into
# world-spanning rectangles under a planar reading - confirmed by checking
# that none of them appear among the areas involved in multi-area point
# assignments.
#
# Note the join legitimately assigns ~6% of points to more than one L4 area
# (shared borders in the TDWG L4 boundaries, e.g. Ukraine, Yunnan,
# Colombia). That is a property of the source data, not of the planar
# fallback, and is unchanged by it.
.s2_was_on <- sf_use_s2()
suppressMessages(sf_use_s2(FALSE))
occ_l4_joined <- suppressWarnings(
  st_join(occ_country_sf, level4_sf, join = st_within, left = FALSE)
)
suppressMessages(sf_use_s2(.s2_was_on))

src_occurrences_l4 <- occ_l4_joined %>%
  st_drop_geometry() %>%
  filter(!is.na(Level4_cod), !is.na(wcfp_name_match)) %>%
  distinct(area_l4 = Level4_cod, area_l4_name = Level_4_Na, area_l3 = Level3_cod,
           taxon = wcfp_name_match, data_source) %>%
  group_by(area_l4, area_l4_name, area_l3, taxon) %>%
  summarise(occurrences_data_source = paste(sort(unique(na.omit(data_source))), collapse = ", "),
            .groups = "drop")

# ---------------------------------------- #
#   L4 Source 2: GRIN (state-matched,      #
#   else smeared across the country)       #
# ---------------------------------------- #
# Rebuilds the same 3-tier taxon-name match as Source 2 above, but keeps
# `state`/`country_code` through the pipeline - the country-level
# src_grin collapses to distinct(country, taxon) too early to reuse here.
grin_with_country_rowid <- grin_with_country %>%
  mutate(.grin_row_id = row_number())

grin_direct_l4 <- grin_with_country_rowid %>%
  filter(name %in% wcfp_plantlist$wcfp_name_match) %>%
  transmute(.grin_row_id, country_code, state, taxon = str_squish(name), grin_status)

grin_via_genus_species_l4 <- grin_with_country_rowid %>%
  filter(!(name %in% wcfp_plantlist$wcfp_name_match)) %>%
  inner_join(taxon_name_lookup, by = c("name" = "taxon_name")) %>%
  transmute(.grin_row_id, country_code, state, taxon = wcfp_name_match, grin_status)

grin_via_hyphen_alias_l4 <- grin_with_country_rowid %>%
  filter(!(name %in% wcfp_plantlist$wcfp_name_match), name %in% names(grin_hyphen_alias)) %>%
  transmute(.grin_row_id, country_code, state, taxon = recode(name, !!!grin_hyphen_alias), grin_status) %>%
  filter(taxon %in% wcfp_plantlist$wcfp_name_match)

grin_taxon_matched <- bind_rows(grin_direct_l4, grin_via_genus_species_l4, grin_via_hyphen_alias_l4)

# ISO3 (GRIN's country_code) -> ISO2 (level4.geojson's ISO_Code)
iso3_to_iso2 <- country_polys %>%
  st_drop_geometry() %>%
  transmute(
    iso3 = if_else(is.na(iso_a3) | iso_a3 == "-99", iso_a3_eh, iso_a3),
    iso2 = if_else(is.na(iso_a2) | iso_a2 == "-99", iso_a2_eh, iso_a2)
  ) %>%
  filter(!is.na(iso3), nzchar(iso3), !is.na(iso2), nzchar(iso2)) %>%
  distinct(iso3, .keep_all = TRUE)

grin_l4_base <- grin_taxon_matched %>%
  left_join(iso3_to_iso2, by = c("country_code" = "iso3")) %>%
  mutate(state_norm = str_squish(state))

# Tier 1: exact state-name match to one specific L4 area within the
# correct country.
grin_l4_state_match <- grin_l4_base %>%
  filter(!is.na(state_norm), nzchar(state_norm), !is.na(iso2)) %>%
  inner_join(level4_lookup_full, by = c("iso2" = "ISO_Code", "state_norm" = "Level_4_Na"))

# Tier 2 (fallback): no usable state match - smeared across every L4 area
# in that country, per explicit user instruction.
grin_l4_smeared <- grin_l4_base %>%
  filter(!is.na(iso2)) %>%
  anti_join(grin_l4_state_match, by = ".grin_row_id") %>%
  inner_join(level4_lookup_full, by = c("iso2" = "ISO_Code"), relationship = "many-to-many")

cat(n_distinct(grin_l4_state_match$.grin_row_id), "of",
    n_distinct(grin_l4_base$.grin_row_id[!is.na(grin_l4_base$iso2)]),
    "GRIN records matched to a specific TDWG level 4 area via `state`; the rest are",
    "smeared across every L4 area in their country.\n")

src_grin_l4 <- bind_rows(grin_l4_state_match, grin_l4_smeared) %>%
  distinct(area_l4 = Level4_cod, area_l4_name = Level_4_Na, area_l3 = Level3_cod,
           taxon, grin_status) %>%
  group_by(area_l4, area_l4_name, area_l3, taxon) %>%
  summarise(grin_status = paste(sort(unique(grin_status)), collapse = ", "), .groups = "drop")

# ---------------------------------------- #
#   L4 Source 3: WCFP (smeared across L4   #
#   children of each TDWG L3 area)         #
# ---------------------------------------- #
src_wcfp_l4 <- wcvp_dist_tdwg3 %>%
  filter(!is.na(area_code_l3), nzchar(area_code_l3)) %>%
  inner_join(level4_lookup_full, by = c("area_code_l3" = "Level3_cod"), relationship = "many-to-many") %>%
  left_join(wcfp_plantlist %>% distinct(WCFP_ID, wcfp_name_match), by = "WCFP_ID") %>%
  filter(!is.na(wcfp_name_match)) %>%
  distinct(area_l4 = Level4_cod, area_l4_name = Level_4_Na, area_l3 = area_code_l3,
           taxon = wcfp_name_match, occurrence_status) %>%
  group_by(area_l4, area_l4_name, area_l3, taxon) %>%
  summarise(occurrence_status = paste(sort(unique(occurrence_status)), collapse = ", "), .groups = "drop")

# ---------------------------------------- #
#   Combine the 3 sources at L4 level      #
# ---------------------------------------- #
area_l4_lookup <- bind_rows(
  src_occurrences_l4 %>% distinct(area_l4, area_l4_name, area_l3),
  src_grin_l4        %>% distinct(area_l4, area_l4_name, area_l3),
  src_wcfp_l4        %>% distinct(area_l4, area_l4_name, area_l3)
) %>%
  distinct(area_l4, .keep_all = TRUE)

combined_long_l4 <- bind_rows(
  src_occurrences_l4 %>% select(-area_l4_name, -area_l3) %>% mutate(data_source_occurrences = "Y"),
  src_grin_l4        %>% select(-area_l4_name, -area_l3) %>% mutate(data_source_GRIN = "Y") %>%
    rename(status_GRIN = grin_status),
  src_wcfp_l4        %>% select(-area_l4_name, -area_l3) %>% mutate(data_source_WCFP = "Y") %>%
    rename(status_WCFP = occurrence_status)
)

combined_l4 <- combined_long_l4 %>%
  group_by(area_l4, taxon) %>%
  summarise(
    data_source_occurrences  = if (any(!is.na(data_source_occurrences))) "Y" else NA_character_,
    data_source_GRIN         = if (any(!is.na(data_source_GRIN)))        "Y" else NA_character_,
    data_source_WCFP         = if (any(!is.na(data_source_WCFP)))        "Y" else NA_character_,
    occurrences_data_source  = paste(sort(unique(na.omit(occurrences_data_source))), collapse = ", "),
    status_GRIN              = paste(sort(unique(na.omit(status_GRIN))), collapse = ", "),
    status_WCFP              = paste(sort(unique(na.omit(status_WCFP))), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    occurrences_data_source = na_if(occurrences_data_source, ""),
    status_GRIN = na_if(status_GRIN, ""),
    status_WCFP = na_if(status_WCFP, ""),
    n_sources = (!is.na(data_source_occurrences)) + (!is.na(data_source_GRIN)) + (!is.na(data_source_WCFP))
  )

combined_confirmed_l4 <- combined_l4 %>%
  filter(n_sources >= 2) %>%
  select(-n_sources) %>%
  left_join(area_l4_lookup, by = "area_l4") %>%
  left_join(wcfp_taxon_lookup, by = c("taxon" = "wcfp_name_match")) %>%
  attach_redlist() %>%
  arrange(area_l4, taxon)

qualifying_l4_areas <- sort(unique(combined_confirmed_l4$area_l4))
n_l3_covered <- n_distinct(combined_confirmed_l4$area_l3)
n_l4_covered <- length(qualifying_l4_areas)
cat(n_l4_covered, "TDWG level 4 areas (across", n_l3_covered,
    "TDWG level 3 areas) have >=1 taxon confirmed by >=2 of the 3 data sources.\n")
cat(nrow(combined_confirmed_l4), "total (L4 area, taxon) rows across all qualifying L4 areas.\n")

# ---------------------------------------- #
#   Write one workbook, one sheet/L4 area  #
# ---------------------------------------- #
wb_l4 <- createWorkbook()
used_sheet_names_l4 <- character(0)

l4_area_summary <- combined_confirmed_l4 %>%
  distinct(area_l4, area_l4_name, area_l3) %>%
  left_join(combined_confirmed_l4 %>% count(area_l4, name = "n_distinct_taxa"), by = "area_l4") %>%
  arrange(area_l4_name)

addWorksheet(wb_l4, "Summary")
writeData(wb_l4, "Summary", l4_area_summary, headerStyle = createStyle(textDecoration = "bold"))
freezePane(wb_l4, "Summary", firstRow = TRUE)
setColWidths(wb_l4, "Summary", cols = 1:4, widths = c(14, 30, 10, 20))
used_sheet_names_l4 <- c(used_sheet_names_l4, "Summary")

for (l4 in qualifying_l4_areas) {
  area_name <- area_l4_lookup$area_l4_name[area_l4_lookup$area_l4 == l4][1]
  sheet_name <- make_sheet_name(coalesce(area_name, l4), used_sheet_names_l4)
  used_sheet_names_l4 <- c(used_sheet_names_l4, sheet_name)

  l4_data <- combined_confirmed_l4 %>%
    filter(area_l4 == l4) %>%
    transmute(
      taxa = coalesce(taxon_name_accepted, taxon),
      authority = taxon_authors_accepted,
      iucn_red_list_category_code = rl_code,
      iucn_red_list_category      = unname(RL_LABEL[rl_code]),
      data_source_occurrences, occurrences_data_source,
      data_source_GRIN, status_GRIN, data_source_WCFP, status_WCFP
    )

  addWorksheet(wb_l4, sheet_name)
  writeData(wb_l4, sheet_name, l4_data, headerStyle = createStyle(textDecoration = "bold"))
  freezePane(wb_l4, sheet_name, firstRow = TRUE)
  setColWidths(wb_l4, sheet_name, cols = 1:10,
               widths = c(38, 22, 16, 24, 22, 22, 14, 20, 14, 20))
}

l4_xlsx_path <- file.path(combined_output_dir, "WCFP_taxa-lists_by_L4area_2plus_sources.xlsx")
saveWorkbook(wb_l4, l4_xlsx_path, overwrite = TRUE)
message("Combined per-L4-area taxa workbook saved to: ", l4_xlsx_path)

# ---------------------------------------- #
#   Popups + map data (combined L4 map)    #
# ---------------------------------------- #
# Keyed on the L4 CODE (area_l4) rather than the name, since level 4 names
# are not unique across countries; the display name is resolved from the
# payload's areaNames lookup at popup-build time.
combined_l4_rows <- combined_confirmed_l4 %>%
  mutate(
    area_key  = area_l4,
    taxa      = coalesce(taxon_name_accepted, taxon),
    authority = coalesce(taxon_authors_accepted, ""),
    info      = vapply(seq_len(n()), function(i) build_data_source_label(
      data_source_occurrences[i], occurrences_data_source[i],
      data_source_GRIN[i], data_source_WCFP[i], status_WCFP[i], status_GRIN[i]
    ), character(1))
  )

l4_area_names <- setNames(
  paste0(area_l4_lookup$area_l4_name, " (L3: ", area_l4_lookup$area_l3, ")"),
  area_l4_lookup$area_l4
)

l4_payload <- list(
  areaNames = as.list(l4_area_names),
  hideLabels = list(COMBINED_LAYER),
  headers   = setNames(list("Distribution data source"), COMBINED_LAYER),
  layers    = setNames(list(build_layer_export(combined_l4_rows)), COMBINED_LAYER)
)

l4_richness <- combined_confirmed_l4 %>% count(area_l4, name = "richness")

map_data_l4 <- level4_sf %>%
  left_join(l4_richness, by = c("Level4_cod" = "area_l4")) %>%
  mutate(
    popup_html = ifelse(
      is.na(richness),
      paste0("<b>", escape_html(Level_4_Na), "</b><br>No taxa confirmed by ≥2 data sources"),
      popup_stub(COMBINED_LAYER, Level4_cod, Level_4_Na, TRUE)
    )
  )

l4_richness_pal <- colorNumeric(
  palette  = "YlGn",
  domain   = map_data_l4$richness,
  # White rather than the original grey90 - reads cleanly against the blue
  # ocean background. The legend no longer carries a "no data" swatch, so
  # this only has to work against the map itself.
  na.color = "#FFFFFF",
  reverse  = FALSE
)

l4_map_bounds <- zoom_in_bbox(st_bbox(map_data_l4))

l4_legend_breaks <- pretty(map_data_l4$richness, n = 6)
l4_legend_breaks <- sort(
  l4_legend_breaks[l4_legend_breaks >= min(map_data_l4$richness, na.rm = TRUE) &
                     l4_legend_breaks <= max(map_data_l4$richness, na.rm = TRUE)],
  decreasing = TRUE
)
l4_legend_html <- build_gradient_legend(l4_richness_pal, l4_legend_breaks)

# Same heading wording as the country map above, so the two read as a pair;
# only the subtitle counts differ (L4/L3 areas rather than countries).
l4_map_title_html <- sprintf(
  paste0(MAP_HEADING,
         "<div style='font-size:13px; color:#000; margin-top:3px;'>",
         "%s taxa across %s TDWG Level-4 botanical areas (%s Level-3 botanical countries)",
         "<br><span style='font-style:italic;'>Distribution of taxa confirmed in &ge;2 data sources</span></div>"),
  scales::comma(n_distinct(combined_confirmed_l4$taxon)),
  scales::comma(n_l4_covered),
  scales::comma(n_l3_covered)
)

# L4 boundary weight (1.3) is a little thicker than the country/continent
# baseline (1) used in the country-level maps above. The L3 outline layer
# (unfilled, non-interactive) gives visual context for how L4 areas group
# into their parent TDWG level 3 area.
combined_l4_leaflet <- leaflet(map_data_l4,
                               options = leafletOptions(zoomSnap = 0.05, zoomDelta = 0.5)) %>%
  # Dedicated pane for the L3 outlines, above Leaflet's overlayPane (z 400)
  # where the L4 polygons are drawn, but below markerPane (600) so labels,
  # tooltips and popups still sit on top. Relying on draw order alone would
  # not hold: the L4 layers use highlightOptions(bringToFront = TRUE), so
  # every hovered L4 area would be lifted above the L3 outlines and stay
  # there, progressively burying them.
  addMapPane("l3_outline", zIndex = 450) %>%
  fitBounds(
    lng1 = as.numeric(l4_map_bounds["xmin"]), lat1 = as.numeric(l4_map_bounds["ymin"]),
    lng2 = as.numeric(l4_map_bounds["xmax"]), lat2 = as.numeric(l4_map_bounds["ymax"])
  ) %>%
  addPolygons(
    fillColor    = ~l4_richness_pal(richness),
    fillOpacity  = 0.8,
    color        = "black",
    weight       = 0.5,
    label        = ~lapply(paste0("<b>", escape_html(Level_4_Na), "</b>: ",
                                  scales::comma(ifelse(is.na(richness), 0, richness)),
                                  " taxa"), htmltools::HTML),
    popup        = ~popup_html,
    popupOptions = WCFP_POPUP_OPTIONS,
    group        = LYR_L4,
    highlightOptions = highlightOptions(weight = 2.5, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data    = level3_sf,
    fill    = FALSE,
    color   = "black",
    # Heavier than the L4 boundaries (1.3) so the parent TDWG level 3 areas
    # read as the dominant division, but not so heavy that the L4 mesh inside
    # them is overwhelmed.
    weight  = 1.1,
    opacity = 1,
    group   = LYR_L3,
    options = pathOptions(interactive = FALSE, pane = "l3_outline")
  ) %>%
  # Both resolutions toggleable and stackable. Both start on, which is how
  # this map has always looked - the control just makes the L3 outlines
  # removable when the L4 mesh needs to be read on its own.
  addLayersControl(
    overlayGroups = c(LYR_L4, LYR_L3),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  addControl(html = l4_legend_html, position = "bottomright") %>%
  addControl(html = map_source_strip_l4, position = "bottomright",
             className = "wcfp-mapsource") %>%
  # Title added before the "Click map to view/download" hint below so that, with
  # both in the same corner, the title stacks above it.
  addControl(html = l4_map_title_html, position = "topleft") %>%
  addControl(
    html = click_hint_l4,
    position = "topleft"
  ) %>%
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(maxWidth = 260, metric = TRUE, imperial = FALSE, updateWhenIdle = TRUE)
  ) %>%
  htmlwidgets::onRender(map_base_onrender_js) %>%
  htmlwidgets::onRender(wcfp_popup_js, data = l4_payload)

l4_html_path <- file.path(combined_output_dir, "WCFP_taxa_distribution_L4area_map_interactive.html")
l4_lib_tmp <- file.path(combined_output_dir, "lib_tmp_l4")
save_widget_selfcontained(
  widget = combined_l4_leaflet,
  file   = l4_html_path,
  libdir = l4_lib_tmp
)
message("Interactive L4-level leaflet map saved to: ", l4_html_path)

# ---------------------------------------- #
#   By-source toggle map at L4 level       #
# ---------------------------------------- #
richness_occurrences_l4 <- src_occurrences_l4 %>% count(area_l4, name = "richness")
richness_grin_l4        <- src_grin_l4        %>% count(area_l4, name = "richness")
richness_wcfp_l4        <- src_wcfp_l4        %>% count(area_l4, name = "richness")

map_data_occurrences_l4 <- level4_sf %>% left_join(richness_occurrences_l4, by = c("Level4_cod" = "area_l4"))
map_data_grin_l4        <- level4_sf %>% left_join(richness_grin_l4,        by = c("Level4_cod" = "area_l4"))
map_data_wcfp_l4        <- level4_sf %>% left_join(richness_wcfp_l4,        by = c("Level4_cod" = "area_l4"))

# Fourth layer: the ">=2 of 3 sources" combined result at L4, reusing the same
# `l4_richness` (derived from `combined_confirmed_l4`) that drives the L4
# richness map above.
map_data_combined_l4    <- level4_sf %>% left_join(l4_richness,             by = c("Level4_cod" = "area_l4"))

# Shared palette/domain across all 4 layers, combined layer included so its
# values can never fall outside the palette range.
by_source_l4_max_richness <- max(
  c(map_data_occurrences_l4$richness, map_data_grin_l4$richness, map_data_wcfp_l4$richness,
    map_data_combined_l4$richness),
  na.rm = TRUE
)
by_source_l4_pal <- colorNumeric(
  palette  = "YlGn",
  domain   = c(0, by_source_l4_max_richness),
  # White rather than the original grey90 - reads cleanly against the blue
  # ocean background. The legend no longer carries a "no data" swatch, so
  # this only has to work against the map itself.
  na.color = "#FFFFFF",
  reverse  = FALSE
)

by_source_l4_legend_breaks <- pretty(c(0, by_source_l4_max_richness), n = 6)
by_source_l4_legend_breaks <- sort(
  by_source_l4_legend_breaks[by_source_l4_legend_breaks <= by_source_l4_max_richness],
  decreasing = TRUE
)
by_source_l4_legend_html <- build_gradient_legend(by_source_l4_pal, by_source_l4_legend_breaks)

n_l3_all <- n_distinct(c(src_occurrences_l4$area_l3, src_grin_l4$area_l3, src_wcfp_l4$area_l3))
n_l4_all <- n_distinct(c(src_occurrences_l4$area_l4, src_grin_l4$area_l4, src_wcfp_l4$area_l4))

# ---------------------------------------- #
#   Per-layer popups + downloadable tables  #
# ---------------------------------------- #
# Same structure as the country by-source map, keyed on the L4 code. Note
# these are the largest payloads in the workflow: WCFP is smeared across
# every L4 child of its L3 area, and ~66% of GRIN records across every L4
# area in their country, so the GRIN and WCFP layers hold considerably
# more rows than the combined >=2-source layer does.
bs4_rows_occ  <- src_occurrences_l4 %>%
  transmute(area_key = area_l4, taxon, info = occurrences_data_source) %>% attach_taxon_names()
bs4_rows_grin <- src_grin_l4 %>%
  transmute(area_key = area_l4, taxon, info = grin_status) %>% attach_taxon_names()
bs4_rows_wcfp <- src_wcfp_l4 %>%
  transmute(area_key = area_l4, taxon, info = occurrence_status) %>% attach_taxon_names()
bs4_rows_comb <- combined_l4_rows %>% select(area_key, taxa, authority, info)

# Display names for every L4 area any source touches, not just the ones
# reaching the >=2-source threshold (area_l4_lookup already spans all three).
bs4_area_names <- setNames(
  paste0(area_l4_lookup$area_l4_name, " (L3: ", area_l4_lookup$area_l3, ")"),
  area_l4_lookup$area_l4
)

by_source_l4_payload <- list(
  areaNames = as.list(bs4_area_names),
  hideLabels = list(COMBINED_LAYER),
  headers   = setNames(as.list(c("Occurrence source", "Detail",
                                 "Native/introduced status", "Distribution data source")),
                       bs_layer_names),
  layers    = setNames(list(build_layer_export(bs4_rows_occ),
                            build_layer_export(bs4_rows_grin),
                            build_layer_export(bs4_rows_wcfp),
                            build_layer_export(bs4_rows_comb)),
                       bs_layer_names)
)

# Same four-line subtitle as the country by-source map, one line per
# toggleable layer. Counts are per level 4 area here rather than per country.
bs4_taxa_counts <- c(
  n_distinct(src_occurrences_l4$taxon),
  n_distinct(src_grin_l4$taxon),
  n_distinct(src_wcfp_l4$taxon),
  n_distinct(combined_confirmed_l4$taxon)
)
bs4_area_counts <- c(
  n_distinct(src_occurrences_l4$area_l4),
  n_distinct(src_grin_l4$area_l4),
  n_distinct(src_wcfp_l4$area_l4),
  n_distinct(combined_confirmed_l4$area_l4)
)

by_source_l4_title_html <- paste0(
  MAP_HEADING,
  "<div style='font-size:13px; color:#000; margin-top:3px;'>",
  paste0(bs_layer_names, ": ", scales::comma(bs4_taxa_counts), " taxa across ",
         scales::comma(bs4_area_counts), " TDWG Level-4 botanical areas", collapse = "<br>"),
  "</div>"
)

by_source_l4_leaflet <- leaflet(options = leafletOptions(zoomSnap = 0.05, zoomDelta = 0.5)) %>%
  # See the combined L4 map above for why the L3 outlines need their own pane
  # rather than just being added last.
  addMapPane("l3_outline", zIndex = 450) %>%
  fitBounds(
    lng1 = as.numeric(l4_map_bounds["xmin"]), lat1 = as.numeric(l4_map_bounds["ymin"]),
    lng2 = as.numeric(l4_map_bounds["xmax"]), lat2 = as.numeric(l4_map_bounds["ymax"])
  ) %>%
  addPolygons(
    data        = map_data_occurrences_l4,
    fillColor   = ~by_source_l4_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 0.5,
    label       = ~lapply(paste0("<b>", escape_html(Level_4_Na), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_OCC, Level4_cod, Level_4_Na, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_OCC,
    highlightOptions = highlightOptions(weight = 2.5, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_grin_l4,
    fillColor   = ~by_source_l4_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 0.5,
    label       = ~lapply(paste0("<b>", escape_html(Level_4_Na), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_GRIN, Level4_cod, Level_4_Na, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_GRIN,
    highlightOptions = highlightOptions(weight = 2.5, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_wcfp_l4,
    fillColor   = ~by_source_l4_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 0.5,
    label       = ~lapply(paste0("<b>", escape_html(Level_4_Na), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_WCFP, Level4_cod, Level_4_Na, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_WCFP,
    highlightOptions = highlightOptions(weight = 2.5, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_combined_l4,
    fillColor   = ~by_source_l4_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 0.5,
    label       = ~lapply(paste0("<b>", escape_html(Level_4_Na), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(COMBINED_LAYER, Level4_cod, Level_4_Na, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = COMBINED_LAYER,
    highlightOptions = highlightOptions(weight = 2.5, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data    = level3_sf,
    fill    = FALSE,
    color   = "black",
    # Heavier than the L4 boundaries (1.3) so the parent TDWG level 3 areas
    # read as the dominant division, but not so heavy that the L4 mesh inside
    # them is overwhelmed.
    weight  = 1.1,
    opacity = 1,
    options = pathOptions(interactive = FALSE, pane = "l3_outline")
  ) %>%
  # Title is added BEFORE the layers control, and both sit in the same
  # corner, so the title stacks above the "Distribution data source" layer toggle.
  addControl(html = by_source_l4_title_html, position = "topleft") %>%
  addLayersControl(
    baseGroups = c("Occurrence records", "GRIN-Global distribution data",
                   "WCFP distribution data", "Confirmed by ≥2 data sources"),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  # Base groups, not overlays: the layers are mutually exclusive, so exactly
  # one source is ever drawn and they cannot be stacked over one another. The
  # combined >=2-source layer is the headline result, so it is the one selected
  # on open; hideGroup() removes the other three, which leaves the radio button
  # for the combined layer as the checked one.
  hideGroup(c(LBL_OCC, LBL_GRIN, LBL_WCFP)) %>%
  # Added after the layers control so it sits below it, keeping the title
  # and the "Distribution data source" toggle adjacent.
  addControl(html = click_hint_l4_bs, position = "topleft") %>%
  addControl(html = by_source_l4_legend_html, position = "bottomright") %>%
  addControl(html = map_source_strip_l4, position = "bottomright",
             className = "wcfp-mapsource") %>%
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(maxWidth = 260, metric = TRUE, imperial = FALSE, updateWhenIdle = TRUE)
  ) %>%
  htmlwidgets::onRender(by_source_onrender_js) %>%
  htmlwidgets::onRender(wcfp_popup_js, data = by_source_l4_payload)

by_source_l4_html_path <- file.path(combined_output_dir, "WCFP_taxa_distribution_by_source_L4area_map_interactive.html")
by_source_l4_lib_tmp <- file.path(combined_output_dir, "lib_tmp_by_source_l4")
save_widget_selfcontained(
  widget = by_source_l4_leaflet,
  file   = by_source_l4_html_path,
  libdir = by_source_l4_lib_tmp
)
message("Interactive by-source L4-level leaflet map saved to: ", by_source_l4_html_path)


## ============================================================================ ##
## TDWG LEVEL 3 CHOROPLETH  (one polygon per level 3 area, shaded by richness)
##
## This is the resolution the WCFP 2026 paper itself reports at, so this map -
## not the level 4 one - is the figure directly comparable to the paper's own
## richness outputs.
##
## IMPORTANT - this is NOT the level 4 map re-shaded. The >=2-source test is
## applied AFTER aggregating to level 3, not before. Those give different
## answers: if occurrences record a taxon in one L4 area and GRIN records it in
## a different L4 area of the SAME parent L3 area, neither L4 area reaches the
## threshold, but at level 3 both sources are saying the same thing - "this
## taxon is in this L3 area" - and it does. Aggregating combined_confirmed_l4
## up to L3 would silently discard those agreements. The count of what the
## native test recovers is printed below.
##
## How each source reaches level 3:
##  - Occurrences: rolled up from src_occurrences_l4 via its area_l3 column.
##    Exactly equivalent to joining the points straight to level3.shp, because
##    L4 areas partition their parent L3 area and the two files carry an
##    identical set of 369 L3 codes (verified: no code in either that is not in
##    the other). Reusing the L4 join also avoids a second expensive spatial
##    join - and avoids a worse s2 problem, since level3.shp has SEVEN
##    s2-invalid features (AND, KRA, MAG, NFL, NGA, NOR, SCI) against
##    level4.geojson's three.
##  - GRIN: rolled up from src_grin_l4 the same way. Note the smearing this
##    inherits is much weaker at L3 than at L4: a GRIN record with no usable
##    `state` is smeared across every L4 area in its country, but collapsing
##    those to their parents just gives the L3 areas that country spans.
##  - WCFP/WCVP: built NATIVELY from wcvp_dist_tdwg3, not rolled up. WCVP is a
##    level 3 dataset to begin with, so at this resolution it involves no
##    smearing and no name matching at all - the cleanest of the three.
## Author: Sarah Gora
## ============================================================================ ##

# Re-collapses an already comma-joined field when several L4 rows fold into one
# L3 row. Splitting first matters: pasting the strings directly would produce
# duplicates like "GBIF, GBIF, GBIF" for a taxon present in three L4 areas.
collapse_csv <- function(x) {
  v <- unlist(strsplit(as.character(na.omit(x)), ",", fixed = TRUE))
  v <- sort(unique(trimws(v)))
  v <- v[nzchar(v)]
  if (!length(v)) NA_character_ else paste(v, collapse = ", ")
}

# Authoritative code -> name lookup for all 369 L3 areas, taken from the
# polygons themselves so every area the map can draw has a display name,
# including those with no data.
area_l3_lookup <- level3_sf %>%
  st_drop_geometry() %>%
  transmute(area_l3 = LEVEL3_COD, area_l3_name = str_squish(LEVEL3_NAM))

# ---------------------------------------- #
#   L3 Source 1: occurrences               #
# ---------------------------------------- #
src_occurrences_l3 <- src_occurrences_l4 %>%
  group_by(area_l3, taxon) %>%
  summarise(occurrences_data_source = collapse_csv(occurrences_data_source),
            .groups = "drop")

# ---------------------------------------- #
#   L3 Source 2: GRIN                      #
# ---------------------------------------- #
src_grin_l3 <- src_grin_l4 %>%
  group_by(area_l3, taxon) %>%
  summarise(grin_status = collapse_csv(grin_status), .groups = "drop")

# ---------------------------------------- #
#   L3 Source 3: WCFP/WCVP (native)        #
# ---------------------------------------- #
# Restricted to codes that level3.shp can actually draw, which also removes
# Antarctica (level3_sf is already filtered). Anything dropped here is
# reported rather than lost silently.
wcvp_l3_codes_all <- wcvp_dist_tdwg3 %>%
  filter(!is.na(area_code_l3), nzchar(area_code_l3)) %>%
  distinct(area_code_l3) %>%
  pull(area_code_l3)
wcvp_l3_codes_unmapped <- setdiff(wcvp_l3_codes_all, level3_sf$LEVEL3_COD)
if (length(wcvp_l3_codes_unmapped) > 0) {
  cat("WCVP level 3 codes with no polygon in level3.shp (dropped from the L3 map): ",
      paste(sort(wcvp_l3_codes_unmapped), collapse = ", "), "\n", sep = "")
}

src_wcfp_l3 <- wcvp_dist_tdwg3 %>%
  filter(!is.na(area_code_l3), nzchar(area_code_l3),
         area_code_l3 %in% level3_sf$LEVEL3_COD) %>%
  left_join(wcfp_plantlist %>% distinct(WCFP_ID, wcfp_name_match), by = "WCFP_ID") %>%
  filter(!is.na(wcfp_name_match)) %>%
  distinct(area_l3 = area_code_l3, taxon = wcfp_name_match, occurrence_status) %>%
  group_by(area_l3, taxon) %>%
  summarise(occurrence_status = collapse_csv(occurrence_status), .groups = "drop")

# ---------------------------------------- #
#   Combine the 3 sources at L3 level      #
# ---------------------------------------- #
combined_long_l3 <- bind_rows(
  src_occurrences_l3 %>% mutate(data_source_occurrences = "Y"),
  src_grin_l3        %>% mutate(data_source_GRIN = "Y") %>% rename(status_GRIN = grin_status),
  src_wcfp_l3        %>% mutate(data_source_WCFP = "Y") %>% rename(status_WCFP = occurrence_status)
)

combined_l3 <- combined_long_l3 %>%
  group_by(area_l3, taxon) %>%
  summarise(
    data_source_occurrences  = if (any(!is.na(data_source_occurrences))) "Y" else NA_character_,
    data_source_GRIN         = if (any(!is.na(data_source_GRIN)))        "Y" else NA_character_,
    data_source_WCFP         = if (any(!is.na(data_source_WCFP)))        "Y" else NA_character_,
    occurrences_data_source  = paste(sort(unique(na.omit(occurrences_data_source))), collapse = ", "),
    status_GRIN              = paste(sort(unique(na.omit(status_GRIN))), collapse = ", "),
    status_WCFP              = paste(sort(unique(na.omit(status_WCFP))), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    occurrences_data_source = na_if(occurrences_data_source, ""),
    status_GRIN = na_if(status_GRIN, ""),
    status_WCFP = na_if(status_WCFP, ""),
    n_sources = (!is.na(data_source_occurrences)) + (!is.na(data_source_GRIN)) + (!is.na(data_source_WCFP))
  )

combined_confirmed_l3 <- combined_l3 %>%
  filter(n_sources >= 2) %>%
  select(-n_sources) %>%
  left_join(area_l3_lookup, by = "area_l3") %>%
  left_join(wcfp_taxon_lookup, by = c("taxon" = "wcfp_name_match")) %>%
  attach_redlist() %>%
  arrange(area_l3, taxon)

qualifying_l3_areas <- sort(unique(combined_confirmed_l3$area_l3))
n_l3_areas_covered  <- length(qualifying_l3_areas)
cat(n_l3_areas_covered, "of", nrow(level3_sf),
    "TDWG level 3 areas have >=1 taxon confirmed by >=2 of the 3 data sources.\n")
cat(nrow(combined_confirmed_l3), "total (L3 area, taxon) rows.\n")

# What the native L3 threshold recovers over simply re-shading the L4 result.
# If this is 0 the two approaches happen to agree; it is reported either way so
# the choice made above is auditable rather than asserted.
l3_rows_from_l4_rollup <- combined_confirmed_l4 %>% distinct(area_l3, taxon)
l3_rows_native         <- combined_confirmed_l3 %>% distinct(area_l3, taxon)
cat("(L3 area, taxon) pairs from rolling the L4 result up  :", nrow(l3_rows_from_l4_rollup), "\n")
cat("(L3 area, taxon) pairs from testing >=2 sources at L3 :", nrow(l3_rows_native), "\n")
cat("  recovered by applying the threshold at L3           :",
    nrow(dplyr::anti_join(l3_rows_native, l3_rows_from_l4_rollup, by = c("area_l3", "taxon"))), "\n")

# ---------------------------------------- #
#   Write one workbook, one sheet/L3 area  #
# ---------------------------------------- #
wb_l3 <- createWorkbook()
used_sheet_names_l3 <- character(0)

l3_area_summary <- combined_confirmed_l3 %>%
  distinct(area_l3, area_l3_name) %>%
  left_join(combined_confirmed_l3 %>% count(area_l3, name = "n_distinct_taxa"), by = "area_l3") %>%
  arrange(area_l3_name)

addWorksheet(wb_l3, "Summary")
writeData(wb_l3, "Summary", l3_area_summary, headerStyle = createStyle(textDecoration = "bold"))
freezePane(wb_l3, "Summary", firstRow = TRUE)
setColWidths(wb_l3, "Summary", cols = 1:3, widths = c(14, 34, 20))
used_sheet_names_l3 <- c(used_sheet_names_l3, "Summary")

for (l3 in qualifying_l3_areas) {
  area_name <- area_l3_lookup$area_l3_name[area_l3_lookup$area_l3 == l3][1]
  sheet_name <- make_sheet_name(coalesce(area_name, l3), used_sheet_names_l3)
  used_sheet_names_l3 <- c(used_sheet_names_l3, sheet_name)

  l3_data <- combined_confirmed_l3 %>%
    filter(area_l3 == l3) %>%
    transmute(
      taxa = coalesce(taxon_name_accepted, taxon),
      authority = taxon_authors_accepted,
      iucn_red_list_category_code = rl_code,
      iucn_red_list_category      = unname(RL_LABEL[rl_code]),
      data_source_occurrences, occurrences_data_source,
      data_source_GRIN, status_GRIN, data_source_WCFP, status_WCFP
    )

  addWorksheet(wb_l3, sheet_name)
  writeData(wb_l3, sheet_name, l3_data, headerStyle = createStyle(textDecoration = "bold"))
  freezePane(wb_l3, sheet_name, firstRow = TRUE)
  setColWidths(wb_l3, sheet_name, cols = 1:10,
               widths = c(38, 22, 16, 24, 22, 22, 14, 20, 14, 20))
}

l3_xlsx_path <- file.path(combined_output_dir, "WCFP_taxa-lists_by_L3area_2plus_sources.xlsx")
saveWorkbook(wb_l3, l3_xlsx_path, overwrite = TRUE)
message("Combined per-L3-area taxa workbook saved to: ", l3_xlsx_path)

# ---------------------------------------- #
#   Popups + map data (combined L3 map)    #
# ---------------------------------------- #
# Keyed on the L3 CODE for consistency with the L4 map, even though level 3
# names happen to be unique - the popup JS resolves the display name from the
# payload's areaNames lookup either way.
combined_l3_rows <- combined_confirmed_l3 %>%
  mutate(
    area_key  = area_l3,
    taxa      = coalesce(taxon_name_accepted, taxon),
    authority = coalesce(taxon_authors_accepted, ""),
    info      = vapply(seq_len(n()), function(i) build_data_source_label(
      data_source_occurrences[i], occurrences_data_source[i],
      data_source_GRIN[i], data_source_WCFP[i], status_WCFP[i], status_GRIN[i]
    ), character(1))
  )

l3_payload <- list(
  areaNames = as.list(setNames(area_l3_lookup$area_l3_name, area_l3_lookup$area_l3)),
  hideLabels = list(COMBINED_LAYER),
  headers   = setNames(list("Distribution data source"), COMBINED_LAYER),
  layers    = setNames(list(build_layer_export(combined_l3_rows)), COMBINED_LAYER)
)

l3_richness <- combined_confirmed_l3 %>% count(area_l3, name = "richness")

map_data_l3 <- level3_sf %>%
  left_join(l3_richness, by = c("LEVEL3_COD" = "area_l3")) %>%
  mutate(
    popup_html = ifelse(
      is.na(richness),
      paste0("<b>", escape_html(LEVEL3_NAM), "</b><br>No taxa confirmed by ≥2 data sources"),
      popup_stub(COMBINED_LAYER, LEVEL3_COD, LEVEL3_NAM, TRUE)
    )
  )

l3_richness_pal <- colorNumeric(
  palette  = "YlGn",
  domain   = map_data_l3$richness,
  na.color = "#FFFFFF",
  reverse  = FALSE
)

l3_map_bounds <- zoom_in_bbox(st_bbox(map_data_l3))

l3_legend_breaks <- pretty(map_data_l3$richness, n = 6)
l3_legend_breaks <- sort(
  l3_legend_breaks[l3_legend_breaks >= min(map_data_l3$richness, na.rm = TRUE) &
                     l3_legend_breaks <= max(map_data_l3$richness, na.rm = TRUE)],
  decreasing = TRUE
)
l3_legend_html <- build_gradient_legend(l3_richness_pal, l3_legend_breaks)

l3_map_title_html <- sprintf(
  paste0(MAP_HEADING,
         "<div style='font-size:13px; color:#000; margin-top:3px;'>",
         "%s taxa across %s TDWG Level-3 botanical countries",
         "<br><span style='font-style:italic;'>Distribution of taxa confirmed in &ge;2 data sources</span></div>"),
  scales::comma(n_distinct(combined_confirmed_l3$taxon)),
  scales::comma(n_l3_areas_covered)
)

click_hint_l3         <- build_click_hint(
  "Click map to view/download food plants taxa list per TDWG Level-3 botanical country")
map_source_strip_l3   <- build_map_source_strip("TDWG WGSRPD Level-3 botanical countries")

# Boundary weight 1.3 matches the L4 polygons on the level 4 maps rather than
# the heavier 2.0 used there for the L3 outlines. On those maps the 2.0 marks
# L3 as the parent division over a finer mesh; here L3 IS the shaded unit, so
# the same weight would just read as heavy.
combined_l3_leaflet <- leaflet(map_data_l3,
                               options = leafletOptions(zoomSnap = 0.05, zoomDelta = 0.5)) %>%
  # Pane for the L4 subdivision outlines, mirroring the l3_outline pane on the
  # level 4 maps: above the overlayPane (z 400) that holds the shaded L3
  # polygons, below markerPane (600) so tooltips and popups stay on top.
  # Needed for the same reason as there - the L3 layer uses
  # highlightOptions(bringToFront = TRUE), so without a dedicated pane every
  # hovered area would be lifted over the outlines and bury them.
  addMapPane("l4_outline", zIndex = 450) %>%
  fitBounds(
    lng1 = as.numeric(l3_map_bounds["xmin"]), lat1 = as.numeric(l3_map_bounds["ymin"]),
    lng2 = as.numeric(l3_map_bounds["xmax"]), lat2 = as.numeric(l3_map_bounds["ymax"])
  ) %>%
  addPolygons(
    fillColor    = ~l3_richness_pal(richness),
    fillOpacity  = 0.8,
    color        = "black",
    weight       = 1.3,
    label        = ~lapply(paste0("<b>", escape_html(LEVEL3_NAM), "</b>: ",
                                  scales::comma(ifelse(is.na(richness), 0, richness)),
                                  " taxa"), htmltools::HTML),
    popup        = ~popup_html,
    popupOptions = WCFP_POPUP_OPTIONS,
    group        = LYR_L3,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  # L4 subdivisions, boundary only. Lighter and thinner than the L3 polygon
  # edges (1.3, black) so the finer mesh reads as subdivision rather than
  # competing with the shaded unit. Non-interactive, so it never intercepts a
  # click meant for the L3 polygon underneath.
  addPolygons(
    data    = level4_sf,
    fill    = FALSE,
    color   = "#444444",
    weight  = 0.6,
    opacity = 0.9,
    group   = LYR_L4,
    options = pathOptions(interactive = FALSE, pane = "l4_outline")
  ) %>%
  # L3 (the shaded layer) starts on, L4 starts off - this map has always shown
  # L3 alone, so the overlay is additive rather than a change of default.
  addLayersControl(
    overlayGroups = c(LYR_L3, LYR_L4),
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  hideGroup(LYR_L4) %>%
  addControl(html = l3_legend_html, position = "bottomright") %>%
  addControl(html = map_source_strip_l3, position = "bottomright",
             className = "wcfp-mapsource") %>%
  addControl(html = l3_map_title_html, position = "topleft") %>%
  addControl(html = click_hint_l3, position = "topleft") %>%
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(maxWidth = 260, metric = TRUE, imperial = FALSE, updateWhenIdle = TRUE)
  ) %>%
  htmlwidgets::onRender(map_base_onrender_js) %>%
  htmlwidgets::onRender(wcfp_popup_js, data = l3_payload)

l3_html_path <- file.path(combined_output_dir, "WCFP_taxa_distribution_L3area_map_interactive.html")
l3_lib_tmp <- file.path(combined_output_dir, "lib_tmp_l3")
save_widget_selfcontained(
  widget = combined_l3_leaflet,
  file   = l3_html_path,
  libdir = l3_lib_tmp
)
message("Interactive L3-level leaflet map saved to: ", l3_html_path)

# ---------------------------------------- #
#   By-source toggle map at L3 level       #
# ---------------------------------------- #
# Same four toggleable layers as the country and L4 by-source maps, so all
# three resolutions offer the same comparison. The per-source L3 tables were
# already built above for the >=2-source test, so this reuses them directly.
richness_occurrences_l3 <- src_occurrences_l3 %>% count(area_l3, name = "richness")
richness_grin_l3        <- src_grin_l3        %>% count(area_l3, name = "richness")
richness_wcfp_l3        <- src_wcfp_l3        %>% count(area_l3, name = "richness")

map_data_occurrences_l3 <- level3_sf %>% left_join(richness_occurrences_l3, by = c("LEVEL3_COD" = "area_l3"))
map_data_grin_l3        <- level3_sf %>% left_join(richness_grin_l3,        by = c("LEVEL3_COD" = "area_l3"))
map_data_wcfp_l3        <- level3_sf %>% left_join(richness_wcfp_l3,        by = c("LEVEL3_COD" = "area_l3"))

# Fourth layer: the ">=2 of 3 sources" result, reusing the same `l3_richness`
# that drives the combined L3 map above.
map_data_combined_l3    <- level3_sf %>% left_join(l3_richness,             by = c("LEVEL3_COD" = "area_l3"))

# Shared palette/domain across all 4 layers, combined layer included so its
# values can never fall outside the palette range.
by_source_l3_max_richness <- max(
  c(map_data_occurrences_l3$richness, map_data_grin_l3$richness, map_data_wcfp_l3$richness,
    map_data_combined_l3$richness),
  na.rm = TRUE
)
by_source_l3_pal <- colorNumeric(
  palette  = "YlGn",
  domain   = c(0, by_source_l3_max_richness),
  na.color = "#FFFFFF",
  reverse  = FALSE
)

by_source_l3_legend_breaks <- pretty(c(0, by_source_l3_max_richness), n = 6)
by_source_l3_legend_breaks <- sort(
  by_source_l3_legend_breaks[by_source_l3_legend_breaks <= by_source_l3_max_richness],
  decreasing = TRUE
)
by_source_l3_legend_html <- build_gradient_legend(by_source_l3_pal, by_source_l3_legend_breaks)

# Per-layer popups. These payloads are much smaller than the L4 equivalents:
# WCFP is native at level 3 rather than smeared across L4 children, and GRIN's
# country-level fallback collapses to just the L3 areas that country spans.
bs3_rows_occ  <- src_occurrences_l3 %>%
  transmute(area_key = area_l3, taxon, info = occurrences_data_source) %>% attach_taxon_names()
bs3_rows_grin <- src_grin_l3 %>%
  transmute(area_key = area_l3, taxon, info = grin_status) %>% attach_taxon_names()
bs3_rows_wcfp <- src_wcfp_l3 %>%
  transmute(area_key = area_l3, taxon, info = occurrence_status) %>% attach_taxon_names()
bs3_rows_comb <- combined_l3_rows %>% select(area_key, taxa, authority, info)

by_source_l3_payload <- list(
  areaNames = as.list(setNames(area_l3_lookup$area_l3_name, area_l3_lookup$area_l3)),
  hideLabels = list(COMBINED_LAYER),
  headers   = setNames(as.list(c("Occurrence source", "Detail",
                                 "Native/introduced status", "Distribution data source")),
                       bs_layer_names),
  layers    = setNames(list(build_layer_export(bs3_rows_occ),
                            build_layer_export(bs3_rows_grin),
                            build_layer_export(bs3_rows_wcfp),
                            build_layer_export(bs3_rows_comb)),
                       bs_layer_names)
)

bs3_taxa_counts <- c(
  n_distinct(src_occurrences_l3$taxon),
  n_distinct(src_grin_l3$taxon),
  n_distinct(src_wcfp_l3$taxon),
  n_distinct(combined_confirmed_l3$taxon)
)
bs3_area_counts <- c(
  n_distinct(src_occurrences_l3$area_l3),
  n_distinct(src_grin_l3$area_l3),
  n_distinct(src_wcfp_l3$area_l3),
  n_distinct(combined_confirmed_l3$area_l3)
)

by_source_l3_title_html <- paste0(
  MAP_HEADING,
  "<div style='font-size:13px; color:#000; margin-top:3px;'>",
  paste0(bs_layer_names, ": ", scales::comma(bs3_taxa_counts), " taxa across ",
         scales::comma(bs3_area_counts), " TDWG Level-3 botanical countries", collapse = "<br>"),
  "</div>"
)

click_hint_l3_bs <- build_click_hint(paste0(
  "Select one distribution data source - layers are exclusive, not stacked<br>",
  "Click map to view/download food plants taxa list per TDWG Level-3 botanical country, for the selected data source"))

by_source_l3_leaflet <- leaflet(options = leafletOptions(zoomSnap = 0.05, zoomDelta = 0.5)) %>%
  fitBounds(
    lng1 = as.numeric(l3_map_bounds["xmin"]), lat1 = as.numeric(l3_map_bounds["ymin"]),
    lng2 = as.numeric(l3_map_bounds["xmax"]), lat2 = as.numeric(l3_map_bounds["ymax"])
  ) %>%
  addPolygons(
    data        = map_data_occurrences_l3,
    fillColor   = ~by_source_l3_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1.3,
    label       = ~lapply(paste0("<b>", escape_html(LEVEL3_NAM), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_OCC, LEVEL3_COD, LEVEL3_NAM, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_OCC,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_grin_l3,
    fillColor   = ~by_source_l3_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1.3,
    label       = ~lapply(paste0("<b>", escape_html(LEVEL3_NAM), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_GRIN, LEVEL3_COD, LEVEL3_NAM, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_GRIN,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_wcfp_l3,
    fillColor   = ~by_source_l3_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1.3,
    label       = ~lapply(paste0("<b>", escape_html(LEVEL3_NAM), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(LBL_WCFP, LEVEL3_COD, LEVEL3_NAM, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = LBL_WCFP,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  addPolygons(
    data        = map_data_combined_l3,
    fillColor   = ~by_source_l3_pal(richness),
    fillOpacity = 0.8,
    color       = "black",
    weight      = 1.3,
    label       = ~lapply(paste0("<b>", escape_html(LEVEL3_NAM), "</b>: ",
                                 scales::comma(ifelse(is.na(richness), 0, richness)),
                                 " taxa"), htmltools::HTML),
    popup       = ~popup_stub(COMBINED_LAYER, LEVEL3_COD, LEVEL3_NAM, !is.na(richness)),
    popupOptions = WCFP_POPUP_OPTIONS,
    group       = COMBINED_LAYER,
    highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                        bringToFront = TRUE)
  ) %>%
  # Title before the layers control, both in the same corner, so the title
  # stacks above the "Distribution data source" toggle.
  addControl(html = by_source_l3_title_html, position = "topleft") %>%
  addLayersControl(
    baseGroups = bs_layer_names,
    position = "topleft",
    options = layersControlOptions(collapsed = FALSE)
  ) %>%
  # Base groups, not overlays: the layers are mutually exclusive, so exactly
  # one source is ever drawn and they cannot be stacked over one another. The
  # combined >=2-source layer is the headline result, so it is the one selected
  # on open; hideGroup() removes the other three, which leaves the radio button
  # for the combined layer as the checked one.
  hideGroup(c(LBL_OCC, LBL_GRIN, LBL_WCFP)) %>%
  addControl(html = click_hint_l3_bs, position = "topleft") %>%
  addControl(html = by_source_l3_legend_html, position = "bottomright") %>%
  addControl(html = map_source_strip_l3, position = "bottomright",
             className = "wcfp-mapsource") %>%
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(maxWidth = 260, metric = TRUE, imperial = FALSE, updateWhenIdle = TRUE)
  ) %>%
  htmlwidgets::onRender(by_source_onrender_js) %>%
  htmlwidgets::onRender(wcfp_popup_js, data = by_source_l3_payload)

by_source_l3_html_path <- file.path(combined_output_dir, "WCFP_taxa_distribution_by_source_L3area_map_interactive.html")
by_source_l3_lib_tmp <- file.path(combined_output_dir, "lib_tmp_by_source_l3")
save_widget_selfcontained(
  widget = by_source_l3_leaflet,
  file   = by_source_l3_html_path,
  libdir = by_source_l3_lib_tmp
)
message("Interactive by-source L3-level leaflet map saved to: ", by_source_l3_html_path)

## ============================================================================ ##
## IUCN RED LIST TOGGLE MAPS AT LEVEL 3
##
## Same shape as the by-source map, but the toggle is the taxon's IUCN Red List
## category instead of the data source. Every taxon shown is one confirmed by
## >=2 of the 3 data sources - these maps subset combined_confirmed_l3, they do
## not re-derive presence.
##
## Three versions are produced from one builder, differing only in how the
## categories are bundled into layers.
##
## SOURCE OF THE CATEGORY: red_list_category_code in the WCFP plantlist. No new
## input file. Checked before use: wcfp_name_match is unique across all 26,632
## rows and no taxon carries two different categories, so the join cannot
## duplicate rows and no severity tie-breaking is actually needed (the code
## below still resolves by severity, as a guard for future data).
##
## COVERAGE - THE HEADLINE CAVEAT: only 10,909 of 26,632 taxa (41%) have any
## assessment, and 8,798 of those are Least Concern. The 15,723 Not Evaluated
## (NE) taxa are given their own layer rather than being dropped, because
## otherwise a reader cannot distinguish "assessed and not threatened" from
## "never assessed" - and on this dataset the second is the larger category
## by far.
## Author: Sarah Gora
## ============================================================================ ##

# RL_RECODE / RL_SEVERITY / RL_LABEL / wcfp_redlist_lookup and the
# unknown-code guard are defined much earlier, alongside wcfp_taxon_lookup,
# because the three taxa-list workbooks carry a Red List column too and are
# written long before this section runs. `rl_code` is therefore already a
# column on combined_confirmed_l3; the alias keeps the rest of this section
# reading as it did when the join happened here.
combined_confirmed_l3_rl <- combined_confirmed_l3

cat("Red List categories across the", n_distinct(combined_confirmed_l3_rl$taxon),
    "taxa confirmed by >=2 sources at level 3:\n")
print(combined_confirmed_l3_rl %>% distinct(taxon, rl_code) %>%
        count(rl_code, sort = TRUE, name = "taxa"))

# The three layer schemes. Each is an ordered list; the FIRST layer is the one
# left visible when the map opens, so the most conservation-relevant layer goes
# first in each.
rl_scheme_grouped <- list(
  "Threatened (CR, EN, VU, EW, EX)" = c("CR", "EN", "VU", "EW", "EX"),
  "Near Threatened (NT)"            = "NT",
  "Least Concern (LC)"              = "LC",
  "Data Deficient (DD)"             = "DD",
  "Not Evaluated (NE)"              = "NE"
)
# EW (4 taxa) and EX (2) are folded into Threatened above rather than given
# their own near-empty layers; they appear separately in the version below.
rl_scheme_categories <- list(
  "Critically Endangered (CR)" = "CR",
  "Endangered (EN)"            = "EN",
  "Vulnerable (VU)"            = "VU",
  "Extinct in the Wild (EW)"   = "EW",
  "Extinct (EX)"               = "EX",
  "Near Threatened (NT)"       = "NT",
  "Least Concern (LC)"         = "LC",
  "Data Deficient (DD)"        = "DD",
  "Not Evaluated (NE)"         = "NE"
)
rl_scheme_threatened <- list(
  "Threatened (CR, EN, VU, EW, EX)" = c("CR", "EN", "VU", "EW", "EX"),
  "Not threatened (NT, LC, DD)"     = c("NT", "LC", "DD"),
  "Not Evaluated (NE)"              = "NE"
)

# Each scheme must cover all nine categories exactly once. A missing category
# would silently drop those taxa from the map; a repeated one would draw them
# on two layers and double-count them in the subtitle. Neither is visible by
# eye on a finished map, so it is checked here.
for (.nm in c("rl_scheme_grouped", "rl_scheme_categories", "rl_scheme_threatened")) {
  .codes <- unlist(get(.nm), use.names = FALSE)
  if (anyDuplicated(.codes) > 0 || !setequal(.codes, RL_SEVERITY)) {
    stop(.nm, " does not partition the Red List categories exactly once.",
         "  missing: ", paste(setdiff(RL_SEVERITY, .codes), collapse = ", "),
         "  | repeated: ", paste(unique(.codes[duplicated(.codes)]), collapse = ", "),
         "  | unknown: ", paste(setdiff(.codes, RL_SEVERITY), collapse = ", "))
  }
}

## ------------------------------------------------------------------ ##
## Per-layer colour scales, and a legend that follows the active layer.
##
## The by-source maps share one palette across their four layers, because
## those layers are of comparable magnitude. Here they are not: Least Concern
## has 8,798 taxa and Extinct has 2. Under a shared scale every threatened
## layer would render as a uniform pale wash - the exact pattern these maps
## exist to show would be the one thing invisible.
##
## So each layer gets its own colour scale, and therefore its own legend. The
## cost is that colours are NOT comparable between layers, which is stated on
## the map itself in the hint box rather than left for the reader to infer.
##
## All the legends are emitted into one control, hidden, and shown/hidden by
## the JS below as layers are toggled. Matching is on the layer name, carried
## in data-layer. Driving it from the checkboxes rather than Leaflet's
## overlayadd/overlayremove events keeps this consistent with the rest of the
## file's DOM-based approach and avoids depending on `this` being the map.
## ------------------------------------------------------------------ ##
build_redlist_legends <- function(layer_names, pals, breaks_list) {
  inner <- vapply(seq_along(layer_names), function(i) {
    paste0("<div class=\"wcfp-rl-legend\" data-layer=\"", escape_attr(layer_names[i]),
           "\" style=\"display:none; margin-bottom:6px;\">",
           build_gradient_legend(pals[[i]], breaks_list[[i]],
                                 subtitle = layer_names[i]),
           "</div>")
  }, character(1))
  paste0("<div class=\"wcfp-rl-legends\">", paste(inner, collapse = ""), "</div>")
}

redlist_onrender_js <- paste0("function(el, x) {", map_onrender_common_js, "
  function addLayersTitle() {
    var lc = el.querySelector('.leaflet-control-layers');
    if (!lc) { return false; }
    if (lc.querySelector('.wcfp-layers-title')) { return true; }
    var h = document.createElement('div');
    h.className = 'wcfp-layers-title';
    h.textContent = 'IUCN Red List category';
    h.style.fontWeight = 'bold';
    h.style.fontSize = '13px';
    h.style.marginBottom = '5px';
    var list = lc.querySelector('.leaflet-control-layers-list') || lc;
    list.insertBefore(h, list.firstChild);
    return true;
  }
  if (!addLayersTitle()) { setTimeout(addLayersTitle, 200); }

  // Show the legend belonging to the selected layer. The layers are base
  // groups, so exactly one is ever on and exactly one legend ever shows -
  // which is what makes the per-layer colour scales safe to read.
  //
  // The input is matched on Leaflet's own class rather than on type=radio:
  // the class is the same whether the control renders radios (base groups) or
  // checkboxes (overlays), so this keeps working if the grouping is ever
  // switched back.
  function syncLegends() {
    var lc  = el.querySelector('.leaflet-control-layers');
    var box = el.querySelector('.wcfp-rl-legends');
    if (!lc || !box) { return false; }
    var active = {};
    var labels = lc.querySelectorAll('label');
    for (var i = 0; i < labels.length; i++) {
      var cb = labels[i].querySelector('input.leaflet-control-layers-selector') ||
               labels[i].querySelector('input');
      var sp = labels[i].querySelector('span');
      if (!cb || !sp) { continue; }
      active[sp.textContent.trim()] = cb.checked;
    }
    var legs = box.querySelectorAll('.wcfp-rl-legend');
    for (var j = 0; j < legs.length; j++) {
      legs[j].style.display = active[legs[j].getAttribute('data-layer')] ? 'block' : 'none';
    }
    return true;
  }
  function initLegendSync() {
    var lc = el.querySelector('.leaflet-control-layers');
    if (!lc) { return false; }
    // Leaflet updates the checkbox state before the change event reaches here,
    // but the deferred call also covers programmatic toggles.
    lc.addEventListener('change', function() { setTimeout(syncLegends, 0); });
    syncLegends();
    return true;
  }
  if (!initLegendSync()) { setTimeout(initLegendSync, 200); }
}")

## Builds one Red List toggle map from a scheme. Everything that differs
## between the three versions is in `layer_defs` and `out_name`.
build_redlist_l3_map <- function(layer_defs, out_name, lib_name) {
  layer_names <- names(layer_defs)

  # One richness surface + palette + legend breaks per layer.
  mdata <- lapply(layer_defs, function(codes) {
    r <- combined_confirmed_l3_rl %>%
      filter(rl_code %in% codes) %>%
      distinct(area_l3, taxon) %>%
      count(area_l3, name = "richness")
    level3_sf %>% left_join(r, by = c("LEVEL3_COD" = "area_l3"))
  })

  pals <- lapply(mdata, function(m) {
    hi <- suppressWarnings(max(m$richness, na.rm = TRUE))
    # A layer can legitimately be empty (no area reaches it); colorNumeric
    # cannot take an all-NA domain, so fall back to a nominal range.
    if (!is.finite(hi)) hi <- 1
    colorNumeric(palette = "YlGn", domain = c(0, hi), na.color = "#FFFFFF", reverse = FALSE)
  })

  brks <- lapply(mdata, function(m) {
    hi <- suppressWarnings(max(m$richness, na.rm = TRUE))
    if (!is.finite(hi)) return(numeric(0))
    b <- pretty(c(0, hi), n = 6)
    b <- b[b >= 0 & b <= hi & b == round(b)]
    sort(unique(b), decreasing = TRUE)
  })

  legend_html <- build_redlist_legends(layer_names, pals, brks)

  payload <- list(
    areaNames = as.list(setNames(area_l3_lookup$area_l3_name, area_l3_lookup$area_l3)),
    hideLabels = list(COMBINED_LAYER),
    # Three detail columns here, not one - the popup JS reads the count from
    # this entry. rep(list(...)) rather than rep(...) so each layer gets the
    # trio as a nested list, which jsonlite renders as a JSON array.
    headers   = setNames(rep(list(as.list(c("IUCN Red List category code",
                                            "IUCN Red List category",
                                            "Distribution data source"))),
                             length(layer_names)), layer_names),
    # Sourced from combined_l3_rows rather than combined_confirmed_l3_rl (the
    # same rows) purely to reuse its `info` column: build_data_source_label()
    # is row-wise, and recomputing it here would repeat ~275k calls once per
    # Red List map for a string that already exists.
    layers    = setNames(lapply(layer_defs, function(codes) {
      build_layer_export(
        combined_l3_rows %>%
          filter(rl_code %in% codes) %>%
          transmute(
            area_key    = area_l3,
            taxa        = coalesce(taxon_name_accepted, taxon),
            authority   = coalesce(taxon_authors_accepted, ""),
            rl_cat_code = rl_code,
            rl_cat_name = unname(RL_LABEL[rl_code]),
            dist_source = info
          ),
        detail_cols = c("rl_cat_code", "rl_cat_name", "dist_source")
      )
    }), layer_names)
  )

  taxa_counts <- vapply(layer_defs, function(codes)
    n_distinct(combined_confirmed_l3_rl$taxon[combined_confirmed_l3_rl$rl_code %in% codes]),
    integer(1))
  area_counts <- vapply(seq_along(layer_defs), function(i)
    sum(!is.na(mdata[[i]]$richness)), integer(1))

  title_html <- paste0(
    MAP_HEADING,
    "<div style='font-size:13px; color:#000; margin-top:3px;'>",
    paste0(layer_names, ": ", scales::comma(taxa_counts), " taxa across ",
           scales::comma(area_counts), " TDWG Level-3 botanical countries", collapse = "<br>"),
    "<br><span style='font-style:italic;'>Distribution of taxa confirmed in &ge;2 data sources</span></div>"
  )

  hint <- build_click_hint(paste0(
    "Select one IUCN Red List category - layers are exclusive, not stacked<br>",
    "<b>Each layer has its own colour scale</b> - shading is not comparable between layers<br>",
    "Click map to view/download food plants taxa list per TDWG Level-3 botanical country, ",
    "for the selected category"))

  mp <- leaflet(options = leafletOptions(zoomSnap = 0.05, zoomDelta = 0.5)) %>%
    fitBounds(
      lng1 = as.numeric(l3_map_bounds["xmin"]), lat1 = as.numeric(l3_map_bounds["ymin"]),
      lng2 = as.numeric(l3_map_bounds["xmax"]), lat2 = as.numeric(l3_map_bounds["ymax"])
    )

  for (i in seq_along(layer_names)) {
    # Colour/label/popup are precomputed as plain columns rather than passed as
    # palette formulas: inside this loop a ~formula would be evaluated later,
    # against whichever `i` the loop had finished on.
    m <- mdata[[i]] %>%
      mutate(
        fill_col = pals[[i]](richness),
        lab_txt  = paste0("<b>", escape_html(LEVEL3_NAM), "</b>: ",
                          scales::comma(ifelse(is.na(richness), 0, richness)), " taxa"),
        pop_html = popup_stub(layer_names[i], LEVEL3_COD, LEVEL3_NAM, !is.na(richness))
      )
    mp <- mp %>% addPolygons(
      data         = m,
      fillColor    = ~fill_col,
      fillOpacity  = 0.8,
      color        = "black",
      weight       = 1.3,
      label        = ~lapply(lab_txt, htmltools::HTML),
      popup        = ~pop_html,
      popupOptions = WCFP_POPUP_OPTIONS,
      group        = layer_names[i],
      highlightOptions = highlightOptions(weight = 4, color = "#000000", opacity = 1,
                                          bringToFront = TRUE)
    )
  }

  mp <- mp %>%
    addControl(html = title_html, position = "topleft") %>%
    addLayersControl(
      baseGroups = layer_names,
      position = "topleft",
      options = layersControlOptions(collapsed = FALSE)
    ) %>%
    # Base groups, so the categories are mutually exclusive. That matters more
    # here than on the by-source maps: each layer has its OWN colour scale, so
    # two drawn at once would be shaded on incomparable ramps. Only the first
    # is selected on open.
    hideGroup(layer_names[-1]) %>%
    addControl(html = hint, position = "topleft") %>%
    addControl(html = legend_html, position = "bottomright") %>%
    addControl(html = map_source_strip_l3, position = "bottomright",
               className = "wcfp-mapsource") %>%
    addScaleBar(
      position = "bottomleft",
      options = scaleBarOptions(maxWidth = 260, metric = TRUE, imperial = FALSE,
                                updateWhenIdle = TRUE)
    ) %>%
    htmlwidgets::onRender(redlist_onrender_js) %>%
    htmlwidgets::onRender(wcfp_popup_js, data = payload)

  path <- file.path(combined_output_dir, out_name)
  save_widget_selfcontained(widget = mp, file = path,
                            libdir = file.path(combined_output_dir, lib_name))
  message("Interactive Red List L3 map saved to: ", path)
  invisible(path)
}

build_redlist_l3_map(
  rl_scheme_grouped,
  "WCFP_taxa_distribution_by_redlist_grouped_L3area_map_interactive.html",
  "lib_tmp_rl_grouped")
build_redlist_l3_map(
  rl_scheme_categories,
  "WCFP_taxa_distribution_by_redlist_categories_L3area_map_interactive.html",
  "lib_tmp_rl_categories")
build_redlist_l3_map(
  rl_scheme_threatened,
  "WCFP_taxa_distribution_by_redlist_threatened_L3area_map_interactive.html",
  "lib_tmp_rl_threatened")

# ---------------------------------------- #
#   Static L3 choropleth (PNG + PDF)       #
# ---------------------------------------- #
# Styled to match the static country map (V1, magma, no graticules) so the two
# publication figures read as a pair. The interactive maps keep YlGn.
map_data_l3_4326 <- st_transform(map_data_l3, 4326)

p_l3 <- ggplot(map_data_l3_4326) +
  geom_sf(aes(fill = richness), color = "gray40", size = 0.15) +
  scale_fill_viridis_c(
    option = "magma",
    direction = -1,
    na.value = "gray95",
    name = "Number of food plant taxa"
  ) +
  theme(
    panel.grid = element_line(color = "transparent"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.03, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.box.margin = margin(0, 0, 0, 0)
  )

l3_png_path <- file.path(output_dir, "WCFP_taxa_distribution_L3area_map_v1_simple.png")
ragg::agg_png(l3_png_path, width = 12, height = 7, units = "in", res = 300)
print(p_l3)
dev.off()
message("Static L3 PNG saved to: ", l3_png_path)

l3_pdf_path <- file.path(output_dir, "WCFP_taxa_distribution_L3area_map_v1_simple.pdf")
ggsave(
  filename = l3_pdf_path,
  plot     = p_l3,
  width    = 12,
  height   = 7,
  device   = cairo_pdf,
  bg       = "white"
)
message("Static L3 PDF saved to: ", l3_pdf_path)


### end
