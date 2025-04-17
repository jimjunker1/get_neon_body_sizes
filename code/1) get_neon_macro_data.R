library(here)
library(neonUtilities)
library(neonstore)
library(tidyverse)
library(lubridate)
library(janitor)

source(here("code/inverts_dw-functions.R"))
source(here("code/update_data_products.R"))
# 1) Download data ---------------------------------------------------------

update <- TRUE
if (update) {
  neonstore::neon_download(
    product = "DP1.20120.001",
    dir = here("data/database-files"),
    .token = Sys.getenv("NEON_TOKEN")
  )
  neonstore::neon_store(
    product = "DP1.20120.001",
    table = "inv_taxonomyProcessed-basic",
    dir = here("data/database-files")
  )
}

# 1.1) check file index
neonstore::neon_index(
  product = "DP1.20120.001",
  dir = here("data/database-files")
) %>%
  slice_sample(n = 1, by = table) %>%
  View()

# 1.2) load data and subset to necessary data
macro_load <- neon_table(
  product = "DP1.20120.001",
  table = "inv_taxonomyProcessed-basic", lazy = TRUE
)

macro <- macro_load %>%
  # select(siteID, collectDate, scientificName, acceptedTaxonID, sizeClass, individualCount, subsamplePercent, estimatedTotalCount, sampleCondition) %>%
  # summarise() %>%
  collect()


# load Length Weight coefficient table (used in part C below)
coeff <- read.csv(here("data/macro_lw_coeffs.csv"))
# neon_token <- source("C:/Users/jfpom/Documents/Wesner/NEON documents/neon_token_source.R")
# stream_sites <- readRDS("data/streams.rds")

# 2) Add LW coefficients, estimate dry weights  ------------------------------------

# add length weight coefficients by taxon
MN.lw <- LW_coef(
  x = macro$inv_taxonomyProcessed,
  lw_coef = coeff,
  percent = TRUE
)

# questionable measurements ####
# filter out individuals that were "damaged" and measurement was affected
# this is a flag which is added by NEON
MN.no.damage <- MN.lw %>%
  filter(!str_detect(
    sampleCondition,
    "measurement"
  )) %>%
  est_dw(fieldData = macro$inv_fieldData)



# 3) filter out NA values in dw
macro_dw <- MN.no.damage %>%
  filter(!is.na(dw), !is.na(no_m2))

nrow(MN.no.damage) / nrow(macro_dw)

saveRDS(macro_dw, file = "data/macro_dw_raw.rds")

macro_dw <- readRDS(file = "data/macro_dw_raw.rds")
# remove taxonomic information and tally density by size class
macro_dw_sizebytaxa <- macro_dw %>%
  group_by(siteID, collectDate, dw, family, genus, acceptedTaxonID) %>%
  reframe(no_m2 = mean(no_m2))

saveRDS(macro_dw_sizebytaxa, file = "data/macro_dw_sizebytaxa.rds")

# remove taxonomic information and tally density by size class
macro_dw_sizeonly <- macro_dw %>%
  group_by(siteID, collectDate, dw) %>%
  reframe(no_m2 = mean(no_m2)) %>%
  group_by(siteID, collectDate, dw) %>%
  reframe(no_m2 = mean(no_m2))

saveRDS(macro_dw_sizeonly, file = "data/macro_dw_sizeonly.rds")

# 4) filter to only ranges that follow a power law
dat_inverts <- macro_dw_sizeonly %>%
  clean_names() %>%
  group_by(site_id, collect_date) %>%
  sample_n(5000, weight = no_m2, replace = T)

dat_inverts_list <- dat_inverts %>%
  group_by(site_id, collect_date) %>%
  group_split()

xmin_inverts_list <- list()

for (i in 1:length(dat_inverts_list)) {
  powerlaw <- conpl$new(dat_inverts_list[[i]]$dw)
  xmin_inverts_list[[i]] <- tibble(
    xmin_clauset = estimate_xmin(powerlaw)$xmin,
    site_id = unique(dat_inverts_list[[i]]$site_id),
    collect_date = unique(dat_inverts_list[[i]]$collect_date)
  )
}

xmins_inverts_clauset <- bind_rows(xmin_inverts_list)

saveRDS(xmins_inverts_clauset, file = "data/xmins_inverts_clauset.rds")

dat_inverts_clauset_xmins <- macro_dw_sizeonly %>%
  clean_names() %>%
  left_join(xmins_inverts_clauset) %>%
  group_by(site_id, collect_date) %>%
  filter(dw >= xmin_clauset) %>%
  mutate(
    xmin = xmin_clauset,
    xmax = max(dw)
  )

saveRDS(dat_inverts_clauset_xmins, file = "data/dat_inverts_clauset_xmins.rds")
