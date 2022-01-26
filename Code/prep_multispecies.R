library(dplR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rjson)
library(readr)
library(stringr)

# Check metadata
meta_files <- list.files("Data/az_nm/metadata/")
meta_list <- list()

n_sites <- length(meta_files)
for(i in 1:n_sites) {
  meta_list[[i]] <- fromJSON(file = paste0("Data/az_nm/metadata/", meta_files[i]))
  names(meta_list)[[i]] <- meta_list[[i]]$studyCode
}

# str(meta_list[[1]])

# dataframe of coordinates and site info for each study site
site_df <- data.frame(study = rep(NA, times = n_sites), 
                      noaa_id = rep(NA, times = n_sites), 
                      sp_code = rep(NA, times = n_sites), 
                      species = rep(NA, times = n_sites), 
                      common = rep(NA, times = n_sites),
                      first_yr = rep(NA, times = n_sites),
                      last_yr = rep(NA, times = n_sites),
                      lat = rep(NA, times = n_sites), 
                      lon = rep(NA, times = n_sites),
                      units = rep(NA, times = n_sites))
for(i in 1:n_sites) {
  tmp <- meta_list[[i]]$site[[1]]$geo$geometry$coordinates
  site_df[i, "study"] <- meta_list[[i]]$studyCode
  site_df[i, "noaa_id"] <- meta_list[[i]]$NOAAStudyId
  site_df[i, "sp_code"] <- meta_list[[i]]$site[[1]]$paleoData[[1]]$species[[1]]$speciesCode
  site_df[i, "species"] <- meta_list[[i]]$site[[1]]$paleoData[[1]]$species[[1]]$scientificName
  site_df[i, "common"] <- meta_list[[i]]$site[[1]]$paleoData[[1]]$species[[1]]$commonName[1]
  site_df[i, "first_yr"] <- meta_list[[i]]$earliestYearCE
  site_df[i, "last_yr"] <- meta_list[[i]]$mostRecentYearCE
  site_df[i, "lat"] <- as.numeric(tmp[1])
  site_df[i, "lon"] <- as.numeric(tmp[2])
  tmpl <- length(meta_list[[i]]$site[[1]]$paleoData[[1]]$dataFile)
  if(tmpl < 3) {
    foo <- try(meta_list[[i]]$site[[1]]$paleoData[[1]]$dataFile[[1]]$variables[[1]]$cvUnit)
  } else {
    foo <- try(meta_list[[i]]$site[[1]]$paleoData[[1]]$dataFile[[3]]$variables[[1]]$cvUnit)
  } # really not sure how to automate this since the units come in at different places in the json input. Maybe try json lite or some way to search the resulting list for units
  if(isTRUE(class(foo)=="try-error")) { 
    site_df[i, "units"] <- NA 
  } else { 
    if(is.null(foo)) {
      site_df[i, "units"] <- NA
    } else {
  site_df[i, "units"] <- foo
    }
  }
}

site_df <- site_df %>%
  dplyr::filter(units == "millimeter")

unique(site_df$species) # need to remove unidentified pine species and bristlecone (PIAR) might be more temperature sensitive than moisture if at elevation creating problems
unique(site_df$common)

# NOAA Climate Divisions
# Shapefiles from: https://www.esrl.noaa.gov/psd/data/usclimdivs/boundaries.html
# Map at: https://www.ncdc.noaa.gov/monitoring-references/maps/images/us-climate-divisions-names.jpg
# New Mexico (29) Northwest Mountains (04)
# Problem - no trans-state-boundary information about common climate
# Need a better and reproducible method in the future

# CONUS NOAA Division Shapefile
library(sp)
library(rgdal)
divs <- readOGR("Data/CONUS_CLIMATE_DIVISIONS.shp/", "GIS.OFFICIAL_CLIM_DIVISIONS")
# plot(divs)

nm <- subset(divs, STATE == "New Mexico")
plot(nm)
str(nm)

# Convert points to spatial points data frame
site_df_sp <- site_df
coordinates(site_df_sp) <- c('lon', 'lat')
proj4string(site_df_sp) <- proj4string(divs)

library(ggplot2)
ggplot(site_df, aes(lon, lat)) + geom_polygon(data = nm, aes(long, lat, group = group), fill = "cadetblue", color = "grey") + geom_point() + coord_fixed(1.3)

# clip (subset) points within Northwest Mountains division (04) of NM
nm04 <- subset(nm, CD_2DIG == "04")
# plot(nm04)

site_df_sp_2904 <- site_df_sp[nm04, ] 

ggplot(as.data.frame(site_df_sp_2904@coords), aes(lon, lat)) + geom_polygon(data = nm04, aes(long, lat, group = group), fill = "cadetblue", color = "grey") + geom_point(color = "red") + coord_fixed(1.3)

site_df_2904 <- as.data.frame(site_df_sp_2904, stringsAsFactors = FALSE)

ggplot(site_df_2904, aes(lon, lat)) + geom_polygon(data = nm, aes(long, lat, group = group), fill = "cadetblue", color = "grey") + geom_point(data = site_df_2904, color = "black") + coord_fixed(1.3)

ggplot(site_df_2904, aes(lon, lat)) + geom_polygon(data = nm, aes(long, lat, group = group), fill = "lightgray", color = "black") + geom_point(data = site_df_2904, aes(color = sp_code)) + coord_fixed(1.3) + xlab("Longitude") + ylab("Latitude") + theme_bw() + labs(color='Species Code') + ggtitle("New Mexico NOAA Climate Divisions")

ggsave("Results/Figures/NM/nm_climate_division_itrdb.tiff")

unique(site_df_2904$common)
unique(site_df_2904$species)

##### import core ring data and subset #####
# import raw data downloaded from ITRDB
sites <- as.character(unique(site_df_2904$study))
tmplog <- list()
tmplog[[1]] <- capture.output(tmp <- read.rwl(paste0("Data/az_nm/data/pub/data/paleo/treering/measurements/northamerica/usa/", sites[1], ".rwl")))
# tmp <- t(tmp)
tmp2 <- tmp
rownames(tmp2) <- c()
raw <- data.frame(year = as.integer(rownames(tmp)),
                  study = site_df_2904$study[1],
                  noaa_id = site_df_2904$noaa_id[1],
                  sp_code = site_df_2904$sp_code[1],
                  tmp2, 
                  stringsAsFactors = FALSE)
raw <- raw %>%
  tidyr::gather(key = core1, value = rwl, -year, -study, -noaa_id, -sp_code)

for(i in 2:length(sites)) {
  tmplog[[i]] <- capture.output(tmp <- try(read.rwl(paste0("Data/az_nm/data/pub/data/paleo/treering/measurements/northamerica/usa/", sites[i], ".rwl")), TRUE))
    if(isTRUE(class(tmp)=="try-error")) { 
      next 
    } else { 
      tmp2 <- tmp
      rownames(tmp2) <- c()
      tmp2 <- data.frame(year = as.integer(rownames(tmp)),
                         study = site_df_2904$study[i],
                         noaa_id = site_df_2904$noaa_id[i],
                         sp_code = site_df_2904$sp_code[i],
                         tmp2, 
                         stringsAsFactors = FALSE)
      tmp2 <- tmp2 %>%
        tidyr::gather(key = core1, value = rwl, -year, -study, -noaa_id, -sp_code)
      raw <- bind_rows(raw, tmp2)
    }
  }

str(raw)
head(raw)

raw <- raw %>%
  unite(core, noaa_id, core1, remove = FALSE)

########### Look for duplicate year-core combos ##########

dups <- raw %>% 
  group_by(year, core) %>% 
  mutate(n = n()) %>% 
  filter(n > 1) %>%
  arrange(year, core)

dups 


##### summaries #####
length(unique(raw$study))
length(unique(raw$year))
length(unique(raw$sp_code))
length(unique(raw$core))
unique(raw$sp_code)
unique(site_df_2904$sp_code)

summary(raw) # problem: rwl up to 22 but all units say millimeters
dev.off()
hist(raw$rwl)

bigs <- raw %>%
  dplyr::filter(rwl > 6)
unique(bigs$study) 
(bad_cores <- unique(bigs$core))
# bigs # many studies, not sure what a cutoff really would be for unrealistic growth. Such a small proportion

if(!dir.exists("Results/Figures/NM")) dir.create("Results/Figures/NM")
g <- ggplot(data = raw, aes(year, rwl)) + geom_line(aes(color = core), alpha = 0.1) + theme(legend.position = "none") + facet_wrap(~sp_code)
ggsave("Results/Figures/NM/raw_series.pdf")

# g <- ggplot(data = raw, aes(year, rwl)) + geom_line() + facet_wrap(~core) + theme(legend.position = "none")
# ggsave("Results/Figures/NM/raw_series_wrap.pdf")

good_cores <- unique(raw$core)[!(unique(raw$core) %in% bad_cores)]

raw <- raw %>%
  group_by(core) %>%
  filter(core %in% good_cores) %>%
  mutate(# rwl = ifelse(rwl > 6, NA, rwl), # filtering by semi-arbitraty cutoff but only losing ~32 of 1440 cores. Fucking ITRDB.
         rwl_norm = (rwl - mean(rwl, na.rm = T)) / sd(rwl, na.rm = T)) %>%
  ungroup()

summary(raw$rwl)
summary(raw$rwl_norm)

g <- ggplot(data = raw, aes(year, rwl_norm)) + geom_line(aes(color = sp_code), alpha = 0.1) + ylab("Normalized RWL") + facet_wrap(~study) + theme(legend.position = "none")
ggsave("Results/Figures/NM/raw_series_norm.pdf")

##### Save raw data #####
save(data = raw, file = "Data/az_nm/raw_2904.RData")


##### Prep Climate Data #####
climate <- read_csv("Data/az_nm/nm_2904_nwmnts_precip_inch.csv")

climate2 <- climate %>%
  dplyr::select(-year) %>%
  dplyr::mutate(precip = rowSums(.)*2.54) %>%
  dplyr::select(precip) %>%
  bind_cols(. , climate[ , "year"])

climate_mo <- climate %>%
  group_by(year) %>%
  gather(month, prec, -year) %>%
  dplyr::mutate(prec = prec * 2.54) %>%
  ungroup()

save(climate2, climate_mo, file = "Data/az_nm/annual_ppt_cm.RData")
