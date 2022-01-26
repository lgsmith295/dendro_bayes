library(dplR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rjson)

# Check metadata
meta_files <- list.files("Data/itrdb_pilo_mount_washington/metadata/")
meta_list <- list()

for(i in 1:length(meta_files)) {
  meta_list[[i]] <- fromJSON(file = paste0("Data/itrdb_pilo_mount_washington/metadata/", meta_files[i]))
  names(meta_list)[[i]] <- meta_list[[i]]$studyCode
}

str(meta_list)

coords <- meta_list$NV520$site[[1]]$geo$geometry$coordinates

# import raw data downloaded from ITRDB
# raw <- read.rwl(paste0("Data/itrdb_pilo_white_mnts/Data/rwi/", names(meta_list)[[1]], ".rwl")) # metadata and filenames inconsistent in capitalization. No bueno.

raw <- read.rwl("Data/norw015.rwl") # changed the first half of MW069x to MW069xc because it looked like might be MW069a, MW069b, and two sets called MW069x so it wouldn't read in as downloaded from the ITRDB

str(raw)
head(raw)
summary(raw)
sort(unique(names(raw)))
length(unique(names(raw)))

raw_df <- data.frame(year = as.integer(rownames(raw)), raw, stringsAsFactors = FALSE)
str(raw_df)

raw_long <- raw_df %>%
  gather(tree, rwi, -year) %>%
  mutate(locs.id = 1,
         loc = names(meta_list)[[1]],
         sp = meta_list$NV520$scienceKeywords[1]) # meta_list$NV520$site[[1]]$paleoData[[1]]$species[[1]]$speciesCode

str(raw_long)

library(treeclim)
norway_prec
precip <- norway_prec %>%
  mutate(year = YEAR)

climate2 <- norway_temp %>%
  mutate(year = YEAR)

load(paste0("Data/prism_", names(meta_list)[[1]], ".RData"))
data <- raw_long %>%
  full_join(climate2)

tree_growth <- raw_long %>%
  group_by(tree) %>%
  summarize(rwi_mean = mean(rwi, na.rm = TRUE),
            rwi_sd = sd(rwi, na.rm = TRUE))

data <- data %>%
  left_join(tree_growth) %>%
  mutate(rwi_norm = (rwi - rwi_mean)/rwi_sd)

str(data)

# ggplot(data = data, aes(year, rwi_norm, color = tree)) + geom_line()

# ggplot(data = data, aes(year, rwi_norm)) + geom_line() + facet_wrap(~tree)

ggplot(data = data, aes(rwi_norm, MAY)) + geom_point() #+ facet_wrap(~month)

ggplot(data = data, aes(rwi_norm, ppt)) + geom_point() + facet_wrap(~month) # + facet_wrap(~type)

ggplot(data = data, aes(JUL, rwi_norm)) + geom_point() + geom_smooth() + facet_wrap(~tree) # + facet_wrap(~month)

ggplot(data = data, aes(ppt, rwi_norm)) + geom_point() + facet_wrap(~month)

cor(data$rwi_norm, data$tmean, use = "complete.obs")
cor(data$rwi_norm, data$ppt, use = "complete.obs")

