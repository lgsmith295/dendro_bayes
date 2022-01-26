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

raw <- read.rwl("Data/itrdb_pilo_mount_washington/rwi/nv520.rwl") # changed the first half of MW069x to MW069xc because it looked like might be MW069a, MW069b, and two sets called MW069x so it wouldn't read in as downloaded from the ITRDB

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


load(paste0("Data/prism_", names(meta_list)[[1]], ".RData"))
tree_climate <- raw_long %>%
  full_join(climate)

tree_growth <- raw_long %>%
  group_by(tree) %>%
  summarize(rwi_mean = mean(rwi, na.rm = TRUE),
            rwi_sd = sd(rwi, na.rm = TRUE))

tree_climate <- tree_climate %>%
  left_join(tree_growth) %>%
  mutate(rwi_norm = (rwi - rwi_mean)/rwi_sd,
         summer = ifelse(month %in% c(5,6,7,8), 1, 0))

str(tree_climate)

save(data = tree_climate, file = "Data/itrdb_pilo_mount_washington/pilo_rwl_climate.RData")

ggplot(data = tree_climate, aes(year, rwi_norm)) + geom_line(alpha = 0.1, color = "blue", aes(group = tree)) + geom_smooth() + theme_bw() + theme(legend.position="none")

ggplot(data = tree_climate, aes(year, rwi)) + geom_line(alpha = 0.1, color = "black", aes(group = tree)) + geom_smooth() +  theme_bw() + theme(legend.position="none")

ggplot(data = tree_climate, aes(year, log(rwi))) + geom_line(alpha = 0.1, color = "black", aes(group = tree)) + geom_smooth() +  theme_bw() + theme(legend.position="none")

ggplot(data = tree_climate, aes(year, rwi_norm)) + geom_line() + facet_wrap(~tree)

ggplot(data = filter(tree_climate, type == "ppt" & summer == 1), aes(rwi_norm, value)) + geom_point() + facet_wrap(~month)

ggplot(data = tree_climate, aes(rwi_norm, value)) + geom_point() + facet_wrap(~type)

ggplot(data = tree_climate, aes(tmean, rwi_norm)) + geom_point() + facet_wrap(~month)

ggplot(data = tree_climate, aes(ppt, rwi_norm)) + geom_point() + facet_wrap(~month)

cor(tree_climate$rwi_norm, tree_climate$tmean, use = "complete.obs")
cor(tree_climate$tree_climate, data$ppt, use = "complete.obs")

