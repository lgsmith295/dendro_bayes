library(dplR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(treeclim)
source("Code/functions.R")

load(file = "Data/az_nm/raw_2904.RData")


##### Look at individual core series for weird ones with potential data entry errors to throw out #####

########### check correlations among cores #############
y_ji <- raw %>%
  dplyr::select(year, core, rwl) %>%
  # mutate(row = row_number()) %>%
  spread(key = core, value = rwl) %>%
  dplyr::select(-year)

core_cors <- corr.rwl.seg(y_ji, make.plot = FALSE)

# 0.4 typical cutoff with overall correlation across full series in east
rho <- as.data.frame(core_cors$overall, stringsAsFactors = FALSE)
rho$core <- rownames((core_cors$overall))
keep_cores <- as.character(rho[which(rho$rho > 0.4), "core"])

raw_cor <- raw %>%
  dplyr::filter(core %in% keep_cores)


##### Detrend with dplR before testing correlations with climate #####
y_ji <- raw_cor %>%
  dplyr::select(year, core, rwl) %>%
  # mutate(row = row_number()) %>%
  spread(key = core, value = rwl) %>%
  dplyr::select(-year)

# y_detrend_negexp <- detrend(as.data.frame(y_ji), method = "ModNegExp") #, return.info = TRUE) # doesn't like negexp for all 
# str(y_detrend_negexp)
# summary(y_detrend_negexp[ , 1])

y_detrend_negexp_info <- detrend(as.data.frame(y_ji), method = "ModNegExp", return.info = TRUE, constrain.nls = "always")

y_detrend_negexp <- detrend(as.data.frame(y_ji), method = "ModNegExp", return.info = FALSE, constrain.nls = "always")

if(FALSE) {
# could run each model and then do the automated detrending based on negative values being produced in each method going down a list of prefernces
y_detrend_info <- detrend(as.data.frame(y_ji), return.info = TRUE)
detrend_method <- rep(NA_character_, times = length(y_detrend_info$model.info))
y_detrend_info$model.info[[1]][[1]]$n.zeros
for(i in 1:length(y_detrend_info$model.info)) {
  detrend_method[i] <- y_detrend_info$model.info[[1]]$ModNegExp$method
}

  names(y_detrend_negexp_info)
  str(y_detrend_negexp_info$series)
  summary(y_detrend_negexp_info$series$`22316_BFM01A`)
  summary(y_detrend_negexp_info$curves$`22316_BFM01A`)
  str(y_detrend_negexp_info$model.info)
  str(y_detrend_negexp_info$data.info)
  
  y_detrend_all <- detrend(as.data.frame(y_ji))
  str(y_detrend_all)
  summary(y_detrend_all$`22316_BFM01A`)

# dplR documentation reads as though it should run a linear or mean detrending if negexp fails (like ARSTAN) but even when it's unconstrained and there are warnings it only returns the negative exponential
detrend_method <- rep(NA_character_, times = length(y_detrend_negexp_info$model.info))
for(i in 1:length(y_detrend_negexp_info$model.info)) {
  detrend_method[i] <- y_detrend_negexp_info$model.info[[1]]$ModNegExp$method
}
detrend_method[which(detrend_method != "NegativeExponential")]

}

##### Make Chronology #####

crn <- chron(y_detrend_negexp)

##### Run treeclim to look for correlated climate #####
# load climate 
load("Data/az_nm/annual_ppt_cm.RData")

months_df <- data.frame(month = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"),
                        mon = 1:12, stringsAsFactors = FALSE)

climate_mo2 <- climate_mo %>%
  left_join(months_df) %>%
  dplyr::select(year, mon, prec) %>%
  rename(month = mon)
  
climate_13mo <- climate_mo2 %>%
  spread(month, prec) %>%
  as.data.frame(. , stringsAsFactors = FALSE)

rownames(crn) <- years

dcc1 <- dcc(chrono = crn, climate = climate_13mo, selection = -12:12)
plot(dcc1)

# Use January through July for correlated climate
climate_jj <- data.frame(year = climate_13mo$year, precip = apply(climate_13mo[ , 2:8], MARGIN = 1, FUN = mean), stringsAsFactors = FALSE)

climate2 <- climate_jj
save(climate2, file = "Data/az_nm/jan_jul_ppt_cm.RData")
  
########### check correlations with climate and remove cores? #############

if(FALSE) {
# make single chronology to check with treeclim using detrended series
chron1 <- chron(y_detrend_negexp_info$series[1950:2010, ])

# Use treeclim to explore
clim <- as.data.frame(climate_mo[which(climate_mo$year %in% (1950:max(raw$year))), ], stringsAsFactors = FALSE)

clim <- clim %>%
  mutate(month = as.integer(as.factor(month)))

dcc1 <- dcc(chron1, climate = clim, moving = TRUE)

plot(dcc1)

dc7 <- dcc(chrono = chron1, climate = clim, 
           selection = .mean(1:12) + .mean(1:9), method = "cor",
           dynamic = TRUE, win_size = 50, sb = FALSE)

plot(dc7)
g_test(dc7, sb = FALSE) 

# sc1 <- seascorr(chron1, clim) # if using prcp and temp

# Not really sure how to use treeclim and don't have time. It's not showing any real correlations.
}

##### Check for correlations manually with each core and annual precip #####
# gather detrended 
detrended <- y_detrend_negexp_info$series %>%
  mutate(year = sort(unique(raw_cor$year))) %>%
  gather(core, detrend_rwl, -year) %>%
  left_join(climate2) %>%
  group_by(core) %>%
  filter(year >= 1895) %>%
  summarise(cor = cor(detrend_rwl, precip, use = "pairwise.complete.obs", method = "spearman"))

summary(detrended)
dim(detrended[which(detrended$cor < 0.4), ])

hist(detrended$cor)

bad_cores <- detrended[which(detrended$cor < 0.4), ]$core # not sure what cutoff to use. A bit arbitrary, but at least get rid of cores with little to no correlation

cores <- unique(raw_cor$core)
keep_cores2 <- keep_cores[which(!(keep_cores %in% bad_cores))]

raw_cor2 <- raw_cor %>%
  dplyr::filter(core %in% keep_cores2)

# cores <- unique(raw_cor$core)


########## check sample depth and remove period with fewer than 5 cores? #####

samp_depth <- raw_cor2 %>%
  ungroup() %>%
  group_by(year, sp_code) %>%
  filter(!is.na(rwl)) %>%
  select(year, sp_code) %>%
  summarise(n = n())

g <- ggplot(samp_depth, aes(year, n)) + geom_line(aes(color = sp_code)) + ylab("Number of cores (sample depth)") + xlab("Year") + theme_bw() + labs(color = "Sp. Code") # consider cutting reconstruction at year 1000 (good depth of 1 species) or 1700 (good depth of 3 species)
ggsave("Results/Figures/NM/sample_depth.pdf")

g + theme_bw_poster()
ggsave("Results/Figures/NM/sample_depth_poster.pdf")

g + geom_hline(yintercept = 10, color = "red")

# ggplot(filter(samp_depth, n >= 10), aes(year, n)) + geom_line(aes(color = sp_code))

##### Save #####

saveRDS(raw_cor2, "Data/az_nm/raw_correlated_2904.RData")


