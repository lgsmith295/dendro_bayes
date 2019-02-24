library(dplR)

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
keep_cores <- as.character(rho[which(rho$rho > 0.5), "core"])
raw_cor <- raw %>%
  dplyr::filter(core %in% keep_cores)


##### Detrend with dplR before testing correlations with climate #####
y_ji <- raw_cor %>%
  dplyr::select(year, core, rwl) %>%
  # mutate(row = row_number()) %>%
  spread(key = core, value = rwl) %>%
  dplyr::select(-year)

y_detrend_negexp <- detrend(as.data.frame(y_ji), method = "ModNegExp") #, return.info = TRUE) # doesn't like negexp for all 
str(y_detrend_negexp)
summary(y_detrend_negexp[ , 1])

y_detrend_negexp_info <- detrend(as.data.frame(y_ji), method = "ModNegExp", return.info = TRUE, constrain.nls = "always")

# could run each model and then do the automated detrending based on negative values being produced in each method going down a list of prefernces
y_detrend_info <- detrend(as.data.frame(y_ji), return.info = TRUE)
detrend_method <- rep(NA_character_, times = length(y_detrend_info$model.info))
y_detrend_info$model.info[[1]][[1]]$n.zeros
for(i in 1:length(y_detrend_info$model.info)) {
  detrend_method[i] <- y_detrend_info$model.info[[1]]$ModNegExp$method
}

if(FALSE) {
  names(y_detrend_negexp_info)
  str(y_detrend_negexp_info$series)
  summary(y_detrend_negexp_info$series$`22316_BFM01A`)
  summary(y_detrend_negexp_info$curves$`22316_BFM01A`)
  str(y_detrend_negexp_info$model.info)
  str(y_detrend_negexp_info$data.info)
  
  y_detrend_all <- detrend(as.data.frame(y_ji))
  str(y_detrend_all)
  summary(y_detrend_all$`22316_BFM01A`)
}

# dplR documentation reads as though it should run a linear or mean detrending if negexp fails (like ARSTAN) but even when it's unconstrained and their are warnings it only returns the negative exponential
detrend_method <- rep(NA_character_, times = length(y_detrend_negexp_info$model.info))
for(i in 1:length(y_detrend_negexp_info$model.info)) {
  detrend_method[i] <- y_detrend_negexp_info$model.info[[1]]$ModNegExp$method
}
detrend_method[which(detrend_method != "NegativeExponential")]

########### check correlations with climate and remove cores? #############

# gather detrended 
detrended <- y_detrend_negexp_info$series %>%
  # add year
  gather(core, detrend_rwl)


cores <- unique(raw_cor$core)




core_cor <- NA
for(i in 1:length(unique(raw$core))) {
  tmp <- raw %>%
    filter(core == cores[i],
           !is.na(rwl))
}

cors <- raw %>%
  group_by(core) %>%
  summarise(core_cor = cor)

########## check sample depth and remove period with fewer than 5 cores? #####

samp_depth <- raw %>%
  ungroup() %>%
  group_by(year, sp_code) %>%
  filter(!is.na(rwl)) %>%
  select(year, sp_code) %>%
  summarise(n = n())

ggplot(samp_depth, aes(year, n)) + geom_line(aes(color = sp_code)) + ylab("Number of cores (sample depth)") + xlab("Year") + theme_bw() + labs(color = "Sp. Code")
ggsave("Results/Figures/NM/sample_depth.pdf")

ggplot(filter(samp_depth, n >= 10), aes(year, n)) + geom_line(aes(color = sp_code))

##### examine with dplR to see what's selected for detrending ####






##### Save #####

save(data = core_climate, file = "Data/itrdb_pilo_mount_washington/pilo_rwl_climate.RData")

ggplot(data = core_climate, aes(year, rwi_norm)) + geom_line(alpha = 0.1, color = "blue", aes(group = core)) + geom_smooth() + theme_bw() + theme(legend.position="none")

ggplot(data = core_climate, aes(year, rwi)) + geom_line(alpha = 0.1, color = "black", aes(group = core)) + geom_smooth() +  theme_bw() + theme(legend.position="none")

ggplot(data = core_climate, aes(year, log(rwi))) + geom_line(alpha = 0.1, color = "black", aes(group = core)) + geom_smooth() +  theme_bw() + theme(legend.position="none")

ggplot(data = core_climate, aes(year, rwi_norm)) + geom_line() + facet_wrap(~core)

ggplot(data = filter(core_climate, type == "ppt" & summer == 1), aes(rwi_norm, value)) + geom_point() + facet_wrap(~month)

ggplot(data = core_climate, aes(rwi_norm, value)) + geom_point() + facet_wrap(~type)

ggplot(data = core_climate, aes(tmean, rwi_norm)) + geom_point() + facet_wrap(~month)

ggplot(data = core_climate, aes(ppt, rwi_norm)) + geom_point() + facet_wrap(~month)

cor(core_climate$rwi_norm, core_climate$tmean, use = "complete.obs")
cor(core_climate$core_climate, data$ppt, use = "complete.obs")

