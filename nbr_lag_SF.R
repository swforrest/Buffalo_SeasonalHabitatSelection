###########################
# Creating lagged variables for the NDWI covariate
# Scott Forrest
###########################

library(tidyverse)
library(TwoStepCLogit) # to make population estimates
library(sjPlot) # to plot model estimates
library(glmmTMB)# for the mixed models 
library(tictoc)
library(beepr)
library(amt)
library(terra)
library(RColorBrewer)

# create a colour scheme
colours <- rev(brewer.pal(11, "RdBu"))

ssf_dat <- read_csv("pheno_start_end_dat_ssf10rs_2024-03-14.csv")

# ssf_dat <- ssf_dat %>% mutate(srtm_start_scaled = scale(srtm_start),
#                               srtm_end_scaled = scale(srtm_end),
#                               water_dist_start_scaled = scale(water_dist_start),
#                               water_dist_end_scaled = scale(water_dist_end),
#                               ndvi_start_scaled = scale(ndvi_start),
#                               ndvi_end_scaled = scale(ndvi_end),
#                               nbr_start_scaled = scale(nbr_start),
#                               nbr_end_scaled = scale(nbr_end),
#                               sl_scaled = scale(sl_),
#                               log_sl_scaled = scale(log_sl_),
#                               cos_ta_scaled = scale(cos_ta_)
# )

nbr <- rast("nbr_raster_amt_cropped_epsg32753.tif")
# nbr <- terra::project(nbr ,"epsg:32753")
# writeRaster(nbr, "nbr_raster_amt_cropped_epsg32753.tif")
terra::time(nbr) <- as.POSIXct(lubridate::ymd("2018-01-01") + months(0:23))

NAFI2018_cropped_resampled <- rast("NAFI2018_cropped_resampled_raster_epsg32753.tif")
NAFI2019_cropped_resampled <- rast("NAFI2019_cropped_resampled_raster_epsg32753.tif")

nbr2018_NAFImasked <- rast("nbr2018_NAFImask_epsg32753.tif")
nbr2019_NAFImasked <- rast("nbr2019_NAFImask_epsg32753.tif")

# NAFI2018 <- rast("2018 firescar image files/FS2018_MTHS.tif")
# NAFI2019 <- rast("2019 firescar image files/FS2019_MTHS.tif")
# # reproject the NAFI raster to match the NBR raster
# NAFI2018 <- terra::project(NAFI2018, crs(nbr))
# NAFI2019 <- terra::project(NAFI2019, crs(nbr))
# NAFI2018_cropped <- crop(NAFI2018, nbr)
# NAFI2019_cropped <- crop(NAFI2019, nbr)
# plot(NAFI2018_cropped)
# plot(NAFI2019_cropped)
# # Resample maskRaster to match the resolution of rasterToMask
# NAFI2018_cropped_resampled <- terra::resample(NAFI2018_cropped, nbr, method="bilinear")
# # plot(NAFI2018_cropped_resampled)
# NAFI2019_cropped_resampled <- terra::resample(NAFI2019_cropped, nbr, method="bilinear")
# # plot(NAFI2018_cropped_resampled)
# writeRaster(NAFI2018_cropped_resampled, "NAFI2018_cropped_resampled_raster_epsg32753.tif")
# writeRaster(NAFI2019_cropped_resampled, "NAFI2019_cropped_resampled_raster_epsg32753.tif")

# mask by the NAFI raster
# nbr2018_NAFImasked <- terra::mask(nbr[[1:12]], NAFI2018_cropped_resampled, maskvalue = 0, updatevalue = 0)
# writeRaster(nbr2018_NAFImasked, "nbr2018_NAFImask_epsg32753.tif")
# nbr2019_NAFImasked <- terra::mask(nbr[[13:24]], NAFI2019_cropped_resampled, maskvalue = 0, updatevalue = 0)
# writeRaster(nbr2019_NAFImasked, "nbr2019_NAFImask_epsg32753.tif")

nbr_NAFImask <- c(nbr2018_NAFImasked, nbr2019_NAFImasked)

for(i in 1:nlyr(nbr)) terra::plot(nbr[[i]], 
                                  col = colours, 
                                  breaks = seq(-1,1, length.out = length(colours) + 1),
                                  main = paste0("dNBR ", time(nbr[[i]])))

for(i in 1:nlyr(nbr_NAFImask)) terra::plot(nbr_NAFImask[[i]], 
                                           col = colours, 
                                           breaks = seq(-1,1, length.out = length(colours) + 1),
                                           main = paste0("dNBR NAFImask ", time(nbr_NAFImask[[i]])))

# Adding time components to NBR -------------------------------------------

nbr_lag_1 <- rep(nbr)
terra::time(nbr_lag_1) <- as.POSIXct(lubridate::ymd("2018-02-01") + months(0:23))

nbr_lag_2 <- rep(nbr)
terra::time(nbr_lag_2) <- as.POSIXct(lubridate::ymd("2018-03-01") + months(0:23))

nbr_lag_3 <- rep(nbr)
terra::time(nbr_lag_3) <- as.POSIXct(lubridate::ymd("2018-04-01") + months(0:23))

# nbr_lag_4 <- rep(nbr)
# terra::time(nbr_lag_4) <- as.POSIXct(lubridate::ymd("2018-05-01") + months(0:23))
# 
# nbr_lag_5 <- rep(nbr)
# terra::time(nbr_lag_5) <- as.POSIXct(lubridate::ymd("2018-06-01") + months(0:23))
# 
# nbr_lag_6 <- rep(nbr)
# terra::time(nbr_lag_6) <- as.POSIXct(lubridate::ymd("2018-07-01") + months(0:23))


nbr_NAFImask_lag_1 <- rep(nbr_NAFImask)
terra::time(nbr_NAFImask_lag_1) <- as.POSIXct(lubridate::ymd("2018-02-01") + months(0:23))

nbr_NAFImask_lag_2 <- rep(nbr_NAFImask)
terra::time(nbr_NAFImask_lag_2) <- as.POSIXct(lubridate::ymd("2018-03-01") + months(0:23))

nbr_NAFImask_lag_3 <- rep(nbr_NAFImask)
terra::time(nbr_NAFImask_lag_3) <- as.POSIXct(lubridate::ymd("2018-04-01") + months(0:23))



# Difference NBR ----------------------------------------------------------

# THE DIFFERENCE HAS ALREADY BEEN CALCULATED

# terra::plot(nbr[[i-1]] - nbr[[i]])

# dNBR_list <- vector(mode = "list", length = nlyr(nbr)-1)
# 
# for(i in 2:nlyr(nbr)){
#   dNBR_list[[i-1]] <- nbr[[i-1]] - nbr[[i]]
# }
# 
# dNBR <- rast(dNBR_list)
# 
# # time starts at that of the second layer (which indicates the difference from the previous layer)
# terra::time(dNBR) <- as.POSIXct(lubridate::ymd("2018-02-01") + months(0:22))
# 
# for(i in 1:nlyr(dNBR)) terra::plot(dNBR[[i]], 
#                                    col = colours, 
#                                    breaks = seq(-1,1, length.out = length(colours) + 1),
#                                    main = paste0("dNBR ", time(dNBR[[i]])))
# 
# writeRaster(dNBR, "dNBR_raster_amt_cropped_epsg32753.tif")

# Sampling values from the covariates -------------------------------------

# buffalo_all <- ssf_dat %>% make_track(id = id,
#                                     x1_,
#                                     y1_, 
#                                     t1_, 
#                                     all_cols = T,
#                                     crs = 3112) 
# 
# steps(buffalo_all, keep_cols = TRUE)

class(ssf_dat) <- c("steps_xyt", class(ssf_dat))
# attr(ssf_dat, "crs") <- 3112
# crs(ssf_dat) <- 3112

buffalo_all <- ssf_dat %>% 
  extract_covariates_var_time(nbr,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_lag_0") %>%
  extract_covariates_var_time(nbr_lag_1,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_lag_1") %>%
  extract_covariates_var_time(nbr_lag_2,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_lag_2") %>%
  extract_covariates_var_time(nbr_lag_3,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_lag_3") %>%
  
  extract_covariates_var_time(nbr_NAFImask,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_NAFImask") %>% 
  extract_covariates_var_time(nbr_NAFImask_lag_1,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_NAFImask_lag_1") %>%
  extract_covariates_var_time(nbr_NAFImask_lag_2,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_NAFImask_lag_2") %>%
  extract_covariates_var_time(nbr_NAFImask_lag_3,
                              where = "both",
                              when = "any",
                              max_time = days(30),
                              name_covar = "nbr_NAFImask_lag_3")
  # extract_covariates_var_time(nbr_lag_4,
  #                             where = "both",
  #                             when = "any",
  #                             max_time = days(30),
  #                             name_covar = "nbr_lag_4") %>%
  # extract_covariates_var_time(nbr_lag_5,
  #                             where = "both",
  #                             when = "any",
  #                             max_time = days(30),
  #                             name_covar = "nbr_lag_5") %>%
  # extract_covariates_var_time(nbr_lag_6,
  #                             where = "both",
  #                             when = "any",
  #                             max_time = days(30),
  #                             name_covar = "nbr_lag_6") %>%
  # extract_covariates_var_time(dNBR,
  #                             where = "both",
  #                             when = "any",
  #                             max_time = days(30),
  #                             name_covar = "dNBR")
  
# buffalo_all <- buffalo_all[, -32:-42]
write_csv(buffalo_all, paste0("pheno_start_end_dat_NBR_lag_NAFImask_ssf10rs_", Sys.Date(), ".csv"))
# "pheno_start_end_dat_nbr_lag_ssf10rs_2024-03-18.csv"

beep(sound = 2)

ggplot() +
  # geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
  #            aes(x = t1_, y = nbr_end), colour = "red", alpha = 0.25) +
  geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
             aes(x = t1_, y = nbr_lag_0_end), colour = "red", alpha = 0.25) +
  geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
             aes(x = t1_, y = nbr_lag_1_end), colour = "orange", alpha = 0.25) +
  geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
             aes(x = t1_, y = nbr_lag_2_end), colour = "purple", alpha = 0.25) +
  geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
             aes(x = t1_, y = nbr_lag_3_end), colour = "pink", alpha = 0.25) +
  # geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
  #            aes(x = t1_, y = nbr_lag_4_end), colour = "green", alpha = 0.25) +
  # geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
  #            aes(x = t1_, y = nbr_lag_5_end), colour = "skyblue", alpha = 0.25) +
  # geom_point(data = buffalo_all %>% filter(id == "2005" & case_ == 1),
  #            aes(x = t1_, y = nbr_lag_6_end), colour = "blue", alpha = 0.25) +
  # geom_smooth(method = "lm") +
  labs(title = "NBR temporal covariate",
       x = "Time",
       y = "NBR") +
  theme_minimal()




