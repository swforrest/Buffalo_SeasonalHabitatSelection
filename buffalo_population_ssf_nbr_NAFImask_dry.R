###########################
#' ## SSF model at the population level
# Created by Kyana Pike as part of a postgraduate internship for CSIRO 
# Feb-Apr 2022
###########################

library(tidyverse)
library(TwoStepCLogit) # to make population estimates
library(sjPlot) # to plot model estimates
library(glmmTMB)# for the mixed models
library(tictoc)
library(beepr)

# Read in annotated available data for SSF modeling
# ssf_dat <- read_csv("ssfdat_filtered_id_2024-01-24.csv")
# ssf_dat <- read_csv("pheno_dat_ssf10rs_2024-02-07.csv")
# ssf_dat <- read_csv("pheno_start_end_dat_nbr_lag_ssf10rs_2024-03-18.csv")
ssf_dat <- read_csv("pheno_start_end_dat_NBR_lag_NAFImask_ssf10rs_2024-03-21.csv")

ssf_dat <- ssf_dat %>% mutate(srtm_start_scaled = scale(srtm_start),
                              srtm_end_scaled = scale(srtm_end),
                              water_dist_start_scaled = scale(water_dist_start),
                              water_dist_end_scaled = scale(water_dist_end),
                              ndvi_start_scaled = scale(ndvi_start),
                              ndvi_end_scaled = scale(ndvi_end),
                              nbr_start_scaled = scale(nbr_start),
                              nbr_end_scaled = scale(nbr_end),
                              nbr_NAFImask_start_scaled = scale(nbr_NAFImask_start),
                              nbr_NAFImask_end_scaled = scale(nbr_NAFImask_end),
                              nbr_NAFImask_lag_1_start_scaled = scale(nbr_NAFImask_lag_1_start),
                              nbr_NAFImask_lag_1_end_scaled = scale(nbr_NAFImask_lag_1_end),
                              nbr_NAFImask_lag_2_start_scaled = scale(nbr_NAFImask_lag_2_start),
                              nbr_NAFImask_lag_2_end_scaled = scale(nbr_NAFImask_lag_2_end),
                              nbr_NAFImask_lag_3_start_scaled = scale(nbr_NAFImask_lag_3_start),
                              nbr_NAFImask_lag_3_end_scaled = scale(nbr_NAFImask_lag_3_end),
                              sl_scaled = scale(sl_),
                              log_sl_scaled = scale(log_sl_),
                              cos_ta_scaled = scale(cos_ta_)
)

# separate the buffalo by season 
ssfdat_dry <- ssf_dat %>% filter(season=="dry") # dry season


################ Creating seasonal mixed effects models ########################

# variables with the '_' at the end indicate they have been scaled and centered, this is true also for sl_1
## for the dry season

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          nbr_end_scaled +
#                          water_dist_end_scaled +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + nbr_end_scaled | id) +
#                          (0 + water_dist_end_scaled | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# buffalo.ssf.dry <- glmmTMB:::fitTMB(buffalo.tmp)
# 
# summary(buffalo.ssf.dry)
# plot_model(buffalo.ssf.dry)
# 
# saveRDS(buffalo.ssf.dry, file = paste0("buffalo_ssf10rs_dry_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry)


################ dry season with quadratics and no NBR ########################

buffalo.tmp <- glmmTMB(case_ ~
                         -1 +
                         ndvi_end_scaled +
                         I(ndvi_end_scaled^2) +
                         water_dist_end_scaled +
                         sl_scaled +
                         log_sl_scaled +
                         cos_ta_scaled +
                         (1 | step_id) +
                         (0 + ndvi_end_scaled | id) +
                         (0 + I(ndvi_end_scaled^2) | id) +
                         (0 + water_dist_end_scaled | id),
                       family = poisson(),
                       data = ssfdat_dry,
                       doFit=FALSE)

# Set variance of random intercept to 10^6
buffalo.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(buffalo.tmp$parameters$theta)
buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
tic()
buffalo.ssf.dry_NOnbr <- glmmTMB:::fitTMB(buffalo.tmp)
toc()

summary(buffalo.ssf.dry_NOnbr)
plot_model(buffalo.ssf.dry_NOnbr)

saveRDS(buffalo.ssf.dry_NOnbr, file = paste0("buffalo_ssf10rs_dry_noNBR_", Sys.Date(), ".rds"))
rm(buffalo.ssf.dry_NOnbr)


################ dry season with quadratics and NBR lag ########################

# load("buffalo_ssf10rs_dry_nbr_NAFImask_lag0_2024-03-18.Rdata")

buffalo.tmp <- glmmTMB(case_ ~
                         -1 +
                         ndvi_end_scaled +
                         I(ndvi_end_scaled^2) +
                         nbr_NAFImask_end_scaled +
                         water_dist_end_scaled +
                         sl_scaled +
                         log_sl_scaled +
                         cos_ta_scaled +
                         (1 | step_id) +
                         (0 + ndvi_end_scaled | id) +
                         (0 + I(ndvi_end_scaled^2) | id) +
                         (0 + nbr_NAFImask_end_scaled | id) +
                         (0 + water_dist_end_scaled | id),
                       family = poisson(),
                       data = ssfdat_dry,
                       doFit=FALSE)

# Set variance of random intercept to 10^6
buffalo.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(buffalo.tmp$parameters$theta)
buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
tic()
buffalo.ssf.dry_nbrNAFI0 <- glmmTMB:::fitTMB(buffalo.tmp)
toc()

summary(buffalo.ssf.dry_nbrNAFI0)
plot_model(buffalo.ssf.dry_nbrNAFI0)

saveRDS(buffalo.ssf.dry_nbrNAFI0, file = paste0("buffalo_ssf10rs_dry_nbrNAFI_lag0_", Sys.Date(), ".rds"))
rm(buffalo.ssf.dry_nbrNAFI0)


################ dry season with quadratics and NBR lag - orthogonal polynomials ########################

# load("buffalo_ssf10rs_dry_OP_nbr_NAFImask_lag0_2024-03-21.Rdata")

# sum(is.na(ssfdat_dry$ndvi_end_scaled))
# ssfdat_dry <- ssfdat_dry %>% filter(!is.na(ndvi_end_scaled))

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          stats::poly(ndvi_end_scaled, 2, raw = FALSE) +
#                          # ndvi_end_scaled +
#                          # I(ndvi_end_scaled^2) +
#                          nbr_NAFImask_lag_0_end_scaled +
#                          water_dist_end_scaled +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + stats::poly(ndvi_end_scaled, 2, raw = FALSE) | id) +
#                          # (0 + ndvi_end_scaled | id) +
#                          # (0 + I(ndvi_end_scaled^2) | id) +
#                          (0 + nbr_NAFImask_lag_0_end_scaled | id) +
#                          (0 + water_dist_end_scaled | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# tic()
# buffalo.ssf.dry_OP_nbr0 <- glmmTMB:::fitTMB(buffalo.tmp)
# toc()
# 
# summary(buffalo.ssf.dry_OP_nbr0)
# plot_model(buffalo.ssf.dry_OP_nbr0)
# 
# saveRDS(buffalo.ssf.dry_OP_nbr0, file = paste0("buffalo_ssf10rs_dry_OP_nbr_NAFImask_lag0_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_OP_nbr0)


################ dry season with quadratics and NBR lag ########################

# load("buffalo_ssf10rs_dry_nbr_NAFImask_lag1_2024-03-18.Rdata")

buffalo.tmp <- glmmTMB(case_ ~
                         -1 +
                         ndvi_end_scaled +
                         I(ndvi_end_scaled^2) +
                         nbr_NAFImask_lag_1_end_scaled +
                         water_dist_end_scaled +
                         sl_scaled +
                         log_sl_scaled +
                         cos_ta_scaled +
                         (1 | step_id) +
                         (0 + ndvi_end_scaled | id) +
                         (0 + I(ndvi_end_scaled^2) | id) +
                         (0 + nbr_NAFImask_lag_1_end_scaled | id) +
                         (0 + water_dist_end_scaled | id),
                       family = poisson(),
                       data = ssfdat_dry,
                       doFit=FALSE)

# Set variance of random intercept to 10^6
buffalo.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(buffalo.tmp$parameters$theta)
buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
tic()
buffalo.ssf.dry_nbrNAFI1 <- glmmTMB:::fitTMB(buffalo.tmp)
toc()

summary(buffalo.ssf.dry_nbrNAFI1)
plot_model(buffalo.ssf.dry_nbrNAFI1)

saveRDS(buffalo.ssf.dry_nbrNAFI1, file = paste0("buffalo_ssf10rs_dry_nbr_NAFImask_lag1_", Sys.Date(), ".rds"))
rm(buffalo.ssf.dry_nbrNAFI1)


################ dry season with quadratics and NBR lag ########################

# load("buffalo_ssf10rs_dry_nbr_NAFImask_lag2_2024-03-18.Rdata")

buffalo.tmp <- glmmTMB(case_ ~
                         -1 +
                         ndvi_end_scaled +
                         I(ndvi_end_scaled^2) +
                         nbr_NAFImask_lag_2_end_scaled +
                         water_dist_end_scaled +
                         sl_scaled +
                         log_sl_scaled +
                         cos_ta_scaled +
                         (1 | step_id) +
                         (0 + ndvi_end_scaled | id) +
                         (0 + I(ndvi_end_scaled^2) | id) +
                         (0 + nbr_NAFImask_lag_2_end_scaled | id) +
                         (0 + water_dist_end_scaled | id),
                       family = poisson(),
                       data = ssfdat_dry,
                       doFit=FALSE)

# Set variance of random intercept to 10^6
buffalo.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(buffalo.tmp$parameters$theta)
buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
tic()
buffalo.ssf.dry_nbrNAFI2 <- glmmTMB:::fitTMB(buffalo.tmp)
toc()

summary(buffalo.ssf.dry_nbrNAFI2)
plot_model(buffalo.ssf.dry_nbrNAFI2)

saveRDS(buffalo.ssf.dry_nbrNAFI2, file = paste0("buffalo_ssf10rs_dry_nbr_NAFImask_lag2_", Sys.Date(), ".rds"))
rm(buffalo.ssf.dry_nbrNAFI2)


################ dry season with quadratics and NBR lag ########################

# load("buffalo_ssf10rs_dry_nbr_NAFImask_lag3_2024-03-18.Rdata")

buffalo.tmp <- glmmTMB(case_ ~
                         -1 +
                         ndvi_end_scaled +
                         I(ndvi_end_scaled^2) +
                         nbr_NAFImask_lag_3_end_scaled +
                         water_dist_end_scaled +
                         sl_scaled +
                         log_sl_scaled +
                         cos_ta_scaled +
                         (1 | step_id) +
                         (0 + ndvi_end_scaled | id) +
                         (0 + I(ndvi_end_scaled^2) | id) +
                         (0 + nbr_NAFImask_lag_3_end_scaled | id) +
                         (0 + water_dist_end_scaled | id),
                       family = poisson(),
                       data = ssfdat_dry,
                       doFit=FALSE)

# Set variance of random intercept to 10^6
buffalo.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(buffalo.tmp$parameters$theta)
buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
tic()
buffalo.ssf.dry_nbrNAFI3 <- glmmTMB:::fitTMB(buffalo.tmp)
toc()

summary(buffalo.ssf.dry_nbrNAFI3)
plot_model(buffalo.ssf.dry_nbrNAFI3)

saveRDS(buffalo.ssf.dry_nbrNAFI3, file = paste0("buffalo_ssf10rs_dry_nbr_NAFImask_lag3_", Sys.Date(), ".rds"))
rm(buffalo.ssf.dry_nbrNAFI3)


################ dry season with quadratics and NBR lag ########################

# load("buffalo_ssf10rs_dry_nbr_NAFImask_lag4_2024-03-21.Rdata")

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          I(ndvi_end_scaled^2) +
#                          nbr_NAFImask_lag_4_end_scaled +
#                          water_dist_end_scaled +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + I(ndvi_end_scaled^2) | id) +
#                          (0 + nbr_NAFImask_lag_4_end_scaled | id) +
#                          (0 + water_dist_end_scaled | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# tic()
# buffalo.ssf.dry_nbr4 <- glmmTMB:::fitTMB(buffalo.tmp)
# toc()
# 
# summary(buffalo.ssf.dry_nbr4)
# plot_model(buffalo.ssf.dry_nbr4)
# 
# saveRDS(buffalo.ssf.dry_nbr4, file = paste0("buffalo_ssf10rs_dry_nbr_NAFImask_lag4_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_nbr4)


################ dry season with quadratics and NBR lag ########################

# load("buffalo_ssf10rs_dry_nbr_NAFImask_lag5_2024-03-21.Rdata")
# 
# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          I(ndvi_end_scaled^2) +
#                          nbr_NAFImask_lag_5_end_scaled +
#                          water_dist_end_scaled +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + I(ndvi_end_scaled^2) | id) +
#                          (0 + nbr_NAFImask_lag_5_end_scaled | id) +
#                          (0 + water_dist_end_scaled | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# tic()
# buffalo.ssf.dry_nbr5 <- glmmTMB:::fitTMB(buffalo.tmp)
# toc()
# 
# summary(buffalo.ssf.dry_nbr5)
# plot_model(buffalo.ssf.dry_nbr5)
# 
# saveRDS(buffalo.ssf.dry_nbr5, file = paste0("buffalo_ssf10rs_dry_nbr_NAFImask_lag5_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_nbr5)


################ dry season with quadratics and NBR lag ########################

# load("buffalo_ssf10rs_dry_nbr_NAFImask_lag6_2024-03-21.Rdata")
# 
# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          I(ndvi_end_scaled^2) +
#                          nbr_NAFImask_lag_6_end_scaled +
#                          water_dist_end_scaled +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + I(ndvi_end_scaled^2) | id) +
#                          (0 + nbr_NAFImask_lag_6_end_scaled | id) +
#                          (0 + water_dist_end_scaled | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# tic()
# buffalo.ssf.dry_nbr6 <- glmmTMB:::fitTMB(buffalo.tmp)
# toc()
# 
# summary(buffalo.ssf.dry_nbr6)
# plot_model(buffalo.ssf.dry_nbr6)
# 
# saveRDS(buffalo.ssf.dry_nbr6, file = paste0("buffalo_ssf10rs_dry_nbr_NAFImask_lag6_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_nbr6)


################ dry season with quadratics ########################

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          I(ndvi_end_scaled^2) +
#                          nbr_end_scaled +
#                          I(nbr_end_scaled^2) +
#                          water_dist_end_scaled +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + I(ndvi_end_scaled^2) | id) +
#                          (0 + nbr_end_scaled | id) +
#                          (0 + I(nbr_end_scaled^2) | id) +
#                          (0 + water_dist_end_scaled | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# buffalo.ssf.dry_ndvi2_nbr2 <- glmmTMB:::fitTMB(buffalo.tmp)
# 
# summary(buffalo.ssf.dry_ndvi2_nbr2)
# plot_model(buffalo.ssf.dry_ndvi2_nbr2)
# 
# saveRDS(buffalo.ssf.dry_ndvi2_nbr2, file = paste0("buffalo_ssf10rs_dry_ndvi2_nbr2_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_ndvi2_nbr2)


################ dry season with quadratics ########################

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          I(ndvi_end_scaled^2) +
#                          nbr_end_scaled +
#                          I(nbr_end_scaled^2) +
#                          water_dist_end_scaled +
#                          I(water_dist_end_scaled^2) +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + I(ndvi_end_scaled^2) | id) +
#                          (0 + nbr_end_scaled | id) +
#                          (0 + I(nbr_end_scaled^2) | id) +
#                          (0 + water_dist_end_scaled | id) +
#                          (0 + I(water_dist_end_scaled^2) | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# buffalo.ssf.dry_ndvi2_nbr2_wdist2 <- glmmTMB:::fitTMB(buffalo.tmp)
# 
# summary(buffalo.ssf.dry_ndvi2_nbr2_wdist2)
# plot_model(buffalo.ssf.dry_ndvi2_nbr2_wdist2)
# 
# saveRDS(buffalo.ssf.dry_ndvi2_nbr2_wdist2, file = paste0("buffalo_ssf10rs_dry_ndvi2_nbr2_wdist2_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_ndvi2_nbr2_wdist2)


################ dry season with quadratics ########################

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          I(ndvi_end_scaled^2) +
#                          nbr_end_scaled +
#                          water_dist_end_scaled +
#                          I(water_dist_end_scaled^2) +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + I(ndvi_end_scaled^2) | id) +
#                          (0 + nbr_end_scaled | id) +
#                          (0 + water_dist_end_scaled | id) +
#                          (0 + I(water_dist_end_scaled^2) | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# buffalo.ssf.dry_ndvi2_wdist2 <- glmmTMB:::fitTMB(buffalo.tmp)
# 
# summary(buffalo.ssf.dry_ndvi2_wdist2)
# plot_model(buffalo.ssf.dry_ndvi2_wdist2)
# 
# saveRDS(buffalo.ssf.dry_ndvi2_wdist2, file = paste0("buffalo_ssf10rs_dry_ndvi2_wdist2_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_ndvi2_wdist2)


################ dry season with quadratics ########################

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          nbr_end_scaled +
#                          I(nbr_end_scaled^2) +
#                          water_dist_end_scaled +
#                          I(water_dist_end_scaled^2) +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + nbr_end_scaled | id) +
#                          (0 + I(nbr_end_scaled^2) | id) +
#                          (0 + water_dist_end_scaled | id) +
#                          (0 + I(water_dist_end_scaled^2) | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# buffalo.ssf.dry_nbr2_wdist2 <- glmmTMB:::fitTMB(buffalo.tmp)
# 
# summary(buffalo.ssf.dry_nbr2_wdist2)
# plot_model(buffalo.ssf.dry_nbr2_wdist2)
# 
# saveRDS(buffalo.ssf.dry_nbr2_wdist2, file = paste0("buffalo_ssf10rs_dry_nbr2_wdist2_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_nbr2_wdist2)


################ dry season with quadratics ########################

# buffalo.tmp <- glmmTMB(case_ ~
#                          -1 +
#                          ndvi_end_scaled +
#                          nbr_end_scaled +
#                          water_dist_end_scaled +
#                          I(water_dist_end_scaled^2) +
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) +
#                          (0 + ndvi_end_scaled | id) +
#                          (0 + nbr_end_scaled | id) +
#                          (0 + water_dist_end_scaled | id) +
#                          (0 + I(water_dist_end_scaled^2) | id),
#                        family = poisson(),
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# buffalo.ssf.dry_wdist2 <- glmmTMB:::fitTMB(buffalo.tmp)
# 
# summary(buffalo.ssf.dry_wdist2)
# plot_model(buffalo.ssf.dry_wdist2)
# 
# saveRDS(buffalo.ssf.dry_wdist2, file = paste0("buffalo_ssf10rs_dry_wdist2_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_wdist2)


################ dry season with quadratics ########################

# buffalo.tmp <- glmmTMB(case_ ~ 
#                          -1 +
#                          ndvi_end_scaled + 
#                          nbr_end_scaled +  
#                          I(nbr_end_scaled^2) + 
#                          water_dist_end_scaled + 
#                          sl_scaled +
#                          log_sl_scaled +
#                          cos_ta_scaled +
#                          (1 | step_id) + 
#                          (0 + ndvi_end_scaled | id) +  
#                          (0 + nbr_end_scaled | id) + 
#                          (0 + I(nbr_end_scaled^2) | id) +
#                          (0 + water_dist_end_scaled | id), 
#                        family = poisson(), 
#                        data = ssfdat_dry,
#                        doFit=FALSE)
# 
# # Set variance of random intercept to 10^6
# buffalo.tmp$parameters$theta[1] <- log(1e3)
# nvarparm<-length(buffalo.tmp$parameters$theta)
# buffalo.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
# buffalo.ssf.dry_nbr2 <- glmmTMB:::fitTMB(buffalo.tmp) 
# 
# summary(buffalo.ssf.dry_nbr2)
# plot_model(buffalo.ssf.dry_nbr2)
# 
# saveRDS(buffalo.ssf.dry_nbr2, file = paste0("buffalo_ssf10rs_nbr2_", Sys.Date(), ".rds"))
# rm(buffalo.ssf.dry_nbr2)