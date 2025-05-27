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

# Read in annotated available data for SSF modeling
# load("~/Desktop/CSIRO-internship/Feral_movement/buffalo-issf/buftrk.Rdata")

# two animals have much smaller samples than the others so will remove them
# remove<- c("2029", # very small sample
#            "2043") # also smaller sample and strange stationary part

# ssfdat <- buftrk %>% filter(!id %in% remove)# this is all the buffalo data but with 2043 and 2029 removed as they have bad/ small samples
# ssfdat <- ssfdat[,c(1:27, 37:39)]
# write_csv(ssfdat, paste0("ssfdat_filtered_id_", Sys.Date(), ".csv"))

ssf_dat <- read_csv("ssfdat_filtered_id_2024-01-24.csv")
# ssf_dat <- read_csv("pheno_dat_ssf10rs_2024-02-07.csv")

ssf_dat <- ssf_dat %>% mutate(srtm_start_scaled = scale(srtm_start),
                            srtm_end_scaled = scale(srtm_end),
                            water_dist_start_scaled = scale(water_dist_start),
                            water_dist_end_scaled = scale(water_dist_end),
                            ndvi_start_scaled = scale(ndvi_start),
                            ndvi_end_scaled = scale(ndvi_end),
                            nbr_start_scaled = scale(nbr_start),
                            nbr_end_scaled = scale(nbr_end),
                            sl_scaled = scale(sl_),
                            log_sl_scaled = scale(log_sl_),
                            cos_ta_scaled = scale(cos_ta_)
)


ggplot(ssf_dat %>% dplyr::filter(case_ != TRUE), aes(x = pheno_end, fill = pheno_end)) + 
  geom_bar() + 
  # facet_wrap(~season) +
  theme_classic() +
  # rotate the x-axis labels 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# remove NAs 
ssf_dat <- ssf_dat %>% drop_na()

### Prepare data for phenology model as not all categories can be in the model 
bad_pheno <- c("no_data", # very small sample
               "saltpan_mudlfat",
                "swamp",
                "water",
                "shrubs",
                "sparse_veg"
                # "forest"
               ) 

# get the locations which have one or more of these categories
pheno <- ssf_dat %>% filter(pheno_end %in% bad_pheno | pheno_start %in% bad_pheno)
pheno %>% count(id, case_) #buffalo that have at least one
pheno_used <- pheno %>% filter(case_==TRUE) # buffalo that have at least one and are used

# only remove full steps that have the USED location in a bad phenology category
bad_phen <- pheno_used$step_id # value that has the strata+id that are associated with all those USED steps

pheno_used <- ssf_dat %>% filter(!step_id %in% bad_phen) #now has filtered the data so that any set of steps that 
# had one or more of the unwanted pheno categories are gone. 
#This is important so all remaining sets have the same # of locations

# check that all have the same number of steps
pheno_used %>% group_by(step_id) %>% count() 

# now filter the AVAILABLE steps that fell into a bad phenology category
pheno_dat <- pheno_used %>% filter(!pheno_end %in% bad_pheno)
pheno_dat %>% count(id, pheno_end) # check that the bad phenology categories are gone

ggplot(pheno_dat %>% filter(case_ == TRUE), aes(x = pheno_end, fill = pheno_end)) + 
  geom_bar() + 
  # facet_wrap(~season) +
  theme_classic() +
  # rotate the x-axis labels 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# check that all have the same number of steps
pheno_dat_counts <- pheno_dat %>% group_by(step_id) %>% count() 
min(pheno_dat_counts$n)
hist(pheno_dat_counts$n)

steps_less_subset_random <- pheno_dat_counts$step_id[which(pheno_dat_counts$n < 11)] # check if there are any single steps with less than 10 random steps (1 used + 10 random = 11)
pheno_dat_ssf <- pheno_dat %>% filter(!step_id %in% steps_less_subset_random)

# the minimum number of random steps should now be 11 (1 used + 10 random = 11)
pheno_dat_counts <- pheno_dat_ssf %>% group_by(step_id) %>% count() 
min(pheno_dat_counts$n)
hist(pheno_dat_counts$n)

# now slice by steps so all have only 10 random steps (1 used + 10 random = 11) - this can be adjust to a different number if needed
pheno_dat_ssf_subset <- pheno_dat_ssf %>% group_by(step_id) %>% slice_head(n = 11)

# check that all have the same number of steps
pheno_dat_ssf_subset_counts <- pheno_dat_ssf_subset %>% group_by(step_id) %>% count() 
min(pheno_dat_ssf_subset_counts$n)
max(pheno_dat_ssf_subset_counts$n)

# YOU MUST REMOVE THE BAD LEVELS OR IT WILL NOT RUN 
# pheno_dat$pheno_end <- droplevels(pheno_dat$pheno_end)
# str(pheno_dat) # double check has the correct nummber of levels 

# there shouldn't be any NAs
for(i in 1:ncol(pheno_dat_ssf_subset)) print(paste(colnames(pheno_dat_ssf_subset[,i]), sum(is.na(pheno_dat_ssf_subset[,i]))), sep = " ")

write.csv(pheno_dat_ssf_subset, paste0("pheno_forest_start_end_dat_ssf10rs_", Sys.Date(), ".csv"))
beep(sound = 2)
