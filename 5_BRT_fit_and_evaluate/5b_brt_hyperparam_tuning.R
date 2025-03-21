### libraries ####
{library(tidyverse)
library(gbm)
library(dismo)
library(here);here <- here::here
library(ggBRT)
library(caret)
library(pROC)}

set.seed(1004)

### load data ####
#CRW data
dat_base <- readRDS(here("data/locs_brts/crw_pas/dat_base.rds")) %>% subset(select = -c(rep))
dat_do <- readRDS(here("data/locs_brts/crw_pas/dat_do.rds")) %>% subset(select = -c(rep))
dat_agi <- readRDS(here("data/locs_brts/crw_pas/dat_agi.rds")) %>% subset(select = -c(rep))

#Bkg data
back_base <- readRDS(here("data/locs_brts/bckg_pas/dat_base_back.rds"))
back_do <- readRDS(here("data/locs_brts/bckg_pas/dat_do_back.rds"))
back_agi <- readRDS(here("data/locs_brts/bckg_pas/dat_agi_back.rds"))

### optimize hyperparameters (CRW PAs) ####
#base
dat_base_nas <- na.omit(dat_base)

#do
dat_do_nas <- na.omit(dat_do)

#agi
dat_agi_nas <- na.omit(dat_agi)

# Set optimization options using caret
fitControl <- trainControl(method = "cv", number = 5) # Can be very slow with high "number"

# Set the range of options for each parameter: I'm varying interaction.depth (tree complexity) and shrinkage (learning rate). Range of values based on BRT paper from Elith et al (2008).
gbmGrid <- expand.grid(interaction.depth = seq(2, 5, by = 1), 
                       n.trees = 8000, #general number of trees models I ran for poster used
                       shrinkage = c(0.001, 0.01, 0.05, 0.1), 
                       n.minobsinnode = 10)

# Now test which combination of parameters works best. For some reason, the presence/absence variable must be a factor. Took ~30 min to run. 
#base model
gbmFit_base <- caret::train(factor(PA) ~ chl_mean + temp_mean + sal_mean + uo_mean + uostr_mean + vo_mean + vostr_mean + ssh_mean + mld_mean + bathy_mean + bathy_sd, 
                        data = dat_base_nas, 
                        method = "gbm", 
                        trControl = fitControl, 
                        verbose = FALSE, 
                        tuneGrid = gbmGrid)
#saveRDS(gbmFit1, here("data/brt/hp_tuning/gbmFit_base.rds"))
plot(gbmFit_base)

#do model
gbmFit_do <- caret::train(factor(PA) ~ o2_mean_0m + o2_mean_60m + o2_mean_250m + chl_mean + temp_mean + sal_mean + uo_mean + uostr_mean + vo_mean + vostr_mean + ssh_mean + mld_mean + bathy_mean + bathy_sd, 
                          data = dat_do_nas, 
                          method = "gbm", 
                          trControl = fitControl, 
                          verbose = FALSE, 
                          tuneGrid = gbmGrid)
#saveRDS(gbmFit_do, here("data/brt/hp_tuning/gbmFit_do.rds"))
plot(gbmFit_do)

#agi model
gbmFit_agi <- caret::train(factor(PA) ~ AGI_0m + AGI_60m + AGI_250m + chl_mean + temp_mean + sal_mean + uo_mean + uostr_mean + vo_mean + vostr_mean + ssh_mean + mld_mean + bathy_mean + bathy_sd, 
                           data = dat_agi_nas, 
                           method = "gbm", 
                           trControl = fitControl, 
                           verbose = FALSE, 
                           tuneGrid = gbmGrid)

#saveRDS(gbmFit_agi, here("data/brt/hp_tuning/gbmFit_agi.rds"))
plot(gbmFit_agi)

# You can plot the results: ideally, a shrinkage value (learning rate) somewhere in the middle of the options you provided will be chosen, otherwise you may need to expand the range
plot(gbmFit1)

# Save the best values for learning rate and tree complexity. I will ultimately use gbm.step to select the best number of trees. 
lr.best_base <- gbmFit1$bestTune$shrinkage #0.05 (~ 0.80)
tc.best_base <- gbmFit1$bestTune$interaction.depth #either 3 or 4 -- within 0.2 accuracy reported by CV

lr.best_do <- gbmFit_do$bestTune$shrinkage #0.01 or 0.05 -- within 0.005 accuracy reported by CV (~ 0.81)
tc.best_do <- gbmFit_do$bestTune$interaction.depth # 3 or 4 trees doesn't make a difference -- within 0.005 accuracy reported by CV

lr.best_agi <- gbmFit_agi$bestTune$shrinkage #0.01 or 0.05, 0.05 somewhat better for 3 trees and 0.01 somewhat better for 4 trees. 
tc.best_agi <- gbmFit_agi$bestTune$interaction.depth #3 and 4 trees look about the same in terms of accuracy as reported by cv. All close to ~ 0.08

# Based on above results, for the CRW PA models, I should use a lr of 0.05, tree complexity of 3, bag fraction of 0.75, and model that selects the optimal number of trees rather than a set number. 

### optimize function ####
#load data 
#plot 75/25
#run caret 
#plot
#save plots 

# Set optimization options using caret
fitControl <- trainControl(method = "cv", number = 5) # Can be very slow with high "number"

# Set the range of options for each parameter: I'm varying interaction.depth (tree complexity) and shrinkage (learning rate). Range of values based on BRT paper from Elith et al (2008).
gbmGrid <- expand.grid(interaction.depth = seq(2, 5, by = 1), 
                       n.trees = 8000, #general number of trees models I ran for poster used
                       shrinkage = c(0.001, 0.01, 0.05, 0.1), 
                       n.minobsinnode = 10)

hyperparam_tune <- function(fitControl, gbmGrid, base_input, do_input, agi_input, hp_file_dest, res = c("ann", "seas"), pa = c("back", "crw")){
  
  #load dat -- need to add subset part back in for crw files
  base_dat <- readRDS(here(base_input)) #%>% subset(select = -c(rep))
  do_dat <- readRDS(here(do_input)) #%>% subset(select = -c(rep))
  agi_dat <- readRDS(here(agi_input)) #%>% subset(select = -c(rep))
  
  #optimize base HPs
  gbmFit_base <- caret::train(factor(PA) ~ chl_mean + temp_mean + sal_mean + uo_mean + uostr_mean + vo_mean + vostr_mean + ssh_mean + mld_mean + bathy_mean + bathy_sd, 
                              data = base_dat, 
                              method = "gbm", 
                              trControl = fitControl, 
                              verbose = FALSE, 
                              tuneGrid = gbmGrid)
  
  #save and plot base file
  saveRDS(gbmFit_base, here(paste0(hp_file_dest,"/gbmFit_base","_",res,"_", pa, ".rds")))
  plot(gbmFit_base)
  
  #optimize do HPs
  gbmFit_do <- caret::train(factor(PA) ~ o2_mean_0m + o2_mean_60m + o2_mean_250m + chl_mean + temp_mean + sal_mean + uo_mean + uostr_mean + vo_mean + vostr_mean + ssh_mean + mld_mean + bathy_mean + bathy_sd, 
                              data = do_dat, 
                              method = "gbm", 
                              trControl = fitControl, 
                              verbose = FALSE, 
                              tuneGrid = gbmGrid)
  
  #save and plot do file
  saveRDS(gbmFit_do, here(paste0(hp_file_dest,"/gbmFit_do","_",res,"_", pa, ".rds")))
  plot(gbmFit_do)
  
  #optimize agi HPs
  gbmFit_agi <- caret::train(factor(PA) ~ AGI_0m + AGI_60m + AGI_250m + chl_mean + temp_mean + sal_mean + uo_mean + uostr_mean + vo_mean + vostr_mean + ssh_mean + mld_mean + bathy_mean + bathy_sd, 
                            data = agi_dat, 
                            method = "gbm", 
                            trControl = fitControl, 
                            verbose = FALSE, 
                            tuneGrid = gbmGrid)
  
  #save and plot agi file
  saveRDS(gbmFit_agi, here(paste0(hp_file_dest,"/gbmFit_agi","_",res,"_", pa, ".rds")))
  plot(gbmFit_agi)

}

#### crw annual ####
hyperparam_tune(fitControl = fitControl, gbmGrid = gbmGrid, 
                base_input = "data/locs_brts/crw_pas_ann/dat_base_ann.rds", 
                do_input = "data/locs_brts/crw_pas_ann/dat_do_ann.rds",
                agi_input = "data/locs_brts/crw_pas_ann/dat_agi_ann.rds",
                hp_file_dest = "data/brt/hp_tuning/annual", 
                res = "ann", pa = "crw")

readRDS(here("data/brt/hp_tuning/annual/gbmFit_base_ann_crw.rds")) %>% plot()
readRDS(here("data/brt/hp_tuning/annual/gbmFit_do_ann_crw.rds")) %>% plot()
readRDS(here("data/brt/hp_tuning/annual/gbmFit_agi_ann_crw.rds")) %>% plot()

#lr notes: 0.05 marginally better than 0.01
#tc notees: 3 or 4

#### crw seasonal ####
hyperparam_tune(fitControl = fitControl, gbmGrid = gbmGrid, 
                base_input = "data/locs_brts/crw_pas_seas/dat_base_seas.rds", 
                do_input = "data/locs_brts/crw_pas_seas/dat_do_seas.rds",
                agi_input = "data/locs_brts/crw_pas_seas/dat_agi_seas.rds",
                hp_file_dest = "data/brt/hp_tuning/seasonal", 
                res = "seas", pa = "crw")

readRDS(here("data/brt/hp_tuning/seasonal/gbmFit_base_seas_crw.rds")) %>% plot()
readRDS(here("data/brt/hp_tuning/seasonal/gbmFit_do_seas_crw.rds")) %>% plot()
readRDS(here("data/brt/hp_tuning/seasonal/gbmFit_agi_seas_crw.rds")) %>% plot()

#lr notes: 0.01 marginally better than 0.05
#tc notees: 3 
