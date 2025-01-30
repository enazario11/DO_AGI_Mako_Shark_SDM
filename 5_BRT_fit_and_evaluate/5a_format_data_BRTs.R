# libraries ####
library(tidyverse)
library(here)

set.seed(1004)

# CRW load data w/ covars ####
#combine all before PA selection and fix PA values before fitting BRTs
dat0 <- readRDS(here("data/locs_w_covar/psat_spot/cmem_locs_covar_AGI_0m.rds"))
dat0$bathy <- replace(dat0$bathy, dat0$bathy == "NaN", NA)
colnames(dat0) <- c("tag", "date", "lon", "lat", "PA", "rep", "dt", "o2_mean_0m", "chl_mean", "temp_mean", "sal_mean", "uo_mean", "uostr_mean", "vo_mean", "vostr_mean", "ssh_mean", "mld_mean", "bathy_mean", "bathy_sd", "pO2_atm_0m", "o2_demand_0m", "AGI_0m")
head(dat0)

dat60 <- readRDS(here("data/locs_w_covar/psat_spot/cmem_locs_covar_AGI_60m.rds"))
colnames(dat60) <- c("tag", "date", "lon", "lat", "PA", "rep", "dt", "o2_mean_60m", "temp_mean", "sal_mean", "pO2_atm_60m", "o2_demand_60m", "AGI_60m")

dat250 <- readRDS(here("data/locs_w_covar/psat_spot/cmem_locs_covar_AGI_250m.rds"))
colnames(dat250) <- c("tag", "date", "lon", "lat", "PA", "rep", "dt", "o2_mean_250m", "temp_mean", "sal_mean", "pO2_atm_250m", "o2_demand_250m", "AGI_250m")

dat_all_temp <- dat0 %>% cbind(dat60$o2_mean_60m,dat250$o2_mean_250m, dat60$AGI_60m, dat250$AGI_250m)

#originally, PA = 0 means a true position. Change so PA = 1 means a true position for fitting the BRT
dat_all_temp$PA <- replace(dat_all_temp$PA, dat_all_temp$PA == 1, 2) #change PAs to temporarily equal 2
dat_all_temp$PA <- replace(dat_all_temp$PA, dat_all_temp$PA == 0, 1) #change true positions to a 1
dat_all_temp$PA <- replace(dat_all_temp$PA, dat_all_temp$PA == 2, 0) #change PA positions to a 0

## randomly select one PA for each tag (CRW dataset only) 
#randomly select 1 PA rep for each tag 
dat_pos <- dat_all_temp %>% filter(PA == 1)
dat_pa <- dat_all_temp %>% filter(PA == 0)

dat_temp <- NULL
for(i in 1:length(unique(dat_pa$tag))){
  #select current id
  curr_ID <- unique(dat_pa$tag)[i]
  temp_df <- dat_pa[dat_pa$tag %in% curr_ID,]
  
  #sample id's randomly
  temp_rep_ID <- sample(unique(temp_df$rep), 1, replace = FALSE)
  
  #narrow your data set
  temp_df2 <- temp_df[temp_df$rep %in% temp_rep_ID, ]
  
  #combine in a single df
  dat_temp <- rbind(dat_temp, temp_df2)
  
}

dat_all_temp2 <- rbind(dat_temp, dat_pos)
dat_all_temp3 <- dat_all_temp2 %>% filter(tag != 96365)
dat_all_temp_fix <- dat_all_temp2 %>% 
  filter(tag == 96365) %>%
  group_by(PA) %>%
  distinct(dt, .keep_all = TRUE) %>%
  ungroup()

dat_all <- rbind(dat_all_temp3, dat_all_temp_fix) %>%
  rename(o2_mean_60m = "dat60$o2_mean_60m", 
         o2_mean_250m = "dat250$o2_mean_250m", 
         AGI_60m = "dat60$AGI_60m", 
         AGI_250m = "dat250$AGI_250m")

#combine depth specific data at 60m and 250m for DO and AGI. Erase for base model dataset. Here, I also change the PA values.
#base model dataset
dat_base <- dat_all %>% 
  subset(select = -c(o2_mean_0m, pO2_atm_0m, o2_demand_0m, AGI_0m, o2_mean_60m, o2_mean_250m, AGI_60m, AGI_250m))
#saveRDS(dat_base, here("data/locs_brts/crw_pas/dat_base.rds"))

#DO model dataset
dat_do <- dat_all %>% 
  subset(select = -c(pO2_atm_0m, o2_demand_0m, AGI_0m, AGI_60m, AGI_250m))
#saveRDS(dat_do, here("data/locs_brts/crw_pas/dat_do.rds"))

#AGI model dataset 
dat_agi <- dat_all %>% 
  subset(select = -c(o2_mean_0m, pO2_atm_0m, o2_demand_0m, o2_mean_60m, o2_mean_250m)) 
#saveRDS(dat_agi, here("data/locs_brts/crw_pas/dat_agi.rds"))

# CRW Annual ####
format_dat_crw_brts <- function(dat0_path, dat60_path, dat250_path, res = c("ann", "seas"), out_path){
  
  set.seed(1004)
  
  #combine all before PA selection and fix PA values before fitting BRTs
  dat0 <- readRDS(here(dat0_path))
  dat0$bathy <- replace(dat0$bathy, dat0$bathy == "NaN", NA)
  colnames(dat0) <- c("tag", "date", "lon", "lat", "PA", "rep", "dt", "dt_ann", "mld_mean", "sal_mean", "ssh_mean", "temp_mean",  "uo_mean", "uostr_mean", "vo_mean", "vostr_mean", "o2_mean_0m", "chl_mean",   "bathy_mean", "bathy_sd", "pO2_atm_0m", "o2_demand_0m", "AGI_0m", "dist_coast")
                      
  dat60 <- readRDS(here(dat60_path))
  colnames(dat60) <- c("tag", "date", "lon", "lat", "PA", "rep", "dt", "dt_ann", "sal_mean", "temp_mean", "o2_mean_60m","pO2_atm_60m", "o2_demand_60m", "AGI_60m" )
  
  dat250 <- readRDS(here(dat250_path))
  colnames(dat250) <- c("tag", "date", "lon", "lat", "PA", "rep", "dt", "dt_ann", "sal_mean", "temp_mean", "o2_mean_250m","pO2_atm_250m", "o2_demand_250m", "AGI_250m" )
  
  dat_all_temp <- dat0 %>% cbind(dat60$o2_mean_60m,dat250$o2_mean_250m, dat60$AGI_60m, dat250$AGI_250m)
  
  #originally, PA = 0 means a true position. Change so PA = 1 means a true position for fitting the BRT
  dat_all_temp$PA <- replace(dat_all_temp$PA, dat_all_temp$PA == 1, 2) #change PAs to temporarily equal 2
  dat_all_temp$PA <- replace(dat_all_temp$PA, dat_all_temp$PA == 0, 1) #change true positions to a 1
  dat_all_temp$PA <- replace(dat_all_temp$PA, dat_all_temp$PA == 2, 0) #change PA positions to a 0
  
  
  
  # randomly select one PA for each tag (CRW dataset only)
  #randomly select 1 PA rep for each tag
  dat_pos <- dat_all_temp %>% filter(PA == 1)
  dat_pa <- dat_all_temp %>% filter(PA == 0)

  dat_temp <- NULL
  for(i in 1:length(unique(dat_pa$tag))){
    #select current id
    curr_ID <- unique(dat_pa$tag)[i]
    temp_df <- dat_pa[dat_pa$tag %in% curr_ID,]

    #sample id's randomly
    temp_rep_ID <- sample(unique(temp_df$rep), 1, replace = FALSE)

    #narrow your data set
    temp_df2 <- temp_df[temp_df$rep %in% temp_rep_ID, ]

    #combine in a single df
    dat_temp <- rbind(dat_temp, temp_df2)
  }
  
  dat_all_temp2 <- rbind(dat_temp, dat_pos)
  dat_all_temp3 <- dat_all_temp2 %>% filter(tag != 96365)
  dat_all_temp_fix <- dat_all_temp2 %>% 
    filter(tag == 96365) %>%
    group_by(PA) %>%
    distinct(dt, .keep_all = TRUE) %>%
    ungroup()
  
  dat_all <- rbind(dat_all_temp3, dat_all_temp_fix) %>%
    rename(o2_mean_60m = "dat60$o2_mean_60m", 
           o2_mean_250m = "dat250$o2_mean_250m", 
           AGI_60m = "dat60$AGI_60m", 
           AGI_250m = "dat250$AGI_250m")
  
  
  #combine depth specific data at 60m and 250m for DO and AGI. Erase for base model dataset. Here, I also change the PA values.
  #base model dataset
  dat_base <- dat_all %>% 
    subset(select = -c(o2_mean_0m, pO2_atm_0m, o2_demand_0m, AGI_0m, o2_mean_60m, o2_mean_250m, AGI_60m, AGI_250m))
  saveRDS(dat_base, here(paste0(out_path,"/dat_base","_", res,".rds")))
  
  #DO model dataset
  dat_do <- dat_all %>% 
    subset(select = -c(pO2_atm_0m, o2_demand_0m, AGI_0m, AGI_60m, AGI_250m))
  saveRDS(dat_do, here(paste0(out_path,"/dat_do","_", res,".rds")))
  
  #AGI model dataset 
  dat_agi <- dat_all %>% 
    subset(select = -c(o2_mean_0m, pO2_atm_0m, o2_demand_0m, o2_mean_60m, o2_mean_250m)) 
  saveRDS(dat_agi, here(paste0(out_path,"/dat_agi","_", res,".rds")))
  
}

dat0_path <- "data/locs_w_covar/psat_spot/annual/cmem_locs_covar_AGI_distcoast_0m_ann.rds"
dat60_path <- "data/locs_w_covar/psat_spot/annual/cmem_locs_covar_AGI_60m_ann.rds"
dat250_path <- "data/locs_w_covar/psat_spot/annual/cmem_locs_covar_AGI_250m_ann.rds"

out_path <- "data/locs_brts/crw_pas_ann"

format_dat_crw_brts(dat0_path = dat0_path, dat60_path = dat60_path, dat250_path = dat250_path, res = "ann", out_path = out_path)

#CRW Seasonal ####
dat0_path <- "data/locs_w_covar/psat_spot/seasonal/cmem_locs_covar_AGI_distcoast_0m_seas.rds"
dat60_path <- "data/locs_w_covar/psat_spot/seasonal/cmem_locs_covar_AGI_60m_seas.rds"
dat250_path <- "data/locs_w_covar/psat_spot/seasonal/cmem_locs_covar_AGI_250m_seas.rds"

out_path <- "data/locs_brts/crw_pas_seas"

format_dat_crw_brts(dat0_path = dat0_path, dat60_path = dat60_path, dat250_path = dat250_path, res = "seas", out_path = out_path)
