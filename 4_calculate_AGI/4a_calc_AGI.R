### load packages ####
library(tidyverse)
library(here)
library(tidyquant)
library(respR)

### read data generated from environmental extraction ####
dat0 <- readRDS(here("data/locs_w_covar/psat_spot/cmem_locs_covar_0m.rds"))
dat60 <- readRDS(here("data/locs_w_covar/psat_spot/cmem_locs_covar_60m.rds"))
dat250 <- readRDS(here("data/locs_w_covar/psat_spot/cmem_locs_covar_250m.rds"))

source(here("functions/oxy_demand_functions.R"))

### convert DO to atm ####
# 0m -- USE AS MODEL
dat0_DOatm <- DO_to_atm(dat0, depth = 0)
thresh0 <- thresh_atm(temp = median(dat0_DOatm$votemper_mean, na.rm = TRUE), so_psu = median(dat0_DOatm$vosaline_mean, na.rm = TRUE), depth = 0) #defualt do value is 2 mL/L from vetter et al., 2008

hist(dat0_DOatm$pO2_0, xlim = c(0, 0.20)) 
abline(v = thresh0, lwd = 2)

# 60m
dat60_DOatm <- DO_to_atm(dat60, depth = 60)
thresh60 <- thresh_atm(temp = median(dat60_DOatm$votemper_mean, na.rm = TRUE), so_psu = median(dat60_DOatm$vosaline_mean, na.rm = TRUE), depth = 60)

hist(dat60_DOatm$pO2_60)
abline(v = thresh60, lwd = 2)

#250m 
dat250_DOatm <- DO_to_atm(dat250, depth = 250)
thresh250 <- thresh_atm(temp = median(dat250_DOatm$votemper_mean, na.rm = TRUE), so_psu = median(dat250_DOatm$vosaline_mean, na.rm = TRUE), depth = 250)

hist(dat250_DOatm$pO2_250)
abline(v = thresh250, lwd = 2)

### static constants ####
# W = 51807.63; average mass in g for juv. makos as estimated by length-weight relationship. Used average FL of 177.7 cm (from study animals)
# d = 0.700; constant, Clarke et al., 2021
# K = 0.070; adult numbers for K (VBGP) for california animals (fishbase.org). No juv info available.
# j2 = 8.000; constant, Clarke et al., 2021
# j1 = 4.500; constant, Clarke et al., 2021
# OxyThresh = 0.05347; 10th percentile from ambient O2 from data
#Tpref = median(dat$sst, na.rm = T)
# Linf = 321; adult numbers for Linf FL for california animals (fishbase.org)
# LwA = 0.01670 #juvenile numbers from fishbase (New Zealand mixed animals-- closest in age class and ocean basin)
# LwB = 2.847 #juvenile numbers from fishbase (New Zealand mixed animals-- closest in age class and ocean basin)

# Tested different values from different age class and populations that had the largest 
# sample sizes (NZ and USA) for LwA and LwB and found that range of oxygen demand values 
# at the surface were 0.75 and 0.88. Because these values were close, we chose the population 
# that was in the Pacific ocean basin and contained a mix of age class individuals. 

### mako specific constants ####
#calculate temp pref
Tpref50 = 16.45201 #50m tpref is 16.452

#run oxygen demand function with mako specific parameters
dat0_DOatm$O2_demand0 <- OxyDemand(Tpref = Tpref50, PO2_thresh = thresh0, T_C = dat0_DOatm$votemper_mean)
dat60_DOatm$O2_demand60 <- OxyDemand(Tpref = Tpref50, PO2_thresh = thresh60, T_C = dat60_DOatm$votemper_mean)
dat250_DOatm$O2_demand250 <- OxyDemand(Tpref = Tpref50, PO2_thresh = thresh250, T_C = dat250_DOatm$votemper_mean)

  #explore outputs
hist(dat0_DOatm$O2_demand0)
hist(dat60_DOatm$O2_demand60)
hist(dat250_DOatm$O2_demand250)
plot(dat0_DOatm$votemper_mean, dat0_DOatm$O2_demand0) #should increase with temp

#calculate AGI
dat0_DOatm$AGI_0m <- dat0_DOatm$pO2_0/dat0_DOatm$O2_demand0
dat60_DOatm$AGI_60m <- dat60_DOatm$pO2_60/dat60_DOatm$O2_demand60
dat250_DOatm$AGI_250m <- dat250_DOatm$pO2_250/dat250_DOatm$O2_demand250

  #explore outputs
hist(dat0_DOatm$AGI_0m)
hist(dat60_DOatm$AGI_60m)
hist(dat250_DOatm$AGI_250m)

plot(dat250_DOatm$votemper_mean, dat250_DOatm$AGI_250m)
plot(dat0_DOatm$pO2_0, dat0_DOatm$AGI_0m)

quantile(dat0_DOatm$AGI_0m, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
mean(dat250_DOatm$AGI_250m, na.rm = T)
min(dat250_DOatm$AGI_250m, na.rm = T)
max(dat250_DOatm$AGI_250m, na.rm = T)
sd(dat250_DOatm$AGI_250m, na.rm = T)

saveRDS(dat0_DOatm, here("data/locs_w_covar/psat_spot/cmem_locs_covar_AGI_0m.rds"))
saveRDS(dat60_DOatm, here("data/locs_w_covar/psat_spot/cmem_locs_covar_AGI_60m.rds"))
saveRDS(dat250_DOatm, here("data/locs_w_covar/psat_spot/cmem_locs_covar_AGI_250m.rds"))
