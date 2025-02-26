#load packages 
library(tidyverse)
library(aniMotum)
library(sf)
library(here);here <- here::here
library(terra)
library(tmvtnorm)
set.seed(1004)

### NOTE ####
#Please replace files paths as needed. This code was developed within an R Project, and thus the referenced file paths originate from the project's root directory. Please feel free to replicate or replace the file paths as needed.

##### load presence data ####
#spot and psat data

#Location data imported here can be found on Dryad as 'psat_spot_io.csv': https://datadryad.org/stash/dataset/doi:10.5061/dryad.31zcrjdxz#readme
all_dat <- readRDS(here("data/presence_locs/psat_spot_domain/psat_spot_data.rds")) #

#format for aniMotum
ssm_dat <- all_dat %>% 
  mutate(posix2 = as.POSIXct(date)) %>% #re-do because was filtering out dates w/o a time stamp
  select(id = "ptt", date ="posix2", lat = "lat",lon = "lon", lc = "lc") %>%
  drop_na(date) #must remove missing lat/lon/date values for ssm to work
 
  #view number of locs by ptt
all_dat %>%
  group_by(ptt) %>%
  summarise(n = n()) %>%
  print(n = 85)

##### CRW SSM ####
#calculate average temporal step length between raw locs (method for determining SSM time step from Maxwell et al., 2019 -- blue sharks spatial segregation)
time_step <- ssm_dat %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(diff = date - lag(date)) %>%
  summarise(med_diff = median(diff, na.rm = T)) %>%
  ungroup() %>%
  summarise(all_mean = mean(med_diff)/3600) #average time step btwn positions for all tracks is 29 hours

ssm_crw <- fit_ssm(ssm_dat,
                   model = "crw",
                   time.step = 29) #average median time step in hours
                    
ssm_crw_r <- route_path(ssm_crw, map_scale = 10, what = "predicted")
aniMotum::map(ssm_crw_r, what = "predicted")|aniMotum::map(ssm_crw_r, what = "rerouted") #plot ssm output for all sharks

#visually inspect residuals and diagnostic plots 
resid_crw <- osar(ssm_crw)

plot(resid_crw, type = "qq", pages = 0)
plot(resid_crw, type = "acf", pages = 0)
plot(resid_crw, type = "ts", pages = 0)

#save RDS file to be used to extract environmental data
saveRDS(ssm_crw_r, file = here("data/presence_locs/psat_spot_domain/processed/animotum_reroute_crw_ssm.rds"))

