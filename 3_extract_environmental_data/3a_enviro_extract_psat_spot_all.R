#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load packages####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(sf)
library(terra)
library(ncdf4)
library(here);here <- here::here
library(raster)

#load ROMSextract, CMEMextract, and pseudo depth functions
source(here("functions/enviro_extract_functions.R"))

set.seed(1004)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load the data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #state space model locs
dat_locs <- readRDS(here("data/presence_locs/psat_spot_domain/processed/psat_spot_animotum.RDS")) %>% mutate(PA = 0, rep = NA)

  #PA locs
pa_locs <- readRDS(here("data/presence_locs/psat_spot_domain/processed/psat_spot_PAs.RDS")) %>% mutate(PA = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#join locs and PAs together####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_locs_comb <- dat_locs %>% 
  subset(select = -c(geometry))
names(dat_locs_comb) <- c("tag", "date", "lon", "lat", "PA", "rep")

pa_locs_comb <- pa_locs %>%
  subset(select = -c(model, x, y, domain))
names(pa_locs_comb) <- c("tag", "rep", "date", "lon", "lat", "PA")
pa_locs_comb <- pa_locs_comb[, c(1, 3, 4, 5, 6, 2)] #reorders columns to match dat_locs_comb DF

all_locs <- rbind(dat_locs_comb, pa_locs_comb)

#set loc df for extraction
    #locs
input_file <- all_locs
input_file$date <- as.factor(as.Date(substr(input_file$date, 1,  10))) #Ensure date format is ok for getvarROMS. 
input_file$dt <- as.POSIXct(strptime(input_file$date, format = "%Y-%m-%d"), tz = "UTC")

    # remove points outside of domain (num of obsv. shouldn't change: 367776)
input_file <- input_file[input_file$lat>=1 & input_file$lat<=49,] 
input_file <- input_file[input_file$lon>=-153 & input_file$lon<=-103,] 

head(input_file)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract bathy and rugosity ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bathy_file <- list.files(here("data/enviro/psat_spot_all/all_processed"), pattern = "bathy_0.25deg2", full.names = TRUE)

all_dat_bathy_cmem <- getBathy(bathy_file, input_file, 'gebco_bathy_0.25deg2', 0.25) #update varid when I delete the other bathy file

#explore outputs
head(all_dat_bathy_cmem)
hist(all_dat_bathy_cmem$bathy, breaks = 30) #bathymetry
hist(all_dat_bathy_cmem$bathy_sd, breaks = 30) #rugosity

ggplot(all_dat_bathy_cmem, aes(bathy)) + geom_histogram(bins = 15, color = "grey") + facet_wrap(~PA, scales = "free") + theme_bw()
ggplot(all_dat_bathy_cmem, aes(bathy_sd)) + geom_histogram(bins = 15, color = "grey") + facet_wrap(~PA, scales = "free") + theme_bw()

all_dat_bathy_cmem %>% 
  group_by(PA) %>% 
  summarise(med_bathy = median(bathy, na.rm = TRUE), 
            med_rug = median(bathy_sd, na.rm = TRUE))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract CMEMS covars ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 0m all dat extract ####
cmem_nc0 <- nc_open(here("data/enviro/psat_spot_all/all_processed/CMEM_DO_CHL_Temp_SO_UO_UOSTR_VO_VOSTR_SSH_MLD_0m_Jan2003_Dec2015_0.25_D.nc"))

#surface extract
xtracto_cmem = function(input_file, nc_file){
  input_file <- getvarCMEM(nc_file, "o2", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "chl", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "votemper", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "vosaline", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "vozocrtx", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "sozotaux", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "vomecrty", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "sometauy", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "sossheig", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "somxlavt", input_file, 0.25, mean, "mean")
}

all_dat_cmem_0m <- xtracto_cmem(input_file, cmem_nc0)

#combine with bathy above
all_cmem_covar_0m <- cbind(all_dat_cmem_0m, all_dat_bathy_cmem$bathy, all_dat_bathy_cmem$bathy_sd)
all_cmem_covar_0m <- all_cmem_covar_0m %>% 
  rename("bathy" = "all_dat_bathy_cmem$bathy", 
         "bathy_sd" = "all_dat_bathy_cmem$bathy_sd")

#saveRDS(all_cmem_covar_0m, here("data/locs_w_covar/psat_spot/cmem_locs_covar_0m.rds"))
head(all_cmem_covar_0m)

    #explore
ggplot(all_cmem_covar_0m, aes(o2_mean)) + geom_histogram(bins = 30, color = "grey") + facet_wrap(~PA, scales = "free") + theme_bw()

cmem_0m_long <- gather(all_cmem_covar_0m, covar, value, o2_mean:bathy_sd) %>% mutate(PA = as.factor(PA))

ggplot(cmem_0m_long, aes(x = value, fill = PA)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~covar, scales = "free") + 
  theme_bw()+
  scale_fill_manual(values = c("dodgerblue4", "darkseagreen4"))

all_cmem_covar_0m %>% 
  group_by(PA) %>% 
  summarise(mean_temp = mean(votemper_mean, na.rm = TRUE), 
            mean_do = mean(o2_mean, na.rm = TRUE),
            mean_sal = mean(vosaline_mean, na.rm = TRUE), 
            mean_mld = mean(ssheig_mean, na.rm = TRUE), 
            mean_ssh = mean(somxlavt_mean, na.rm = TRUE), 
            mean_uo = mean(vozocrtx_mean, na.rm = TRUE), 
            mean_uostr = mean(sozotaux, na.rm = TRUE),
            mean_vo = mean(vomecrty_mean, na.rm = TRUE), 
            mean_vostr = mean(sometauy, na.rm = TRUE),
            mean_chl = mean(chl_mean, na.rm = TRUE))

#60m extract ####
# dat extract ###
cmem_nc60 <- nc_open(here("data/enviro/psat_spot_all/all_processed/CMEM_DO_Temp_SO_60m_Jan2003_Dec2015_0.25_D.nc"))

xtracto_cmem_depth = function(input_file, nc_file){
  input_file <- getvarCMEM(nc_file, "o2", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "votemper", input_file, 0.25, mean, "mean")
  input_file <- getvarCMEM(nc_file, "vosaline", input_file, 0.25, mean, "mean")
}

all_dat_cmem_60m <- xtracto_cmem_depth(input_file, cmem_nc60)

#saveRDS(all_dat_cmem_60m, here("data/locs_w_covar/psat_spot/cmem_locs_covar_60m.rds"))
head(all_dat_cmem_60m)

#explore
ggplot(all_dat_cmem_60m, aes(o2_mean)) + geom_histogram(bins = 30, color = "grey") + facet_wrap(~PA, scales = "free") + theme_bw()

cmem_60m_long <- gather(all_dat_cmem_60m, covar, value, o2_mean:vosaline_mean) %>% mutate(PA = as.factor(PA))

ggplot(cmem_60m_long, aes(x = value, fill = PA)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~covar, scales = "free") + 
  theme_bw()+
  scale_fill_manual(values = c("dodgerblue4", "darkseagreen4"))

all_dat_cmem_60m %>% 
  group_by(PA) %>% 
  summarise(mean_temp = mean(votemper_mean, na.rm = TRUE), 
            mean_do = mean(o2_mean, na.rm = TRUE),
            mean_sal = mean(vosaline_mean, na.rm = TRUE))

#250m extract ####
# dat extract ###
cmem_nc250 <- nc_open(here("data/enviro/psat_spot_all/all_processed/CMEM_DO_Temp_SO_250m_Jan2003_Dec2015_0.25_D.nc"))

all_dat_cmem_250m <- xtracto_cmem_depth(input_file, cmem_nc250)

#saveRDS(all_dat_cmem_250m, here("data/locs_w_covar/psat_spot/cmem_locs_covar_250m.rds"))
head(all_dat_cmem_250m)

#explore
ggplot(all_dat_cmem_250m, aes(o2_mean)) + geom_histogram(bins = 30, color = "grey") + facet_wrap(~PA, scales = "free") + theme_bw()

cmem_250m_long <- gather(all_dat_cmem_250m, covar, value, o2_mean:vosaline_mean) %>% mutate(PA = as.factor(PA))

ggplot(cmem_250m_long, aes(x = value, fill = PA)) + 
  geom_density(alpha = 0.5) + 
  facet_wrap(~covar, scales = "free") + 
  theme_bw()+
  scale_fill_manual(values = c("dodgerblue4", "darkseagreen4"))

all_dat_cmem_250m %>% 
  group_by(PA) %>% 
  summarise(mean_temp = mean(votemper_mean, na.rm = TRUE), 
            mean_do = mean(o2_mean, na.rm = TRUE),
            mean_sal = mean(vosaline_mean, na.rm = TRUE))
