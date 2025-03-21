## note: new averaged environmental data files are created throughout this doc that have been commented out
## to prevent accidental overwriting. Please uncomment lines as needed.

### NOTE ####
#Please replace files paths as needed. This code was developed within an R Project, and thus the referenced file paths originate from the project's root directory. Please feel free to replicate or replace the file paths as needed.


### libraries ####
{ library(tidyverse)
  library(terra)
  library(here)
  library(tidyterra)
}

### load raster data ####
dat0_phys <- list.files(here("data/enviro/psat_spot_all/phys_merc/0m/processed"), full.names = TRUE)
dat60_phys <- list.files(here("data/enviro/psat_spot_all/phys_merc/60m/processed"), full.names = TRUE)
dat250_phys <- list.files(here("data/enviro/psat_spot_all/phys_merc/250m/processed"), full.names = TRUE)

dat0_bgc <- list.files(here("data/enviro/psat_spot_all/biogeo_cmem/0m"), full.names = TRUE)
dat60_bgc <- list.files(here("data/enviro/psat_spot_all/biogeo_cmem/60m"), full.names = TRUE)
dat250_bgc <- list.files(here("data/enviro/psat_spot_all/biogeo_cmem/250m"), full.names = TRUE)

template_rast <- rast(
  crs = 'EPSG:4326',
  extent = ext(-153, -103, 1, 49), #study domain (+/- 2 deg of min and max lat/lons from observed and PA locs)
  resolution = 0.25 #coarsest spatial resolution
)


### average by year ####
yr_avg_merge <- function(Phys_Input, BGC_Input, template_rast = template_rast, out_rast = NULL, all_names, longnames_input, units_input, phys_layers, bgc_layers){
  for(i in 1:length(Phys_Input)){
    phys_rast <- rast(Phys_Input[i])
    phys_name <- phys_layers[i]
    
    phys_rast_h <- resample(phys_rast, template_rast)
    phys_yr <- tapp(phys_rast_h, "years", mean)
    time(phys_yr) <- seq(as.Date("2003-01-01"), by = "year", length = 13)
      
    assign(phys_name, phys_yr)
    all_names <- c(all_names, phys_name)
    
  }
  
  for(i in 1:length(BGC_Input)){
    bgc_rast <- rast(BGC_Input[i])
    bgc_name <- bgc_layers[i]
    
    bgc_rast_h <- resample(bgc_rast, template_rast)
    bgc_yr <- tapp(bgc_rast_h, "years", mean)
    time(bgc_yr) <- seq(as.Date("2003-01-01"), by = "year", length = 13)
    
    assign(bgc_name, bgc_yr)
    all_names <- c(all_names, bgc_name)
  }
  
  out_rast <- sds(mget(all_names))
  names(out_rast) <- all_names
  longnames(out_rast) <- longnames_input
  units(out_rast) <- units_input
  
  return(out_rast)
}

#### 0m ####
all_names = NULL
dat0_yr <- yr_avg_merge(Phys_Input = dat0_phys, BGC_Input = dat0_bgc,
                        template_rast = template_rast,
                        all_names = NULL,
                        longnames_input = c("mixed layer depth", "salinity", "sea surface height", "sea water temperature", "eastward velocity", "eastward wind stress", "northward velocity", "northward wind stress", "dissolved oxygen", "chlorophyll"), 
                        units_input = c("m", "PSU", "m", "C", "m/s", "Pa", "m/s", "Pa", "m", "m", "mmol/m^3", "mg/m^3"), 
                        phys_layers = c("somxlavt", "vosaline", "sossheig", "votemper", "vozocrtx", "sozotaux", "vomecrty", "sometauy"), 
                        bgc_layers = c("chl", "o2"))

#writeCDF(dat0_yr, here("data/enviro/psat_spot_all/all_processed/annual_res/dat_0m_annual3.nc"), overwrite = TRUE)

#### 60m ####
all_names = NULL
dat60_yr <- yr_avg_merge(Phys_Input = dat60_phys, BGC_Input = dat60_bgc[2],
                        template_rast = template_rast,
                        all_names = NULL,
                        longnames_input = c("salinity", "sea water temperature", "dissolved oxygen"), 
                        units_input = c("PSU", "C", "mmol/m^3"), 
                        phys_layers = c("vosaline", "votemper"), 
                        bgc_layers = c("o2"))

#writeCDF(dat60_yr, here("data/enviro/psat_spot_all/all_processed/annual_res/dat_60m_annual.nc"), overwrite = TRUE)

#### 250m #####
all_names = NULL
dat250_yr <- yr_avg_merge(Phys_Input = dat250_phys, BGC_Input = dat250_bgc[2],
                        template_rast = template_rast,
                        all_names = NULL,
                        longnames_input = c("salinity", "sea water temperature", "dissolved oxygen"), 
                        units_input = c("PSU", "C", "mmol/m^3"), 
                        phys_layers = c("vosaline", "votemper"), 
                        bgc_layers = c( "o2"))

#writeCDF(dat250_yr, here("data/enviro/psat_spot_all/all_processed/annual_res/dat_250m_annual.nc"), overwrite = TRUE)

### average by season (static) ####
seas_avg_merge <- function(Phys_Input, BGC_Input, template_rast = template_rast, out_rast = NULL, all_names, longnames_input, units_input, phys_layers, bgc_layers){
  for(i in 1:length(Phys_Input)){
    phys_rast <- rast(Phys_Input[i])
    phys_name <- phys_layers[i]
    
    phys_rast_h <- resample(phys_rast, template_rast)
    phys_mo <- tapp(phys_rast_h, "months", mean)
    
    phys_wn <- mean(phys_mo[[c(12, 1, 2)]], na.rm = TRUE)
    phys_sp <- mean(phys_mo[[c(3, 4, 5)]], na.rm = TRUE)
    phys_su <- mean(phys_mo[[c(6, 7, 8)]], na.rm = TRUE)
    phys_fa <- mean(phys_mo[[c(9, 10, 11)]], na.rm = TRUE)
    
    phys_seas <- c(phys_wn, phys_sp, phys_su, phys_fa)
    
    time(phys_seas) <- as.Date(c("2003-12-01", "2003-03-01", "2003-06-01", "2003-09-01"), format = "%Y-%m-%d")
    
    assign(phys_name, phys_seas)
    all_names <- c(all_names, phys_name)
    
  }
  
  for(i in 1:length(BGC_Input)){
    bgc_rast <- rast(BGC_Input[i])
    bgc_name <- bgc_layers[i]
    
    bgc_rast_h <- resample(bgc_rast, template_rast)
    bgc_mo <- tapp(bgc_rast_h, "months", mean)
    
    bgc_wn <- mean(bgc_mo[[c(12, 1, 2)]], na.rm = TRUE)
    bgc_sp <- mean(bgc_mo[[c(3, 4, 5)]], na.rm = TRUE)
    bgc_su <- mean(bgc_mo[[c(6, 7, 8)]], na.rm = TRUE)
    bgc_fa <- mean(bgc_mo[[c(9, 10, 11)]], na.rm = TRUE)
    
    bgc_seas <- c(bgc_wn, bgc_sp, bgc_su, bgc_fa)
    
    time(bgc_seas) <- as.Date(c("2003-12-01", "2003-03-01", "2003-06-01", "2003-09-01"), format = "%Y-%m-%d")
    
    assign(bgc_name, bgc_seas)
    all_names <- c(all_names, bgc_name)
  }
  
  out_rast <- sds(mget(all_names))
  names(out_rast) <- all_names
  longnames(out_rast) <- longnames_input
  units(out_rast) <- units_input
  
  return(out_rast)
}

#### 0m ####
all_names = NULL
dat0_seas <- seas_avg_merge(Phys_Input = dat0_phys, BGC_Input = dat0_bgc,
                        template_rast = template_rast,
                        all_names = NULL,
                        longnames_input = c("mixed layer depth", "salinity", "sea surface height", "sea water temperature", "eastward velocity", "eastward wind stress", "northward velocity", "northward wind stress", "dissolved oxygen", "chlorophyll"), 
                        units_input = c("m", "PSU", "m", "C", "m/s", "Pa", "m/s", "Pa", "m", "m", "mmol/m^3", "mg/m^3"), 
                        phys_layers = c("somxlavt", "vosaline", "sossheig", "votemper", "vozocrtx", "sozotaux", "vomecrty", "sometauy"), 
                        bgc_layers = c("chl", "o2"))

#writeCDF(dat0_seas, here("data/enviro/psat_spot_all/all_processed/season_res/dat_0m_season.nc"), overwrite = TRUE)

#### 60m ####
all_names = NULL
dat60_seas <- seas_avg_merge(Phys_Input = dat60_phys, BGC_Input = dat60_bgc[2],
                            template_rast = template_rast,
                            all_names = NULL,
                            longnames_input = c("salinity", "sea water temperature", "dissolved oxygen"), 
                            units_input = c("PSU", "C", "mmol/m^3"), 
                            phys_layers = c("vosaline", "votemper"), 
                            bgc_layers = c( "o2"))

#writeCDF(dat60_seas, here("data/enviro/psat_spot_all/all_processed/season_res/dat_60m_season.nc"), overwrite = TRUE)

#### 250m ####
all_names = NULL
dat250_seas <- seas_avg_merge(Phys_Input = dat250_phys, BGC_Input = dat250_bgc[2],
                             template_rast = template_rast,
                             all_names = NULL,
                             longnames_input = c("salinity", "sea water temperature", "dissolved oxygen"), 
                             units_input = c("PSU", "C", "mmol/m^3"), 
                             phys_layers = c("vosaline", "votemper"), 
                             bgc_layers = c( "o2"))

#writeCDF(dat250_seas, here("data/enviro/psat_spot_all/all_processed/season_res/dat_250m_season.nc"), overwrite = TRUE)

