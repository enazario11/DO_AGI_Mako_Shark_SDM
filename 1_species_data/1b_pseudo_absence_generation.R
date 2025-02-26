#load libraries ####
library(tidyverse)
library(aniMotum)
library(sf)
library(here);here <- here::here
library(terra)
library(tmvtnorm)

### NOTE ####
#Please replace files paths as needed. This code was developed within an R Project, and thus the referenced file paths originate from the project's root directory. Please feel free to replicate or replace the file paths as needed.

#load data ####
#checking min and max lat/lon for describing PA domain
ssm_crw_r <- readRDS(here("data/presence_locs/psat_spot_domain/processed/animotum_reroute_crw_ssm.RDS"))

dat_pred <- ssm_crw_r %>% 
  rowwise() %>% 
  mutate(routed = list(st_transform(ssm$rerouted, 4326)),
         routed = list(subset(routed, select = c(geometry, date)))) %>%
  ungroup()

#unnest prediction column
ssm_locs <- dat_pred %>%
  unnest(routed) %>%
  subset(select = c(id, geometry, date))

#add coordinates as new lat lon cols - also add WC indicator
ssm_locs$lon_p<-st_coordinates(ssm_locs$geometry)[,1] # get coordinates
ssm_locs$lat_p<-st_coordinates(ssm_locs$geometry)[,2] # get coordinates

#to be used for pseudoabsence generation below and for environmental data extraction (folder 3)
saveRDS(ssm_locs, file = here("data/presence_locs/psat_spot_domain/processed/psat_spot_animotum.RDS"))

#observed min and max
min(ssm_locs$lon_p) #-151.08
max(ssm_locs$lon_p) #-105.41
min(ssm_locs$lat_p) #2.99
max(ssm_locs$lat_p) #47.37

#buffered min and max for CMEMS data download
min(ssm_locs$lon_p) -2 #-153.41
max(ssm_locs$lon_p) +2 #-103.41
min(ssm_locs$lat_p) -2 #0.98
max(ssm_locs$lat_p) +2 #49.37

# aniMotum CRW PA generation ####
# create spatVect for CMEMS domain used to generate gradient
#get the bounding box of the two x & y coordintates, make sfc -- min and max of interpolated data + or - 2
ylims <- c(2.99, 47.37)
xlims <- c(-151.08, -105.41)
box_coords <- tibble(x = xlims, y = ylims) %>% 
  st_as_sf(coords = c("x", "y")) %>% 
  st_set_crs(st_crs(4326))

bounding_box <- st_bbox(box_coords) %>% st_as_sfc()

#reproject to CRS that aligns with SSM output
bb_merc <- st_transform(bounding_box, crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")

#read in continents polygon
land_vect <- read_sf(here("data/enviro/continents_shp/World_Continents.shp"))

land_merc = land_vect %>%
  st_as_sfc() %>%
  st_transform(crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs")

land_subset <- st_intersection(land_merc, bb_merc)
plot(land_subset)

#crop where ROMS domain polygon and continents polygons intersect to get a final polygon of the CMEMS domain
grad_poly <- st_difference(bb_merc, land_subset)
plot(grad_poly)

df <- data.frame(id = seq(length(grad_poly)))
df$geometry <- grad_poly
grad_sf <- st_as_sf(df)

grad_spatVect <- vect(grad_poly)

# calculate gradient between CMEMS polygon (that fits within CMEMS dataset) and template that is 3x the CMEMS domain
#domain that is x3 area of ROMS domain
CMEMS_large_rast <- rast(
  crs = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs",
  extent = ext(-18924.31, -9448.548, -2016.141, 10116.68), 
  resolution =  32.01272
)

#create 2D gradient
x <- rasterize(grad_spatVect, CMEMS_large_rast, fun = "mean") 

## generate gradient rasters
dist <- distance(x)
x1 <- terrain(dist, v = "slope", unit = "radians")
y1 <- terrain(dist, v = "aspect", unit = "radians")
grad.x <- -1 * x1 * cos(0.5 * pi - y1)
grad.y <- -1 * x1 * sin(0.5 * pi - y1)
grad <- c(grad.x, grad.y)

plot(grad)

#develop pseudo crw locs using the predicted times (constrained to have same number of locs)
#load 04/08/24 updated sim_fit() function
source(here("functions/sim_fit.R"))

#selected a beta value of -375 as this was the smallest absolute value number that resulted in min of 90% of the PAs remaining in the study area
pa_crw <- sim_fit(ssm_crw_r, grad = grad, beta = c(-350, -350), what = "predicted", reps = 100);filter_pa <- sim_filter(pa_crw, keep = 0.30, var = c("lon", "lat"), FUN = "mean");routed_pa <- route_path(filter_pa, centroids = TRUE)
plot(routed_pa[4,])

saveRDS(routed_pa, file = here("data/presence_locs/psat_spot_domain/processed/PA_routed_350beta_30perc.RDS"))

#IDing how many tracks are outside of the study domain 
#unnest prediction column
dat_pa <- routed_pa %>%
  unnest(sims)

#ID by simulation rep how many reps stay in study area
dat_pa2 <- dat_pa %>% 
  mutate(domain = ifelse(lon <= -151.08 |
                           lon >= -105.41 |
                           lat <= 2.99 | 
                           lat >= 47.37, "omit", "keep"))

dat_pa_filt <- dat_pa2 %>%
  group_by(id, rep) %>%
  filter(!any(domain == "omit")) %>%
  ungroup()

test <- dat_pa_filt %>% 
  group_by(id) %>% 
  summarise(n_reps = n_distinct(rep)) %>% #min number of reps inside domain is XXX, sampling rest down to that count
  print(n = 73) #min is 39 tracks -- sample to that

mean(test$n_reps/30*100) #mean is 101% stay in

#randomly sample reps so that all have same number
dat_PA_22 <- NULL
for(i in 1:length(unique(dat_pa_filt$id))){
  #select current id
  curr_ID <- unique(dat_pa_filt$id)[i]
  print(curr_ID)
  temp_df <- dat_pa_filt[dat_pa_filt$id %in% curr_ID,]
  
  #sample 52 id's randomly
  temp_rep_ID <- sample(unique(temp_df$rep), 22, replace = FALSE)
  
  #narrow your data set
  temp_df2 <- temp_df[temp_df$rep %in% temp_rep_ID, ]
  
  #combine in a single df
  dat_PA_22 <- rbind(dat_PA_24, temp_df2)
  
}

#check if it worked 
dat_PA_22 %>%
  group_by(id) %>%
  summarise(n_reps = n_distinct(rep)) %>%
  print(n = 73)

#save pseudoabsence location dataframe
saveRDS(dat_PA_22, file = here("data/presence_locs/psat_spot_domain/processed/psat_spot_PAs.rds"))


