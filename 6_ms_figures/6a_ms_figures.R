### libraries ####
{library(here);here <- here::here #plyr's here function masks here::here
  library(MetBrewer)
  library(terra)
  library(ggBRT)
  library(patchwork)
  library(ggrepel)
  library(tidyverse)
  library(tidyr)
  set.seed(1004)
  source(here("functions/hsi_rast_functions.R"))
  source(here("functions/BRT_evaluation_functions.R"))
  source(here("functions/oxy_demand_functions_rev.R"))
  source(here("functions/avg_functions.R"))
  source(here("functions/partial_plot.R"))
}

### saved custom themes ####
theme_ms_map <- function(){ 
  font <- "Arial"   #assign font family up front
  
  theme_minimal() %+replace%    #replace elements we want to change
    
    theme(
      
      plot.title = element_text(             #axis titles
        family = font,            #font family
        color = "black",
        size = 16), 
      
      #grid elements
      panel.grid.major = element_blank(),    #strip major gridlines
      panel.grid.minor = element_blank(),    #strip minor gridlines
      axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      axis.title = element_text(             #axis titles
        family = font,            #font family
        color = "black",
        size = 16),               #font size
      
      axis.text = element_text(              #axis text
        family = font,            #axis famuly
        color = "black",
        size = 16),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10)),
      
      legend.position = "none"
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}

### Figure 1: Conceptual fig of hypotheses ###
#coastline data
north_map = map_data("world") %>% group_by(group)
shore     = north_map[north_map$region=="Canada" 
                      | north_map$region=="USA"
                      | north_map$region=="Mexico",]

ggplot(shore, aes(long, lat)) +
  #plot coastline
  coord_map("mercator", xlim = c(-153, -103), ylim = c(1, 49))+
  geom_polygon(aes(group=group), fill="grey75",lwd=1)+
  theme_ms_map()

ggsave(here("figs/ms/fig1_concept_hyp/coast.png"), height = 7 , width = 5, units = c("in"))

### Figure 2: locs over study area bathymetry ####
#### presence location data ####
#shark data
ani_locs <- readRDS(here("data/presence_locs/psat_spot_domain/processed/psat_spot_animotum.RDS")) %>% 
  mutate(PA = 0, rep = NA) %>% 
  subset(select = -c(geometry))

names(ani_locs) <- c("tag", "date", "lon", "lat", "PA", "rep")

#bathy data
r = rast(here("data/enviro/psat_spot_all/bathy_gebco/processed/gebco_bathy_0.25deg2.nc"))
bathy_df <- as.data.frame(r,xy = TRUE)
bathy_df <- bathy_df %>%
  filter(gebco_bathy_0.25deg2 <= 0)

mycolors <- c("#08306B","#023858", "#034B76", "#0B559F", "#045D92","#0469A6","#1379B4","#2F8BBD","#6BAED6", "#509AC6","#74A9CF","#88BEDC", "#90B4D5","#ACBFDC","#C4CBE2" )

ggplot(shore, aes(long, lat)) +
  #plot coastline
  coord_map("mercator", xlim = c(-153, -103), ylim = c(1, 49))+
  geom_polygon(aes(group=group), fill="grey75",lwd=1) +
  #plot bathymetry map
  geom_contour_filled(data=bathy_df, 
                      aes(x,y,z=gebco_bathy_0.25deg2)) +#breaks=seq(0,1,by=0.2)
  scale_fill_manual(values = mycolors)+
  #plot shark locations
  geom_path(data=test, aes(lon, lat, colour=as.factor(tag)), size = 0.5) +
  scale_color_manual(values = met.brewer("OKeeffe2", 73)) +
theme_ms_map() +
  xlim(-153, -103)+
  ylim(1, 49) +
  theme(legend.position = "none", 
        axis.title = element_blank()) + 
  guides(fill = guide_legend(title = "Bathymetry (m)", reverse = TRUE), 
         color = FALSE)

ggsave(here("figs/ms/fig2_tracks_bathy/tracks_bathy.png"), height = 7, width = 5, units = c("in"))

### Figure 3: AGI maps ####
#neutral year - 2013/2014
hsi_rast_gen(date_start = c("2013-09-01"), date_end = c("2014-01-31"), season = "FW", output_name = "neut_FW_Sept2013_Jan2014")

#La Niña - 2007/2008
hsi_rast_gen(date_start = c("2007-11-01"), date_end = c("2008-01-31"), season = "FW", output_name = "LN_FW_2007_2008")

#La Niña - 2010
hsi_rast_gen(date_start = c("2010-09-01"), date_end = c("2010-11-30"), season = "F", output_name = "LN_F_2010")

#EL Niño - 2014/2015
hsi_rast_gen(date_start = c("2014-11-01"), date_end = c("2015-01-31"), season = "FW", output_name = "EN_FW_Nov2014_Jan2015")

agi_250m_layered <- agi_maps_layerd(rast_folder_base = here("data/enviro/psat_spot_all/hsi_rasts/agi_rasts/neut_FW_Sept2013_Jan2014"), 
                                    rast_folder_LN = here("data/enviro/psat_spot_all/hsi_rasts/agi_rasts/LN_F_2010"), 
                                    rast_folder_EN = here("data/enviro/psat_spot_all/hsi_rasts/agi_rasts/EN_FW_Nov2014_Jan2015"))

ggsave(here("figs/ms/fig3_agi/agi_250m_layered.png"), agi_250m_layered, height = 8, width = 9, units = c("in"))

### Figure 4: performance metrics overall ####
#entire domain and study period
mod_metric_files <- list.files(here("data/brt/mod_outputs/revised"), pattern = ".rds", full.names = TRUE)

base_file <- readRDS(mod_metric_files[2])
do_file <- readRDS(mod_metric_files[3])
agi_file <- readRDS(mod_metric_files[1])

base_file$mod_type <- "Base model"
agi_file$mod_type <- "AGI model"
do_file$mod_type <- "DO model"

mod_metrics <- rbind(base_file, agi_file, do_file)
mod_metrics <- mod_metrics %>% mutate(st_id = "Overall")

#combine datasets
mod_metrics <- mod_metrics %>% mutate(dev_exp = dev_exp*100)
all_sum <- mod_metrics %>%
  group_by(mod_type) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE), 
            sd_auc = sd(AUC, na.rm = TRUE), 
            mean_tss = mean(TSS, na.rm = TRUE), 
            sd_tss = sd(TSS, na.rm = TRUE), 
            mean_dev = mean(dev_exp, na.rm = TRUE),
            sd_dev = sd(dev_exp, na.rm = TRUE)) %>%
  ungroup() 

#overall plots
TSS_overall <- all_sum %>% mutate(mod_type = as.factor(mod_type), 
                                  mod_type = fct_relevel(mod_type, c("Base model", "AGI model", "DO model"))) %>%
  arrange(desc(mean_tss)) %>%
  ggplot(aes(x = mod_type, y=mean_tss)) +
  geom_errorbar(aes(ymin = mean_tss - 2*sd_tss, ymax = mean_tss + 2*sd_tss), color = "black", size =  1, width = 0, linewidth = 1)+
  geom_point(color = "black", size = 4)+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("TSS") + 
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16)) 

AUC_overall <- all_sum %>% mutate(mod_type = as.factor(mod_type), 
                                  mod_type = fct_relevel(mod_type, c("Base model", "AGI model", "DO model"))) %>%
  arrange(desc(mean_auc)) %>%
  ggplot(aes(x = mod_type, y=mean_auc)) +
  geom_errorbar(aes(ymin = mean_auc - 2*sd_auc, ymax = mean_auc + 2*sd_auc), color = "black", size =  1, width = 0, linewidth = 1)+
  geom_point(color = "black", size = 4)+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("AUC") + 
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16)) 

dev_overall <- all_sum %>% mutate(mod_type = as.factor(mod_type), 
                                                 mod_type = fct_relevel(mod_type, c("Base model", "AGI model", "DO model"))) %>%
  arrange(desc(mean_dev)) %>%
  ggplot(aes(x = mod_type, y=mean_dev)) +
  geom_errorbar(aes(ymin = mean_dev - 2*sd_dev, ymax = mean_dev + 2*sd_dev), color = "black", size =  1, width = 0, linewidth = 1)+
  geom_point(color = "black", size = 4)+
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("% Deviance explained") + 
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16)) 

overall_metric_plots <- TSS_overall/AUC_overall/dev_overall
ggsave(here("figs/ms/fig4_perf_metric/overall_metrics.png"), overall_metric_plots, height = 13, width = 5, units = c("in"))

# Figure 5: predictor relative importance ####
avg_inf_df <- avg_pred(base_mods = list.files(here("data/brt/mod_outputs/perf_metric_iters_revised/base"), full.names = TRUE),
                       do_mods = list.files(here("data/brt/mod_outputs/perf_metric_iters_revised/do"), full.names = TRUE), 
                       agi_mods = list.files(here("data/brt/mod_outputs/perf_metric_iters_revised/agi"), full.names = TRUE)) 

avg_inf_sum <- avg_inf_df %>% 
  group_by(model, var) %>%
  summarise(inf_mean = mean(inf_val), 
            inf_sd = sd(inf_val))

for(i in 1:nrow(avg_inf_sum)){
  if(str_detect(avg_inf_sum$var[i], "AGI")){
    avg_inf_sum$var_type[i] = "AGI predictor"
  }
  
  if(str_detect(avg_inf_sum$var[i], "DO")){
    avg_inf_sum$var_type[i] = "DO predictor"
  }
  
  if(!str_detect(avg_inf_sum$var[i], "DO")&!str_detect(avg_inf_sum$var[i], "AGI")){
    avg_inf_sum$var_type[i] = "Environmental predictor"
  }
}

avg_inf_sum <- avg_inf_sum %>% 
  mutate(model = as.factor(model), 
         model= fct_relevel(model, c("Base", "DO", "AGI")), 
         var = as.factor(var), 
         var = fct_relevel(var, c("AGI, annual, 250m", "AGI, daily, 0m", "AGI, seasonal, 0m", "AGI, seasonal, 250m","AGI, annual, 0m", "AGI, daily, 250m")))

avg_inf_sum$model <- factor(avg_inf_sum$model, levels=c('Base', 'DO', 'AGI'))

#facet by model
pred_fig <- avg_inf_sum %>%
  ungroup() %>%
  mutate(var = as.factor(var), 
         var = fct_reorder(var, inf_mean)) %>%
  ggplot(aes(x = inf_mean, y = var, color = model, fill = model))+
  geom_errorbarh(aes(xmax = inf_mean + inf_sd, xmin = inf_mean - inf_sd, height = 0), linewidth = 1)+
  geom_point(size = 3)+
  scale_color_manual(values = c("#224B5E","#6A8D80", "#ABB98B")) +
  xlab("Relative importance (%)")+
  ylab("Predictor variable")+
  guides(fill="none", color = "none") +
  tidyquant::theme_tq()+
  theme(axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14), 
        legend.position = "top", 
        legend.justification = "left") +
facet_wrap(~model, scales = "free_y"); pred_fig

ggsave(here("figs/ms/fig5_pred/point_pred_facet.png"), pred_fig, width = 11, height = 5, units = c("in"))

### Figure 6: HSI maps study period ####
#generate raster
hsi_rast_gen(date_start = c("2003-01-01"), date_end = c("2015-12-31"), season = "SuFWSp", output_name = "Jan2003_Dec2015")

#plot HSI
all_maps_avg <- hsi_maps_avg(rast_folder = "data/enviro/psat_spot_all/hsi_rasts/Jan2003_Dec2015",
                             mod_folder = "data/brt/mod_outputs/revised",
                             fig_folder = "figs/ms/fig6_hsi_all",
                             ms = "Y", 
                             iter = 20)

ggsave(here("figs/ms/fig6_hsi_all/all_maps_avg_20.svg"), all_maps_avg, height = 7, width = 10, units = c("in"))

### Figure 7: ENSO HSI maps ####
#have to save using export button otherwise adds border, using height of 750 and width 500 (LN width 300)
#base year
enso_base <- hsi_maps_enso_avg(rast_folder = "data/enviro/psat_spot_all/hsi_rasts/neut_FW_Sept2013_Jan2014", enso = "diff", test_type = "revised", iter = 20)

#LN year - 2008
enso_LN_08 <- hsi_maps_difference_enso_avg(enso_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/LN_FW_2007_2008", test_type = "revised", neut_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/neut_FW_Sept2013_Jan2014", enso = "LN", iter = 20)

#LN year - 2010
enso_LN <- hsi_maps_difference_enso_avg(enso_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/LN_F_2010", test_type = "revised", neut_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/neut_FW_Sept2013_Jan2014", enso = "LN", iter = 20)

#EN year
enso_EN <- hsi_maps_difference_enso_avg(enso_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/EN_FW_Nov2014_Jan2015", test_type = "revised", neut_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/neut_FW_Sept2013_Jan2014", enso = "EN", iter = 20)

enso_all <- enso_base | enso_LN | enso_LN_08 | enso_EN
ggsave(here("figs/fig7_enso/enso_maps.png"), enso_all, height = 10, width = 10, units = c("in"))

#calculate % area HSI > 0.25 in NEC in strong LN year between models
calc_perc_area <- function(mod_rast, mod_type, test_type, iter){
  
  if(mod_type == "base"){  
  names(mod_rast) <- c("bathy_mean", "temp_mean", "sal_mean", "chl_mean", "ssh_mean", "bathy_sd", "mld_mean")
  }
  if(mod_type == "do"){
  names(mod_rast) <- c("o2_mean_0m", "o2_mean_250m_ann", "o2_mean_0m_seas", "temp_mean", "o2_mean_250m_seas", "bathy_mean", "sal_mean", "chl_mean", "o2_mean_0m_ann", "o2_mean_250m", "ssh_mean", "mld_mean", "bathy_sd")
  }
  if(mod_type == "agi"){
  names(mod_rast) <- c("temp_mean", "AGI_250m_ann", "AGI_0m", "bathy_mean", "AGI_0m_seas", "sal_mean", "AGI_250m_seas", "AGI_0m_ann", "chl_mean", "AGI_250m", "bathy_sd", "mld_mean", "ssh_mean")
  }
  
  #creating map dfs -------------------------------------------------------------------------------------------------
  mod_folder <- list.files(here(paste0("data/brt/mod_outputs/revised/", test_type,"/", mod_type)), full.names = TRUE)
  
  bbox <- ext(-153, -103, 1, 49)

  map_list <- rast()
  #for loop to create raster for each model iteration
  for(i in 1:iter){
    #creating map dfs -------------------------------------------------------------------------------------------------
    print(i)
    
    #predict
    mod_file <- readRDS(mod_folder[i])
    map_pred <- predict(mod_rast, mod_file, type = "response", n.trees = mod_file$gbm.call$best.trees, na.rm = FALSE)
    map_pred <- crop(map_pred, bbox)
    map_list <- c(map_list, map_pred)
  }
  
  #take average of rasters produced from each model
  pred_avg <- mean(map_list)

  # mask out land
  land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  land <- vect(land)

  land <- crop(land, bbox)
  pred_avg <- mask(pred_avg, land, inverse = TRUE)
  
  #filter for just NEC and calculate percent area w hsi > 0.75
  rast_nec_filt <- pred_avg %>% filter(y <= 15)

  hsi_nec <- raster::clamp(rast_nec_filt, lower = 0.25, values = FALSE)

  hsi_area_map <- expanse(hsi_nec)
  rast_area_map <- expanse(rast_nec_filt)

  perc_area <- (hsi_area_map$area/rast_area_map$area)*100
  
  return(perc_area)
}

#Base
base_ln_rast_10 <- rast(here("data/enviro/psat_spot_all/hsi_rasts/LN_F_2010/LN_F_2010_base_rast.nc"))
base_ln_rast_10 <- base_ln_rast_10 %>% filter(y > 2) #filter out weird artefact of averaging rasters
nec_base_10 <- calc_perc_area(mod_rast = base_ln_rast_10, mod_type = "base", test_type = "revised", iter = 20)

base_ln_rast_08 <- rast(here("data/enviro/psat_spot_all/hsi_rasts/LN_FW_2007_2008/LN_FW_2007_2008_base_rast.nc"))
base_ln_rast_08 <- base_ln_rast_08 %>% filter(y > 2) #filter out weird artefact of averaging rasters
nec_base_08 <- calc_perc_area(mod_rast = base_ln_rast_08, mod_type = "base", test_type = "revised", iter = 20)

#DO
do_ln_rast_10 <- rast(here("data/enviro/psat_spot_all/hsi_rasts/LN_F_2010/LN_F_2010_do_rast.nc"))
do_ln_rast_10 <- do_ln_rast_10 %>% filter(y > 2) #filter out weird artefact of averaging rasters
nec_do_10 <- calc_perc_area(mod_rast = do_ln_rast_10, mod_type = "do", test_type = "revised", iter = 20)

do_ln_rast_08 <- rast(here("data/enviro/psat_spot_all/hsi_rasts/LN_FW_2007_2008/LN_FW_2007_2008_do_rast.nc"))
do_ln_rast_08 <- do_ln_rast_08 %>% filter(y > 2) #filter out weird artefact of averaging rasters
nec_do_08 <- calc_perc_area(mod_rast = do_ln_rast_08, mod_type = "do", test_type = "revised", iter = 20)

#AGI
agi_ln_rast_10 <- rast(here("data/enviro/psat_spot_all/hsi_rasts/LN_F_2010/LN_F_2010_agi_rast.nc"))
agi_ln_rast_10 <- agi_ln_rast_10 %>% filter(y > 2) #filter out weird artefact of averaging rasters
nec_agi_10 <- calc_perc_area(mod_rast = agi_ln_rast_10, mod_type = "agi", test_type = "revised", iter = 20)

agi_ln_rast_08 <- rast(here("data/enviro/psat_spot_all/hsi_rasts/LN_FW_2007_2008/LN_FW_2007_2008_agi_rast.nc"))
agi_ln_rast_08 <- agi_ln_rast_08 %>% filter(y > 2) #filter out weird artefact of averaging rasters
nec_agi_08 <- calc_perc_area(mod_rast = agi_ln_rast_08, mod_type = "agi", test_type = "revised", iter = 20)

df <- data.frame(area = rep(c("NEC_LN_2010", "NEC_LC_2008"), times = 3),
            model_type = c("base", "do", "agi", "base", "do", "agi"),
            perc_area  = c(nec_base_10, nec_base_08, nec_do_10, nec_do_08, nec_agi_10, nec_agi_08)
); df