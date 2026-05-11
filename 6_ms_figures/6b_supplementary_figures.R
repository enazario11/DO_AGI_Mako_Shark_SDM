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
  source(here("functions/oxy_demand_functions.R"))
  source(here("functions/avg_functions.R"))
  #source(here("functions/partial_plot.R"))
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

#### ST1: Shark metadata ####
md_df <- readRDS(here("data/presence_locs/psat_spot_domain/psat_spot_data.rds"))
spot <- read.csv(here("data/presence_locs/mako_spot_filtered_1_step_per_day.csv"))
spot_ptts <- unique(spot$PTT)

md_df2 <- md_df %>% 
  group_by(ptt) %>%
  arrange(date) %>%
  mutate(deploy_dur = difftime(tail(date, n = 1), head(date, n = 1), units = "days"), 
         tag_type = ifelse(ptt %in% spot_ptts, "SPOT", "SPOT + PSAT")) %>%
  distinct(ptt, .keep_all = TRUE) %>%
  subset(select = c(ptt, tag_type, FL, sex, deploy_dur, date, lat, lon))

#### SF1: FL histos
ggplot(md_df2, aes(FL)) +  
  geom_histogram(color = "black", fill = "red4", bins = 15) + 
  facet_wrap(~sex) + 
  theme_tq() + 
  xlab("Fork length (cm)") + 
  theme(axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14))

#### SF2: Model performance metrics scatter ####
output_sum <- read.csv(here::here("data/brt/mod_outputs/brt_supp_table.csv"))
output_sum$deviance_exp <- output_sum$deviance_exp*100

mod_metrics <- ggplot(output_sum, aes(AUC, TSS, color = deviance_exp, label = model)) +
  geom_point(size = 5) +
  xlab('AUC') +
  ylab('TSS') +
  labs(color = "Deviance explained (%)")+
  scale_color_gradientn(colors = MetBrewer::met.brewer("Greek")) +
  ggrepel::geom_label_repel(aes(label = model),
                            box.padding   = 0.35,
                            point.padding = 1,
                            segment.color = 'grey50',
                            max.overlaps = 40,
                            label.size = 0.5)+
  theme_minimal()+
  theme(legend.position = "right", 
        axis.title = element_text(color = "black"), 
        axis.text = element_text(color = "black"))
ggsave(here::here("figs/ms/supp_figs/explore_mod_metrics.png"), mod_metrics, width = 10, height = 8, units = c("in"))

#### SF3: DO, 0m, 2003-2015 ####
rast_0m <- rast(here("data/enviro/psat_spot_all/all_processed/CMEM_DO_CHL_Temp_SO_UO_UOSTR_VO_VOSTR_SSH_MLD_0m_Jan2003_Dec2015_0.25_D.nc"))
rast_do_0m <- rast_0m["o2"]
rast_do_0m_mean <- mean(rast_do_0m)

do_0m_map <- ggplot() +
  geom_spatraster(data = rast_do_0m_mean) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))
  
ggsave(here("figs/ms/supp_figs/do_2003_2015_0m.png"), do_0m_map, height = 3, width = 6, units = c("in"))

#### SF4: DO, 250m, 2003-2015 ####
rast_250m <- rast(here("data/enviro/psat_spot_all/all_processed/CMEM_DO_Temp_SO_250m_Jan2003_Dec2015_0.25_D.nc"))
rast_do_250m <- rast_250m["o2"]
rast_do_250m_mean <- mean(rast_do_250m)

do_250m_map <- ggplot() +
  geom_spatraster(data = rast_do_250m_mean) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
theme_map() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

ggsave(here("figs/ms/supp_figs/do_2003_2015_250m.png"), do_250m_map, height = 3, width = 6, units = c("in"))

#### SF5: DO, 0m, seasonal averages ####
rast_0m_seas <- rast(here("data/enviro/psat_spot_all/all_processed/season_res/dat_0m_season.nc"))

  #winter Do 0m
rast_seas_0m_subW <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 12)
rast_seas_0m_do_W <- rast_seas_0m_subW["o2_1"]

do_0m_W <- ggplot() +
  geom_spatraster(data = rast_seas_0m_do_W) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(150, 350)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Winter, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

  #spring DO 0m 
rast_seas_0m_subSp <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 3)
rast_seas_0m_do_Sp <- rast_seas_0m_subSp["o2_2"]

do_0m_Sp <- ggplot() +
  geom_spatraster(data = rast_seas_0m_do_Sp) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(150, 350)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Spring, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

  #summer DO 0m
rast_seas_0m_subSu <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 6)
rast_seas_0m_do_Su <- rast_seas_0m_subSu["o2_3"]

do_0m_Su <- ggplot() +
  geom_spatraster(data = rast_seas_0m_do_Su) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(150, 350)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Summer, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

  #fall DO 0m 
rast_seas_0m_subF <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 9)
rast_seas_0m_do_F <- rast_seas_0m_subF["o2_4"]

do_0m_F <- ggplot() +
  geom_spatraster(data = rast_seas_0m_do_F) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(150, 350)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Fall, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))


all_do_0m <- (do_0m_W|do_0m_Sp)/(do_0m_Su|do_0m_F)+
  plot_layout(guides = "collect") & theme(legend.position = 'right', legend.title = element_text(size = 16), legend.text = element_text(size = 14)) & labs(fill = "Dissolved oxygen (mmol/m^3)")

ggsave(here("figs/ms/supp_figs/do_0m_seasonal.png"), all_do_0m, height = 8, width = 10, units = c("in"))

#### SF6: DO, 250m, seasonal averages ####
rast_250m_seas <- rast(here("data/enviro/psat_spot_all/all_processed/season_res/dat_250m_season.nc"))

#winter Do 250m
rast_seas_250m_subW <- subset(rast_250m_seas, month(time(rast_250m_seas)) == 12)
rast_seas_250m_do_W <- rast_seas_250m_subW["o2_1"]

do_250m_W <- ggplot() +
  geom_spatraster(data = rast_seas_250m_do_W) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(0, 250)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Winter, 250m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

#spring DO 250m 
rast_seas_250m_subSp <- subset(rast_250m_seas, month(time(rast_250m_seas)) == 3)
rast_seas_250m_do_Sp <- rast_seas_250m_subSp["o2_2"]

do_250m_Sp <- ggplot() +
  geom_spatraster(data = rast_seas_250m_do_Sp) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(0, 250)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Spring, 250m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

#summer DO 250m
rast_seas_250m_subSu <- subset(rast_250m_seas, month(time(rast_250m_seas)) == 6)
rast_seas_250m_do_Su <- rast_seas_250m_subSu["o2_3"]

do_250m_Su <- ggplot() +
  geom_spatraster(data = rast_seas_250m_do_Su) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(0, 250)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Summer, 250m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

#fall DO 250m 
rast_seas_250m_subF <- subset(rast_250m_seas, month(time(rast_250m_seas)) == 9)
rast_seas_250m_do_F <- rast_seas_250m_subF["o2_4"]

do_250m_F <- ggplot() +
  geom_spatraster(data = rast_seas_250m_do_F) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1, limits = c(0, 250)) +
  labs(fill = "Dissolved oxygen (mmol/m^3)")+
  ggtitle("Fall, 250m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))


all_do_250m <- (do_250m_W|do_250m_Sp)/(do_250m_Su|do_250m_F)+
  plot_layout(guides = "collect") & theme(legend.position = 'right', legend.title = element_text(size = 16), legend.text = element_text(size = 14)) & labs(fill = "Dissolved oxygen (mmol/m^3)")

ggsave(here("figs/ms/supp_figs/do_250m_seasonal.png"), all_do_250m, height = 8, width = 10, units = c("in"))

#### SF7: temperature, 0m, 2003-2015 ####
rast_temp_0m <- rast_0m["votemper"]
rast_temp_0m_mean <- mean(rast_temp_0m)

#plot
temp_0m_map <- ggplot() +
  geom_spatraster(data = rast_temp_0m_mean) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = -1) +
  labs(fill = "Temperature (C)")+
  theme_ms_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5), 
        legend.position = "right", 
        legend.title =element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"))

ggsave(here("figs/ms/supp_figs/temp_2003_2015_0m.png"), temp_0m_map, height = 3, width = 6, units = c("in"))

#### SF8: temperature, 0m, seasonal averages ####
rast_0m_seas <- rast(here("data/enviro/psat_spot_all/all_processed/season_res/dat_0m_season.nc"))

#winter temp 0m
rast_seas_0m_subW <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 12)
rast_seas_0m_temp_W <- rast_seas_0m_subW["votemper_1"]

temp_0m_W <- ggplot() +
  geom_spatraster(data = rast_seas_0m_temp_W) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = 1, limits = c(0, 35)) +
  labs(fill = "Temperature (C)")+
  ggtitle("Winter, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

#spring temp 0m 
rast_seas_0m_subSp <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 3)
rast_seas_0m_temp_Sp <- rast_seas_0m_subSp["votemper_2"]

temp_0m_Sp <- ggplot() +
  geom_spatraster(data = rast_seas_0m_temp_Sp) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = 1, limits = c(0, 35)) +
  labs(fill = "Temperature (C)")+
  ggtitle("Spring, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

#summer temp 0m
rast_seas_0m_subSu <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 6)
rast_seas_0m_temp_Su <- rast_seas_0m_subSu["votemper_3"]

temp_0m_Su <- ggplot() +
  geom_spatraster(data = rast_seas_0m_temp_Su) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = 1, limits = c(0, 35)) +
  labs(fill = "Temperature")+
  ggtitle("Summer, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))

#fall temp 0m 
rast_seas_0m_subF <- subset(rast_0m_seas, month(time(rast_0m_seas)) == 9)
rast_seas_0m_temp_F <- rast_seas_0m_subF["votemper_4"]

temp_0m_F <- ggplot() +
  geom_spatraster(data = rast_seas_0m_temp_F) +
  geom_map(data = testt,map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", direction = 1, limits = c(0, 35)) +
  labs(fill = "Temperature (C)")+
  ggtitle("Fall, 0m")+
  theme_map()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5))


all_temp_0m <- (temp_0m_W|temp_0m_Sp)/(temp_0m_Su|temp_0m_F)+
  plot_layout(guides = "collect") & theme(legend.position = 'right', legend.title = element_text(size = 16), legend.text = element_text(size = 14)) & labs(fill = "Temperature (C)")

ggsave(here("figs/ms/supp_figs/temp_0m_seasonal.png"), all_temp_0m, height = 8, width = 10, units = c("in"))

#### SF9: spatiotemporal performance metrics ####
base_dat_r <- readRDS(here("data/brt/mod_outputs/revised/brts_st/region/base_metrics.rds")) %>% drop_na() %>% mutate(mod_type = "Base model")
do_dat_r <- readRDS(here("data/brt/mod_outputs/revised/brts_st/region/do_metrics.rds")) %>% drop_na() %>% mutate(mod_type = "DO model")
agi_dat_r <- readRDS(here("data/brt/mod_outputs/revised/brts_st/region/agi_metrics.rds")) %>% drop_na() %>% mutate(mod_type = "AGI model")

all_reg <- rbind(base_dat_r, do_dat_r, agi_dat_r) %>%
  mutate(region = NA, 
         dev_exp = dev_exp*100)

for(i in 1:nrow(all_reg)){
  if(grepl('ccs',all_reg$st_id[i])){
    all_reg$region[i] = "CCS"
  } 
  if(grepl('nec',all_reg$st_id[i])){
    all_reg$region[i] = "NEC"
  } 
  if(grepl('nep', all_reg$st_id[i])){
    all_reg$region[i] = "NEP"
  }
}

sum_reg <- all_reg %>%
  group_by(region, mod_type) %>%
  summarise(mean_tss = mean(TSS), 
            sd_tss = sd(TSS), 
            mean_auc = mean(AUC), 
            sd_auc = sd(AUC), 
            mean_dev = mean(dev_exp), 
            sd_dev = sd(dev_exp))

#enso boxplots
TSS_neut_e <- sum_enso %>% mutate(mod_type = as.factor(mod_type), 
                              mod_type = fct_relevel(mod_type, c("Base model", "DO model", "AGI model")), 
                              ENSO = as.factor(ENSO), 
                              ENSO = fct_relevel(ENSO, c("Neutral", "El Niño", "La Niña"))) %>%
  filter(ENSO == "Neutral") %>%
  ggplot(aes(x = ENSO, y=mean_tss)) +
  geom_point(aes(fill = mod_type, color = mod_type), position = position_dodge(width = 1), shape = 22, size = 5) +
  geom_linerange(aes(ymin = mean_tss-sd_tss, ymax = mean_tss+sd_tss, color = mod_type), position = position_dodge(width = 1), size = 1) +
  theme_tq() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("TSS") + 
  scale_fill_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  scale_color_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  labs(fill = "")+
  guides(color = "none")+
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top", 
        legend.justification = "left", 
        strip.text.x = element_text(size = 14), 
        rect = element_rect(fill = "transparent") ) 


TSS_stress_e <- sum_enso %>% mutate(mod_type = as.factor(mod_type), 
                                  mod_type = fct_relevel(mod_type, c("Base model", "DO model", "AGI model")), 
                                  ENSO = as.factor(ENSO), 
                                  ENSO = fct_relevel(ENSO, c("Neutral", "El Niño", "La Niña"))) %>%
  filter(ENSO != "Neutral") %>%
  ggplot(aes(x = ENSO, y=mean_tss)) +
  geom_point(aes(fill = mod_type, color = mod_type), position = position_dodge(width = 1), shape = 22, size = 5) +
  geom_linerange(aes(ymin = mean_tss-sd_tss, ymax = mean_tss+sd_tss, color = mod_type), position = position_dodge(width = 1), size = 1) +
  theme_tq() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("TSS") + 
  scale_fill_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  scale_color_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  labs(fill = "")+
  guides(color = "none")+
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top", 
        legend.justification = "left", 
        strip.text.x = element_text(size = 14)) 

TSS_neut_r <- sum_reg %>% mutate(mod_type = as.factor(mod_type), 
                              mod_type = fct_relevel(mod_type, c("Base model", "DO model", "AGI model")), 
                              region = as.factor(region), 
                              region = fct_relevel(region, c("NEP", "CCS", "NEC"))) %>%
  filter(region != "NEC") %>%
  ggplot(aes(x = region, y=mean_tss)) +
  geom_point(aes(fill = mod_type, color = mod_type), position = position_dodge(width = 1), shape = 22, size = 5) +
  geom_linerange(aes(ymin = mean_tss-sd_tss, ymax = mean_tss+sd_tss, color = mod_type), position = position_dodge(width = 1), size = 1) +
  theme_tq() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("TSS") + 
  scale_fill_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  scale_color_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  labs(fill = "")+
  guides(color = "none")+
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top", 
        legend.justification = "left", 
        strip.text.x = element_text(size = 14),
        rect = element_rect(fill = "transparent"))

TSS_stress_r <- sum_reg %>% mutate(mod_type = as.factor(mod_type), 
                                 mod_type = fct_relevel(mod_type, c("Base model", "DO model", "AGI model")), 
                                 region = as.factor(region), 
                                 region = fct_relevel(region, c("NEP", "CCS", "NEC"))) %>%
  filter(region == "NEC") %>%
  ggplot(aes(x = region, y=mean_tss)) +
  geom_point(aes(fill = mod_type, color = mod_type), position = position_dodge(width = 1), shape = 22, size = 5) +
  geom_linerange(aes(ymin = mean_tss-sd_tss, ymax = mean_tss+sd_tss, color = mod_type), position = position_dodge(width = 1), size = 1) +  
  theme_tq() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  ylab("TSS") + 
  scale_fill_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  scale_color_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  labs(fill = "")+
  guides(color = "none")+
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top", 
        legend.justification = "left", 
        strip.text.x = element_text(size = 14),
        rect = element_rect(fill = "transparent") ) 

#### SF10: base model partial plots ####
##### function for the partial plots #####
ggPD_boot <- function (gbm.object, predictor = NULL, n.plots = length(pred.names), 
                       list.4.preds = NULL, booted.preds = NULL, nrow = NULL, ncol = NULL, 
                       col.line = "darkorange", cex.line = 0.5, type.ci = "lines", 
                       col.ci = "grey80", cex.ci = 0.3, lty.ci = 2, alpha.ci = 0.5, 
                       smooth = FALSE, col.smooth = "blue", cex.smooth = 0.3, span = 0.3, 
                       rug = FALSE, rug.pos = "t", common.scale = TRUE, type = NULL, cis = c(0.025, 
                                                                                0.975), y.label = "Fitted function", x.label = paste(var.name, "  (", 
                                                                                                                                     round(gbm.object$contributions[predictor, 2], 
                                                                                                                                           1), "%)", sep = ""), 
                       ...) 
{
  gbm.call <- gbm.object$gbm.call
  pred.names <- gbm.call$predictor.names
  ggPD_boot.plots <- function(gbm.object) {
    if (!requireNamespace("gbm")) {
      stop("you need to install the gbm package to run this function")
    }
    if (is.null(booted.preds)) {
      stop("you need to set booted.preds as the array from the bootstrap run\n           (eg testboot$function.preds using testboot<-gbm.bootstrap.functions())")
    }
    if (is.null(list.4.preds)) {
      stop("you need to set list.4.preds as the result of plot.gbm.4list()")
    }
    requireNamespace("splines")
    gbm.x <- gbm.call$gbm.x
    response.name <- gbm.call$response.name
    nt <- gbm.object$n.trees
    data <- gbm.call$dataframe
    max.vars <- length(gbm.object$contributions$var)
    if (n.plots > max.vars) {
      n.plots <- max.vars
      warning("reducing no of plotted predictors to maximum available (", 
              max.vars, ")")
    }
    predictors <- list(rep(NA, n.plots))
    responses <- list(rep(NA, n.plots))
    responses.lower <- list(rep(NA, n.plots))
    responses.upper <- list(rep(NA, n.plots))
    for (j in c(1:max.vars)) {
      k <- match(gbm.object$contributions$var[j], pred.names)
      if (is.null(x.label)) {
        var.name <- gbm.call$predictor.names[k]
      }
      else {
        var.name <- x.label
      }
      pred.data <- data[, gbm.call$gbm.x[k]]
      response.matrix <- gbm::plot.gbm(gbm.object, i.var = k, 
                                       n.trees = nt, return.grid = TRUE, ...)
      predictors[[j]] <- response.matrix[, 1]
      if (is.factor(data[, gbm.call$gbm.x[k]])) {
        predictors[[j]] <- factor(predictors[[j]], levels = levels(data[, 
                                                                        gbm.call$gbm.x[k]]))
      }
      responses[[j]] <- response.matrix[, 2] - mean(response.matrix[, 
                                                                    2])
      num.values <- nrow(response.matrix)
      temp <- apply(booted.preds[, k, ] - mean(booted.preds[, 
                                                            k, ]), 1, function(x) {
                                                              quantile(x, cis[1], na.rm = T)
                                                            })
      responses.lower[[j]] <- temp[1:num.values]
      temp <- apply(booted.preds[, k, ] - mean(booted.preds[, 
                                                            k, ]), 1, function(x) {
                                                              quantile(x, cis[2], na.rm = T)
                                                            })
      responses.upper[[j]] <- temp[1:num.values]
      if (j == 1) {
        ymin = min(responses.lower[[j]])
        ymax = max(responses.upper[[j]])
        dat <- data.frame(pred.data)
      }
      else {
        ymin = min(ymin, min(responses.lower[[j]]))
        ymax = max(ymax, max(responses.upper[[j]]))
        dat <- data.frame(dat, pred.data)
      }
    }
    if (is.null(predictor)) {
      fittedFunc <- list()
      fittedFunc.lower <- list()
      fittedFunc.upper <- list()
      fittedVal <- list()
      ribbon <- list()
      ggPD <- list()
      for (i in 1:n.plots) {
        k <- match(gbm.object$contributions$var[i], pred.names)
        var.name <- gbm.call$predictor.names[k]
        fittedFunc[[i]] <- data.frame(predictors[i], 
                                      responses[i])
        colnames(fittedFunc[[i]]) <- c("x", "y")
        fittedFunc.lower[[i]] <- data.frame(predictors[i], 
                                            responses.lower[i])
        colnames(fittedFunc.lower[[i]]) <- c("x", "y")
        fittedFunc.upper[[i]] <- data.frame(predictors[i], 
                                            responses.upper[i])
        colnames(fittedFunc.upper[[i]]) <- c("x", "y")
        fittedVal[[i]] <- data.frame(gbm.object$fitted, 
                                     dat[i])
        colnames(fittedVal[[i]]) <- c("y", "x")
        ribbon[[i]] <- data.frame(x = fittedFunc.lower[[i]]$x, 
                                  ylow = fittedFunc.lower[[i]]$y, yup = fittedFunc.upper[[i]]$y)
        if (is.factor(fittedFunc[[i]]$x)) {
          ggPD[[i]] <- ggplot(fittedFunc[[i]], aes(x = x, 
                                                   y = y)) + geom_boxplot(color = col.line, 
                                                                          size = cex.line) + geom_boxplot(data = fittedFunc.lower[[i]], 
                                                                                                          aes(x = x, y = y), color = col.ci) + geom_boxplot(data = fittedFunc.upper[[i]], 
                                                                                                                                                            aes(x = x, y = y), color = col.ci) + ylab(y.label) + 
            xlab(paste(var.name, "  (", round(gbm.object$contributions[i, 
                                                                       2], 1), "%)", sep = "")) + theme_bw() + 
            theme(panel.grid.minor = element_line(linetype = "blank"), 
                  panel.grid.major = element_line(linetype = "blank"), 
                  axis.text = element_text(size = 6), axis.title = element_text(size = 10), 
                  axis.line.y = element_line(size = 0.1), 
                  axis.line.x = element_line(size = 0.1))
          if (common.scale == T) {
            ggPD[[i]] <- ggPD[[i]] + ylim(c(ymin, ymax))
          }
        }
        if (type.ci == "lines") {
          ggPD[[i]] <- ggplot(fittedFunc[[i]], aes(x = x, 
                                                   y = y)) + geom_line(color = col.line, size = cex.line) + 
            geom_line(data = fittedFunc.lower[[i]], aes(x = x, 
                                                        y = y), size = cex.ci, color = col.ci, 
                      linetype = lty.ci) + geom_line(data = fittedFunc.upper[[i]], 
                                                     aes(x = x, y = y), size = cex.ci, color = col.ci, 
                                                     linetype = lty.ci) + ylab(y.label) + xlab(paste(var.name, 
                                                                                                     "  (", round(gbm.object$contributions[i, 
                                                                                                                                           2], 1), "%)", sep = "")) + theme_bw() + 
            theme(panel.grid.minor = element_line(linetype = "blank"), 
                  panel.grid.major = element_line(linetype = "blank"), 
                  axis.title.x = element_text(size = 10), 
                  axis.line.y = element_line(size = 0.1), 
                  axis.line.x = element_line(size = 0.1))
          if (smooth == T) {
            ggPD[[i]] <- ggPD[[i]] + geom_smooth(span = span, 
                                                 size = 0.3, color = col.smooth, se = F, 
                                                 linetype = 2)
          }
          if (rug == T) {
            ggPD[[i]] <- ggPD[[i]] + geom_rug(data = fittedVal[[i]], 
                                              aes(x = x, y = y), sides = rug.pos, position = "jitter", 
                                              color = "#EBEBEB")
          }
          if (common.scale == T) {
            ggPD[[i]] <- ggPD[[i]] + ylim(c(ymin, ymax))
          }
        }
        if (type.ci == "ribbon") {
          ggPD[[i]] <- ggplot() + geom_ribbon(data = ribbon[[i]], 
                                              aes(x = x, ymin = ylow, ymax = yup), fill = col.ci, 
                                              alpha = alpha.ci) + geom_line(data = fittedFunc[[i]], 
                                                                            aes(x = x, y = y), color = col.line, size = cex.line) + 
            ylab(y.label) + xlab(paste(var.name, "  (", 
                                       round(gbm.object$contributions[i, 2], 1), 
                                       "%)", sep = "")) + theme_bw() + theme(panel.grid.minor = element_line(linetype = "blank"), 
                                                                             panel.grid.major = element_line(linetype = "blank"), 
                                                                             axis.title = element_text(size = 10), axis.line.y = element_line(size = 0.1), 
                                                                             axis.line.x = element_line(size = 0.1))
          if (smooth == T) {
            ggPD[[i]] <- ggPD[[i]] + geom_smooth(data = fittedFunc[[i]], 
                                                 aes(x = x, y = y), span = span, size = 0.3, 
                                                 color = col.smooth, se = F, linetype = 2)
          }
          if (rug == T) {
            ggPD[[i]] <- ggPD[[i]] + geom_rug(data = fittedVal[[i]], 
                                              aes(x = x, y = y), sides = rug.pos, position = "jitter", 
                                              color = "#EBEBEB")
          }
          if (common.scale == T) {
            ggPD[[i]] <- ggPD[[i]] + ylim(c(ymin, ymax))
          }
        }
      }
      list(ggPD = ggPD)
    }
    else {
      if (is.character(predictor)) {
        predictor <- match(predictor, gbm.object$contributions$var)
      }
      k <- match(gbm.object$contributions$var[predictor], 
                 pred.names)
      var.name <- gbm.call$predictor.names[k]
      fittedFunc <- data.frame(predictors[predictor], responses[predictor])
      colnames(fittedFunc) <- c("x", "y")
      fittedFunc.lower <- data.frame(predictors[predictor], 
                                     responses.lower[predictor])
      colnames(fittedFunc.lower) <- c("x", "y")
      fittedFunc.upper <- data.frame(predictors[predictor], 
                                     responses.upper[predictor])
      colnames(fittedFunc.upper) <- c("x", "y")
      ribbon <- data.frame(x = fittedFunc.lower$x, ylow = fittedFunc.lower$y, 
                           yup = fittedFunc.upper$y)
      fittedVal <- data.frame(gbm.object$fitted, dat[predictor])
      colnames(fittedVal) <- c("y", "x")
      if (is.factor(fittedFunc$x)) {
        ggPD <- ggplot(fittedFunc, aes(x = x, y = y)) + 
          geom_boxplot(color = col.line, size = cex.line) + 
          geom_boxplot(data = fittedFunc.lower, aes(x = x, 
                                                    y = y), color = col.ci) + geom_boxplot(data = fittedFunc.upper, 
                                                                                           aes(x = x, y = y), color = col.ci) + ylab(y.label) + 
          xlab(paste(var.name, "  (", round(gbm.object$contributions[predictor, 
                                                                     2], 1), "%)", sep = "")) + theme_bw() + theme(panel.grid.minor = element_line(linetype = "blank"), 
                                                                                                                   panel.grid.major = element_line(linetype = "blank"), 
                                                                                                                   axis.text = element_text(size = 6), axis.title = element_text(size = 10), 
                                                                                                                   axis.line.y = element_line(size = 0.1), axis.line.x = element_line(size = 0.1))
        if (common.scale == T) {
          ggPD <- ggPD + ylim(c(ymin, ymax))
        }
      }
      if (type.ci == "lines") {
        ggPD <- ggplot(fittedFunc, aes(x = x, y = y)) + 
          geom_line(color = col.line, size = cex.line) + 
          geom_line(data = fittedFunc.lower, aes(x = x, 
                                                 y = y), size = cex.ci, color = col.ci, linetype = lty.ci) + 
          geom_line(data = fittedFunc.upper, aes(x = x, 
                                                 y = y), size = cex.ci, color = col.ci, linetype = lty.ci) + 
          ylab(y.label) + xlab(paste(var.name, "  (", 
                                     round(gbm.object$contributions[predictor, 2], 
                                           1), "%)", sep = "")) + theme_bw() + theme(panel.grid.minor = element_line(linetype = "blank"), 
                                                                                     panel.grid.major = element_line(linetype = "blank"), 
                                                                                     axis.title = element_text(size = 10), axis.line.y = element_line(size = 0.1), 
                                                                                     axis.line.x = element_line(size = 0.1))
        if (smooth == T) {
          ggPD <- ggPD + geom_smooth(span = span, size = 0.3, 
                                     color = col.smooth, se = F, linetype = 2)
        }
        if (rug == T) {
          ggPD <- ggPD + geom_rug(data = fittedVal, aes(x = x, 
                                                        y = y), sides = rug.pos, position = "jitter", 
                                  color = "#EBEBEB")
        }
        if (common.scale == T) {
          ggPD <- ggPD + ylim(c(ymin, ymax))
        }
      }
      if (type.ci == "ribbon" & type == "base") {
        ggPD <- ggplot() + geom_ribbon(data = ribbon, 
                                       aes(x = x, ymin = ylow, ymax = yup), fill = col.ci, 
                                       alpha = alpha.ci) + geom_line(data = fittedFunc, 
                                                                     aes(x = x, y = y), color = col.line, size = cex.line) + 
          ylab(y.label) + xlab(x.label) + theme_minimal() + theme(panel.grid.minor = element_line(linetype = "blank"), 
                                                                                     panel.grid.major = element_line(linetype = "blank"), 
                                                                                     axis.title = element_text(size = 14), axis.line.y = element_line(size = 0.1), 
                                                                                     axis.line.x = element_line(size = 0.1), 
                                                                  axis.text = element_text(size = 14, color = "black"))
        if (smooth == T) {
          ggPD <- ggPD + geom_smooth(data = fittedFunc, 
                                     aes(x = x, y = y), span = span, size = 0.3, 
                                     color = col.smooth, se = F, linetype = 2)
        }
        if (rug == T) {
          ggPD <- ggPD + geom_rug(data = fittedVal, aes(x = x, 
                                                        y = y), sides = rug.pos, position = "jitter", 
                                  color = "#EBEBEB")
        }
        if (common.scale == T) {
          ggPD <- ggPD + ylim(c(ymin, ymax))
        }
      }
      
      if (type.ci == "ribbon" & type == "do") {
        ggPD <- ggplot() + geom_ribbon(data = ribbon, 
                                       aes(x = x, ymin = ylow, ymax = yup), fill = col.ci, 
                                       alpha = alpha.ci) + geom_line(data = fittedFunc, 
                                                                     aes(x = x, y = y), color = col.line, size = cex.line) + 
          ylab(y.label) + xlab(x.label) + theme_minimal() + theme(panel.grid.minor = element_line(linetype = "blank"), 
                                                             panel.grid.major = element_line(linetype = "blank"), 
                                                             axis.title = element_text(size = 14), axis.line.y = element_blank(), 
                                                             axis.line.x = element_line(size = 0.1), 
                                                             axis.text = element_text(size = 14, color = "black"))
        if (smooth == T) {
          ggPD <- ggPD + geom_smooth(data = fittedFunc, 
                                     aes(x = x, y = y), span = span, size = 0.3, 
                                     color = col.smooth, se = F, linetype = 2)
        }
        if (rug == T) {
          ggPD <- ggPD + geom_rug(data = fittedVal, aes(x = x, 
                                                        y = y), sides = rug.pos, position = "jitter", 
                                  color = "#EBEBEB")
        }
        if (common.scale == T) {
          ggPD <- ggPD + ylim(c(ymin, ymax))
        }
      }
      
      if (type.ci == "ribbon" & type == "agi") {
        ggPD <- ggplot() + geom_ribbon(data = ribbon, 
                                       aes(x = x, ymin = ylow, ymax = yup), fill = col.ci, 
                                       alpha = alpha.ci) + geom_line(data = fittedFunc, 
                                                                     aes(x = x, y = y), color = col.line, size = cex.line) + 
          ylab(y.label) + xlab(x.label) + theme_minimal() + theme(panel.grid.minor = element_line(linetype = "blank"), 
                                                                  panel.grid.major = element_line(linetype = "blank"), 
                                                                  axis.title = element_text(size = 14), axis.line.y = element_blank(), 
                                                                  axis.line.x = element_line(size = 0.1), 
                                                                  axis.text = element_text(size = 14, color = "black"))
        if (smooth == T) {
          ggPD <- ggPD + geom_smooth(data = fittedFunc, 
                                     aes(x = x, y = y), span = span, size = 0.3, 
                                     color = col.smooth, se = F, linetype = 2)
        }
        if (rug == T) {
          ggPD <- ggPD + geom_rug(data = fittedVal, aes(x = x, 
                                                        y = y), sides = rug.pos, position = "jitter", 
                                  color = "#EBEBEB")
        }
        if (common.scale == T) {
          ggPD <- ggPD + ylim(c(ymin, ymax))
        }
      }
      
      list(ggPD = ggPD)
    }
  }
  plot <- ggPD_boot.plots(gbm.object)
  if (is.null(predictor)) {
    do.call(grid.arrange, c(plot$ggPD, list(nrow = nrow, 
                                            ncol = ncol)))
  }
  else grid.draw(plot$ggPD)
}

##### Base model plot ####
# Boostrap the BRT 1000 times to build confidence intervals
base_mod <- readRDS(here("data/brt/mod_outputs/revised/base/base_1.rds"))
brt1.prerun_base<- plot.gbm.4list(base_mod)
base_boot <- gbm.bootstrap.functions(base_mod, list.predictors=brt1.prerun_base, n.reps=20)

# Draw partial dependency plots a given predictor (i.e., "Complexity")
#base model
plot_list <- list()
base_names <- c("z", "temp", "sal", "chl-a", "SSH", "z_sd", "MLD")
for(i in 1:nrow(base_mod$contributions)){

  var.name = base_names[i]

  plot_temp <- ggPD_boot(base_mod, 
                         predictor = base_mod$contributions[i, 1], 
                         list.4.preds = brt1.prerun_base, 
                         booted.preds = base_boot$function.preds, 
                         type.ci = "ribbon",
                         type = "base",
                         col.line = "#92351e",
                         cex.line = 1.5,
                         rug = T, 
                         alpha.ci = 0.75, 
                         y.label = "Probability of presence", 
                         var.name = var.name)
                        #  x.label = paste(base_names[i], "  (", 
                        #                  round(base_mod$contributions[i, 2], 
                        #                        1), "%)", sep = ""))
  
  plot_list[[i]] <- plot_temp
}

base_plots <- do.call(grid.arrange, c(plot_list, ncol = 5))

#### SF11: DO model partial plots ####
do_mod_fin <- readRDS(here("data/brt/mod_outputs/revised/do/do_1.rds"))
brt1.prerun_do<- plot.gbm.4list(do_mod_fin)
do_boot <- gbm.bootstrap.functions(do_mod_fin, list.predictors=brt1.prerun_do, n.reps=20)

plot_list <- list()
do_names <- c("DO, annual, 250m", "DO, daily, 0m", "DO, seasonal, 0m", "temp", "sal", "DO, seasonal, 250m", "z", "chl-a", "DO, daily, 250m", "DO, annual, 0m",  "SSH", "MLD", "z_sd")
for(i in 1:nrow(do_mod_fin$contributions)){

  var.name = do_names[i]

  plot_temp <- ggPD_boot(do_mod_fin, 
                         predictor = do_mod_fin$contributions[i, 1], 
                         list.4.preds = brt1.prerun_do, 
                         booted.preds = do_boot$function.preds, 
                         type.ci = "ribbon",
                         col.line = "#92351e", 
                         cex.line = 1.5,
                         type = "do",
                         rug = T, 
                         alpha.ci = 0.75,
                         y.label = "Probability of presence", 
                         var.name = var.name)
  
  plot_list[[i]] <- plot_temp
}

do_plots <- do.call(grid.arrange, c(plot_list, ncol = 5))

#### SF12: AGI model partial plots ####
agi_mod_fin <- readRDS(here("data/brt/mod_outputs/revised/agi/agi_1.rds"))
brt1.prerun_agi<- plot.gbm.4list(agi_mod_fin)
agi_boot <- gbm.bootstrap.functions(agi_mod_fin, list.predictors=brt1.prerun_agi, n.reps=20)

plot_list <- list()
agi_names <- c("AGI, annual, 250m", "temp", "AGI, daily, 0m", "z", "AGI, seasonal, 0m", "sal", "AGI, seasonal, 250m", "SSH", "AGI, annual, 0m", "AGI, daily, 250m", "chl-a", "z_sd", "MLD")
for(i in 1:nrow(agi_mod_fin$contributions)){
  plot_temp <- ggPD_boot(agi_mod_fin, 
                         predictor = agi_mod_fin$contributions[i, 1], 
                         list.4.preds = brt1.prerun_agi, 
                         booted.preds = agi_boot$function.preds, 
                         type.ci = "ribbon",
                         rug = T, 
                         type = "agi",
                         col.line = "#92351e", 
                         cex.line = 1.5,
                         alpha.ci = 0.75, 
                         y.label = "Probability of presence",
                         x.label = "")
  
  plot_list[[i]] <- plot_temp
}

agi_plots <- do.call(grid.arrange, c(plot_list, ncol = 5))

#### SF13: SHAP values ####
### target point
target_loc <- vect(cbind(-141,12), crs="EPSG:4326")

### load models
brt_do <- readRDS(here("data/brt/mod_outputs/revised/brts_st/enso/do/do_1.rds"))
brt_agi <- readRDS(here("data/brt/mod_outputs/revised/brts_st/enso/agi/agi_1.rds"))

pred_fun <- function(model, newdata) {
  predict(model, newdata = newdata, type = "response")
}

hsi_rast_gen(date_start = c("2010-10-26"), date_end = c("2010-10-26"), season = "F", output_name = "LN_F_Oct26_2010")

### DO 
### Training data
test_do <- readRDS(here("data/brt/mod_outputs/revised/brts_st/enso/do/test/do_test1.rds")) 
train_do <- filter(dat_do_st, !(row_id %in% test_do$row_id)) %>% 
  dplyr::select(-c(tag, date, lon, lat, rep, dt, uo_mean, uostr_mean, vo_mean, vostr_mean, soi, st_id))

# Extract target location's env
do_rast <- rast("data/enviro/psat_spot_all/hsi_rasts/LN_F_Oct26_2010/LN_F_Oct26_2010_do_rast.nc")
names(do_rast) <- c("o2_mean_0m", "o2_mean_250m_ann", "o2_mean_0m_seas", "temp_mean", "o2_mean_250m_seas", "bathy_mean", "sal_mean", "chl_mean", "o2_mean_0m_ann", "o2_mean_250m", "ssh_mean", "mld_mean", "bathy_sd")
target_env_do <- terra::extract(do_rast, target_loc) %>% 
  as_tibble() %>% 
  select(-ID)

#plot
x = train_do[, c("chl_mean", "temp_mean", "sal_mean", "ssh_mean", "mld_mean", "bathy_mean", "bathy_sd", "o2_mean_0m", "o2_mean_250m", "o2_mean_0m_seas", "o2_mean_0m_ann", "o2_mean_250m_seas", "o2_mean_250m_ann")]
predictor <- Predictor$new(brt_do, data = x, y = train_do$PA, predict.function = pred_fun)
shapley <- Shapley$new(predictor, x.interest = target_env_do)

print(shapley$results)
do_shap <- plot(shapley) + tidyquant::theme_tq()
ggsave(here("figs/ms/supp_figs/shap/do_shap.png"), do_shap, width = 6, height = 6, units = c("in"))

# AGI 
# Training data
test_agi <- readRDS(here("data/brt/mod_outputs/revised/brts_st/enso/agi/test/agi_test1.rds")) 
train_agi <- filter(dat_agi_st, !(row_id %in% test_agi$row_id)) %>% 
  dplyr::select(-c(tag, date, lon, lat, rep, dt, uo_mean, uostr_mean, vo_mean, vostr_mean, soi, st_id))

# Extract target location's env
agi_rast <- rast("data/enviro/psat_spot_all/hsi_rasts/LN_F_Oct26_2010/LN_F_Oct26_2010_agi_rast.nc")
names(agi_rast) <- c("temp_mean", "AGI_250m_ann", "AGI_0m", "bathy_mean", "AGI_0m_seas", "sal_mean", "AGI_250m_seas", "AGI_0m_ann", "chl_mean", "AGI_250m", "bathy_sd", "mld_mean", "ssh_mean")
target_env_agi <- terra::extract(agi_rast, target_loc) %>% 
  as_tibble() %>% 
  select(-ID)

#plot
x = train_agi[,c("chl_mean", "temp_mean", "sal_mean", "ssh_mean", "mld_mean", "bathy_mean", "bathy_sd", "AGI_0m", "AGI_250m", "AGI_0m_seas", "AGI_0m_ann", "AGI_250m_seas", "AGI_250m_ann")]
predictor <- Predictor$new(brt_agi, data = x, y = train_agi$PA, predict.function = pred_fun)
shapley <- Shapley$new(predictor, x.interest = target_env_agi)

print(shapley$results)
agi_shap <- plot(shapley) + tidyquant::theme_tq()
ggsave(here("figs/ms/supp_figs/shap/agi_shap.png"), agi_shap, width = 6, height = 6, units = c("in"))

#### SF14: DO difference plots, 0 and 250m, neutral and LN year ####
map.world = map_data(map="world")
testt = map.world %>% filter(long <= 180)

#neutral rast: 2013
test_n <- rast(here("data/enviro/psat_spot_all/hsi_rasts/Jan13_Dec13/Jan13_Dec13_do_rast.nc"))
names(test_n) <- c("o2_mean_0m", "o2_mean_250m_ann", "o2_mean_0m_seas", "temp_mean", "o2_mean_250m_seas", "bathy_mean", "sal_mean", "chl_mean", "o2_mean_0m_ann", "o2_mean_250m", "ssh_mean", "mld_mean", "bathy_sd")

#LN rast: 2010
test_l <- rast(here("data/enviro/psat_spot_all/hsi_rasts/LN_F_2010/LN_F_2010_do_rast.nc"))
names(test_l) <- c("o2_mean_0m", "o2_mean_250m_ann", "o2_mean_0m_seas", "temp_mean", "o2_mean_250m_seas", "bathy_mean", "sal_mean", "chl_mean", "o2_mean_0m_ann", "o2_mean_250m", "ssh_mean", "mld_mean", "bathy_sd")

#plot
plot(test_n$o2_mean_0m)
plot(test_l$o2_mean_0m)
diff_0m <- diff(c(test_n$o2_mean_0m, test_l$o2_mean_0m))

plot(test_n$o2_mean_250m)
plot(test_l$o2_mean_250m)
diff_250m <- diff(c(test_n$o2_mean_250m, test_l$o2_mean_250m))

#map
diff_map_0m <- ggplot() +
  geom_spatraster(data = diff_0m) +
  geom_map(data = testt, map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", limits = c(-20, 20), direction = -1) +
  theme_ms_map() +
  ggtitle("Difference in DO at 0m")+
  theme(legend.position = "right")

diff_map_250m <- ggplot() +
  geom_spatraster(data = diff_250m) +
  geom_map(data = testt, map = testt,aes(map_id = region, x = long, y = lat), fill = "grey75", color = "black") +
  scale_x_continuous(expand =c(0,0),limits = c(-153,-103)) +
  scale_y_continuous(expand=c(0,0),limits = c(1,49)) +
  scale_fill_whitebox_c(palette = "muted", limits = c(-100, 100), direction = -1) +
  theme_ms_map() +
  ggtitle("Difference in DO at 250m")+
  theme(legend.position = "right")

all_diff <- diff_map_0m|diff_map_250m
ggsave(here("figs/ms/supp_figs/change_do.png"), all_diff, width = 12, height = 5, units = c("in"))