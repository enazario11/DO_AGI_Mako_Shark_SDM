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
  source(here("scripts/7a_diet_data.R"))
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


### Figure 2: AGI map ####
#neutral year
hsi_rast_gen(date_start = c("2013-09-01"), date_end = c("2014-01-31"), season = "FW", output_name = "neut_FW_Sept2013_Jan2014")

#La Niña
hsi_rast_gen(date_start = c("2010-09-01"), date_end = c("2010-11-30"), season = "F", output_name = "LN_F_2010")


#EL Niño
hsi_rast_gen(date_start = c("2014-11-01"), date_end = c("2015-01-31"), season = "FW", output_name = "EN_FW_Nov2014_Jan2015")


agi_250m_layered <- agi_maps_layerd(rast_folder_base = here("data/enviro/psat_spot_all/hsi_rasts/agi_rasts/Jan13_Dec13"), 
                                    rast_folder_LN = here("data/enviro/psat_spot_all/hsi_rasts/agi_rasts/LN_F_2010"), 
                                    rast_folder_EN = here("data/enviro/psat_spot_all/hsi_rasts/agi_rasts/EN_FW_Nov2014_Jan2015"))

ggsave(here("figs/ms/fig2_agi/agi_250m_layered.png"), agi_250m_layered, height = 7, width = 8, units = c("in"))

### Figure 3a: performance metrics overall ####
#entire domain and study period
mod_metric_files <- list.files(here("data/brt/mod_outputs/perf_metric_iters"), pattern = ".rds", full.names = TRUE)

base_file <- readRDS(mod_metric_files[2])
do_file <- readRDS(mod_metric_files[4])
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
  summarise(mean_auc = mean(AUC), 
            sd_auc = sd(AUC), 
            mean_tss = mean(TSS), 
            sd_tss = sd(TSS), 
            mean_dev = mean(dev_exp),
            sd_dev = sd(dev_exp)) %>%
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

ggsave(here("figs/ms/fig5_metrics_all/overall_metrics.png"), overall_metric_plots, height = 13, width = 5, units = c("in"))

#Figure 3b: performance metrics st ####
#spatiotemporal analysis
base_st <- readRDS(here("data/brt/mod_outputs/perf_metric_iters/brts_st/base_metrics.rds")) %>% mutate(mod_type = "Base model")
do_st <- readRDS(here("data/brt/mod_outputs/perf_metric_iters/brts_st/do_metrics.rds")) %>% mutate(mod_type = "DO model")
agi_st <- readRDS(here("data/brt/mod_outputs/perf_metric_iters/brts_st/agi_metrics.rds")) %>% mutate(mod_type = "AGI model")

all_st <- rbind(base_st, do_st, agi_st) %>%
  mutate(region = NA, 
         ENSO = NA, 
         dev_exp = dev_exp*100)

for(i in 1:nrow(all_st)){
  if(grepl('ccs',all_st$st_id[i])){
    all_st$region[i] = "CCS"
  } 
  if(grepl('nec',all_st$st_id[i])){
    all_st$region[i] = "NEC"
  } 
  if(grepl('nep', all_st$st_id[i])){
    all_st$region[i] = "NEP"
  }
  if(grepl('Overall',all_st$st_id[i])){
    all_st$region[i] = "Overall"
  }  
  if(grepl('en',all_st$st_id[i])){
    all_st$ENSO[i] = "El Niño"
  } 
  if(grepl('ln',all_st$st_id[i])){
    all_st$ENSO[i] = "La Niña"
  } 
  if(grepl('neut', all_st$st_id[i])){
    all_st$ENSO[i] = "Neutral"
  }
  if(grepl('Overall',all_st$st_id[i])){
    all_st$ENSO[i] = "Overall"}
}

#spatiotemporal plots
TSS_plot <- all_st %>% mutate(mod_type = as.factor(mod_type), 
                               mod_type = fct_relevel(mod_type, c("Base model", "AGI model", "DO model")), 
                               st_id = as.factor(st_id), 
                               region = as.factor(region), 
                               region = fct_relevel(region, c("NEP", "CCS", "NEC")), 
                               ENSO = as.factor(ENSO), 
                               ENSO = fct_relevel(ENSO, c("Neutral", "El Niño", "La Niña"))) %>%
  ggplot(aes(x = mod_type, y=TSS)) +
  geom_boxplot(aes(fill = region, color = region), position = position_dodge(width = 1), lwd = 1, alpha = 0.55)+
  theme_tq() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(~ENSO, scales = "free_y")+
  xlab("") +
  ylab("TSS") + 
  scale_fill_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  scale_color_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  labs(fill = "Region:")+
  guides(color = "none")+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "top", 
        legend.justification = "left", 
        strip.text.x = element_text(size = 14)) 

AUC_plot <- all_st %>% mutate(mod_type = as.factor(mod_type), 
                               mod_type = fct_relevel(mod_type, c("Base model", "AGI model", "DO model")), 
                               st_id = as.factor(st_id), 
                               region = as.factor(region), 
                               region = fct_relevel(region, c("NEP", "CCS", "NEC")), 
                               ENSO = as.factor(ENSO), 
                               ENSO = fct_relevel(ENSO, c("Neutral", "El Niño", "La Niña"))) %>%
  ggplot(aes(x = mod_type, y=AUC)) +
  geom_boxplot(aes(fill = region, color = region), position = position_dodge(width = 1), lwd = 1, alpha = 0.55)+
  theme_tq() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(~ENSO, scales = "free_y")+
  xlab("") +
  ylab("AUC") + 
  scale_color_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  scale_fill_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  #coord_cartesian(ylim = c(0.4, 0.65))+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "none", 
        strip.text.x = element_text(size = 14)) 

perc_exp_plot <- all_st %>% mutate(mod_type = as.factor(mod_type), 
                                    mod_type = fct_relevel(mod_type, c("Base model", "AGI model", "DO model")), 
                                    st_id = as.factor(st_id), 
                                    region = as.factor(region), 
                                    region = fct_relevel(region, c("NEP", "CCS", "NEC")), 
                                    ENSO = as.factor(ENSO), 
                                    ENSO = fct_relevel(ENSO, c("Neutral", "El Niño", "La Niña"))) %>%
  ggplot(aes(x = mod_type, y=dev_exp)) +
  geom_boxplot(aes(fill = region, color = region), position = position_dodge(width = 1), lwd = 1, alpha = 0.55)+
  theme_tq() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(~ENSO, scales = "free_y")+
  xlab("") +
  ylab("% Deviance explained") + 
  scale_color_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  scale_fill_manual(values = c("#224B5E", "#527875", "#83A58C"))+
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 16), 
        legend.title = element_text(size = 16, color = "black"), 
        legend.text = element_text(size = 14, color = "black"),
        legend.position = "none", 
        strip.text.x = element_text(size = 14)) 

all_metric_plots <- TSS_plot/AUC_plot/perc_exp_plot

ggsave(here("figs/ms/fig6_metrics_st/Figure_6_Metrics_st.png"), all_metric_plots, height = 10, width = 13, units = c("in"))

# Figure 4: predictor relative importance ####
#list models
base_mod <- readRDS(here("data/brt/mod_outputs/final_mods/brt_base_0m_dail_no_wind.rds"))
do_mod_fin <- readRDS(here("data/brt/mod_outputs/final_mods/brt_do_0m_250m_dail_seas_ann.rds"))
agi_mod_fin <- readRDS(here("data/brt/mod_outputs/final_mods/brt_agi_0m_250m_dail_seas_ann.rds"))
do_agi_comb <- readRDS(here("data/brt/mod_outputs/final_mods/brt_agi_250_DO_0_dail_seas_ann.rds"))

#base model 
base_inf <- as.data.frame(ggBRT::ggInfluence(base_mod, plot = FALSE)) %>% rownames_to_column()
colnames(base_inf) <- c("Predictor_variable", "relative_influence")
base_inf$Predictor_variable <- gsub("bathy_mean", "z", base_inf$Predictor_variable)
base_inf$Predictor_variable <- gsub("temp_mean", "temp", base_inf$Predictor_variable)
base_inf$Predictor_variable <- gsub("sal_mean", "sal", base_inf$Predictor_variable)
base_inf$Predictor_variable <- gsub("chl_mean", "chl-a", base_inf$Predictor_variable)
base_inf$Predictor_variable <- gsub("bathy_sd", "z_sd", base_inf$Predictor_variable)
base_inf$Predictor_variable <- gsub("ssh_mean", "SSH", base_inf$Predictor_variable)
base_inf$Predictor_variable <- gsub("mld_mean", "MLD", base_inf$Predictor_variable)

inf_df <- data.frame(model = c("Base", "Base", "Base", "Base", "Base", "Base", "Base"), 
                     var = base_inf$Predictor_variable, 
                     inf_val = base_inf$relative_influence)
#do model
do_inf <- as.data.frame(ggBRT::ggInfluence(do_mod_fin, plot = FALSE)) %>% rownames_to_column()
colnames(do_inf) <- c("Predictor_variable", "relative_influence")
do_inf$Predictor_variable <- gsub("\\<o2_mean_0m\\>", "DO, daily, 0m", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("o2_mean_250m_ann", "DO, annual, 250m", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("o2_mean_0m_seas", "DO, seasonal, 0m", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("o2_mean_250m_seas", "DO, seasonal, 250m", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("temp_mean", "temp", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("sal_mean", "sal", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("bathy_mean", "z", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("chl_mean", "chl-a", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("o2_mean_0m_ann", "DO, annual, 0m", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("o2_mean_250m", "DO, daily, 250m", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("ssh_mean", "SSH", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("mld_mean", "MLD", do_inf$Predictor_variable)
do_inf$Predictor_variable <- gsub("bathy_sd", "z_sd", do_inf$Predictor_variable)

do_inf <- do_inf %>% mutate(model = "DO") 
colnames(do_inf) <- c("var", "inf_val", "model")
inf_df <- rbind(inf_df, do_inf)

#agi model
agi_inf <- as.data.frame(ggBRT::ggInfluence(agi_mod_fin, plot = FALSE)) %>% rownames_to_column()
colnames(agi_inf) <- c("Predictor_variable", "relative_influence")
agi_inf$Predictor_variable <- gsub("\\<AGI_0m\\>", "AGI, daily, 0m", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("AGI_250m_ann", "AGI, annual, 250m", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("AGI_0m_seas", "AGI, seasonal, 0m", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("AGI_250m_seas", "AGI, seasonal, 250m", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("temp_mean", "temp", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("sal_mean", "sal", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("bathy_mean", "z", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("chl_mean", "chl-a", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("AGI_0m_ann", "AGI, annual, 0m", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("AGI_250m", "AGI, daily, 250m", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("ssh_mean", "SSH", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("mld_mean", "MLD", agi_inf$Predictor_variable)
agi_inf$Predictor_variable <- gsub("bathy_sd", "z_sd", agi_inf$Predictor_variable)

agi_inf <- agi_inf %>% mutate(model = "AGI") 
colnames(agi_inf) <- c("var", "inf_val", "model")
inf_df <- rbind(inf_df, agi_inf)

avg_inf_df <- avg_pred(base_mods = list.files(here("data/brt/mod_outputs/perf_metric_iters/base"), full.names = TRUE),
                       do_mods = list.files(here("data/brt/mod_outputs/perf_metric_iters/do"), full.names = TRUE), 
                       agi_mods = list.files(here("data/brt/mod_outputs/perf_metric_iters/agi"), full.names = TRUE)) 

##### revised plot  #####
avg_inf_sum <- avg_inf_df %>% 
  group_by(model, var) %>%
  summarise(inf_val = mean(inf_val))

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

do_agi <- avg_inf_sum %>% 
  filter(var_type != "Environmental predictor") %>% 
  mutate(var = as.factor(var)) %>%
  mutate(var = as.factor(var), var = fct_reorder(var, -inf_val)) %>%
  ggplot(aes(x = var, y = inf_val))+
  geom_bar(aes(fill = model), stat = "identity", position = position_dodge2(width = 0.8, preserve = "single"), color = "white")+
  scale_fill_manual(values = met.brewer("Hokusai1", direction = -1, n = 15)[2:4])+
  xlab("")+
  ylab("Relative importance (%)")+
  labs(fill = "Model")+
  facet_wrap(~var_type, scales = "free_x")+
 tidyquant::theme_tq()+
  theme(axis.text.y = element_text(size = 14, color = "black"), 
        axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
        axis.title = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14),
        legend.position = "none")

enviro <- avg_inf_sum %>%
  filter(var_type == "Environmental predictor") %>% 
  mutate(var = as.factor(var), var = fct_reorder(var, -inf_val)) %>%
  ggplot(aes(x = var, y = inf_val))+
  geom_bar(aes(fill = model), stat = "identity", position = "dodge", color = "white", width = 0.8)+
  scale_fill_manual(values = met.brewer("Hokusai1", direction = -1, n = 15))+
  xlab("")+
  ylab("Relative importance (%)")+
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 30, by = 10)) +
  labs(fill = "Model")+
  annotate(geom = "rect", fill = "#2c3e50", color = "#2c3e50", xmin = -Inf, xmax = +Inf, ymin = 34, ymax = +Inf)+
  annotate(label = "Environmental predictor", geom = "text", color = "white", x = 4, y = 38, size = 6)+
  tidyquant::theme_tq()+
  theme(axis.text = element_text(size = 14, color = "black"), 
        axis.title = element_text(size = 16, color = "black"),
        strip.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14), 
        legend.position = "top", 
        legend.justification = "left")

pred_fig <- enviro/do_agi
pred_fig

ggsave(here("figs/ms/fig6_pred/bar_pred.png"), pred_fig, width = 11, height = 8, units = c("in"))
  
#ggsave(here("figs/ms/rel_inf_pred.png"), all_pred,  height = 15, width = 15, units = c("in"))
ggsave(here("figs/ms/fig3_pred/base_pred.png"), base_pred, height = 7, width = 7, units = c("in"))
ggsave(here("figs/ms/fig3_pred/do_pred.png"), do_pred, height = 7, width = 7, units = c("in"))
ggsave(here("figs/ms/fig3_pred/agi_pred.png"), agi_pred, height = 7, width = 7, units = c("in"))
ggsave(here("figs/ms/fig3_pred/do_agi_pred.png"), do_agi_pred, height = 7, width = 7, units = c("in"))

### Figure 5: HSI maps study period ####
all_maps_avg <- hsi_maps_avg(rast_folder = "data/enviro/psat_spot_all/hsi_rasts/Jan03_Dec15", ms = "Y")
ggsave(here("figs/ms/fig7_hsi_all/all_maps_avg_20.png"), all_maps_avg, height = 5, width = 8, units = c("in"))

### Figure 6: ENSO HSI maps ####
#have to save using export button otherwise adds border, using height of 750 and width 500 (LN width 300)
#base year
enso_base <- hsi_maps_enso_avg(rast_folder = "data/enviro/psat_spot_all/hsi_rasts/Jan13_Dec13", enso = "diff")

#LN year 
enso_LN <- hsi_maps_difference_enso_avg(enso_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/LN_F_2010", neut_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/Jan13_Dec13", enso = "LN")

#EN year
enso_EN <- hsi_maps_difference_enso_avg(enso_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/EN_FW_Nov2014_Jan2015", neut_rast_folder = "data/enviro/psat_spot_all/hsi_rasts/Jan13_Dec13", enso = "EN")

#diet data 
#neutral year 
diet_neutral <- rmpq_prey_year2 %>% ungroup() %>%
  filter(perc_GII >= 1 & Year == 2013) %>% 
  top_n(3) %>%
  ggplot(aes(x = reorder(Common_Name, -perc_GII), y = perc_GII)) +
  geom_bar(stat = "identity", fill = "#224B5E", alpha = 0.85) +
  theme_minimal() + 
  theme(axis.text.y=element_text(size=20, color = "black"), 
        axis.title=element_text(size=22, color = "black"), 
        axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  xlab('')+
  ylab('% GII')+
  scale_x_discrete(position = "top") 
ggsave(here("figs/ms/fig7_enso_diet/diet_neutral.png"), diet_neutral, height = 4, width = 6, units = c("in"))


#LN year 
diet_LN <- rmpq_prey_year2 %>% ungroup() %>%
  filter(perc_GII >= 1 & Year == 2010) %>% 
  top_n(3) %>%
  ggplot(aes(x = reorder(Common_Name, -perc_GII), y = perc_GII)) +
  geom_bar(stat = "identity", fill = "#224B5E", alpha = 0.85) +
  theme_minimal() + 
  theme(axis.text=element_blank(), 
        axis.title=element_blank(), 
        #axis.text.y = element_blank(), 
        #axis.text.x = element_blank(),
        panel.grid = element_blank()) +
  xlab('')+
  ylab('% GII')+
  scale_x_discrete(position = "top") 
ggsave(here("figs/ms/fig7_enso_diet/diet_LN.png"), diet_LN, height = 4, width = 6, units = c("in"))


#EN year 
diet_EN <- rmpq_prey_year2 %>%
  filter(perc_GII >= 1 & Year == 2014) %>% 
  top_n(3) %>%
  ggplot(aes(x = reorder(Common_Name, -perc_GII), y = perc_GII)) +
  geom_bar(stat = "identity", fill = "#224B5E", alpha = 0.85) +
  theme_minimal() + 
  theme(axis.text=element_blank(), 
        axis.title=element_blank(), 
        #axis.text.y = element_blank(), 
        #axis.text.x = element_text(),
        panel.grid = element_blank()) +
  xlab('')+
  ylab('% GII')+
  scale_x_discrete(position = "top") 

ggsave(here("figs/ms/fig7_enso_diet/diet_EN.png"), diet_EN, height = 4, width = 6, units = c("in"))
