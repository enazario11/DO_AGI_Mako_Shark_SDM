
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

#### SF2: Model performance metrics scatter

#### SF3: DO, 0m, 2003-2015


#### SF4: DO, 250m, 2003-2015


#### SF5: DO, 0m, seasonal averages


#### SF6: DO, 250m, seasonal averages



#### SF7: temperature, 0m, seasonal averages


#### SF8: temperature, 0m, 2003-2015


#### SF9: spatiotemporal performance metrics



#### SF10: base model partial plots


#### SF11: DO model partial plots



#### SF12: AGI model partial plots



#### SF13: SHAP values


#### SF14: DO difference plots, 0 and 250m, neutral and LN year