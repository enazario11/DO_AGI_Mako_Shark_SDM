#load packages
library(tidyverse)
library(RColorBrewer)
library(tidyquant)
library(here)

options(scipen=999) #scientific notation turned off

##### Data loading and tidying #############
#load diet and location data files
dat_coll <- read.csv("data/diet/collection.csv")
dat_pred <- read.csv("data/diet/predator.csv")
dat_prey <- read.csv("data/diet/prey_comp.csv")
dat_prey_size <- read.csv("data/diet/prey_size.csv")

  #limits diet data to SCB Ecoregion
scbE_coll <- dat_coll %>% 
  filter(Latitude < 47.33 & Latitude > 2.76 & Longitude > -150.94 & Longitude < -105.63)

  #filters predator DF based on IDs from SCB Ecoregion described above. Filters out FL that weren't in the loc DF
dat_pred <- dat_pred %>% #filters out measurement entries that are blank
  filter(Predator_Measurement_2 == "FL") 

scbE_coll_id <- scbE_coll$Collection_ID #478 IDs
scbE_pred <- dat_pred[dat_pred$Collection_ID %in% scbE_coll_id, ]

FL_femaleD <- scbE_pred %>%
  filter(Predator_Sex == "F") %>%
  filter(Predator_Length_2 < 249 & Predator_Length_2 > 100)

FL_maleD <- scbE_pred %>%
  filter(Predator_Sex == "M") %>%
  filter(Predator_Length_2 < 180 & Predator_Length_2 > 100)

scbE_pred2 <- rbind(FL_femaleD, FL_maleD)

  #filters prey DF based on IDs from SCB and FL filtered predator DF 
scbE_pred_id <- scbE_pred2$Predator_ID #371 IDs
scbE_prey <- dat_prey[dat_prey$Predator_ID %in% scbE_pred_id, ]

#### RMPQ and GII Calculations -- all years#####
scbE_prey <- scbE_prey %>%
  mutate(freq_metric = 1)

rmpq_prey <- scbE_prey %>%
  group_by(Common_Name) %>%
  summarise(N = sum(Prey_N), 
            W = sum(Prey_Wt), 
            Freq = sum(freq_metric), 
            perc_N = (N/1353)*100, 
            perc_W = (W/90099.18)*100, 
            perc_F = (Freq/264)*100, 
            GII = (perc_N + perc_W + perc_F)/sqrt(3), 
            perc_GII = (perc_N + perc_W + perc_F)/3)

##### Annual diet data summaries #######
#open scbE_coll, scbE_pred2, and scbE_prey DFs
#seasonal comparisons not possible (diet data collected Aug - Jan). Complete annual changes in GII (all species). Then do same but for top 5-10 species
#Annual comparisons all -- GII (geometric index of importance)

#Match dfs by IDs and bring over the collection year to the prey DF
yr_coll <- scbE_coll %>%
  subset(select = c("Collection_ID", "Year"))

scbE_pred3 <- merge(scbE_pred2, yr_coll, by = "Collection_ID", all = T)
yr_pred <- scbE_pred3 %>%
  subset(select = c("Year", "Predator_ID"))
scbE_prey2 <- merge(scbE_prey, yr_pred, by = "Predator_ID", all = T) %>%
  drop_na(Prey_ID)

#calculate the RMPQ values by species and year
rmpq_prey_year <- scbE_prey2 %>%
  group_by(Year) %>%
  mutate(pred_metric = length(unique(Predator_ID))) %>%
  ungroup() %>%
  group_by(Year, Common_Name) %>%
  summarise(N = sum(Prey_N), 
            W = sum(Prey_Wt), 
            Freq = sum(freq_metric), 
            pred_metric = unique(pred_metric)) %>%
  ungroup() %>%
  filter(Year != 1998)

#calculate the GII by year using the by species data
rmpq_prey_year2 <- rmpq_prey_year %>%
  group_by(Year) %>%
  mutate(perc_N = (N/sum(N))*100, 
         perc_W = (W/sum(W))*100, 
         perc_F = (Freq/pred_metric)*100,
         GII = (perc_N + perc_W + perc_F)/sqrt(3), 
         perc_GII = (perc_N + perc_W + perc_F)/3)

#plot all prey data across years in my study
#per year
by_year <- rmpq_prey_year2 %>%
  filter(GII > 4.06 & Year <= 2014 & Year >= 2003) %>% #25% quantile value
  ggplot(aes(x = reorder(Common_Name, -GII), y = GII)) +
  geom_bar(stat = "identity", fill = "#92351e") +
  facet_wrap(~Year, scales = "free_x")+
  theme_tq() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95)) +
  xlab('') 

#pooled years
all_years <- rmpq_prey %>%
  filter(perc_GII >= 1) %>% 
  ggplot(aes(x = reorder(Common_Name, -perc_GII), y = perc_GII)) +
  geom_bar(stat = "identity", fill = "#92351e") +
  theme_tq() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2, hjust = 0.95), 
        axis.text=element_text(size=14), 
        axis.title=element_text(size=14)) +
  xlab('')+
  ylab('% GII')

