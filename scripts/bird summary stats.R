library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(RAM)
library(vegan)

bird_data <- read.csv("data/bird data.csv")
bird_data <- janitor::clean_names(bird_data)
bird_data <- filter(bird_data, site %in% c("Conference","Levera","Mt. Hartman","Westerhall"))
    # 5 observations missing site and date information. Either from Westerhall or Conference

bird_data$site <- as.factor(bird_data$site)
bird_data$date <- dmy(bird_data$date)
bird_data$species <- as.factor(bird_data$species)
bird_data$transect <- as.factor(bird_data$transect)
bird_data$plot <- as.factor(bird_data$plot)
bird_data <- mutate(bird_data, plotID = str_c(site, transect, plot, sep = "0"))

total_richness <- n_distinct(bird_data$species)
total_abundance <- sum(bird_data$abundance)

summary <- bird_data %>% group_by(site) %>%
  summarize(n_counts = n_distinct(date),
            richness = n_distinct(species),
            abundance = sum(abundance))

rank <- bird_data %>% group_by(species, site) %>%
  summarize(abundance = sum(abundance))

rank %>% ggplot()+
  geom_col(aes(x = site, y = abundance)) +
  facet_wrap(. ~ species, scales = "free_y")

plot <- bird_data %>% group_by(plot) %>%
  summarize(richness = n_distinct(species), 
            abundance = sum(abundance))
