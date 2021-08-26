library(janitor)
library(tidyverse)
library(lubridate)
library(ggplot2)

veg <- read.csv("data/veg data for ivi updated.csv")
veg <- janitor::clean_names(veg)
veg$date <- dmy(veg$date)
veg$start_time <-hm(veg$start_time)
veg$end_time <-hm(veg$end_time)
veg$site <- as.factor(veg$site)
veg$transect <- as.factor(veg$transect)
veg$plot <- as.factor(veg$plot)
veg <- mutate(veg, plotID = str_c(site, transect, plot, sep = "0"))
veg$plotID <- as.factor(veg$plotID)
veg <- mutate(veg, position = factor(position, levels = c("seaward", "interior", "landward")))
veg$redox_depth_cm <- as.numeric(veg$redox_depth_cm)

veg <- select(veg, -red_seedlings, - black_seedlings, -white_seedlings, -notes)

#size classes
a <- 2.01e-05
b <- 0.000246888
c <- 0.003148181
d <- 0.018132763
e <- 0.0322394

veg_long <- pivot_longer(veg, c(-site, -date, -start_time, -end_time, -transect, 
                                -plot, -plotID, -position, -redox_depth_cm), 
                         names_to = "size", values_to = "count")

size_a <- filter(veg_long, size %in% c("red_a", "white_a", "black_a"))
size_a <- mutate(size_a, basal = count*a)

size_b <- filter(veg_long, size %in% c("red_b", "white_b", "black_b"))
size_b <- mutate(size_b, basal = count*b)

size_c <- filter(veg_long, size %in% c("red_c", "white_c", "black_c"))
size_c <- mutate(size_c, basal = count*c)

size_d <- filter(veg_long, size %in% c("red_d", "white_d", "black_d"))
size_d <- mutate(size_d, basal = count*d)

size_e <- filter(veg_long, size %in% c("red_e", "white_e", "black_e"))
size_e <- mutate(size_e, basal = count*e)

veg_size <- bind_rows(size_a, size_b, size_c, size_d, size_e)

spp_red <- filter(veg_size, size %in% c("red_a", "red_b", "red_c", "red_d", "red_e"))
spp_red <- mutate(spp_red, species = "red")

spp_black <- filter(veg_size, size %in% c("black_a", "black_b", "black_c", "black_d", "black_e"))
spp_black <- mutate(spp_black, species = "black")

spp_white <- filter(veg_size, size %in% c("white_a", "white_b", "white_c", "white_d", "white_e"))
spp_white <- mutate(spp_white, species = "white")

veg_spp <- bind_rows(spp_red, spp_black, spp_white)
veg_spp$species <- as.factor(veg_spp$species)

veg_all <- filter(veg_spp, count != 0)

#ivi = rel density (abundance) + rel dominance (area) + rel frequency (occurrence)
# - Density in % = (no. of individuals of a species / total no. of individuals) x 100
total_no <- sum(veg_all$count)

# - Dominance in % = (total basal area of a species / basal area of all species) x 100
total_basal <- sum(veg_all$basal)

# - Freqr in % = (frequency of a species / sum-frequency of all species) x 100
total_plots <- n_distinct(veg_all$plotID)

#overall ivi
summary <- veg_all %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/total_plots)
sum_freq <- sum(summary$freq)
summary <- veg_all %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/total_no)*100, 
            dominance = sum(basal), reldom = (dominance/total_basal)*100, 
            freq = n_distinct(plotID)/total_plots, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)

plot <- summary %>% ggplot()+
  geom_col(aes(x = species, y = ivi, fill = species))+
  theme_bw()+
  labs(x = "Mangrove species", y = "Importance Value Index")
print(plot)

ggsave("figures/updated overall plot.png", plot, width = 160, height = 100, unit = "mm")

#by site
sites <- veg_all %>% group_by(site) %>%
  summarize(total_no = sum(count), 
            total_basal = sum(basal), 
            total_plots = n_distinct(plotID))

Conf <- filter(veg_all, site %in% "Conference")
confsummary <- Conf %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/8)
sum_freq = sum(confsummary$freq)
confsummary <- Conf %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/162)*100, 
            dominance = sum(basal), reldom = (dominance/1.070891)*100, 
            freq = n_distinct(plotID)/8, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)
confsummary <- mutate(confsummary, site = "Conference")

Lev <- filter(veg_all, site %in% "Levera")
levsummary <- Lev %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/21)
sum_freq = sum(levsummary$freq)
levsummary <- Lev %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/229)*100, 
            dominance = sum(basal), reldom = (dominance/3.8554640)*100, 
            freq = n_distinct(plotID)/21, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)
levsummary <- mutate(levsummary, site = "Levera")

Hart <- filter(veg_all, site %in% "Mt. Hartman")
hartsummary <- Hart %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/6)
sum_freq = sum(hartsummary$freq)
hartsummary <- Hart %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/109)*100, 
            dominance = sum(basal), reldom = (dominance/0.403547)*100, 
            freq = n_distinct(plotID)/6, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)
hartsummary <- mutate(hartsummary, site = "Mt. Hartman")

West <- filter(veg_all, site %in% "Westerhall")
westsummary <- West %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/12)
sum_freq = sum(westsummary$freq)
westsummary <- West %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/284)*100, 
            dominance = sum(basal), reldom = (dominance/1.851027)*100, 
            freq = n_distinct(plotID)/12, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)
westsummary <- mutate(westsummary, site = "Westerhall")

sitesummary <- bind_rows(confsummary, levsummary, hartsummary, westsummary)

iviplot <- sitesummary %>% ggplot()+
  geom_col(aes(x = species, y = ivi, fill = species))+
  facet_wrap(.~site)+
  theme_bw()+
  labs(x = "Mangrove species", y = "Importance Value Index")
print(iviplot)

ggsave("figures/updated ivi plot by site.png", iviplot, width = 160, height = 100, unit = "mm")

#by plot/zone
plots <- veg_all %>% group_by(position) %>%
  summarize(total_no = sum(count), 
            total_basal = sum(basal), 
            total_plots = n_distinct(plotID))

plot1 <- filter(veg_all, position %in% "seaward")
plot1summary <- plot1 %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/15)
sum_freq <- sum(plot1summary$freq)
plot1summary <- plot1 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/322)*100, 
            dominance = sum(basal), reldom = (dominance/2.688635)*100, 
            freq = n_distinct(plotID)/15, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)

plot2 <- filter(veg_all, position %in% "interior")
plot2summary <- plot2 %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/17)
sum_freq <- sum(plot2summary$freq)
plot2summary <- plot2 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/207)*100, 
            dominance = sum(basal), reldom = (dominance/2.312809)*100, 
            freq = n_distinct(plotID)/17, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)     

plot3 <- filter(veg_all, position %in% "landward")
plot3summary <- plot3 %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/15)
sum_freq <- sum(plot3summary$freq)
plot3summary <- plot3 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/255)*100, 
            dominance = sum(basal), reldom = (dominance/2.179485)*100, 
            freq = n_distinct(plotID)/15, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)  

plotivis <- bind_rows(plot1summary, plot2summary, plot3summary, .id = "zone")

zoneplot <- plotivis %>% ggplot()+
  geom_col(aes(x = species, y = ivi, fill = species))+
  facet_wrap(.~zone)+
  theme_bw()+
  labs(x = "Mangrove species", y = "Importance Value Index")
print(zoneplot)

ggsave("figures/updated ivi plots by zone.png", zoneplot, width = 160, height = 100, unit = "mm")

