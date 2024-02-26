library(janitor)
library(tidyverse)
library(lubridate)
# library(BiodiversityR)
library(ggplot2)

veg <- read.csv("data/Veg data for ivi updated.csv")
veg <- janitor::clean_names(veg)
veg$date <- dmy(veg$date)
veg$site <- as.factor(veg$site)
veg$transect <- as.factor(veg$transect)
veg$plot <- as.factor(veg$plot)
veg <- mutate(veg, plotID = str_c(site, transect, plot, sep = "0"))
veg$plotID <- as.factor(veg$plotID)

veg <- select(veg, -red_seedlings, - black_seedlings, -white_seedlings, -notes)

#basal area by size classes, calculated separately
a <- 2.01e-05
b <- 0.000246888
c <- 0.003148181
d <- 0.018132763
e <- 0.0322394

veg_long <- pivot_longer(veg, c(-site, -date, -start_time, -end_time, -transect, 
                                -plot, -position, -plotID, -redox_depth_cm), 
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

veg_all <- dplyr::select(veg_spp, forest = site, plotID, species, count, basal)

# importancevalue(veg_all, site = 'plotID', species = 'species', count = 'count', 
               # basal = 'basal', factor = 'forest', level = 'Levera')

# importancevalue.comp(veg_all, site="plotID", species="species",
                     # count="count", basal="basal",
                     # factor="forest")
#ivi code consistently returning errors: "Error in if (xi > xj) 1L else -1L : missing value where TRUE/FALSE needed
# In addition: Warning message:
#   In Ops.factor(xi, xj) : '>' not meaningful for factors"

#doing it manually
veg_all <- filter(veg_spp, count != 0)

# - Density in % = (no. of individuals of a species / total no. of individuals) x 100
total_no <- sum(veg_all$count)

# - Dominance in % = (total basal area of a species / basal area of all species) x 100
total_basal <- sum(veg_all$basal)

# - Freqr in % = (frequency of a species / sum-frequency of all species) x 100
total_plots <- n_distinct(veg_all$plotID)
summary <- veg_all %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/total_plots)
sum_freq <- sum(summary$freq)

summary <- veg_all %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/total_no)*100, 
            dominance = sum(basal), reldom = (dominance/total_basal)*100, 
            freq = n_distinct(plotID)/total_plots, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)
write.csv(summary, "data/overall ivi results.csv")


plot <- summary %>% ggplot()+
  geom_col(aes(x = species, y = ivi, fill = species), fill = c("black", "red", "gray"))+
  theme_bw()+
  labs(x = "Mangrove species", y = "Importance Value Index")
print(plot)

ggsave("figures/updated overall ivi plot.png", plot, width = 160, height = 100, unit = "mm")

sites <- veg_all %>% group_by(site) %>%
  summarize(total_no = sum(count), 
            total_basal = sum(basal), 
            total_plots = n_distinct(plotID))

Conf <- filter(veg_all, site %in% "Conference")
confsummary <- Conf %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/6)
sum_freq = sum(confsummary$freq)
confsummary <- Conf %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/42)*100, 
            dominance = sum(basal), reldom = (dominance/0.689957)*100, 
            freq = n_distinct(plotID)/6, relfreq = (freq/sum_freq)*100, 
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

write.csv(sitesummary, "data/ivi results by site.csv")

iviplot <- sitesummary %>% ggplot()+
  geom_col(aes(x = species, y = ivi, fill = species), 
           fill = c("black", "red", "gray", "black", "red", "gray", "black", "red", "gray", "black", "red", "gray"))+
  facet_wrap(.~site)+
  theme_bw()+
  labs(x = "Mangrove species", y = "Importance Value Index")
print(iviplot)

ggsave("figures/updated ivi plot by site.png", iviplot, width = 160, height = 100, unit = "mm")
