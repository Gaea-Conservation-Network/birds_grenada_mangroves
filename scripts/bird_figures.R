
# libraries ---------------------------------------------------------------

library(tidyverse)

library(lubridate)

library(vegan)

# for figures

library(patchwork)
library(viridis)

# data --------------------------------------------------------------------

bird_data <- read.csv("data/bird data.csv")
bird_data <- janitor::clean_names(bird_data)
bird_data <- filter(bird_data, site %in% c("Conference","Levera","Mt. Hartman","Westerhall"))


# select relevant columns & sum plots for each visit/date

data <- bird_data %>% 
  select(site,date,plot:abundance) %>% 
  group_by(site, date, species) %>% 
  summarize_at(vars(abundance), sum) 

write.csv(data, "data/bird_summed.csv")

matrix <- read.csv("data/bird_matrix.csv")

# Univariate data ---------------------------------------------------------

matrix$rich <- specnumber(matrix[,3:36]) # species richness in vegan
matrix$abundance <- rowSums(matrix[,3:36]) # abundance
matrix$H <- diversity(matrix[,3:36]) # Shannon Weiner
matrix$D1 <- diversity(matrix[,3:36], index = "simpson") # 1 - D
matrix$J <- matrix$H/log(specnumber(matrix[,3:36]))


univariate <- matrix %>% select(site,date,rich:J)

write.csv(univariate , "data/bird_univariate.csv")


# boxplots --------------------------------------------------------------
uni <- read.csv("data/bird_univariate.csv", row.names = 1)

ab.box <- ggplot(uni, aes(x = site, y = abundance)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 4,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Abundance",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none") 

s.box <- ggplot(uni, 
                aes(x = site, y = richness)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 4,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Species Richness",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  ylim(0,30)

sw.box <- ggplot(uni, 
                aes(x = site, y = H)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 4,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Shannon-Weiner (H)",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  ylim(0,3)

J.box <- ggplot(uni, 
                 aes(x = site, y = J)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 4,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Pielou's Evenness (J)",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  ylim(0,1)



panel <- ab.box + s.box + sw.box + J.box +
  plot_annotation(tag_levels = "A")


ggsave("figures/bird_abundance_boxplots.jpeg",
       panel)



# community comp plot -----------------------------------------------------

unique(data$site)

spp <- ggplot(data, aes(x = site, y = abundance + 1)) +
  geom_boxplot(aes(fill = site), alpha = 0.8,
             shape = 21) +
  scale_y_log10() +
  labs(y = "Species count",
       x = ' ') +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none")

spp.facet <- ggplot(data, aes(x = site, y = abundance)) +
  geom_point(aes(fill = site),
             size = 6, alpha = 0.8,
             shape = 21) +
  facet_wrap(~species) +
  labs(y = "Abundance",
       x = ' ') +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none")


ggsave("figures/bird_spp_facet.jpeg",
       spp.facet,
       width = 30,
       height = 20)

