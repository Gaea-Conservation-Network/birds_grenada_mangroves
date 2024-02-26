library(tidyverse)
library(viridis)

species.data <- read.csv("data/bird_spp_sum.csv")
species.data <- species.data[2:36]

species.long <- species.data %>% pivot_longer(!site, names_to = "species",
                              values_to = "count")



species.long.abundant <- species.long %>% filter(count >= 5)

write.csv(species.long.abundant, "data/species_long_abundance.csv")

species.long.ab <- read.csv("data/species_long_abundance.csv")

species.gg <-ggplot(species.long.ab, aes(fill = label, y = count, x = site)) + 
  geom_bar(stat = "identity") +
  theme_classic(14) +
  labs(x = " ",
       y = "Count of most abundant birds (>= 5)") +
  theme(panel.border = element_rect(fill = NA),
        legend.position = "none") +
  theme(axis.text = element_text(size = 16)) +
  scale_fill_viridis(discrete = TRUE,
                     alpha = 0.5) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 5)
ggsave("figures/species_abundant.jpeg")
            



spp <- species.long %>% 
  group_by(species) %>% 
  summarize(total = sum(count))
