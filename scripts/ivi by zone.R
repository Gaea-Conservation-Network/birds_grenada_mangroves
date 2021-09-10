plots <- veg_all %>% group_by(plot) %>%
  summarize(total_no = sum(count), 
            total_basal = sum(basal), 
            total_plots = n_distinct(plotID))

plot1 <- filter(veg_all, plot %in% 1)
plot1summary <- plot1 %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/15)
sum_freq <- sum(plot1summary$freq)
plot1summary <- plot1 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/218)*100, 
            dominance = sum(basal), reldom = (dominance/2.566978)*100, 
            freq = n_distinct(plotID)/15, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)

plot2 <- filter(veg_all, plot %in% 2)
plot2summary <- plot2 %>% group_by(species) %>%
  summarize(freq = n_distinct(plotID)/15)
sum_freq <- sum(plot2summary$freq)
plot2summary <- plot2 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/191)*100, 
            dominance = sum(basal), reldom = (dominance/2.053532)*100, 
            freq = n_distinct(plotID)/15, relfreq = (freq/sum_freq)*100, 
            ivi = reldensity + reldom + relfreq)     

plot3 <- filter(veg_all, plot %in% 3)
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
  geom_col(aes(x = species, y = ivi, fill = species), 
           fill = c("black", "red", "gray", "black", "red", "gray", "black", "red", "gray"))+
  facet_wrap(.~zone)+
  theme_bw()+
  labs(x = "Mangrove species", y = "Importance Value Index")

print(zoneplot)

ggsave("figures/updated ivi plots by zone.png", zoneplot, width = 160, height = 100, unit = "mm")
