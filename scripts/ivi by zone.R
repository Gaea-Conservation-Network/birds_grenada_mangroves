plots <- veg_all %>% group_by(plot) %>%
  summarize(total_no = sum(count), 
            total_basal = sum(basal), 
            total_plots = n_distinct(plotID))

plot1 <- filter(veg_all, plot %in% 1)
plot1summary <- plot1 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/218)*100, 
            dominance = sum(basal), reldom = (dominance/2.566978)*100, 
            freq = n_distinct(plotID)/15, relfreq = (freq/1.6001)*100, 
            ivi = reldensity + reldom + relfreq)

plot2 <- filter(veg_all, plot %in% 2)
plot2summary <- plot2 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/192)*100, 
            dominance = sum(basal), reldom = (dominance/2.191676)*100, 
            freq = n_distinct(plotID)/15, relfreq = (freq/1.8)*100, 
            ivi = reldensity + reldom + relfreq)     

plot3 <- filter(veg_all, plot %in% 3)
plot3summary <- plot3 %>% group_by(species) %>%
  summarize(density = sum(count), reldensity = (density/199)*100, 
            dominance = sum(basal), reldom = (dominance/1.849592)*100, 
            freq = n_distinct(plotID)/15, relfreq = (freq/1.733)*100, 
            ivi = reldensity + reldom + relfreq)  

plotivis <- bind_rows(plot1summary, plot2summary, plot3summary, .id = "zone")

zoneplot <- plotivis %>% ggplot()+
  geom_col(aes(x = species, y = ivi, fill = species))+
  facet_wrap(.~zone)+
  theme_bw()+
  labs(x = "Mangrove species", y = "Importance Value Index")
print(zoneplot)

ggsave("outputs/ivi plots by zone.png", zoneplot, width = 160, height = 100, unit = "mm")
