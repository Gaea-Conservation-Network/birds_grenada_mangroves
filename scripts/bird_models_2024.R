library(tidyverse)
library(performance) # check model assumptions (check_model)
library(car) #ANOVA
library(patchwork)

# tree height based on procrustes
data <- read.csv("data/FunctionalDiversity.csv")
data <- data %>% select(SiteName:WhiteBasalArea)

# figure design
colours_man = c("Black mangrove" = "#22223b",
                "Red mangrove" = "#9e2a2b",
                "White mangrove" = "#dda15e")

shape = c("Black mangrove" = 21 ,
          "Red mangrove" = 22,
          "White mangrove" = 23)

# check data --------------------------------------------------------------
str(data)

hist(data$Species.Abundance)
hist(data$Species.Richness)



data_long <- data %>% pivot_longer(cols = BlackBasalArea:WhiteBasalArea,
                          names_to = "type",
                          values_to = "area")

data_long$type <- recode_factor(data_long$type,
                                "BlackBasalArea" = "Black mangrove",
                                "RedBasalArea" = "Red mangrove",
                                "WhiteBasalArea" = "White mangrove")
                                
data_long$logrichness <- log10(data_long$Species.Richness)
data_long$logabundance <- log10(data_long$Species.Abundance)


write.csv(data_long, "data/basal_models_2024/functionaldiversity_long_April2024.csv")

# ANCOVA Models ------------------------------------------------------------------
## bird diversity metric x basal area
data_long <- read.csv("data/basal_models_2024/functionaldiversity_long_April2024.csv")

##### bird species richness ####
colnames(data_long)

# interaction 
rich_basal <- lm(Species.Richness ~ type*area, data = data_long)
summary(rich_basal)

write.csv(Anova(rich_basal, type = "3"), "data/basal_models_2024/richness_basal_output.csv")

# Anova Table (Type III tests)
# 
# Response: Species.Richness
#              Sum Sq Df F value    Pr(>F)    
# (Intercept) 1431.46  1 30.6307 2.293e-06 ***
# type           8.69  2  0.0930    0.9114    
# area           1.17  1  0.0250    0.8751    
# type:area      8.02  2  0.0858    0.9180    
# Residuals   1822.58 39                      
# ---

# no interaction
rich_basal2 <- lm(logrichness ~ type + area, data = data_long)
summary(rich_basal2)

Anova(rich_basal2, type = "3")

# Anova Table (Type III tests)
# 
# Response: logrichness
#              Sum Sq Df  F value Pr(>F)    
# (Intercept) 15.8214  1 321.2629 <2e-16 ***
# type         0.0002  2   0.0024 0.9976    
# area         0.0014  1   0.0286 0.8664    
# Residuals    2.0192 41  

check_model(rich_basal2)


## Species Richness figure
richness_basal_fig <- data_long %>% ggplot(aes(y = Species.Richness, 
                                               x = area,
                         shape = type, fill = type)) +
  geom_point(size = 5, stroke = 1.5) +
  stat_smooth(method = lm, level = 0.95, aes(colour = type)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  scale_x_continuous(breaks = round(seq(min(data_long$area), max(data_long$area),
                                        by = 0.5),1)) +
  scale_fill_manual(values = colours_man) +
  scale_colour_manual(values = colours_man) +
  scale_shape_manual(values = shape) +
  labs(x = "Basal Area (m)",
       y = "Bird Species Richness")

##### abundance ####

ab_basal <- lm(Species.Abundance ~ type*area, data = data_long)
summary(ab_basal)

write.csv(Anova(ab_basal, type = "3"), "data/basal_models_2024/ab_basal_output.csv")

# Anova Table (Type III tests)
# 
# Response: Species.Abundance
# Sum Sq Df F value   Pr(>F)   
# (Intercept)  12222  1  9.7706 0.003343 **
#   type            31  2  0.0123 0.987766   
# area            68  1  0.0547 0.816345   
# type:area       15  2  0.0060 0.993984   
# Residuals    48784 39

# no interaction
ab_basal2 <- lm(logabundance ~ type + area, data = data_long)
summary(ab_basal2)

Anova(ab_basal2, type = "3")

# Anova Table (Type III tests)
# 
# Response: logabundance
# Sum Sq Df  F value Pr(>F)    
# (Intercept) 29.5550  1 352.3860 <2e-16 ***
#   type         0.0004  2   0.0024 0.9976    
# area         0.0024  1   0.0287 0.8662    
# Residuals    3.4387 41 

check_model(ab_basal2)


## Species Richness figure
ab_basal_fig <- data_long %>% ggplot(aes(y = Species.Abundance, x = area,
                                           shape = type, fill = type)) +
  geom_point(size = 5, stroke = 1.5) +
  stat_smooth(method = lm, level = 0.95, aes(colour = type)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = c(0.8, 0.9)) +
  scale_x_continuous(breaks = round(seq(min(data_long$area), max(data_long$area),
                                        by = 0.5),1)) +
  scale_fill_manual(values = colours_man) +
  scale_colour_manual(values = colours_man) +
  scale_shape_manual(values = shape) +
  labs(x = "Basal Area (m)",
       y = "Bird Abundance")





##### Shannon-Weiner ####
colnames(data_long)

sw_basal <- lm(Shannon.Weiner ~ type*area, data = data_long)
summary(sw_basal)

write.csv(Anova(sw_basal, type = "3"), "data/basal_models_2024/sw_basal_output.csv")
# 
# Anova Table (Type III tests)
# 
# Response: Shannon.Weiner
# Sum Sq Df  F value Pr(>F)    
# (Intercept) 45.410  1 214.0864 <2e-16 ***
#   type         0.064  2   0.1509 0.8605    
# area         0.030  1   0.1427 0.7076    
# type:area    0.110  2   0.2594 0.7729    
# Residuals    8.272 39 

# no interaction
sw_basal2 <- lm(Shannon.Weiner ~ type + area, data = data_long)
summary(sw_basal2)

Anova(sw_basal2, type = "3")

# Anova Table (Type III tests)
# 
# Response: Shannon.Weiner
# Sum Sq Df  F value Pr(>F)    
# (Intercept) 70.719  1 345.9078 <2e-16 ***
# type         0.004  2   0.0108 0.9892    
# area         0.026  1   0.1271 0.7232    
# Residuals    8.382 41

check_model(sw_basal2)


## Shannon.Weiner figure
sw_basal_fig <- data_long %>% ggplot(aes(y = Shannon.Weiner, x = area,
                                         shape = type, fill = type)) +
  geom_point(size = 5, stroke = 1.5) +
  stat_smooth(method = lm, level = 0.95, aes(colour = type)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  scale_x_continuous(breaks = round(seq(min(data_long$area), max(data_long$area),
                                        by = 0.5),1)) +
  scale_fill_manual(values = colours_man) +
  scale_colour_manual(values = colours_man) +
  scale_shape_manual(values = shape) +
  labs(x = "Basal Area (m)",
       y = "Shannon Weiner")


##### Func. Dispersion ####
colnames(data_long)

fd_basal <- lm(FunctionalDispersion ~ type*area, data = data_long)
summary(fd_basal)

write.csv(Anova(fd_basal, type = "3"), "data/basal_models_2024/fd_basal_output.csv")

# Anova Table (Type III tests)
# 
# Response: FunctionalDispersion
# Sum Sq Df  F value    Pr(>F)    
# (Intercept) 6290.3  1 145.6474 9.648e-15 ***
#   type          39.2  2   0.4541    0.6383    
# area          92.8  1   2.1480    0.1508    
# type:area     96.2  2   1.1142    0.3384    
# Residuals   1684.4 39


# no interaction
fd_basal2 <- lm(FunctionalDispersion ~ type + area, data = data_long)
summary(fd_basal2)

Anova(fd_basal2, type = "3")

# Anova Table (Type III tests)
# 
# Response: FunctionalDispersion
#             Sum Sq Df  F value  Pr(>F)    
# (Intercept) 10417.0  1 239.8618 < 2e-16 ***
# type           29.4  2   0.3390 0.71447    
# area          173.1  1   3.9863 0.05254 .  
# Residuals    1780.6 41 

check_model(fd_basal2) # kinda fucky 


## Functional dispersion figure
fd_basal_fig <- data_long %>% ggplot(aes(y = FunctionalDispersion, x = area,
                                         shape = type, fill = type)) +
  geom_point(size = 5, stroke = 1.5) +
  stat_smooth(method = lm, level = 0.95, aes(colour = type)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none") +
  scale_x_continuous(breaks = round(seq(min(data_long$area), max(data_long$area),
                                        by = 0.5),1)) +
  scale_fill_manual(values = colours_man) +
  scale_colour_manual(values = colours_man) +
  scale_shape_manual(values = shape) +
  labs(x = "Basal Area (m)",
       y = "Functional Dispersion")




##### Func. Divergence ####
colnames(data_long)

fdiv_basal <- lm(FunctionalDivergence ~ type*area, data = data_long)
summary(fdiv_basal)

write.csv(Anova(fdiv_basal, type = "3"), "data/basal_models_2024/fdiv_basal_output.csv")

# Anova Table (Type III tests)
# 
# Response: FunctionalDivergence
# Sum Sq Df   F value Pr(>F)    
# (Intercept) 4.5132  1 1834.8203 <2e-16 ***
# type        0.0025  2    0.4989 0.6110    
# area        0.0031  1    1.2776 0.2653    
# type:area   0.0052  2    1.0567 0.3573    
# Residuals   0.0959 39 

# no interaction
fdiv_basal2 <- lm(FunctionalDivergence ~ type + area, data = data_long)
summary(fd_basal2)

Anova(fdiv_basal2, type = "3")

# Anova Table (Type III tests)
# 
# Response: FunctionalDivergence
# Sum Sq Df   F value Pr(>F)    
# (Intercept) 7.0107  1 2842.2727 <2e-16 ***
# type        0.0001  2    0.0201 0.9801    
# area        0.0006  1    0.2366 0.6293    
# Residuals   0.1011 41

check_model(fdiv_basal2) 


## Functional dispersion figure
fdiv_basal_fig <- data_long %>% ggplot(aes(y = FunctionalDivergence, x = area,
                                         shape = type, fill = type)) +
  geom_point(size = 5, stroke = 1.5) +
  stat_smooth(method = lm, level = 0.95, aes(colour = type)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = c(0.8,0.9)) +
  scale_x_continuous(breaks = round(seq(min(data_long$area), max(data_long$area),
                                        by = 0.5),1)) +
  scale_fill_manual(values = colours_man) +
  scale_colour_manual(values = colours_man) +
  scale_shape_manual(values = shape) +
  labs(x = "Basal Area (m)",
       y = "Functional Divergence")



# Panel Figure ------------------------------------------------------------

div <- ab_basal_fig + richness_basal_fig / sw_basal_fig + plot_annotation(tag_levels = "A")

ggsave("figures/ab_s_sw_basal.png",
       div,
       height = 9.06,
       width = 14.6)


fun <- fd_basal_fig + fdiv_basal_fig + plot_annotation(tag_levels = "A")


ggsave("figures/fun_div_dis_basal.png",
       fun,
       height = 7.73,
       width = 13.8)
