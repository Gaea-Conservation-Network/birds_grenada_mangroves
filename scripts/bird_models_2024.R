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


write.csv(data_long, "data/functionaldiversity_long_April2024.csv")

# ANCOVA Models ------------------------------------------------------------------
## bird diversity metric x basal area
data_long <- read.csv("data/functionaldiversity_long_April2024.csv")

##### bird species richness ####
colnames(data_long)

# interaction 
rich_basal <- lm(Species.Richness ~ type*area, data = data_long)
summary(rich_basal)

Anova(rich_basal, type = "3")

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
        legend.position = c(0.8, 0.9)) +
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

Anova(ab_basal, type = "3")

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

sw_basal <- lm(Species.Abundance ~ type*area, data = data_long)
summary(ab_basal)

Anova(ab_basal, type = "3")

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