
# libraries ---------------------------------------------------------------

library(tidyverse)

# glm 
library(car)
library(agricolae)
library(lmPerm)

library(vegan)

library(performance)
library(see)

# perMANOVA pair-wise comparison
library(devtools)
install_github("GuillemSalazar/EcolUtils")

library(EcolUtils)

# mv glm

library(mvabund)

# for figures

library(ggrepel)
library(patchwork)
library(viridis)


# load data ---------------------------------------------------------------
bird_data <- read.csv("data/bird_data_independent reps.csv")
bird_data <- janitor::clean_names(bird_data)

str(bird_data)
unique(bird_data$site)
unique(bird_data$site_id_ind)

bird_data <- filter(bird_data, site %in% c("Conference","Levera","Mt. Hartman","Westerhall"))


# make data frames  ------------------------------------------------------

colnames(bird_data)

matrix1 <- bird_data %>% 
  mutate(date = dmy(date)) %>% 
  group_by(site, site_id_ind, date, species) %>% 
  summarize_at(vars(abundance), sum) %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

write.csv(matrix1, "data/bird_matrix.csv")

## summing across visits

matrix_sum <- bird_data %>% 
  group_by(site, site_id_ind, species) %>% 
  summarize_at(vars(abundance), sum) %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

write.csv(matrix_sum, "data/bird_matrix_sum.csv")


spp_sum <- bird_data %>% 
  group_by(site, species) %>% 
  summarize_at(vars(abundance), sum) %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

write.csv(spp_sum, "data/bird_spp_sum.csv")



matrix_sum_uni <- matrix_sum %>% 
  mutate(abundance = rowSums(across(ANCH:SHCO)),
         richness = rowSums(across(ANCH:SHCO) > 0),
         H = diversity(across(ANCH:SHCO), index = "shannon"),
         J = H/log(specnumber(across(ANCH:SHCO))))
  
sum_univariate <- matrix_sum_uni %>% 
  select(site, site_id_ind, abundance:J)

write.csv(sum_univariate, "data/bird_univariate_sum.csv")

## keeping in the dates

matrix_uni_time <- matrix1 %>% 
  mutate(abundance = rowSums(across(BANA:SHCO)),
         richness = rowSums(across(BANA:SHCO) > 0),
         H = diversity(across(BANA:SHCO), index = "shannon"),
         J = H/log(specnumber(across(BANA:SHCO))))

univariate_time <- matrix_uni_time %>% 
  select(site:date,abundance:J)

write.csv(univariate_time, "data/bird_univariate_time.csv")
write.csv(matrix1, "data/bird_matrix_time.csv")

# import data -------------------------------------------------------------

matrix <- read.csv("data/bird_matrix_sum.csv", row.names = 1)
uni <- read.csv("data/bird_univariate_sum.csv", row.names = 1)


# Univariate analyses -----------------------------------------------------

histogram(uni$rich)
histogram(uni$abundance)
histogram(uni$H)
histogram(uni$J)

ggplot(uni, aes(x = site, y = abundance)) +
  geom_boxplot() +
  geom_jitter()

ggplot(uni, aes(x = site, y = richness)) +
  geom_boxplot() +
  geom_jitter()

# general linear model comparing among sites

# abundance

ab.mod <- lm(abundance ~ site, data = uni)

Anova(ab.mod, type = 3)

#Anova Table (Type III tests)
#
#Response: abundance
#             Sum Sq Df F value Pr(>F)  
#(Intercept) 3605.3  1  7.6262 0.0185 *
#site         162.4  3  0.1145 0.9498  
#Residuals   5200.3 11     
                
check_model(ab.mod) # residuals not super normal

# permutational lm

ab.lmp <- lmp(abundance ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(ab.lmp)

#Anova Table (Type II tests)
#
#Response: abundance
#          Sum Sq Df F value Pr(>F)
#site1      162.4  3  0.1145 0.9498
#Residuals 5200.3 11 



## species richness

s.mod <- lm(richness ~ site, data = uni)
Anova(s.mod, type = 3)

#Anova Table (Type III tests)
#
#Response: richness
#             Sum Sq Df F value  Pr(>F)  
#(Intercept) 243.00  1  7.4926 0.01933 *
#site         18.98  3  0.1951 0.89754  
#Residuals   356.75 11

check_model(s.mod) # also a mess

# permutation lm

s.lmp <- lmp(richness ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(s.lmp)

#Anova Table (Type II tests)
#
#Response: richness
#           Sum Sq Df F value Pr(>F)
#site1      18.98  3  0.1951 0.8975
#Residuals 356.75 11 


## Shannon-weiner
# GLM

h.mod <- lm(H ~ site, data = uni)
Anova(h.mod, type = 3)

#Anova Table (Type III tests)
#
#Response: H
#Sum Sq Df F value    Pr(>F)    
#(Intercept) 9.7546  1 40.1095 5.574e-05 ***
#  site        0.5884  3  0.8065    0.5161    
#Residuals   2.6752 11 

check_model(h.mod)

# permutation lm

h.lmp <- lmp(H ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(h.lmp)

#Anova Table (Type II tests)
#
#Response: H
#Sum Sq Df F value Pr(>F)
#site1     0.58842  3  0.8065 0.5161
#Residuals 2.67518 11 

# Pielou's evenness

j.mod <- lm(J ~ site, data = uni)
Anova(j.mod, type = 3)

#Anova Table (Type III tests)
#
#Response: J
#Sum Sq Df  F value    Pr(>F)    
#(Intercept) 2.13973  1 630.3949 4.603e-11 ***
#  site        0.01148  3   1.1273    0.3802    
#Residuals   0.03734 11 

check_model(j.mod)

# permutational lm

j.lmp <- lmp(J ~ site, data = uni, 
    perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(j.lmp)

#Anova Table (Type II tests)
#
#Response: J
#Sum Sq Df F value Pr(>F)
#site1     0.011479  3  1.1273 0.3802
#Residuals 0.037337 11  



# Univariate figures ------------------------------------------------------

ab.box <- ggplot(uni, aes(x = site, y = abundance)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 5,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Abundance",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") 


s.box <- ggplot(uni, 
                aes(x = site, y = richness)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 5,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Species Richness",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") 

sw.box <- ggplot(uni, 
                 aes(x = site, y = H)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 5,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Shannon-Weiner (H)",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") 

J.box <- ggplot(uni, 
                aes(x = site, y = J)) +
  geom_boxplot(lwd = 1) +
  geom_jitter(aes(fill = site),
              size = 5,
              shape = 21,
              stroke = 1.5) +
  labs(y = "Pielou's Evenness (J)",
       x = " ") +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") 


panel <- ab.box + s.box + sw.box + J.box 

ggsave("figures/bird_abundance_boxplots.jpeg",
       panel)



# Multivariate analyses ---------------------------------------------------

env <- matrix %>% select(site, site_id_ind)
spp <- matrix %>% select(ANCH:SHCO)

spp.rel <- decostand(spp, "max", 2, na.rm = NULL) # rel by column max


# betadisper --------------------------------------------------------------

spp.b <- vegdist(spp.rel, method = "bray")

groups <- factor(matrix$site)

(dispersion <- betadisper(spp.b, groups))

#Average distance to median:

#Conference      Levera Mt. Hartman  Westerhall 
#0.3646      0.5428      0.4265      0.3196 


anova(dispersion)

#Analysis of Variance Table
#
#Response: Distances
#          Df  Sum Sq  Mean Sq F value  Pr(>F)  
#Groups     3 0.13785 0.045951  3.6696 0.04719 *
#Residuals 11 0.13774 0.012522

boxplot(dispersion)
plot(dispersion)


# perMANOVA ---------------------------------------------------------------

spp.pmv <- adonis2(spp.rel ~ site,
                   data = env,
                   method = "bray")

#adonis2(formula = spp.rel ~ site, data = env, method = "bray")
#          Df SumOfSqs      R2     F Pr(>F)  
#site      3   1.2722 0.29537 1.537  0.059 .
#Residual 11   3.0348 0.70463               
#Total    14   4.3070 1.00000  

# mv glm ------------------------------------------------------------

# dates are all over the place so this doesn't really make sense, but this is the code/results
env <- matrix %>% select(site, site_id)
env$site <- as.factor(env$site)

spp <- matrix %>% select(ANCH:SHCO)


spp.mv <- mvabund(spp)

plot(spp.mv ~ env$site) # check out dist

spp.mod <- manyglm(spp.mv ~ site,
                     data = env, family = "negative.binomial")

plot(spp.mod)
output.sppmod <- anova(spp.mod, p.uni = "adjusted")

#Multivariate test:
#            Res.Df Df.diff   Dev Pr(>Dev)    
#(Intercept)     19                           
#site            16       3 187.1    0.001 ***

p.uni <- as.data.frame(output.sppmod$uni.p) %>% t
p.uni


write.csv(p.uni, "data/mvabund_species_puni.csv")


long <- matrix %>% 
  pivot_longer(ANCH:SHCO,
               names_to = "species",
               values_to = "count")



raw.ab <- ggplot(long, aes(x = site, y = count, 
                              fill = site,
                              shape = site)) +
  geom_boxplot() +
  geom_jitter(size = 3, alpha = 0.8) +
  labs(y = "Bird Abundance",
       x = ' ') +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shapes) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") 

ggsave("figures/raw_bird_abundances.jpeg")

# community x time --------------------------------------------------------


# this is by time, rather than summed across visits 

# I set Feb 15 as time "0" and calculated duration manually bc I am impatient
# could probably do it with lubridate and diffdate with a set starting point

matrix.time <- read.csv("data/bird_matrix_time.csv",row.names = 1)

spp.time <- matrix.time %>% select(BANA:SHCO)
env.time <- matrix.time %>% 
  select(site, site_id, date, time)

env.time$site <- as.factor(env.time$site)

meanvar.plot(spp.time)


# negative binomial
sp.mv <- mvabund(spp.time)

plot(sp.mv ~ env.time$site) # check out dist

str(env.time)


sp.mvmod2 <- manyglm(sp.mv ~ site * time,
                     data = env.time, family = "poisson")

plot(sp.mvmod2)
output.mvmod <- anova(sp.mvmod, p.uni = "adjusted")

#Multivariate test:
#  Res.Df Df.diff    Dev Pr(>Dev)    
#(Intercept)     38                            
#site            35       3 189.92    0.001 ***
#  time            34       1  50.65    0.072 .  
#site:time       31       3 101.13    0.001 ***
#  ---


sp.mvmod <- manyglm(sp.mv ~ site * time,
                    data = env.time, family = "negative.binomial")

plot(sp.mvmod) # go with this one
output.mvmod <- anova(sp.mvmod, p.uni = "adjusted")

#Multivariate test:
#            Res.Df Df.diff    Dev Pr(>Dev)    
#(Intercept)     38                            
#site            35       3 189.92    0.001 ***
#time            34       1  50.65    0.077 .  
#site:time       31       3 101.13    0.001 ***


### figure 

time.long <- matrix.time %>% 
  mutate(date = ymd(date)) %>% 
  pivot_longer(ANCH:SHCO,
               names_to = "species",
               values_to = "count")

shapes = c("Conference" = 21,
           "Levera" = 22,
           "Mt. Hartman" = 23,
           "Westerhall" = 24)


time <- ggplot(time.long, aes(x = time, y = count, 
                    fill = site,
                    shape = site)) +
  geom_jitter(size = 3, alpha = 0.8) +
  facet_wrap(~site) +
  labs(y = "Bird Abundance",
       x = "Days since first survey") +
  theme_bw() +
  geom_smooth(aes(colour = site)) +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shapes) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") 

ggsave("figures/significant_bird_abundances.jpeg")



byspecies <- ggplot(time.long, aes(x = time, y = count, 
                              fill = site,
                              shape = site)) +
  geom_jitter(size = 3, alpha = 0.8) +
  facet_wrap(~species) +
  labs(y = "Bird Abundance",
       x = 'Days since first survey') +
  theme_bw() +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shapes) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11)) 

# NMDS ordination ---------------------------------------------------------

# figure out number of dimension
k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
spp.nms <- metaMDSdist(spp.rel)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(spp.nms, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D looks good

#### NMDS analysis 

set.seed(120) 

nms <- metaMDS(spp.rel, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms

#metaMDS(comm = spp.rel, distance = "bray", k = 3, trymax = 1000,      autotransform = FALSE) 
#
#global Multidimensional Scaling using monoMDS
#
#Data:     spp.rel 
#Distance: bray 
#
#Dimensions: 3 
#Stress:     0.07515997 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

layout(matrix(1:2, ncol = 2))
plot(nms, main = "Bird NMDS plot"); stressplot(nms, main = "Shepard plot")
layout(1)

ordiplot(nms, type = "n")
orditorp(nms, display = "species")
orditorp(nms, display = "sites")

# how many iterations of the NMDS
nms$iters # 82

(g <- goodness(nms)) 
sum(g^2)
nms$stress^2  # 0.005649022

1-nms$stress^2 # 0.994351 #analogous to square correlation coefficient


# NMDS plotting -----------------------------------------------------------

## extract the scores for plotting 
scr <- as.data.frame(scores(nms, display = "sites")) # extract NMDS scores

# adding categorical info to scores
env$NMDS1 <- scr$NMDS1
env$NMDS2 <- scr$NMDS2
env$NMDS3 <- scr$NMDS3

scores <- env

write.csv(scores,"data/bird_nmds/NMDS_scores.csv") # save this as a csv

## species correlated with axis 1 & 2

alltaxa12 <- envfit(nms, spp.rel,
                    choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa12$vectors)$arrows,
                          (alltaxa12$vectors)$r,
                          (alltaxa12$vectors)$pvals) #take list and make into data frame


corr.sp12 <- all.taxa.df %>% 
  filter(X.alltaxa12.vectors..r > 0.4) %>% 
  rownames_to_column("species")

target12 <- corr.sp12$species # string of the Family names

axis12.vectors <- spp.rel %>% select(all_of(target12)) # make a matrix of just those

(nmds.vectors.12 <- envfit(nms$points, axis12.vectors,
                           permutations = 999, choices = c(1,2)))                        

corr.vectors.12 <- as.data.frame(nmds.vectors.12$vectors$arrows*sqrt(nmds.vectors.12$vectors$r)) #scaling vectors
corr.vectors.12$species <- rownames(corr.vectors.12) # add Family as a column

write.csv(corr.vectors.12, "data/bird_nmds/NMDS_correlatedvectors_axis12.csv")

## Species correlated with axis 1 & 3

alltaxa13 <- envfit(nms, spp.rel,
                    choices = c(1,3)) #produces a list with r3, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa13$vectors)$arrows,
                          (alltaxa13$vectors)$r,
                          (alltaxa13$vectors)$pvals) #take list and make into data frame


corr.sp13 <- all.taxa.df %>% 
  filter(X.alltaxa13.vectors..r > 0.4) %>% 
  rownames_to_column("species")

target13 <- corr.sp13$species # string of the Family names

axis13.vectors <- spp.rel %>% select(all_of(target13)) # make a matrix of just those

(nmds.vectors.13 <- envfit(nms$points, axis13.vectors,
                           permutations = 999, choices = c(1,3)))                        

corr.vectors.13 <- as.data.frame(nmds.vectors.13$vectors$arrows*sqrt(nmds.vectors.13$vectors$r)) #scaling vectors
corr.vectors.13$species <- rownames(corr.vectors.13) # add Family as a column

write.csv(corr.vectors.13, "data/bird_nmds/NMDS_correlatedvectors_axis13.csv")


# NMDS figure -------------------------------------------------------------

scores <- read.csv("data/bird_nmds/NMDS_scores.csv")
axis12 <- read.csv("data/bird_nmds/NMDS_correlatedvectors_axis12.csv")
axis13 <- read.csv("data/bird_nmds/NMDS_correlatedvectors_axis13.csv")

unique(scores$site)


shapes = c("Conference" = 21,
           "Levera" = 22,
           "Mt. Hartman" = 23,
           "Westerhall" = 24)

axis12.p <- ggplot(data = scores,
       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = scores, 
             aes(x = NMDS1, y = NMDS2, 
                 shape = site, fill = site),
             size = 6, stroke = 1.5) +
  geom_segment(data = axis12 , 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "black") +
  geom_label_repel(data = axis12, 
                   aes(x = MDS1, y = MDS2, label = species),
                   color= "black",
                   size = 5,
                   force = 3,
                   box.padding = 1) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shapes) +
  theme(legend.position = "none") +
  ylim(-1.75, 1.75) +
  xlim(-1.75, 1.75)


axis13.p <- ggplot(data = scores,
                   aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = scores, 
             aes(x = NMDS1, y = NMDS3, 
                 shape = site, fill = site),
             size = 6, stroke = 1.5) +
  geom_segment(data = axis13 , 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "black") +
  geom_label_repel(data = axis13, 
                   aes(x = MDS1, y = MDS3, label = species),
                   color="black",
                   size = 5,
                   force = 2,
                   box.padding = 1) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shapes) +
  theme(legend.title = element_blank()) +
  ylim(-1.75, 1.75) +
  xlim(-1.75, 1.75)


nms <- axis12.p + axis13.p
nms

ggsave("figures/birds_nmds.jpeg",
       nms)




# traits nmds -------------------------------------------------------------

traits <- read.csv("data/traits_matrix1.csv")



