
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

bird_data <- filter(bird_data, site %in% c("Conference","Levera","Mt. Hartman","Westerhall"))


# make data frames  ------------------------------------------------------

colnames(bird_data)

matrix1 <- bird_data %>% 
  mutate(date = dmy(date)) %>% 
  group_by(site, site_id, date, species) %>% 
  summarize_at(vars(abundance), sum) %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

write.csv(matrix1, "data/bird_matrix.csv")

## summing across visits

matrix_sum <- bird_data %>% 
  group_by(site, site_id, species) %>% 
  summarize_at(vars(abundance), sum) %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

write.csv(matrix_sum, "data/bird_matrix_sum.csv")


matrix_sum_uni <- matrix_sum %>% 
  mutate(abundance = rowSums(across(ANCH:SHCO)),
         richness = rowSums(across(ANCH:SHCO) > 0),
         H = diversity(across(ANCH:SHCO), index = "shannon"),
         J = H/log(specnumber(across(ANCH:SHCO))))
  
sum_univariate <- matrix_sum_uni %>% 
  select(site, site_id, abundance:J)

write.csv(sum_univariate, "data/bird_univariate_sum.csv")

## keeping in the dates

matrix_uni_time <- matrix1 %>% 
  mutate(abundance = rowSums(across(BANA:SHCO)),
         richness = rowSums(across(BANA:SHCO) > 0),
         H = diversity(across(BANA:SHCO), index = "shannon"),
         J = H/log(specnumber(across(BANA:SHCO))))

univariate_time <- matrix_uni %>% 
  select(site:visit,abundance:J)

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
#              Sum Sq Df F value   Pr(>F)    
#  (Intercept) 6962.0  1 47.6154 3.57e-06 ***
#  site        2607.7  3  5.9449 0.006351 ** 
#  Residuals   2339.4 16                 

check_model(ab.mod) # residuals not super normal

# permutational lm

ab.lmp <- lmp(abundance ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(ab.lmp)

#Anova Table (Type II tests)
#
#Response: abundance
#Sum Sq Df F value   Pr(>F)   
#site1     2607.7  3  5.9449 0.006351 **
#Residuals 2339.4 1

site.ab <- HSD.test(ab.lmp, "site")

# abundance groups
# Conference   59.00000      a
# Westerhall   35.25000     ab
# Levera       25.06000      b
# Mt. Hartman  16.33333      b


## species richness

s.mod <- lm(richness ~ site, data = uni)
Anova(s.mod, type = 3)

#Response: richness
#              Sum Sq Df F value    Pr(>F)    
#(Intercept) 288.000  1 29.6884 5.358e-05 ***
#site         89.738  3  3.0835    0.0572 .  
#Residuals   155.212 16   

check_model(s.mod) # also a mess

# permutation lm

s.lmp <- lmp(richness ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(s.lmp)

# Response: richness
# Sum Sq Df F value Pr(>F)  
# site1      89.738  3  3.0835 0.0572 .
# Residuals 155.212 16  


## Shannon-weiner
# GLM

h.mod <- lm(H ~ site, data = uni)
Anova(h.mod, type = 3)

#Response: H
#              Sum Sq Df F value    Pr(>F)    
#  (Intercept) 8.1745  1 93.1971 4.474e-08 ***
#  site        1.2651  3  4.8076   0.01424 *  
#  Residuals   1.4034 16

check_model(h.mod)

# permutation lm

h.lmp <- lmp(H ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(h.lmp)

#Response: H
#            Sum Sq Df F value  Pr(>F)  
#  site1     1.2651  3  4.8076 0.01424 *
#  Residuals 1.4034 16  

h.hsd <- HSD.test(h.lmp, "site")

# H groups
# Westerhall  2.114178      a
# Mt. Hartman 2.038855     ab
# Conference  2.021697     ab
# Levera      1.566060      b


# Pielou's evenness

j.mod <- lm(J ~ site, data = uni)
Anova(j.mod, type = 3)

#Response: J
#             Sum Sq Df  F value    Pr(>F)    
#(Intercept) 1.33158  1 441.8948 4.434e-13 ***
#site        0.01670  3   1.8474    0.1793    
#Residuals   0.04821 16 

check_model(j.mod)

# permutational lm

j.lmp <- lmp(J ~ site, data = uni, 
    perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(j.lmp)

# Response: J
# Sum Sq Df F value Pr(>F)
# site1     0.016701  3  1.8474 0.1793
# Residuals 0.048213 16   


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
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  annotate("text", x = c(1:4), y = c(70, 50, 55, 55),
           label = c("a","b","b","ab"),
           size = 6)


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
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
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
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  annotate("text", x = c(1:4), y = 2.6,
           label = c("ab","b","ab","a"),
           size = 6)

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
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),
        legend.position = "none") 


panel <- ab.box + s.box + sw.box + J.box +
  plot_annotation(tag_levels = "A")


ggsave("figures/bird_abundance_boxplots.jpeg",
       panel)



# Multivariate analyses ---------------------------------------------------

env <- matrix %>% select(site, site_id)
spp <- matrix %>% select(ANCH:SHCO)

spp.rel <- decostand(spp, "max", 2, na.rm = NULL) # rel by column max


# betadisper --------------------------------------------------------------

spp.b <- vegdist(spp.rel, method = "bray")

groups <- factor(matrix$site)

(dispersion <- betadisper(spp.b, groups))

#Average distance to median:
#Conference   Levera Mt. Hartman  Westerhall 
#0.2328      0.5249      0.4718      0.3417 

anova(dispersion)

#Response: Distances
#          Df  Sum Sq  Mean Sq F value   Pr(>F)   
#Groups     3 0.20439 0.068129   7.236 0.002767 **
#Residuals 16 0.15064 0.009415  

boxplot(dispersion)
plot(dispersion)


# perMANOVA ---------------------------------------------------------------

spp.pmv <- adonis2(spp.rel ~ site,
                   data = env,
                   method = "bray")

#adonis2(formula = spp.rel ~ site, data = env, method = "bray")
#         Df SumOfSqs      R2      F Pr(>F)  
#site      3   1.3667 0.23676 1.6544  0.015 *
#Residual 16   4.4058 0.76324                
#Total    19   5.7725 1.00000  

(adonis.pair(spp.b, groups,
             nper = 1000,
             corr.method = "bonferroni"))

# combination SumsOfSqs   MeanSqs  F.Model        R2    P.value P.value.corrected
# 1      Conference <-> Levera 0.3019544 0.3019544 1.017578 0.0846741 0.43556444         1.0000000
# 2 Conference <-> Mt. Hartman 0.4510855 0.4510855 1.733594 0.3662321 0.10000000         0.6000000
# 3  Conference <-> Westerhall 0.4237486 0.4237486 2.933754 0.4231119 0.06666667         0.4000000
# 4     Levera <-> Mt. Hartman 0.5644450 0.5644450 1.769403 0.1285025 0.02997003         0.1798202
# 5      Levera <-> Westerhall 0.5611739 0.5611739 2.012383 0.1340482 0.01898102         0.1138861
# 6 Mt. Hartman <-> Westerhall 0.3271773 0.3271773 1.432890 0.2227444 0.13686314         0.8211788

# mv glm ------------------------------------------------------------

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

shapes = c("Conference" = 21,
           "Levera" = 22,
           "Mt. Hartman" = 23,
           "Westerhall" = 24)


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
       x = ' ') +
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

#Dimensions: 3 
#Stress:     0.1706756 
#Stress type 1, weak ties
#Two convergent solutions found after 49 tries

layout(matrix(1:2, ncol = 2))
plot(nms, main = "Bird NMDS plot"); stressplot(nms, main = "Shepard plot")
layout(1)

ordiplot(nms, type = "n")
orditorp(nms, display = "species")
orditorp(nms, display = "sites")

# how many iterations of the NMDS
nms$iters # 145

(g <- goodness(nms)) 
sum(g^2)
nms$stress^2  # 0.02913016

1-nms$stress^2 # 0.9708698 #analogous to square correlation coefficient


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
  filter(X.alltaxa12.vectors..r > 0.2) %>% 
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
  filter(X.alltaxa13.vectors..r > 0.2) %>% 
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
                   aes(x = MDS1, y = MDS2, label = spp),
                   color="black",
                   size = 5) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shapes) +
  theme(legend.position = "none") +
  ylim(-1, 1) +
  xlim(-1.5, 1.5)


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
                   aes(x = MDS1, y = MDS3, label = spp),
                   color="black",
                   size = 5) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shapes) +
  theme(legend.title = element_blank()) +
  ylim(-1, 1) +
  xlim(-1.5, 1.5)


nms <- axis12.p + axis13.p

ggsave("figures/birds_nmds.jpeg",
       nms)

