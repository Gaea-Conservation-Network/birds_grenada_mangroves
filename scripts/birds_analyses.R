
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
bird_data <- read.csv("data/bird data.csv")
bird_data <- janitor::clean_names(bird_data)
bird_data <- filter(bird_data, site %in% c("Conference","Levera","Mt. Hartman","Westerhall"))


# make into a matrix ------------------------------------------------------

matrix <- data %>% 
  pivot_wider(names_from = "species",
              values_from = "abundance",
              values_fill = 0)


write.csv(matrix, "data/bird_matrix.csv")



uni <- read.csv("data/bird_univariate.csv")


# Univariate analyses -----------------------------------------------------

histogram(uni$rich)
histogram(uni$abundance)
histogram(uni$H)
histogram(uni$J)

# general linear model comparing among sites

# abundance

ab.mod <- lm(abundance ~ site, data = uni)

Anova(ab.mod, type = 3)

#Anova Table (Type III tests)

#Response: abundance
#            Sum Sq Df F value  Pr(>F)  
#(Intercept)  4641.3  1  4.6324 0.08400 .
#site        21228.3  3  7.0625 0.03015 *
#Residuals    5009.7  5

ab.hsd <- HSD.test(ab.mod, "site")

#            abundance groups
#Levera      154.50000      a
#Westerhall   65.00000     ab
#Conference   39.33333      b
#Mt. Hartman  24.50000      b

check_model(ab.mod) # residuals not super normal

# permutational lm

ab.lmp <- lmp(abundance ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(ab.lmp)

#Response: abundance
#          Sum Sq Df F value  Pr(>F)  
#site1     21228.3  3  7.0625 0.03015 *
#Residuals  5009.7  5 


## species richness

s.mod <- lm(rich ~ site, data = uni)
Anova(s.mod, type = 3)

#Response: rich
#Sum Sq Df F value  Pr(>F)  
#(Intercept) 261.333  1 13.4477 0.01449 *
#site         91.056  3  1.5618 0.30879  
#Residuals    97.167  5  

check_model(s.mod) # also a mess

# permutation lm

s.lmp <- lmp(rich ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(s.lmp)

#Response: rich
#Sum Sq Df F value Pr(>F)
#site1     91.056  3  1.5618 0.3088
#Residuals 97.167  5 


## Shannon-weiner
# GLM

h.mod <- lm(H ~ site, data = uni)
Anova(h.mod, type = 3)

#Response: H
#          Sum Sq Df F value    Pr(>F)    
#(Intercept) 9.8604  1 94.7316 0.0001946 ***
#site        0.3258  3  1.0434 0.4495016    
#Residuals   0.5204  5  

check_model(h.mod)

# permutation lm

h.lmp <- lmp(H ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(h.lmp)

#Response: H
#Sum Sq Df F value Pr(>F)
#site1     0.32583  3  1.0434 0.4495
#Residuals 0.52044  5  


# Pielou's evenness

j.mod <- lm(J ~ site, data = uni)
Anova(j.mod, type = 3)

#Response: J
#Sum Sq Df  F value    Pr(>F)    
#(Intercept) 1.98855  1 386.5701 6.285e-06 ***
#  site        0.04236  3   2.7446    0.1525    
#Residuals   0.02572  5  

check_model(j.mod)

# permutational lm

j.lmp <- lmp(J ~ site, data = uni, 
    perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(j.lmp)

#Response: J
#Sum Sq Df F value Pr(>F)
#site1     0.042355  3  2.7446 0.1525
#Residuals 0.025720  5 


# Multivariate analyses ---------------------------------------------------
env <- matrix %>% select(site, date)
spp <- matrix %>% select(BANA:SHCO)

spp.rel <- decostand(spp, "max", 2, na.rm = NULL) # rel by column max


# betadisper --------------------------------------------------------------

spp.b <- vegdist(spp.rel, method = "bray")

groups <- factor(matrix$site)

(dispersion <- betadisper(spp.b, groups))

#Average distance to median:
#Conference      Levera Mt. Hartman  Westerhall 
#0.3165      0.2867      0.4006      0.3543 

anova(dispersion)

#Response: Distances
#          Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     3 0.015004 0.0050012  0.7957 0.5467
#Residuals  5 0.031427 0.0062854  

boxplot(dispersion)
plot(dispersion)



# perMANOVA ---------------------------------------------------------------


spp.pmv <- adonis2(spp.rel ~ site,
                   data = env,
                   method = "bray")

#adonis2(formula = spp.rel ~ site, data = env, method = "bray")
#          Df SumOfSqs      R2      F Pr(>F)   
#site      3   1.2074 0.53423 1.9117  0.003 **
#Residual  5   1.0527 0.46577                 
#Total     8   2.2601 1.00000  

(adonis.pair(spp.b, groups,
             nper = 1000,
             corr.method = "bonferroni"))

#                 combination SumsOfSqs   MeanSqs  F.Model        R2   P.value P.value.corrected
#1      Conference <-> Levera 0.4301938 0.4301938 3.735109 0.5545729 0.1000000               0.6
#2 Conference <-> Mt. Hartman 0.3731016 0.3731016 3.114176 0.5093370 0.1000000               0.6
#3  Conference <-> Westerhall 0.2693335 0.2693335 1.646667 0.3543759 0.1000000               0.6
#4     Levera <-> Mt. Hartman 0.5594699 0.5594699 4.540916 0.6942324 0.3333333               1.0
#5      Levera <-> Westerhall 0.2618101 0.2618101 1.386417 0.4094052 0.6666667               1.0
#6 Mt. Hartman <-> Westerhall 0.2561766 0.2561766 1.308445 0.3954864 0.3333333               1.0



# mv glm ------------------------------------------------------------

sp.mv <- mvabund(spp)

sp.mvmod <- manyglm(sp.mv ~ site,
                    data = env, family = "negativebinomial")

output.mvmod <- anova(sp.mvmod, p.uni = "adjusted")

#Multivariate test:
#              Res.Df Df.diff   Dev Pr(>Dev)  
#(Intercept)      8                         
#site             5       3 199.5    0.039 *

p.uni <- as.data.frame(output.mvmod$uni.p) %>% t
p.uni

write.csv(p.uni, "data/mvabund_species_puni.csv")

# no specific species respond strongly to site but there are comm. differences among sites



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

#Data:     spp.rel 
#Distance: bray 
#
#Dimensions: 3 
#Stress:     0.03284544 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

layout(matrix(1:2, ncol = 2))
plot(nms, main = "Bird NMDS plot"); stressplot(nms, main = "Shepard plot")
layout(1)

ordiplot(nms, type = "n")
orditorp(nms, display = "species")
orditorp(nms, display = "sites")

# how many iterations of the NMDS
nms$iters # 73

(g <- goodness(nms)) 
sum(g^2)
nms$stress^2  # 0.001078823

1-nms$stress^2 # 0.9989212 #analogous to square correlation coefficient


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
  filter(X.alltaxa12.vectors..r > 0.3) %>% 
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
  filter(X.alltaxa13.vectors..r > 0.3) %>% 
  rownames_to_column("species")

target13 <- corr.sp13$species # string of the Family names

axis13.vectors <- spp.rel %>% select(all_of(target13)) # make a matrix of just those

(nmds.vectors.13 <- envfit(nms$points, axis13.vectors,
                           permutations = 999, choices = c(1,3)))                        

corr.vectors.13 <- as.data.frame(nmds.vectors.13$vectors$arrows*sqrt(nmds.vectors.13$vectors$r)) #scaling vectors
corr.vectors.13$species <- rownames(corr.vectors.13) # add Family as a column

write.csv(corr.vectors.13, "data/bird_nmds/NMDS_correlatedvectors_axis13.csv")


# NMDS figure -------------------------------------------------------------


shape = c("Conference" = 21,
          "Levera" = 22,
          "Mt. Hartman" = 23,
          "Westerhall" =24)

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
                   aes(x = MDS1, y = MDS2, label = species),
                   color="black",
                   size = 5) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shape) +
  theme(legend.position = "none")


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
                   size = 5) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  scale_colour_viridis(discrete = TRUE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = shape) +
  theme(legend.title = element_blank())


nms <- axis12.p + axis13.p

ggsave("figures/birds_nmds.jpeg",
       nms)

