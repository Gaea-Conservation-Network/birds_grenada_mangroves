
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
bird_data <- filter(bird_data, site %in% c("Conference","Levera","Mt. Hartman","Westerhall"))


# make data frames  ------------------------------------------------------

colnames(bird_data)

matrix <- bird_data %>% 
  group_by(site, site_id, date, visit, species) %>% 
  summarize_at(vars(abundance), sum) %>% 
  pivot_wider(names_from = species,
              values_from = abundance,
              values_fill = 0)

write.csv(matrix, "data/bird_matrix.csv")


matrix_uni <- matrix %>% 
  mutate(abundance = rowSums(across(BANA:SHCO)),
         richness = rowSums(across(BANA:SHCO) > 0),
         H = diversity(across(BANA:SHCO), index = "shannon"),
         J = H/log(specnumber(across(BANA:SHCO))))

univariate <- matrix_uni %>% 
  select(site:visit,abundance:J)

write.csv(univariate, "data/bird_univariate.csv")


# import data -------------------------------------------------------------

matrix <- read.csv("data/bird_matrix.csv")
uni <- read.csv("data/bird_univariate.csv", row.names = 1)


# Univariate analyses -----------------------------------------------------

histogram(uni$rich)
histogram(uni$abundance)
histogram(uni$H)
histogram(uni$J)

ggplot(uni, aes(x = site, y = abundance)) +
  geom_violin(trim = FALSE)

ggplot(uni, aes(x = site, y = richness)) +
  geom_violin(trim = FALSE)

# general linear model comparing among sites

# abundance

ab.mod <- lm(abundance ~ site, data = uni)

Anova(ab.mod, type = 3)

#Anova Table (Type III tests)
#
#Response: abundance
#            Sum Sq Df F value    Pr(>F)    
#(Intercept) 2320.7  1 21.5308 4.734e-05 ***
#site         541.3  3  1.6739    0.1904    
#Residuals   3772.4 35                      

check_model(ab.mod) # residuals not super normal

# permutational lm

ab.lmp <- lmp(abundance ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(ab.lmp)

#Anova Table (Type II tests)

#Response: abundance
#Sum Sq Df F value Pr(>F)
#site1      541.3  3  1.6739 0.1904
#Residuals 3772.4 35


## species richness

s.mod <- lm(richness ~ site, data = uni)
Anova(s.mod, type = 3)

#Response: richness
#              Sum Sq Df F value    Pr(>F)    
#(Intercept) 240.667  1 44.1110 1.114e-07 ***
#site         64.632  3  3.9487   0.01586 *  
#Residuals   190.958 35 

check_model(s.mod) # also a mess

# permutation lm

s.lmp <- lmp(richness ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(s.lmp)

#Response: richness
#          Sum Sq Df F value  Pr(>F)  
#site1      64.632  3  3.9487 0.01586 *
#Residuals 190.958 35 

## Shannon-weiner
# GLM

h.mod <- lm(H ~ site, data = uni)
Anova(h.mod, type = 3)

#Response: H
#                Sum Sq Df F value    Pr(>F)    
#  (Intercept) 13.5652  1 93.5024 2.027e-11 ***
#  site         2.5101  3  5.7672  0.002583 ** 
#  Residuals    5.0778 35 

h.hsd <- HSD.test(h.mod, "site")

#$groups
#H groups
#Westerhall  1.800013      a
#Mt. Hartman 1.700604      a
#Conference  1.503619     ab
#Levera      1.183062      b

check_model(h.mod)

# permutation lm

h.lmp <- lmp(H ~ site, data = uni, 
             perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(h.lmp)

#Response: H
#          Sum Sq Df F value   Pr(>F)   
#site1     2.5101  3  5.7672 0.002583 **
#Residuals 5.0778 35  

# Pielou's evenness

j.mod <- lm(J ~ site, data = uni)
Anova(j.mod, type = 3)

#Response: J
#Sum Sq Df  F value Pr(>F)    
#(Intercept) 4.1379  1 653.8733 <2e-16 ***
#  site        0.0398  3   2.0981 0.1188    
#Residuals   0.2152 34 

check_model(j.mod)

# permutational lm

j.lmp <- lmp(J ~ site, data = uni, 
    perm = "Prob", Ca = 0.0001, maxIter = 999)

Anova(j.lmp)


#Response: J
#Sum Sq Df F value Pr(>F)
#site1     0.039832  3  2.0981 0.1188
#Residuals 0.215163 34  

# Multivariate analyses ---------------------------------------------------
env <- matrix %>% select(site, site_id, date)
spp <- matrix %>% select(BANA:SHCO)

spp.rel <- decostand(spp, "max", 2, na.rm = NULL) # rel by column max


# betadisper --------------------------------------------------------------

spp.b <- vegdist(spp.rel, method = "bray")

groups <- factor(matrix$site)

(dispersion <- betadisper(spp.b, groups))

#Average distance to median:
#Conference  Levera   Mt. Hartman    Westerhall 
#0.4424      0.5819      0.5510      0.4578 

anova(dispersion)

#Response: Distances
#Df  Sum Sq  Mean Sq F value   Pr(>F)   
#Groups     3 0.13658 0.045525  5.0195 0.005342 **
#  Residuals 35 0.31744 0.009070 

boxplot(dispersion)
plot(dispersion)



# perMANOVA ---------------------------------------------------------------


spp.pmv <- adonis2(spp.rel ~ site,
                   data = env,
                   method = "bray")

#adonis2(formula = spp.rel ~ site, data = env, method = "bray")
#          Df SumOfSqs     R2      F Pr(>F)    
#site      3    2.008 0.1468 2.0073  0.001 ***
#Residual 35   11.671 0.8532                  
#Total    38   13.679 1.0000  

(adonis.pair(spp.b, groups,
             nper = 1000,
             corr.method = "bonferroni"))
#combination SumsOfSqs   MeanSqs  F.Model         R2     P.value P.value.corrected
#1      Conference <-> Levera 0.6877700 0.6877700 2.023639 0.07221187 0.035964036        0.21578422
#2 Conference <-> Mt. Hartman 0.5734069 0.5734069 1.857361 0.17106929 0.024975025        0.14985015
#3  Conference <-> Westerhall 0.8030135 0.8030135 3.197026 0.24225352 0.006993007        0.04195804
#4     Levera <-> Mt. Hartman 0.6191196 0.6191196 1.689874 0.06331518 0.056943057        0.34165834
#5      Levera <-> Westerhall 0.7904896 0.7904896 2.311239 0.08163681 0.006993007        0.04195804
#6 Mt. Hartman <-> Westerhall 0.3874040 0.3874040 1.230096 0.12024285 0.252747253        1.00000000

# mv glm ------------------------------------------------------------

sp.mv <- mvabund(spp)

sp.mvmod <- manyglm(sp.mv ~ site,
                    data = env, family = "negativebinomial")

output.mvmod <- anova(sp.mvmod, p.uni = "adjusted")


#Multivariate test:
#             Res.Df  Df.diff  Dev   Pr(>Dev)    
#(Intercept)     38                         
#site            35       3    175    0.001 ***

p.uni <- as.data.frame(output.mvmod$uni.p) %>% t
p.uni


write.csv(p.uni, "data/mvabund_species_puni.csv")

#LBHE, EADO, SPTH are significantly related to a certain site

shapes = c("Conference" = 21,
           "Levera" = 22,
           "Mt. Hartman" = 23,
           "Westerhall" = 24)

sig.birds <- bird_data %>% 
  filter(species %in% c("LBHE", "EADO", "SPTH"))

ab <- ggplot(sig.birds, aes(x = site, y = abundance)) +
  geom_jitter(aes(shape = site,
             fill = site),
             size = 4) +
  facet_wrap(~species, ncol =4) +
  scale_shape_manual(values = shapes) +
  labs(x = " ", y = "Abundance") +
  theme(legend.position = "none")

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

