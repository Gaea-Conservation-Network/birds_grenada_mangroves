
# packages ----------------------------------------------------------------
library(tidyverse)
library(vegan)

library(ggrepel)
library(patchwork)

# data --------------------------------------------------------------------

matrix <- read.csv("data/bird_matrix.csv", row.names = 1)

traits <- read.csv("data/traits_matrix.csv", row.names = 1)

trt.rel <- decostand(traits, "max", 2, na.rm = NULL) # rel by column max
trt.rel$sites <- matrix$site

trts <- trt.rel[,1:34]


# figure out number of dimension
k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
spp.nms <- metaMDSdist(trts)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(spp.nms, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 2D looks good

#### NMDS analysis 

set.seed(120) 

nms <- metaMDS(trts, distance = "bray", # species data, bray-curtis dissimilarity
               autotransform = FALSE,  # NMDS will do autotransformations for you
               k = 2, trymax = 1000)   # k = number of axes
nms

#Data:     trts 
#Distance: bray 
#
#Dimensions: 2 
#Stress:     0.06967904 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

layout(matrix(1:2, ncol = 2))
plot(nms, main = "Bird Trait NMDS plot"); stressplot(nms, main = "Shepard plot")
layout(1)

ordiplot(nms, type = "n")
orditorp(nms, display = "species")
orditorp(nms, display = "sites")

# how many iterations of the NMDS
nms$iters # 38

(g <- goodness(nms)) 
sum(g^2)
nms$stress^2  # 0.004855168

1-nms$stress^2 # 0.9951448 #analogous to square correlation coefficient


# NMDS plotting -----------------------------------------------------------

## extract the scores for plotting 
scr <- as.data.frame(scores(nms, display = "sites")) # extract NMDS scores

# adding categorical info to scores
scr$site <- matrix$site

scores <- scr

write.csv(scores,"data/trait_nmds/NMDS_scores.csv") # save this as a csv

## species correlated with axis 1 & 2

alltaxa12 <- envfit(nms, trts,
                    choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa12$vectors)$arrows,
                          (alltaxa12$vectors)$r,
                          (alltaxa12$vectors)$pvals) #take list and make into data frame


corr.sp12 <- all.taxa.df %>% 
  filter(X.alltaxa12.vectors..r > 0.3) %>% 
  rownames_to_column("species")

target12 <- corr.sp12$species # string of the Family names

axis12.vectors <- trts %>% select(all_of(target12)) # make a matrix of just those

(nmds.vectors.12 <- envfit(nms$points, axis12.vectors,
                           permutations = 999, choices = c(1,2)))                        

corr.vectors.12 <- as.data.frame(nmds.vectors.12$vectors$arrows*sqrt(nmds.vectors.12$vectors$r)) #scaling vectors
corr.vectors.12$species <- rownames(corr.vectors.12) # add Family as a column

write.csv(corr.vectors.12, "data/trait_nmds/NMDS_correlatedvectors_axis12.csv")

# NMDS figure -------------------------------------------------------------


shape = c("Conference" = 21,
          "Levera" = 22,
          "Mt. Hartman" = 23,
          "Westerhall" =24)

scores <- read.csv("data/trait_nmds/NMDS_scores.csv")
axis12 <- read.csv("data/trait_nmds/NMDS_correlatedvectors_axis12.csv")

unique(scores$site)


axis12.plot <- ggplot(data = scores,
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
  scale_shape_manual(values = shape) 

ggsave("figures/traits_nmds.jpeg",
       axis12.plot)

