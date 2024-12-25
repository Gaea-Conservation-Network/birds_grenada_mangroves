######################################################################################################
###################################### NMDS: Number of Axis ##########################################
###################################################################################################### 


# Function to determine the optimal number of axes for an NMDS #

# We use this function to determine what the optimal number of axis are in a NMDS
# This function saves the stress from each NMDS. From axis 1 to 10, we repeat the NMDS noRuns number of times
NMDS.ScreePoints<-function(data,noRuns) { #where data is the name of the data frame, noRuns is the number of runs
  runr = matrix(NA, nrow = 10, ncol = noRuns)
  for (i in 1:10) {
    runr[i,] = replicate(noRuns, metaMDS(data,autotransform=T,trymax=100,k=i)$stress)
    
  }
  return(runr)
}

######################################################################################################
######################################################################################################
############################# Biological Interactions  ###############################################
######################################################################################################
######################################################################################################

CongruenceCrossTaxon <- function(DataVeg, DataBird, ColTraitVeg){
  Veg <- vegdist(DataVeg %>% arrange(Independent)%>% ungroup() %>%select(contains(ColTraitVeg)))
  Bird <- vegdist(DataBird %>% arrange(Independent)%>% ungroup() %>%select(ANCH:ZEND))
  
 MantelOverall <- mantel(vegdist(wisconsin(Veg)),
                         vegdist(wisconsin(Bird)),method="spearman")
  
  
  ProtestOverall <- protest(decostand(Veg,method = "hellinger"),
                            decostand(Bird,method = "hellinger"))
  Result <- tibble(Group = "Overall",
                   Test = c("Mantel", "Procustes"),
                   Strength = c(MantelOverall$statistic, ProtestOverall$t0),
                   `P-Value` = c(MantelOverall$signif, ProtestOverall$signif))
  Result2 <-
    map_dfr(.x = unique(DataVeg$SiteName), ~{
      
      Veg <- vegdist(DataVeg %>% arrange(Independent)%>% ungroup() %>% filter(SiteName == .x)%>%
                       select(contains(ColTraitVeg)))
      Bird <- vegdist(DataBird %>% arrange(Independent)%>% ungroup() %>% filter(SiteName == .x)%>%
                        select(ANCH:ZEND))
      
      MantelOverall <- mantel(vegdist(wisconsin(Veg)),vegdist(wisconsin(Veg)),method="spearman")
      
      ProtestOverall <- protest(decostand(Veg,method = "hellinger"),
                                decostand(Bird,method = "hellinger"))
      Result <- tibble(Group = .x,
                       Test = c("Mantel", "Procustes"),
                       Strength = c(MantelOverall$statistic, ProtestOverall$t0),
                       `P-Value` = c(MantelOverall$signif, ProtestOverall$signif))
      
    })
  
  Result3 <-
    map_dfr(.x = unique(DataVeg$SiteType), ~{
      
      Veg <- vegdist(DataVeg %>% arrange(Independent)%>% ungroup() %>% filter(SiteType == .x)%>%
                       select(contains(ColTraitVeg)))
      Bird <- vegdist(DataBird %>% arrange(Independent)%>% ungroup() %>% filter(SiteType == .x)%>%
                        select(ANCH:ZEND))
      
      MantelOverall <- mantel(vegdist(wisconsin(Veg)),vegdist(wisconsin(Veg)),method="spearman")
      
      ProtestOverall <- protest(decostand(Veg,method = "hellinger"),
                                decostand(Bird,method = "hellinger"))
      Result <- tibble(Group = .x,
                       Test = c("Mantel", "Procustes"),
                       Strength = c(MantelOverall$statistic, ProtestOverall$t0),
                       `P-Value` = c(MantelOverall$signif, ProtestOverall$signif))
      
    })
  
  Result <- Result %>% bind_rows(Result2) %>% bind_rows(Result3)
  
  return(Result)
}


#######################################################################################################


######################################################################################################
######################################################################################################
############################# Biological Interactions  ###############################################
######################################################################################################
######################################################################################################

DispersionPlots <- function(DataX, DataType, groups){
  
  if(DataType == "Birds"){
    
    DataX <-
      DataX%>% arrange(Independent)%>% ungroup()
    
    MatriX <- 
      DataX %>%
      select(ANCH:ZEND)
    
  }
  
  else if(DataType == "Traits"){
    DataX <-
      DataX%>% arrange(SiteID)%>% ungroup() 
   
     # if any of the rows sum to 0, i need to drop them
    KeepRows <- DataX %>%
      select(`Aquatic invertebrates`:`Resident`) %>% 
      mutate(SumRows = rowSums(.)) %>% 
      bind_cols(DataX %>% select(SiteID)) %>% filter(SumRows>0) %>% select(SiteID) %>% reduce(c)
    
    DataX <-
      DataX%>% filter(SiteID%in%KeepRows) 
    
    MatriX <- 
      DataX %>%
      select(`Aquatic invertebrates`:`Resident`)
      
    
  }
  else{
    DataX <-
      DataX %>%
      arrange(SiteID)%>% ungroup() 
    
    if(DataType == "BasalArea")
    {
      MatriX <- 
        decostand(DataX %>%
        select(contains(c("BasalArea"))), "max", 2)
    }
    else if(DataType == "CanopyWidth"){
      
      MatriX <- 
        decostand(DataX %>%
                    select(contains(c("CanopyWidth"))), "max", 2)
    }
    else{
      MatriX <- 
        decostand(DataX %>%
                    select(contains(c("TreeHeight"))), "max", 2)
      
    }
  }
 
  Plot <-
    map(1:length(groups), ~{
     
    
      
    Group <- DataX %>% select(!!sym(groups[.x])) %>% ungroup() %>% reduce(c)
    
    if(any(DataType%in% c("CanopyWidth", "TreeHeight", "BasalArea"))){
      
      dispersion <- betadisper(vegdist(decostand(MatriX, method = "hellinger")), Group)
      
    }
  
  else{
    dispersion <- betadisper(vegdist(decostand(MatriX, "max", 2)), Group)
  }
    
    OrdiPlot <-
      DataX %>% select(!!sym(groups[.x]))%>%
      rename(Grouping = groups[.x])%>%
      bind_cols(dispersion$vectors%>% as_data_frame())
    
    Centriods <- dispersion$centroids %>% as_data_frame()%>% 
      mutate(Grouping = unique(Group))  %>% 
      rename(PCoA1C = PCoA1) %>% rename(PCoA2C = PCoA2)%>% select(Grouping, PCoA1C, PCoA2C)
    
    OrdiPlot <-
      OrdiPlot %>% left_join(Centriods, by = "Grouping")
    
    OrdiPlot %>%
      ggplot()+
      geom_mark_hull(aes(x=PCoA1, y=PCoA2, 
                         colour = stage(Grouping, after_scale = alpha(color, 0.75)),
                         fill = stage(Grouping, after_scale = alpha(fill, 0.2))),
                     expand = unit(0, "mm"),
                     radius = unit(1, "mm"),
                     concavity = 10, show.legend=FALSE) +
      geom_segment(aes(x=PCoA1C, y=PCoA2C, xend=PCoA1, yend=PCoA2, colour=Grouping), 
                   show.legend=FALSE) +
      geom_point(aes(x=PCoA1, y=PCoA2, colour=Grouping))+
      scale_fill_viridis_d(option = "magma")+
      scale_color_viridis_d(option = "magma")+
      theme_bw()+
      theme(
        # X-axis title: no angle, all the way to the right
        axis.title.x = element_text(angle = 0, hjust = 1, vjust = 0),
        
        # Y-axis title: normal vertical orientation, all the way to the left
        axis.title.y = element_text(angle = 90, hjust = 1, vjust = 0),
        
        # Adjust plot margins to accommodate titles
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        legend.position = "bottom",
      )+
      geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
      geom_hline(yintercept = c(0), color = "grey70", linetype = 2) + 
      labs(color = " ", fill = " ")+
      labs(x = paste0("PCoA 1 (", dispersion$eig[1]%>% round(2), "%)"), 
           y = paste0("PCoA 2 (", dispersion$eig[2] %>% round(2), "%)"))
    
    
  })
  
  png(paste0(here::here(), "/figures/", paste(c(DataType, groups), collapse = ""), ".png"),
      width = 6, height = 6, res = 600, units = "in", bg = "white")
  gridExtra::grid.arrange(grobs = Plot)
  
  dev.off()
  Plots <- gridExtra::grid.arrange(grobs = Plot)
  
  
  
  return(Plots)
}




DispersionResults <- function(DataX, DataType, groups){
  
  if(DataType == "Birds"){
    
    DataX <-
      DataX%>% arrange(Independent)%>% ungroup()
    
    MatriX <- 
      DataX %>%
      select(ANCH:ZEND)
    
    
  }
  else if(DataType == "Traits"){
    DataX <-
      DataX%>% arrange(SiteID)%>% ungroup() 
    
    # if any of the rows sum to 0, i need to drop them
    KeepRows <- DataX %>%
      select(`Aquatic invertebrates`:`Resident`) %>% 
      mutate(SumRows = rowSums(.)) %>% 
      bind_cols(DataX %>% select(SiteID)) %>% filter(SumRows>0) %>% select(SiteID) %>% reduce(c)
    
    DataX <-
      DataX%>% filter(SiteID%in%KeepRows) 
    
    MatriX <- 
      DataX %>%
      select(`Aquatic invertebrates`:`Resident`)
    
  }
  else{
    DataX <-
      DataX %>%
      arrange(SiteID)%>% ungroup() 
    
    if(DataType == "BasalArea")
    {
      MatriX <- 
        decostand(DataX %>%
                    select(contains(c("BasalArea"))), "max", 2)
    }
    else if(DataType == "CanopyWidth"){
      
      MatriX <- 
        decostand(DataX %>%
                    select(contains(c("CanopyWidth"))), "max", 2)
    }
    else{
      MatriX <- 
        decostand(DataX %>%
                    select(contains(c("TreeHeight"))), "max", 2)
      
    }
  }
  
  Results <-
    map_dfr(1:length(groups), ~{
      Group <- DataX %>% select(!!sym(groups[.x])) %>% ungroup() %>% reduce(c)
      
      if(any(DataType%in% c("CanopyWidth", "TreeHeight", "BasalArea"))){
      
      dispersion <- betadisper(vegdist(wisconsin(MatriX)), Group)
      
    }
      
      else{
        dispersion <- betadisper(vegdist(decostand(MatriX, "max", 2)), Group)
      }
      
      result <- anova(dispersion)
      
      save <- tibble(Group = groups[.x],
                     Df = result$Df[1], `F value` = result$`F value`[1], `Pr(>F)` = result$`Pr(>F)`[1])
    })
  
  Result2 <-
    map_dfr(1:length(groups), ~{
      Group <- DataX %>% select(!!sym(groups[.x])) %>% ungroup() %>% reduce(c)
      
      if(any(DataType%in% c("CanopyWidth", "TreeHeight", "BasalArea"))){
      
      dispersion <- betadisper(vegdist(decostand(MatriX, method = "hellinger")), Group)
      
    }
  
  else{
    dispersion <- betadisper(vegdist(decostand(MatriX, "max", 2)), Group)
  }
      
      result <- TukeyHSD(dispersion)
      
      save <- tibble(Group = result$group %>% as.data.frame() %>% row.names(),
                     `P-Value (diff)` = result$group[,"diff"]%>% unname(), 
                     `P-Value (lwr)` = result$group[,"lwr"]%>% unname(), 
                     `P-Value (upr)` = result$group[,"upr"]%>% unname(), 
                     `P-Value (adj)` = result$group[,"p adj"]%>% unname())
    })
  Final <- list(ANOVA = Results,
                TukeyHSD = Result2)
  
  return(Final)
}


######################################################################################################
######################################################################################################
############################# Biological Interactions  ###############################################
######################################################################################################
######################################################################################################


NMDSPlotsBirdsPlants <- function(NMDSList,
                                 DataX,
                                 DataY, 
                                 Groups,
                                 Type){
  
  FileNames <- names(NMDSList)
  Files<- 
    map(1:length(NMDSList), function(x)
    {
    
    DataY <-
      map(Groups, function(y)
        { # Species	Food	Behaviour	Nesting	Habitat	Type
        
        if(y == "Species"){
          DataY <-
            suppressWarnings(suppressMessages(DataY[[1]] %>%
            arrange(Independent)%>% ungroup() %>%
            select(ANCH:ZEND)))
        }
        
        else if(y == "Food"){
          DataY <-
            suppressWarnings(suppressMessages(DataY[[2]] %>%
            arrange(Independent)%>% ungroup() %>%
            select(Independent, 
                   ANCH:ZEND) %>%
            pivot_longer(!Independent, names_to = "Code",
                         values_to = "Abundance")%>%
            left_join(read_csv(paste0(here::here(),
                                      "/data/bird_traits_raw.csv")) %>%
                        select(Code, Food))%>%
            group_by(Independent,
                     Food)%>%
            summarise(Abundance = round(mean(Abundance)))%>%
            filter(!is.na(Food))%>%
            pivot_wider(names_from = "Food",
                        values_from = "Abundance")%>%
            arrange(Independent) %>% ungroup() %>%
            select(-Independent)))
        }
        
        else if(y == "Behaviour"){
          DataY <-
            suppressWarnings(suppressMessages(DataY[[2]] %>%
            arrange(Independent)%>% ungroup() %>%
            select(Independent, 
                   ANCH:ZEND) %>%
            pivot_longer(!Independent, names_to = "Code",
                         values_to = "Abundance")%>%
            left_join(read_csv(paste0(here::here(),
                                      "/data/bird_traits_raw.csv")) %>%
                        select(Code, Behaviour))%>%
            group_by(Independent,
                     Behaviour)%>%
            summarise(Abundance = round(mean(Abundance)))%>%
            filter(!is.na(Behaviour))%>%
            pivot_wider(names_from = "Behaviour",
                        values_from = "Abundance")%>%
            arrange(Independent) %>% ungroup() %>%
            select(-Independent)))
        }
        
        
        else if(y == "Nesting"){
          DataY <-
            suppressWarnings(suppressMessages(DataY[[2]] %>%
            arrange(Independent)%>% ungroup() %>%
            select(Independent, 
                   ANCH:ZEND) %>%
            pivot_longer(!Independent, names_to = "Code",
                         values_to = "Abundance")%>%
            left_join(read_csv(paste0(here::here(),
                                      "/data/bird_traits_raw.csv")) %>%
                        select(Code, Nesting))%>%
            group_by(Independent,
                     Nesting)%>%
            summarise(Abundance = round(mean(Abundance)))%>%
            filter(!is.na(Nesting))%>%
            pivot_wider(names_from = "Nesting",
                        values_from = "Abundance")%>%
            arrange(Independent) %>% ungroup() %>%
            select(-Independent)))
        }
        
        
        else if(y == "Habitat"){
          DataY <-
            suppressWarnings(suppressMessages(DataY[[2]] %>%
            arrange(Independent)%>% ungroup() %>%
            select(Independent, 
                   ANCH:ZEND) %>%
            pivot_longer(!Independent, names_to = "Code",
                         values_to = "Abundance")%>%
            left_join(read_csv(paste0(here::here(),
                                      "/data/bird_traits_raw.csv")) %>%
                        select(Code, Habitat))%>%
            group_by(Independent,
                     Habitat)%>%
            summarise(Abundance = round(mean(Abundance)))%>%
            filter(!is.na(Habitat))%>%
            pivot_wider(names_from = "Habitat",
                        values_from = "Abundance")%>%
            arrange(Independent) %>% ungroup() %>%
            select(-Independent)))
        }
        
        else if(y == "Type"){
          DataY <-
            suppressWarnings(suppressMessages(DataY[[2]] %>%
            arrange(Independent)%>% ungroup() %>%
            select(Independent, 
                   ANCH:ZEND) %>%
            pivot_longer(!Independent, names_to = "Code",
                         values_to = "Abundance")%>%
            left_join(read_csv(paste0(here::here(),
                                      "/data/bird_traits_raw.csv")) %>%
                        select(Code, Type))%>%
            group_by(Independent,
                     Type)%>%
            summarise(Abundance = round(mean(Abundance)))%>%
            filter(!is.na(Type))%>%
            pivot_wider(names_from = "Type",
                        values_from = "Abundance")%>%
            arrange(Independent) %>% ungroup() %>%
            select(-Independent)))
        }
        
        NMS.Vectors <- data.frame(scores(envfit(NMDSList[[x]], # extract species scores and fit data
                                                DataY, 
                                                permutations = 999), # no more than these permutations
                                         "vectors")) # exract vectors and not scores
        # need to make a column that will these vectors to the orignal data
        if(y != "Species"){
          NMS.Vectors <- 
            NMS.Vectors %>%
            mutate(Label = row.names(NMS.Vectors))}
        
        else{
          NMS.Vectors <- 
            suppressWarnings(suppressMessages(NMS.Vectors %>%
            mutate(Code = row.names(NMS.Vectors))%>%
            left_join(read_csv(paste0(here::here(), "/data/bird_traits_raw.csv")))%>%
            rename(Label = Species)))
        }
        
        # save the vector correlations
        NMS.Correlation <- data.frame(envfit(NMDSList[[x]],
                                             DataY, 
                                             permutations = 999)$vectors[2])
        
        NMS.Correlation <-
          NMS.Correlation %>%
          mutate(Code = row.names(NMS.Correlation))
        
        if(y != "Species"){
          AxisCorrelation <-
            suppressWarnings(suppressMessages(NMS.Vectors %>% left_join(NMS.Vectors)))
          
          NMS.Scores <- 
            DataX %>% arrange(Independent)%>% 
            ungroup() %>% select(SiteName, SiteType, Independent)%>%
            bind_cols(scores(NMDSList[[x]])$sites)
        }
        
        else{
          AxisCorrelation <- suppressWarnings(suppressMessages(NMS.Correlation %>% 
            left_join(NMS.Vectors)  %>% filter(r>0.2) ))
          NMS.Scores <- 
            DataX %>% arrange(Independent)%>% 
            ungroup() %>% select(SiteName, SiteType, Independent)%>%
            bind_cols(scores(NMDSList[[x]])$sites)
        }
        if(Type == "Site"){
          ggplot()+
            geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                       color = SiteName),size = 3)+
            geom_text_repel(data = AxisCorrelation,
                            aes(x = NMDS1, y = NMDS2,label= Label),
                            box.padding = unit(0.75, "lines"),
                            point.padding = unit(0.15, "lines"))+                                      
            geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                         arrow = arrow(length = unit(0.35, "cm"),
                                       type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
            scale_color_brewer(palette = "BrBG")+
            scale_fill_brewer(palette = "BrBG")+
            theme_gaea()+
            labs(color = " ", fill = " ", 
                 caption = paste0(FileNames[x],
                                  " as Proxy for Community Structure"))+
            theme(legend.position = "bottom")
          
        }
        else{
          ggplot()+
            geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                       color = SiteType),size = 3)+
            geom_text_repel(data = AxisCorrelation,
                            aes(x = NMDS1, y = NMDS2,label= Label),
                            box.padding = unit(0.75, "lines"),
                            point.padding = unit(0.15, "lines"))+                                      
            geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                         arrow = arrow(length = unit(0.35, "cm"),
                                       type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
            scale_color_brewer(palette = "BrBG")+
            scale_fill_brewer(palette = "BrBG")+
            theme_gaea()+
            labs(color = " ", fill = " ", 
                 caption = paste0(FileNames[x],
                                  " as Proxy for Community Structure"))+
            theme(legend.position = "bottom")
        }
        }
        )
  
    
  }
  )
  
  files_save <- 
    map(1:length(NMDSList), function(x){
      map(1:length(Groups), function(y){
        paste0(here::here(), "/figures/PlantsBirds", Type, str_remove_all(FileNames[x], " "),
               Groups[y], "NMDS.png")
      })})
  
  
  map(1:length(NMDSList), function(x){
    map(1:length(Groups), function(y){
      png(paste0(here::here(), "/figures/PlantsBirds",Type, str_remove_all(FileNames[x], " "),
                 Groups[y], "NMDS.png"), width = 7, height = 5,units = "in", res = 600, bg = "white")
      plot(Files[[x]][[y]])
      dev.off()
    })})
    
  return(Files)

}




NMDSPlotsBirds <- function(NMDSList,
                           FileNames,
                           DataX,
                           Choices,
                           Groups,
                           Type)
  {
  
  Files <-
        map(Groups, function(y)
        { # Species	Food	Behaviour	Nesting	Habitat	Type
          
          if(y == "Species"){
            DataY <-
              suppressWarnings(suppressMessages(DataX %>%ungroup()%>%
                                                  select(ANCH:ZEND)))
          }
          
          else if(y == "Food"){
            DataY <-
              suppressWarnings(suppressMessages(DataX%>% ungroup()%>%
                                                  select(Independent, 
                                                         ANCH:ZEND) %>%
                                                  pivot_longer(!Independent, names_to = "Code",
                                                               values_to = "Abundance")%>%
                                                  left_join(read_csv(paste0(here::here(),
                                                                            "/data/bird_traits_raw.csv")) %>%
                                                              select(Code, Food))%>%
                                                  group_by(Independent,
                                                           Food)%>%
                                                  summarise(Abundance = round(mean(Abundance)))%>%
                                                  filter(!is.na(Food))%>%
                                                  pivot_wider(names_from = "Food",
                                                              values_from = "Abundance")%>%
                                                  arrange(Independent) %>% ungroup() %>%
                                                  select(-Independent)))
          }
          
          else if(y == "Behaviour"){
            DataY <-
              suppressWarnings(suppressMessages(DataX %>%ungroup()%>% 
                                                  select(Independent, 
                                                         ANCH:ZEND) %>%
                                                  pivot_longer(!Independent, names_to = "Code",
                                                               values_to = "Abundance")%>%
                                                  left_join(read_csv(paste0(here::here(),
                                                                            "/data/bird_traits_raw.csv")) %>%
                                                              select(Code, Behaviour))%>%
                                                  group_by(Independent,
                                                           Behaviour)%>%
                                                  summarise(Abundance = round(mean(Abundance)))%>%
                                                  filter(!is.na(Behaviour))%>%
                                                  pivot_wider(names_from = "Behaviour",
                                                              values_from = "Abundance")%>%
                                                  arrange(Independent) %>% ungroup() %>%
                                                  select(-Independent)))
          }
          
          
          else if(y == "Nesting"){
            DataY <-
              suppressWarnings(suppressMessages(DataX %>% ungroup() %>%
                                                  select(Independent, 
                                                         ANCH:ZEND) %>%
                                                  pivot_longer(!Independent, names_to = "Code",
                                                               values_to = "Abundance")%>%
                                                  left_join(read_csv(paste0(here::here(),
                                                                            "/data/bird_traits_raw.csv")) %>%
                                                              select(Code, Nesting))%>%
                                                  group_by(Independent,
                                                           Nesting)%>%
                                                  summarise(Abundance = round(mean(Abundance)))%>%
                                                  filter(!is.na(Nesting))%>%
                                                  pivot_wider(names_from = "Nesting",
                                                              values_from = "Abundance")%>%
                                                  arrange(Independent) %>% ungroup() %>%
                                                  select(-Independent)))
          }
          
          
          else if(y == "Habitat"){
            DataY <-
              suppressWarnings(suppressMessages(DataX %>% ungroup() %>%
                                                  select(Independent, 
                                                         ANCH:ZEND) %>%
                                                  pivot_longer(!Independent, names_to = "Code",
                                                               values_to = "Abundance")%>%
                                                  left_join(read_csv(paste0(here::here(),
                                                                            "/data/bird_traits_raw.csv")) %>%
                                                              select(Code, Habitat))%>%
                                                  group_by(Independent,
                                                           Habitat)%>%
                                                  summarise(Abundance = round(mean(Abundance)))%>%
                                                  filter(!is.na(Habitat))%>%
                                                  pivot_wider(names_from = "Habitat",
                                                              values_from = "Abundance")%>%
                                                  arrange(Independent) %>% ungroup() %>%
                                                  select(-Independent)))
          }
          
          else if(y == "Type"){
            DataY <-
              suppressWarnings(suppressMessages(DataX %>% ungroup() %>%
                                                  select(Independent, 
                                                         ANCH:ZEND) %>%
                                                  pivot_longer(!Independent, names_to = "Code",
                                                               values_to = "Abundance")%>%
                                                  left_join(read_csv(paste0(here::here(),
                                                                            "/data/bird_traits_raw.csv")) %>%
                                                              select(Code, Type))%>%
                                                  group_by(Independent,
                                                           Type)%>%
                                                  summarise(Abundance = round(mean(Abundance)))%>%
                                                  filter(!is.na(Type))%>%
                                                  pivot_wider(names_from = "Type",
                                                              values_from = "Abundance")%>%
                                                  arrange(Independent) %>% ungroup() %>%
                                                  select(-Independent)))
          }
          
          NMS.Vectors <- data.frame(scores(envfit(NMDSList, # extract species scores and fit data
                                                  DataY, 
                                                  choices= Choices,
                                                  permutations = 999), # no more than these permutations
                                           "vectors")) # exract vectors and not scores
          # need to make a column that will these vectors to the orignal data
          if(y != "Species"){
            NMS.Vectors <- 
              NMS.Vectors %>%
              mutate(Label = row.names(NMS.Vectors))}
          
          else{
            NMS.Vectors <- 
              suppressWarnings(suppressMessages(NMS.Vectors %>%
                                                  mutate(Code = row.names(NMS.Vectors))%>%
                                                  left_join(read_csv(paste0(here::here(), "/data/bird_traits_raw.csv")))%>%
                                                  rename(Label = Species)))
          }
          
          # save the vector correlations
          NMS.Correlation <- data.frame(envfit(NMDSList,
                                               DataY, 
                                               permutations = 999)$vectors[2])
          
          NMS.Correlation <-
            NMS.Correlation %>%
            mutate(Code = row.names(NMS.Correlation))
          
          if(y != "Species"){
            AxisCorrelation <-
              suppressWarnings(suppressMessages(NMS.Vectors %>% left_join(NMS.Vectors)))
            
            NMS.Scores <- 
              DataX %>% 
              ungroup() %>% select(SiteName, SiteType, Independent)%>%
              bind_cols(scores(NMDSList)$sites)
          }
          
          else{
            AxisCorrelation <- suppressWarnings(suppressMessages(NMS.Correlation %>% 
                                                                   left_join(NMS.Vectors)  %>% filter(r>0.2) ))
            NMS.Scores <- 
              DataX %>% 
              ungroup() %>% select(SiteName, SiteType, Independent)%>%
              bind_cols(scores(NMDSList)$sites)
          }
          if(Type == "Site"){
           
            if(any(str_detect(colnames(NMS.Scores), "NMDS3"))){
              
              Plot1 <- ggplot()+
                geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                           color = SiteName, fill = SiteName),size = 3)+
                stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2,
                             data=NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                                 fill = SiteName))+
                geom_text_repel(data = AxisCorrelation,
                                aes(x = NMDS1, y = NMDS2,label= Label),
                                box.padding = unit(0.75, "lines"),
                                point.padding = unit(0.15, "lines"))+      
                geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                             arrow = arrow(length = unit(0.35, "cm"),
                                           type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
                scale_color_brewer(palette = "BrBG")+
                scale_fill_brewer(palette = "BrBG")+
                theme_gaea()+
                labs(color = " ", fill = " ")+
                theme(legend.position = "none")
              
              Plot2 <- 
                ggplot()+
                geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS3, 
                                                           color = SiteName, fill = SiteName),size = 3)+
                stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2,
                             data=NMS.Scores,aes(x = NMDS1, y = NMDS3,
                                                 fill = SiteName))+
                geom_text_repel(data = AxisCorrelation,
                                aes(x = NMDS1, y = NMDS3,label= Label),
                                box.padding = unit(0.75, "lines"),
                                point.padding = unit(0.15, "lines"))+     
                geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS3),
                             arrow = arrow(length = unit(0.35, "cm"),
                                           type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
                scale_color_brewer(palette = "BrBG")+
                scale_fill_brewer(palette = "BrBG")+
                theme_gaea()+
                labs(color = " ", fill = " ")+
                theme(legend.position = "none")
              Plot3 <- 
                lemon::g_legend(ggplot()+
                                  geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS3,    color = SiteName, fill = SiteName),size = 3)+
                                  stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2, data=NMS.Scores,aes(x = NMDS1, y = NMDS3,
                                                                                             fill = SiteName))+
                                  geom_text_repel(data = AxisCorrelation,  aes(x = NMDS1, y = NMDS3,label= Label),
                                                                            box.padding = unit(0.75, "lines"),
                                                                            point.padding = unit(0.15, "lines"))+ 
                                  geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS3),
                                               arrow = arrow(length = unit(0.35, "cm"),
                                                             type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
                                  scale_color_brewer(palette = "BrBG")+
                                  scale_fill_brewer(palette = "BrBG")+
                                  theme_gaea()+
                                  labs(color = " ", fill = " ")+
                                  theme(legend.position = "bottom"))
              Plots <- list(Plot1, Plot2, Plot3)
              
            }
            else{
              ggplot()+
                geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                           color = SiteName, fill = SiteName),size = 3)+
                                            stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2,
                                                         data=NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                                                             fill = SiteName))+
                                            geom_text_repel(data = AxisCorrelation,
                                                            aes(x = NMDS1, y = NMDS2,label= Label),
                                                            box.padding = unit(0.75, "lines"),
                                                            point.padding = unit(0.15, "lines"))+                                      
                                            geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                                                         arrow = arrow(length = unit(0.35, "cm"),
                                                                       type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
                scale_color_brewer(palette = "BrBG")+
                scale_fill_brewer(palette = "BrBG")+
                theme_gaea()+
                                            labs(color = " ", fill = " ")+
                                            theme(legend.position = "bottom")
            }
                                          
          
                                         
            }
            
            
          else{
            if(any(str_detect(colnames(NMS.Scores), "NMDS3"))){
              
            Plot1 <- 
              ggplot()+
              geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                         color = SiteType, fill = SiteType),size = 3)+
              stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2,
                           data=NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                               fill = SiteType))+
              geom_text_repel(data = AxisCorrelation,
                              aes(x = NMDS1, y = NMDS2,label= Label),
                              box.padding = unit(0.75, "lines"),
                              point.padding = unit(0.15, "lines"))+                                      
              geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                           arrow = arrow(length = unit(0.35, "cm"),
                                         type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
              scale_color_brewer(palette = "BrBG")+
              scale_fill_brewer(palette = "BrBG")+
              theme_gaea()+
              labs(color = " ", fill = " ")+
              theme(legend.position = "none")
            Plot2 <- 
              ggplot()+
              geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS3,
                                                         color = SiteType, fill = SiteType),size = 3)+
              stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2,
                           data=NMS.Scores,aes(x = NMDS1, y = NMDS3,
                                               fill = SiteType))+
              geom_text_repel(data = AxisCorrelation,
                              aes(x = NMDS1, y = NMDS3,label= Label),
                              box.padding = unit(0.75, "lines"),
                              point.padding = unit(0.15, "lines"))+                                      
              geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS3),
                           arrow = arrow(length = unit(0.35, "cm"),
                                         type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
              scale_color_brewer(palette = "BrBG")+
              scale_fill_brewer(palette = "BrBG")+
              theme_gaea()+
              labs(color = " ", fill = " ")+
              theme(legend.position = "none")
            Plot3 <- 
              lemon::g_legend(ggplot()+
                       geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS3,
                                                                  color = SiteType, fill = SiteType),size = 3)+
                       stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2,
                                    data=NMS.Scores,aes(x = NMDS1, y = NMDS3,
                                                        fill = SiteType))+
                       geom_text_repel(data = AxisCorrelation,
                                       aes(x = NMDS1, y = NMDS3,label= Label),
                                       box.padding = unit(0.75, "lines"),
                                       point.padding = unit(0.15, "lines"))+                                      
                       geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS3),
                                    arrow = arrow(length = unit(0.35, "cm"),
                                                  type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
                         scale_color_brewer(palette = "BrBG")+
                         scale_fill_brewer(palette = "BrBG")+
                         theme_gaea()+
                       labs(color = " ", fill = " ")+
                       theme(legend.position = "bottom"))
            
            Plots <- list(Plot1, Plot2, Plot3)
            }
            
            else{
              ggplot()+
                geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                           color = SiteType, fill = SiteType),size = 3)+
                stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/2,
                             data=NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                                 fill = SiteType))+
                geom_text_repel(data = AxisCorrelation,
                                aes(x = NMDS1, y = NMDS2,label= Label),
                                box.padding = unit(0.75, "lines"),
                                point.padding = unit(0.15, "lines"))+                                      
                geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                             arrow = arrow(length = unit(0.35, "cm"),
                                           type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
                scale_color_brewer(palette = "BrBG")+
                scale_fill_brewer(palette = "BrBG")+
                theme_gaea()+
                labs(color = " ", fill = " ")+
                theme(legend.position = "bottom")
            }
          }
        
        }
        )
  
  files_save <- 
    map(1:length(Groups), function(y){
        paste0(here::here(), "/figures/Birds", Type, str_remove_all(FileNames, " "),
               Groups[y], "NMDS.png")
      })
  
  
  map(1:length(Groups), function(y){
    
    if(length(Files[[y]])== 1){
      
      png(paste0(here::here(), "/figures/Birds",Type, str_remove_all(FileNames, " "),
                 Groups[y], "NMDS.png"), width = 7, height = 5,units = "in", res = 600, bg = "white")
      cowplot::ggdraw(Files[[y]])+
        theme(plot.background = element_rect(fill="white", color = NA))
      dev.off()
    }
    else{
      
      png(paste0(here::here(), "/figures/Birds",Type, str_remove_all(FileNames, " "),
                 Groups[y], "NMDS.png"), width = 7, height = 5,units = "in", res = 600, bg = "white")
     cowplot::ggdraw(gridExtra::grid.arrange(grobs = list(Files[[y]][[1]],
                                                Files[[y]][[2]],
                                                Files[[y]][[3]]),
                                   layout_matrix = rbind(c(1,1,1,1,2,2,2,2),
                                                         c(1,1,1,1,2,2,2,2),
                                                         c(1,1,1,1,2,2,2,2),
                                                         c(1,1,1,1,2,2,2,2),
                                                         c(3))))+
       theme(plot.background = element_rect(fill="white", color = NA))
                        
      dev.off()
      
    }
      
    })
  
  Files <-  map(1:length(Groups), function(y){
    
    if(length(Files[[y]])== 1){
     cowplot::ggdraw(Files[[y]])+
        theme(plot.background = element_rect(fill="white", color = NA))
    }
    else{
      cowplot::ggdraw(gridExtra::grid.arrange(grobs = list(Files[[y]][[1]],
                                                Files[[y]][[2]],
                                                Files[[y]][[3]]),
                                   layout_matrix = rbind(c(1,1,1,1,2,2,2,2),
                                                         c(1,1,1,1,2,2,2,2),
                                                         c(1,1,1,1,2,2,2,2),
                                                         c(1,1,1,1,2,2,2,2),
                                                         c(3)))) +
                        theme(plot.background = element_rect(fill="white", color = NA))
      
    }
    })
    
     
  
  return(Files)
  
}




NMDSPlotsPlants <- function(NMDSList,
                       DataX, Type){
  
  FileNames <- names(NMDSList)
  Files<- 
    map(1:length(NMDSList), function(x)
    {
      
      DataY <-
        suppressWarnings(suppressMessages(DataX %>%
            arrange(SiteID)%>% ungroup() %>%
            select(contains(c("TreeHeight", "BasalArea", "CanopyWidth")))))

          
          NMS.Vectors <- data.frame(scores(envfit(NMDSList[[x]], # extract species scores and fit data
                                                  DataY, 
                                                  permutations = 999), # no more than these permutations
                                           "vectors")) # exract vectors and not scores
          # need to make a column that will these vectors to the orignal data
          NMS.Vectors <- 
              NMS.Vectors %>%
              mutate(Code = row.names(NMS.Vectors),
                     Label = case_when(str_detect(Code, "TreeHeight") ~ "Tree Height",
                                       str_detect(Code, "BasalArea") ~ "Basal Area",
                                       str_detect(Code, "CanopyWidth") ~ "Canopy Width"),
                     Label = case_when(str_detect(Code, "Black") ~ paste0(Label, "(", "Black", ")"),
                                       str_detect(Code, "Red") ~  paste0(Label, "(", "Red", ")"),
                                       str_detect(Code, "White") ~  paste0(Label, "(", "White", ")")))
          
          # save the vector correlations
          NMS.Correlation <- data.frame(envfit(NMDSList[[x]],
                                               DataY, 
                                               permutations = 999)$vectors[2])
          
          NMS.Correlation <-
            NMS.Correlation %>%
            mutate(Code = row.names(NMS.Correlation))
          
            NMS.Scores <- 
              DataX %>% arrange(SiteID)%>% 
              ungroup() %>% select(SiteName, SiteType, SiteID)%>%
              bind_cols(scores(NMDSList[[x]])$sites)
          
            AxisCorrelation <- suppressWarnings(suppressMessages(NMS.Correlation %>% left_join(NMS.Vectors)
                                                %>% filter(r>0.2)) )
          
            
          if(Type == "Site"){
            
            
            ggplot()+
              geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                         color = SiteName),size = 3)+
              stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/8,
                           data=NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                               fill = SiteName))+
              stat_ellipse(level = 0.90,data = NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                                              color = SiteName))+
              geom_text_repel(data = AxisCorrelation,
                              aes(x = NMDS1, y = NMDS2,label= Label),
                              box.padding = unit(0.75, "lines"),
                              point.padding = unit(0.15, "lines"))+                                      
              geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                           arrow = arrow(length = unit(0.35, "cm"),
                                         type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
              scale_color_brewer(palette = "BrBG")+
              scale_fill_brewer(palette = "BrBG")+
              theme_gaea()+
              labs(color = " ", fill = " ", 
                   caption = paste0(FileNames[x],
                                    " as Proxy for Community Structure"))+
              theme(legend.position = "bottom")
            
          }
            else{
              
              ggplot()+
                geom_point(data=data.frame(NMS.Scores),aes(x = NMDS1, y = NMDS2,
                                                           color = SiteType),size = 3)+
                stat_ellipse(level = 0.90,geom = "polygon", alpha = 1/8,
                             data=NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                                 fill = SiteType))+
                stat_ellipse(level = 0.90,data = NMS.Scores,aes(x = NMDS1, y = NMDS2,
                                                                color = SiteType))+
                geom_text_repel(data = AxisCorrelation,
                                aes(x = NMDS1, y = NMDS2,label= Label),
                                box.padding = unit(0.75, "lines"),
                                point.padding = unit(0.15, "lines"))+                                      
                geom_segment(data= AxisCorrelation,aes(x=0,xend = NMDS1, y=0, yend = NMDS2),
                             arrow = arrow(length = unit(0.35, "cm"),
                                           type="closed"),size=1,color = "grey15",inherit.aes=TRUE)+
                scale_color_brewer(palette = "BrBG")+
                scale_fill_brewer(palette = "BrBG")+
                theme_gaea()+
                labs(color = " ", fill = " ", 
                     caption = paste0(FileNames[x],
                                      " as Proxy for Community Structure"))+
                theme(legend.position = "bottom")}
          
          }
    )
  
  files_save <- 
    map(1:length(NMDSList), function(x){
        paste0(here::here(), "/figures/Plants",Type, str_remove_all(FileNames[x], " "), "NMDS.png")})
  
  
  map(1:length(NMDSList), function(x){
    png(paste0(here::here(), "/figures/Plants", Type, str_remove_all(FileNames[x], " "), "NMDS.png"), 
        width = 7, height = 5,units = "in", res = 600, bg = "white")
      plot(Files[[x]])
      dev.off()
      })
  
  return(Files)
  
}

######################################################################################################
########################################### Theme: Gaea #############################################
###################################################################################################### 


theme_gaea <-  function(){ 
  font <- "Helvetica"   
  
  
  theme_bw() %+replace%   
    theme(text = element_text(family = font),
          # X-axis title: no angle, all the way to the right
          axis.title.x = element_text(angle = 0, hjust = 1, vjust = 1),
          
          # Y-axis title: normal vertical orientation, all the way to the left
          axis.title.y = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.position = "bottom",
          panel.border = element_rect(fill = NA, colour = "grey60"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.10, colour = "grey70"),
          plot.subtitle = element_text(hjust = 1, face = "bold"))
  }

######################################################################################################
################################ Models - Plants & Birds #############################################
###################################################################################################### 

DiversityAlpha <- function(data,
                           SiteID,
                           ColStart,
                           ColEnd){
  library(vegetarian)
  library(tidyverse)
  A <- data %>% ungroup()%>% select(!!sym(ColStart):!!sym(ColEnd))
  
  Diversity <- 
    data %>% ungroup()%>% 
    select(!!sym(SiteID))%>%
    mutate(Abundance = rowSums(A))%>%
    bind_cols( map_dfc(0:2, function(x){
      
      map(1:nrow(A), function(y){
        
        H(A[y,], lev = "alpha", q = x)
      })%>% reduce(c)
      
    }))
 
colnames(Diversity) <- c(SiteID, "Species Abundance", "Species Richness", "Shannon-Weiner", "Simpson")
  return(Diversity)
}
