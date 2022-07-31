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
  Bird <- vegdist(DataBird %>% arrange(Independent)%>% ungroup() %>%select(ANCH:SHCO))
  
  MantelOverall <- mantel(vegdist(decostand(Veg,  "max", 2)),
                          vegdist(decostand(Bird,  "max", 2)),method="spearman")
  
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
                        select(ANCH:SHCO))
      
      MantelOverall <- mantel(vegdist(Veg),vegdist(Bird),method="spearman")
      
      ProtestOverall <- protest(decostand(Veg,method = "hellinger"),
                                decostand(Bird,method = "hellinger"))
      Result <- tibble(Group = .x,
                       Test = c("Mantel", "Procustes"),
                       Strength = c(MantelOverall$statistic, ProtestOverall$t0),
                       `P-Value` = c(MantelOverall$signif, ProtestOverall$signif))
      
    })
  
  Result <- Result %>% bind_rows(Result2)
  
  return(Result)
}



######################################################################################################
######################################################################################################
############################# Biological Interactions  ###############################################
######################################################################################################
######################################################################################################


NMDSPlotsBirdsPlants <- function(NMDSList,
                                 DataX,
                                 DataY, 
                                 Groups){
  
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
            select(ANCH:SHCO)))
        }
        
        else if(y == "Food"){
          DataY <-
            suppressWarnings(suppressMessages(DataY[[2]] %>%
            arrange(Independent)%>% ungroup() %>%
            select(Independent, 
                   ANCH:SHCO) %>%
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
                   ANCH:SHCO) %>%
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
                   ANCH:SHCO) %>%
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
                   ANCH:SHCO) %>%
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
                   ANCH:SHCO) %>%
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
            ungroup() %>% select(SiteName, Independent)%>%
            bind_cols(scores(NMDSList[[x]])$sites)
        }
        
        else{
          AxisCorrelation <- suppressWarnings(suppressMessages(NMS.Correlation %>% 
            left_join(NMS.Vectors)  %>% filter(r>0.2) ))
          NMS.Scores <- 
            DataX %>% arrange(Independent)%>% 
            ungroup() %>% select(SiteName, Independent)%>%
            bind_cols(scores(NMDSList[[x]])$sites)
        }
        
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
          theme_classic()+
          scale_colour_viridis(discrete = TRUE) +
          scale_fill_viridis(discrete = TRUE) +
          labs(color = " ", fill = " ", 
               caption = paste0(FileNames[x],
                                " as Proxy for Community Structure"))+
          theme(legend.position = "bottom")
        }
        )
  
    
  }
  )
  
  files_save <- 
    map(1:length(NMDSList), function(x){
      map(1:length(Groups), function(y){
        paste0(here::here(), "/outputs/PlantsBirds", str_remove_all(FileNames[x], " "),
               Groups[y], "NMDS.png")
      })})
  
  
  map(1:length(NMDSList), function(x){
    map(1:length(Groups), function(y){
      png(paste0(here::here(), "/outputs/PlantsBirds", str_remove_all(FileNames[x], " "),
                 Groups[y], "NMDS.png"), width = 7, height = 5,units = "in", res = 600)
      plot(Files[[x]][[y]])
      dev.off()
    })})
    
  return(Files)

}




NMDSPlotsPlants <- function(NMDSList,
                       DataX){
  
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
              ungroup() %>% select(SiteName, SiteID)%>%
              bind_cols(scores(NMDSList[[x]])$sites)
          
            AxisCorrelation <- suppressWarnings(suppressMessages(NMS.Correlation %>% left_join(NMS.Vectors)
                                                %>% filter(r>0.2)) )
            
          
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
            theme_classic()+
            scale_colour_viridis(discrete = TRUE) +
            scale_fill_viridis(discrete = TRUE) +
            labs(color = " ", fill = " ", 
                 caption = paste0(FileNames[x],
                                  " as Proxy for Community Structure"))+
            theme(legend.position = "bottom")
          }
    )
  
  files_save <- 
    map(1:length(NMDSList), function(x){
        paste0(here::here(), "/outputs/Plants", str_remove_all(FileNames[x], " "), "NMDS.png")})
  
  
  map(1:length(NMDSList), function(x){
    png(paste0(here::here(), "/outputs/Plants", str_remove_all(FileNames[x], " "), "NMDS.png"), 
        width = 7, height = 5,units = "in", res = 600)
      plot(Files[[x]])
      dev.off()
      })
  
  return(Files)
  
}
