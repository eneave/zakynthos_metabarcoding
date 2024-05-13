#############################################
## Bubbleplot for Greece data
#############################################

## open foranalysis_greece_motutable_final.RData

## Create long version of data frame

# keep relavant data for the long dataframe
motus4 <- motus3[c(1,19,21:50)]

# Remove Parablennius yatabei
motus4 <- motus4[!motus4$scientific_name_final== "Parablennius yatabei",]

library(tidyverse)

long.df <- motus4 %>%
  gather(key = "id_seq", value = "reads", 
         -id,
         -nonnative,
         -best_identity_final,
         -assignment_method,
         -rank_final,
         -scientific_name_final,
         -total_reads)

# combine metadata and long.df

motus5 <- merge(long.df, meta, by=c("id_seq"), all.x = TRUE) 

# calculate proportional read counts

motus5$PRC <- (motus5$reads/motus5$reads_persample)*100

## calculate the log reads
motus5$log_reads <- log(motus5$reads)


## bubbleplot 
library(ggplot2)

## each sample separated, proportional read counts; native and non-native species
windows()
  ggplot(motus5, 
         aes(x = id_seq, y = scientific_name_final, 
             fill = sample_type, size = PRC)) + 
    scale_y_discrete(limits = rev) +
    geom_point(pch = 21) +
    scale_size_continuous(name = "Proportional \nread counts (%)",
                          range = c(1, 6),
                          limits = c(0.01,100),
                          breaks = c(0.01, 1, 10, 30, 50, 70),
                          labels = c("< 1%","1-10%","10-30%","30-50%","50%-70%", "> 70%")) + 
    facet_grid(nonnative ~ location, scales = "free", space = "free") +
    theme_light() +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical") 
  

## each sample separated, proportional read counts; eDNA and sponges separated
  windows()
  ggplot(motus5, 
         aes(x = id_seq, y = scientific_name_final, 
             fill = nonnative, size = PRC)) + 
    scale_y_discrete(limits = rev) +
    geom_point(pch = 21) +
    scale_size_continuous(name = "Proportional \nread counts (%)",
                          range = c(1, 6),
                          limits = c(0.01,100),
                          breaks = c(0.01, 1, 10, 30, 50, 70),
                          labels = c("< 1%","1-10%","10-30%","30-50%","50%-70%", "> 70%")) + 
    facet_grid(sample_type ~ location, scales = "free", space = "free") +
    theme_light() +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical") 
  
  

## each sample separated, log read counts
  windows()
  ggplot(motus5, 
         aes(x = id_seq, y = scientific_name_final, 
             fill = sample_type, size = log_reads)) + 
    scale_y_discrete(limits = rev) +
    geom_point(pch = 21) +
    scale_size_continuous(name = "Log read counts",
                          range = c(1, 6),
                          limits = c(0,10),
                          breaks = c(0, 2, 4, 6, 8, 10),
                          labels = c("0","2","4","6","8", "10")) + 
    facet_grid(nonnative ~ location, scales = "free", space = "free") +
    theme_light() +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical") 


## combining sample types (sort of)
  windows()
  ggplot(motus5, 
         aes(x = sponge_type, y = scientific_name_final, 
             fill = sample_type, size = PRC)) + 
    scale_y_discrete(limits = rev) +
    geom_point(pch = 21) +
    scale_size_continuous(name = "Proportional \nread counts (%)",
                          range = c(1, 6),
                          limits = c(0.01,100),
                          breaks = c(0.01, 1, 10, 30, 50, 70),
                          labels = c("< 1%","1-10%","10-30%","30-50%","50%-70%", "> 70%")) + 
    facet_grid(nonnative ~ location, scales = "free", space = "free") +
    theme_light() +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical") 


#####
## Repeat but with visual census data added in
#####

motus3.5 <- read.csv("motus_contam_removed_over100reads_edit_UVCadded.csv") 

# keep relavant data for the long dataframe
  motus4.5 <- motus3.5[c(1,19,21:75)]
  
# Remove Parablennius yatabei
motus4.5 <- motus4.5[!motus4.5$scientific_name_final== "Parablennius yatabei",]
  
  
  long.df.2 <- motus4.5 %>%
    gather(key = "id_seq", value = "counts", 
           -id,
           -nonnative,
           -best_identity_final,
           -assignment_method,
           -rank_final,
           -scientific_name_final,
           -total_UVC,
           -total_reads)
  
# combine metadata and long.df
  
  motus6 <- merge(long.df.2, meta, by=c("id_seq"), all.x = TRUE) 
  
# calculate proportional read/visual counts
  
  motus6$prop <- (motus6$counts/motus6$x_persample)*100
  

## each sample separated, proportional read counts; native and non-native species
windows()
  ggplot(motus6, 
         aes(x = id_seq, y = scientific_name_final, 
             fill = sample_type, size = prop)) + 
    scale_y_discrete(limits = rev) +
    geom_point(pch = 21) +
    scale_size_continuous(name = "Relative Proportion (%)",
                          range = c(1, 6),
                          limits = c(0.01,100),
                          breaks = c(0.01, 1, 10, 30, 50, 70),
                          labels = c("< 1%","1-10%","10-30%","30-50%","50%-70%", "> 70%")) + 
    facet_grid(nonnative ~ location, scales = "free", space = "free") +
    theme_light() +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical")  
  
## Presence-absence plot
  
motus6$pa <- ifelse(motus6$counts==0,0,1)

## need to collapse the motus by sample type

windows()
  ggplot(motus6, 
         aes(x = sample_type, y = scientific_name_final, 
             fill = factor(pa))) + 
    scale_fill_manual(values = c("white", "black")) +
    scale_y_discrete(limits = rev) +
    geom_tile(colour="grey") +
    facet_grid(nonnative ~ location, scales = "free", space = "free") +
    labs(x = "", y ="") +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "white"),
          panel.grid.minor = element_line(colour = "white"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical")  

#####
# Prepare data for geom_tile plot to do gradient relative proportions
#####
  
## remove sponge samples
## calculate relative proportions by groups, UVC and eDNA from two locations
  
motus7 <- motus6 %>% filter(sample_type=='UVC' | sample_type=='eDNA') %>%
    group_by(scientific_name_final, sample_type, location) %>%
    summarise(sum_counts=sum(counts),
              sum_pergroup= sum(x_persample),
              .groups = 'drop') %>%
    as.data.frame()  

## Add proportions
motus7$prop <- (motus7$sum_counts/motus7$sum_pergroup)*100 

motus7$prop2 <- ifelse(motus7$prop== 0, NA, motus7$prop)

## Add nonnative
motus7$nis <- ifelse(motus7$scientific_name_final== "Atherinomorus forskalii", "NIS",
                    ifelse(motus7$scientific_name_final== "Fistularia commersonii", "NIS",
                          ifelse(motus7$scientific_name_final== "Lutjanus sp.", "NIS",
                                 ifelse(motus7$scientific_name_final== "Odonus niger", "NIS",
                                        #ifelse(motus7$scientific_name_final== "Parablennius yatabei", "NIS",
                                               ifelse(motus7$scientific_name_final== "Siganus sp.", "NIS",
                                                      ifelse(motus7$scientific_name_final== "Siganus luridus", "NIS",
                                                             ifelse(motus7$scientific_name_final== "Siganus rivulatus", "NIS",
                                                                    ifelse(motus7$scientific_name_final== "Tylosurus crocodilus", "NIS",
                                                                           ifelse(motus7$scientific_name_final== "Tylosurus sp.", "NIS", "Native")))))))))#)


## version 1 native vs NIS; other Siganus species aren't labeled NIS properly for some reason

windows()
ggplot(motus7, aes(x = sample_type, y = scientific_name_final, fill = prop2)) + 
  scale_y_discrete(limits = rev) +
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1) +
  scale_fill_gradient(low = "#1E88E5", high = "#D81B60", na.value = "white",
                      name = "Relative \nProportion (%)",
                      limits = c(0,50),
                      breaks = c(0.0001, 10, 20, 30, 40, 50),
                      labels = c("> 0","10","20","30","40", "50")) +
  facet_grid(nis ~ location, scales = "free", space = "free") +
  labs(x = "", y ="") +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"), 
        axis.ticks = element_line(colour = "#000000"),
        strip.text.y = element_text(colour = "#000000", angle = 360),
        axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "#000000", face = "italic"),
        text = element_text(size = 14),
        legend.direction = "vertical", 
        legend.box = "vertical") 

## version 2 family sorted

meta2 <- read.csv("meta2.csv")
motus8 <- merge(motus7, meta2, by=c("scientific_name_final"), all.x = TRUE) 

heatmap1 <-
ggplot(motus8, aes(x = sample_type, y = scientific_name_final, fill = prop2)) + 
  scale_y_discrete(limits = rev) +
  geom_tile(color = "black",
            lwd = 0.3,
            linetype = 1) +
  scale_fill_gradient(low = "#1E88E5", high = "#D81B60", na.value = "white",
                      name = "Relative \nProportion (%)",
                      limits = c(0,50),
                      breaks = c(0.0001, 10, 20, 30, 40, 50),
                      labels = c("> 0","10","20","30","40", "50")) +
  facet_grid(family ~ location, scales = "free", space = "free") +
  labs(x = "", y ="") +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"), 
        axis.ticks = element_line(colour = "#000000"),
        strip.text.y = element_text(colour = "#000000", angle = 360),
        axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "#000000", face = "italic"),
        text = element_text(size = 14),
        legend.direction = "vertical", 
        legend.box = "vertical") 

ggsave(filename="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/Objective2_Erika/sequences_bubbleplot/heatmap1.jpg", 
       plot = heatmap1, width = 8, height = 12, dpi = 300, units = "in")

## modify heatmap1

#fix nis column
motus8$nis <- ifelse(motus8$nis2== "NIS", "NIS", "Native")

#heatmap2 <-
  ggplot(motus8, aes(x = sample_type, y = scientific_name_final)) + 
  scale_y_discrete(limits = rev) +
  geom_point(aes(shape = nis, fill = prop2), size=5) +
  #geom_tile(color = "black",
  #          lwd = 0.3,
   #         linetype = 1) +
  scale_shape_manual(values=c(21,23))+
  scale_fill_gradient(low = "#1E88E5", high = "#D81B60", na.value = "white",
                      name = "Relative \nProportion (%)",
                      limits = c(0,50),
                      breaks = c(0.0001, 10, 20, 30, 40, 50),
                      labels = c("> 0","10","20","30","40", "50")) +
  facet_grid(family ~ location, scales = "free", space = "free") +
  labs(x = "", y ="") +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"), 
        axis.ticks = element_line(colour = "#000000"),
        strip.text.y = element_text(colour = "#000000", angle = 360),
        axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "#000000", face = "italic"),
        text = element_text(size = 14),
        legend.direction = "vertical", 
        legend.box = "vertical") 


heatmap3 <-
    ggplot(motus8, aes(x = sample_type, y = scientific_name_final2, fill = prop2)) + 
    scale_y_discrete(limits = rev) +
    geom_tile(color = "black",
              lwd = 0.3,
              linetype = 1) +
    scale_fill_gradient(low = "#1E88E5", high = "#D81B60", na.value = "white",
                        name = "Relative \nProportion (%)",
                        limits = c(0,50),
                        breaks = c(0.0001, 10, 20, 30, 40, 50),
                        labels = c("> 0","10","20","30","40", "50")) +
    facet_grid(family ~ location, scales = "free", space = "free") +
    labs(x = "", y ="") +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "lightgrey"),
          panel.grid.minor = element_line(colour = "lightgrey"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical") 
  
  ggsave(filename="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/Objective2_Erika/sequences_bubbleplot/heatmap_nis.jpg", 
         plot = heatmap3, width = 8, height = 12, dpi = 300, units = "in")

## repeat heatmap 3 but add sponges

motus7.5 <- motus6 %>% 
    group_by(scientific_name_final, sample_type, location) %>%
    summarise(sum_counts=sum(counts),
              sum_pergroup= sum(x_persample),
              .groups = 'drop') %>%
    as.data.frame()  
  
  ## Add proportions
  motus7.5$prop <- (motus7.5$sum_counts/motus7.5$sum_pergroup)*100 
  
  motus7.5$prop2 <- ifelse(motus7.5$prop== 0, NA, motus7.5$prop)
  
  max(motus7.5$prop)
  
  ## Add metadata

  motus8.5 <- merge(motus7.5, meta2, by=c("scientific_name_final"), all.x = TRUE) 
  
heatmap4 <-
    ggplot(motus8.5, aes(x = sample_type, y = scientific_name_final2, fill = prop2)) + 
    scale_y_discrete(limits = rev) +
    geom_tile(color = "black",
              lwd = 0.3,
              linetype = 1) +
    scale_fill_gradient(low = "#1E88E5", high = "#D81B60", na.value = "white",
                        name = "Relative \nProportion (%)",
                        limits = c(0,63),
                        breaks = c(0.0001, 10, 20, 30, 40, 50, 63),
                        labels = c("> 0","10","20","30","40", "50", "60")) +
    facet_grid(family ~ location, scales = "free", space = "free") +
    labs(x = "", y ="") +
    theme(strip.background.y = element_rect(fill= "lightgrey"),
          strip.background.x = element_rect(fill= "lightgrey"),
          strip.text.x= element_text(colour = "#000000"),
          panel.grid.major = element_line(colour = "lightgrey"),
          panel.grid.minor = element_line(colour = "lightgrey"), 
          axis.ticks = element_line(colour = "#000000"),
          strip.text.y = element_text(colour = "#000000", angle = 360),
          axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
          axis.text.y = element_text(colour = "#000000", face = "italic"),
          text = element_text(size = 14),
          legend.direction = "vertical", 
          legend.box = "vertical") 
  
ggsave(filename="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/Objective2_Erika/sequences_bubbleplot/heatmap_inclsponge_nis.jpg", 
         plot = heatmap4, width = 8, height = 12, dpi = 300, units = "in")

#####
# Prepare data for bubble plot comparing eDNA and sponges
#####

## need to make a dataframe that has the eDNA samples combined, the sponge DNA samples combined, 
## and the sponge samples separated out
  
## remove UVC samples
## calculate relative proportions by groups, sponge and eDNA from two locations
  
motus9 <- motus6 %>% filter(sample_type=='sponge' | sample_type=='eDNA') %>%
    group_by(scientific_name_final, sample_type, location) %>%
    summarise(counts=sum(counts),
              reads_persample= sum(reads_persample),
              .groups = 'drop') %>%
    as.data.frame()  

## subset just the sponge samples and add them to motus 9
  
motus6.5 <- subset(motus6, sample_type== "sponge")  
motus6.5 <- motus6.5[c(7,11,13,10,17,1)]
motus9$id_seq <- motus9$sample_type

##
motus10 <- rbind(motus6.5, motus9)

motus10$sqrt_reads <- sqrt(motus10$counts)
min(motus10$sqrt_reads) #0
max(motus10$sqrt_reads) #150
motus10$sqrt_reads <- ifelse(motus10$sqrt== 0, NA, motus10$sqrt_reads)


motus10 <- merge(motus10, meta2, by=c("scientific_name_final"), all.x = TRUE) 
motus10$nis <- ifelse(motus10$nis2== "NIS", "NIS", "Native")

motus11 <- na.omit(motus10)

meta3 <- read.csv("meta3.csv")
motus11 <- merge(motus11, meta3, by=c("id_seq"), all.x = TRUE) 

motus11$id1 <- factor(motus11$id1,
                               levels = c("All eDNA","All sponges","A1",     
                                          "A2", "A3", "A4",    
                                          "B1", "B2", "D1", "D2", 
                                          "D3", "D4", 
                                          "D5", "D6",  
                                          "D7", "D8"))

## prepare bubbleplot

#windows()
bubbleplot1 <-
ggplot(motus11, aes(x = id1, y = scientific_name_final, fill = sample_type, size = sqrt_reads)) + 
  scale_y_discrete(limits = rev) +
  geom_point(pch = 21) +
  scale_fill_manual(name = "Sample Type", values =c("#004D40","#FFC107")) +
  scale_size_continuous(name = expression(sqrt(Reads)),
                        range = c(1, 6),
                        limits = c(1,150),
                        breaks = c(1, 25, 50, 75, 100, 150),
                        labels = c("1","25","50","75","100", "150")) + 
  facet_grid(nis ~ location, scales = "free", space = "free") +
  labs(x = "", y ="") +
  theme_light() +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"), 
        axis.ticks = element_line(colour = "#000000"),
        strip.text.y = element_text(colour = "#000000", angle = 360),
        axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "#000000", face = "italic"),
        text = element_text(size = 14),
        legend.direction = "vertical", 
        legend.box = "vertical")  

ggsave(filename="C:/Users/beseneav/OneDrive - Liverpool John Moores University/PhD/Objective2_Erika/sequences_bubbleplot/bubblplot_molecular_nis.jpg", 
       plot = bubbleplot1, width = 10, height = 8, dpi = 300, units = "in")


