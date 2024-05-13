###################################
## NMDS and PERMANOVA for Greece ##
###################################

library(vegan)
library(tidyverse)

## create second dataframe for presence absence
samp2 <- ifelse(samp == 0, 0, 1)

## calculate richness
mdata$Richness <- colSums(samp2)

#####
## Create distance matrix and NMDS
#####

## Calculate jaccard dissimilariy matrix
jaccard.dist <- vegdist(t(samp2), method = "jaccard", binary =  TRUE)

## NMDS with jaccard distances

set.seed(12)

nmds.jaccard.ord <- metaMDS(jaccard.dist, distance = "jaccard", trymax = 1000)
plot(nmds.jaccard.ord)
 
nmds.jaccard.ord.k3 <- metaMDS(jaccard.dist, distance = "jaccard", k=3, trymax = 1000)
nmds.jaccard.ord.k3

## extract site scores from NMDS and add metadata

data.scores.k2 <- as.data.frame(scores(nmds.jaccard.ord))  
data.scores.k2$id_seq <- colnames(samp)
nmds_pa <- merge(data.scores.k2, meta, by.x = "id_seq")


#####
## Create distance matrix and NMDS for UVC data
#####

samp_uvc <- motus3.5[motus3.5$total_UVC != 0, ]  

samp_uvc_2 <- samp_uvc[c(25:48)]
samp_uvc_pa <- ifelse(samp_uvc_2 == 0, 0, 1)

## Calculate jaccard dissimilariy matrix
jaccard.dist.uvc <- vegdist(t(samp_uvc_pa), method = "jaccard", binary =  TRUE)

## Hellinger's standardization
uvc.hell <- decostand(t(samp_uvc_2), method = "hellinger")

## calculate Bray-Curtis dissimilarity matrix
bray.dist.uvc <- vegdist(uvc.hell, method = "bray")


## nmds
set.seed(13)

nmds.jaccard.ord.uvc <- metaMDS(jaccard.dist.uvc, distance = "jaccard", trymax = 1000)
nmds.jaccard.ord.uvc

plot(nmds.jaccard.ord.uvc, type = "t")


nmds.bray.ord.uvc <- metaMDS(bray.dist.uvc, distance = "bray", trymax = 1000)
nmds.bray.ord.uvc

plot(nmds.bray.ord.uvc, type = "t")

#####
## Repeat NMDS for UVC data but add observers observations together
#####

write.csv(samp_uvc_pa, "samp_uvc_pa.csv")
samp_uvc_pa2 <- read.csv("samp_uvc_pa_observersadd.csv")

jaccard.dist.uvc2 <- vegdist(t(samp_uvc_pa2), method = "jaccard", binary =  TRUE)

set.seed(24)

nmds.jaccard.ord.2 <- metaMDS(jaccard.dist.uvc2, distance = "jaccard", trymax = 1000)
plot(nmds.jaccard.ord.2, type ="t")

write.csv(samp_uvc_2, "samp_uvc.csv")
samp_uvc_add <- read.csv("samp_uvc_observeradd.csv")

## Hellinger's standardization
uvc.hell <- decostand(t(samp_uvc_add), method = "hellinger")

## calculate Bray-Curtis dissimilarity matrix
bray.dist.uvc <- vegdist(uvc.hell, method = "bray")

nmds.bray.ord.uvc <- metaMDS(bray.dist.uvc, distance = "bray", trymax = 1000)
nmds.bray.ord.uvc

plot(nmds.bray.ord.uvc, type = "t")
#####
## Create distance matrix and NMDS for eDNA data
#####

samp_edna <- motus3.5[-c(25:49,51,53:56,59,61,62,66:68,70,71,74,75)]
samp_edna$total_reads <- rowSums(samp_edna[c(25:35)])

samp_edna1 <- samp_edna[samp_edna$total_reads != 0, ]  

samp_edna_2 <- samp_edna1[c(25:35)]
samp_edna_pa <- ifelse(samp_edna_2 == 0, 0, 1)

## Calculate jaccard dissimilariy matrix
jaccard.dist.edna <- vegdist(t(samp_edna_pa), method = "jaccard", binary =  TRUE)

## Hellinger's standardization
edna.hell <- decostand(t(samp_edna_2), method = "hellinger")

## calculate Bray-Curtis dissimilarity matrix
bray.dist.edna <- vegdist(edna.hell, method = "bray")


## nmds
set.seed(14)

nmds.jaccard.ord.edna <- metaMDS(jaccard.dist.edna, distance = "jaccard", trymax = 1000)
nmds.jaccard.ord.edna

plot(nmds.jaccard.ord.edna, type = "t")
stressplot(nmds.jaccard.ord.edna)

nmds.bray.ord.edna <- metaMDS(bray.dist.edna, distance = "bray", trymax = 1000)
nmds.bray.ord.edna

##### Saved plot as the outlier beta-diversity
plot(nmds.bray.ord.edna, type = "t") # This plot shows outliers, maybe remove ?
stressplot(nmds.bray.ord.edna)

#####
### Repeat but remove eDNA data with 1-3 detections 
#####

samp_ednarm <- motus3.5[-c(25:49,51,53:56,59,61,62:64,66:68,70,71,74,75)]
samp_ednarm$total_reads <- rowSums(samp_ednarm[c(25:33)])

samp_edna1rm <- samp_ednarm[samp_ednarm$total_reads != 0, ]  

samp_edna_2rm <- samp_edna1rm[c(25:33)]
samp_edna_parm <- ifelse(samp_edna_2rm == 0, 0, 1)

## Calculate jaccard dissimilariy matrix
jaccard.dist.ednarm <- vegdist(t(samp_edna_parm), method = "jaccard", binary =  TRUE)

## Hellinger's standardization
edna.hellrm <- decostand(t(samp_edna_2rm), method = "hellinger")
##sqrt read transformation
edna.sqrt <- sqrt(t(samp_edna_2rm))

## calculate Bray-Curtis dissimilarity matrix
bray.dist.ednarm <- vegdist(edna.hellrm, method = "bray")
bray.sr.dist.ednarm <- vegdist(edna.sqrt, method = "bray")

## nmds
set.seed(15)

nmds.jaccard.ord.ednarm <- metaMDS(jaccard.dist.ednarm, distance = "jaccard", trymax = 1000)
nmds.jaccard.ord.ednarm

plot(nmds.jaccard.ord.ednarm, type = "t")


nmds.bray.ord.ednarm <- metaMDS(bray.dist.ednarm, distance = "bray", trymax = 1000)
nmds.bray.ord.ednarm

plot(nmds.bray.ord.ednarm, type = "t")


nmds.bray.sr.ord.ednarm <- metaMDS(bray.sr.dist.ednarm, distance = "bray", trymax = 1000)

plot(nmds.bray.sr.ord.ednarm, type = "t")

#####
### Repeat but remove eDNA data from three samples that were outliers
#####

## Need to edit below code

samp_ednarm <- motus3.5[-c(25:49,51,53:56,59,61,62:64,66:68,70,71,74,75)]
samp_ednarm_v2 <- samp_ednarm[-c(29)] #also remove KE_B_2eDNA
samp_ednarm_v2$total_reads <- rowSums(samp_ednarm_v2[c(25:32)])

samp_edna1rm_v2 <- samp_ednarm_v2[samp_ednarm_v2$total_reads != 0, ]  

## Remove Parablennius yatabei - only 93% to genus and species level...
samp_edna1rm_v2.2 <- samp_edna1rm_v2[!samp_edna1rm_v2$scientific_name_final== "Parablennius yatabei",]


samp_edna_nooutlier <- samp_edna1rm_v2.2[c(25:32)]
samp_edna_pa_nooutlier <- ifelse(samp_edna_nooutlier == 0, 0, 1)

## Calculate jaccard dissimilariy matrix
jaccard.dist.nooutlier <- vegdist(t(samp_edna_pa_nooutlier), method = "jaccard", binary =  TRUE)

## Hellinger's standardization
edna.hell.nooutlier <- decostand(t(samp_edna_nooutlier), method = "hellinger")

## calculate Bray-Curtis dissimilarity matrix
bray.dist.edna.nooutlier <- vegdist(edna.hell.nooutlier, method = "bray")

## nmds
set.seed(34)

nmds.jaccard.ord.edna.nooutlier <- metaMDS(jaccard.dist.nooutlier, distance = "jaccard", trymax = 1000)

plot(nmds.jaccard.ord.edna.nooutlier, type = "t")
stressplot(nmds.jaccard.ord.edna.nooutlier)


nmds.bray.ord.nooutlier <- metaMDS(bray.dist.edna.nooutlier, distance = "bray", trymax = 1000)
nmds.bray.ord.nooutlier

plot(nmds.bray.ord.nooutlier, type = "t")
stressplot(nmds.bray.ord.nooutlier)



#####
## Plot jaccard NMDS plots in ggplot2
#####

library(ggplot2)
library(ggtext)

##### UVC

## extract site scores from NMDS and add metadata
data.scores.uvc <- as.data.frame(scores(nmds.jaccard.ord.uvc, display ="sites"))  
data.scores.uvc$id_seq <- colnames(samp_uvc_pa)
nmds_uvc <- merge(data.scores.uvc, meta, by.x = "id_seq")
nmds_uvc$Richness <- colSums(samp_uvc_pa)

## colour by location

ggplot() + 
  geom_point(data = nmds_uvc, aes(x = NMDS1, y = NMDS2, fill = location), 
             colour = c("black"), shape = 21, alpha=0.8, size = 9)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  geom_text(data = nmds_uvc, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="black")+
  annotate(geom = "text", x = 0.5, y= -0.7, label = "stress = 0.1607",size = 5,colour ="black")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"),
        #legend.position = c(0.1,0.5),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.key=element_blank(), #hides grey background in legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", colour = "Depth range (m)", y = "NMDS2", size = "Richness")  + 
  labs(x = "NMDS1", y = "NMDS2", fill="Location") +
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) +
  xlim(-0.8,0.8) +
  ylim(-0.8,0.8) 
  
## colour by observer

ggplot() + 
  geom_point(data = nmds_uvc, aes(x = NMDS1, y = NMDS2, fill = observer), 
             colour = c("black"), shape = 21, alpha=0.8, size = 9)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  geom_text(data = nmds_uvc, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="black")+
  annotate(geom = "text", x = 0.5, y= -0.7, label = "stress = 0.1607",size = 5,colour ="black")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"),
        #legend.position = c(0.1,0.5),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.key=element_blank(), #hides grey background in legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", colour = "Depth range (m)", y = "NMDS2", size = "Richness")  + 
  labs(x = "NMDS1", y = "NMDS2", fill="Location") +
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) +
  xlim(-0.8,0.8) +
  ylim(-0.8,0.8)

###### 
##UVC OBSERVERS combined
#####
## extract site scores from NMDS and add metadata
data.scores.uvc2 <- as.data.frame(scores(nmds.jaccard.ord.2, display ="sites"))  
data.scores.uvc2$id_seq <- colnames(samp_uvc_pa2)


id_seq <- c("DM_A", "DM_B", "DM_C", "DE_A", "DE_B", "DE_C", "KM_A", "KM_B", "KM_C", "KE_A", "KE_B", "KE_C")
Location <- c("Dafni Beach", "Dafni Beach", "Dafni Beach",
              "Dafni Beach", "Dafni Beach", "Dafni Beach",
              "Korakonissi", "Korakonissi", "Korakonissi",
              "Korakonissi", "Korakonissi", "Korakonissi")
mdata2 <- data.frame(id_seq,Location)
rm(id_seq,Location)

nmds_uvc2 <- merge(data.scores.uvc2, mdata2, by.x = "id_seq")
nmds_uvc2$Richness <- colSums(samp_uvc_pa2)

uvc_nmds <-
ggplot() + 
  geom_point(data = nmds_uvc2, aes(x = NMDS1, y = NMDS2, fill = Location), 
             colour = c("black"), shape = 21, alpha=0.8, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  geom_text(data = nmds_uvc2, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 0.5, y= -0.7, label = "stress = 0.1607",size = 5,colour ="black")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"),
        #legend.position = c(0.1,0.5),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.key=element_blank(), #hides grey background in legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", colour = "Depth range (m)", y = "NMDS2", size = "Richness")  + 
  labs(x = "NMDS1", y = "NMDS2", fill="Location") +
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) +
  xlim(-0.7,0.7) +
  ylim(-0.7,0.7) 


uvc_nmds_no_legend <-
  ggplot() + 
  geom_point(data = nmds_uvc2, aes(x = NMDS1, y = NMDS2, fill = Location), 
             colour = c("black"), shape = 21, alpha=0.8, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  geom_text(data = nmds_uvc2, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 0.5, y= -0.7, label = "stress = 0.1607",size = 5,colour ="black")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        #legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        #legend.text = element_text(size = 11, face ="bold", colour ="black"),
        legend.position = "none",
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        #legend.key=element_blank(), #hides grey background in legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", colour = "Depth range (m)", y = "NMDS2", size = "Richness")  + 
  labs(x = "NMDS1", y = "NMDS2", fill="Location") +
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) +
  xlim(-0.7,0.7) +
  ylim(-0.7,0.7) 

##### 
##eDNA
#####

## extract site scores from NMDS and add metadata

data.scores.edna.nooutlier.j <- as.data.frame(scores(nmds.jaccard.ord.edna.nooutlier, display ="sites"))  
data.scores.edna.nooutlier.j$id_seq <- rownames(data.scores.edna.nooutlier.j)
nmds_edna <- merge(data.scores.edna.nooutlier.j, mdata, by.x = "id_seq")
#data.scores.edna <- as.data.frame(scores(nmds.jaccard.ord.ednarm, display ="sites"))  
#data.scores.edna$id_seq <- colnames(samp_edna_parm)
#nmds_edna <- merge(data.scores.edna, mdata, by.x = "id_seq")


## colour by location

edna_nmds <-
ggplot() + 
  geom_point(data = nmds_edna, aes(x = NMDS1, y = NMDS2, fill = location), 
             colour = c("black"), shape = 21, alpha=0.8, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  geom_text(data = nmds_edna, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 1, y= -1.5, label = "stress = 0.0373",size = 5,colour ="black")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"),
        legend.direction = "vertical", 
        legend.box = "vertical",
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) +
  labs(x = "NMDS1", colour = "Depth range (m)", y = "NMDS2", size = "Richness")  + 
  labs(x = "NMDS1", y = "NMDS2", fill="Location") +
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) +
  xlim(-1.5,1.5) +
  ylim(-1.5,1.5) 

#####
## Test for differences in beta-diversity
#####

library("vegan")

## Calculate beta-dispersion
uvc_bd <- betadisper(jaccard.dist.uvc2, nmds_uvc2$Location)

## Test for differences in group dispersion
anova(uvc_bd)
TukeyHSD(uvc_bd)

## Test whether the groups have different compositions
uvc_ad <- adonis2(jaccard.dist.uvc2~nmds_uvc2$Location, permutations = 999) ## significant p = 0.003

data.scores.uvc <- as.data.frame(scores(nmds.bray.ord.uvc, display ="sites")) 
data.scores.uvc$id_seq <- rownames(data.scores.uvc)
nmds_uvc_b <- merge(data.scores.uvc, mdata2, by.x = "id_seq")

uvc_ad2 <- adonis2(bray.dist.uvc~nmds_uvc_b$Location, permutations = 999) ## significant p = 0.004

## repeat for eDNA

## Calculate beta-dispersion
edna_bd <- betadisper(jaccard.dist.nooutlier, nmds_edna_nooutlier$Location3)

## Test for differences in group dispersion
anova(edna_bd)
TukeyHSD(edna_bd)

## Test whether the groups have different compositions
edna_ad <- adonis2(jaccard.dist.ednarm~nmds_edna$location, permutations = 999) ## not significant

## 2 samples removed, 1 detection and 3 detections
edna_ad2 <- adonis2(bray.dist.ednarm~nmds_edna$location, permutations = 999) ##significant with bray-curtis

## all samples
data.scores.edna.orig <- as.data.frame(scores(nmds.bray.ord.edna, display ="sites"))  
Location2 <- c("Dafni Beach", "Dafni Beach", "Dafni Beach",
              "Dafni Beach", 
              "Korakonissi", "Korakonissi", "Korakonissi",
              "Korakonissi", "Korakonissi", "Korakonissi",
              "Korakonissi")
nmds_edna_orig <- data.frame(data.scores.edna.orig,Location2)

edna_ad3 <- adonis2(bray.dist.edna~nmds_edna_orig$Location2, permutations = 999) ##not significant

## what if you remove all outliers 
data.scores.edna.nooutlier <- as.data.frame(scores(nmds.bray.ord.nooutlier, display ="sites"))  
data.scores.edna.nooutlier.j <- as.data.frame(scores(nmds.jaccard.ord.edna.nooutlier, display ="sites"))  
Location3 <- c("Dafni Beach", "Dafni Beach", "Dafni Beach",
               "Dafni Beach", 
               "Korakonissi", "Korakonissi", "Korakonissi",
               "Korakonissi")
nmds_edna_nooutlier <- data.frame(data.scores.edna.nooutlier,Location3)
nmds_edna_nooutlier_j <- data.frame(data.scores.edna.nooutlier.j,Location3)


edna_ad4 <- adonis2(bray.dist.edna.nooutlier~nmds_edna_nooutlier$Location3, permutations = 999) ## significant p = 0.027

edna_sd5 <- adonis2(jaccard.dist.nooutlier~nmds_edna_nooutlier$Location3, permutations = 100) ## significant

############################
## Bray-Curtis NMDS plots ##
############################

## UVC Bray-Curtis NMDS

## extract site scores from NMDS and add metadata
data.scores.uvc3 <- as.data.frame(scores(nmds.bray.ord.uvc, display ="sites"))  
data.scores.uvc3$id_seq <- colnames(samp_uvc_add)


id_seq <- c("DM_A", "DM_B", "DM_C", "DE_A", "DE_B", "DE_C", "KM_A", "KM_B", "KM_C", "KE_A", "KE_B", "KE_C")
Location <- c("Dafni Beach", "Dafni Beach", "Dafni Beach",
              "Dafni Beach", "Dafni Beach", "Dafni Beach",
              "Korakonissi", "Korakonissi", "Korakonissi",
              "Korakonissi", "Korakonissi", "Korakonissi")
mdata2 <- data.frame(id_seq,Location)
rm(id_seq,Location)

nmds_uvc3 <- merge(data.scores.uvc3, mdata2, by.x = "id_seq")
nmds_uvc3$Richness <- colSums(samp_uvc_pa2)


uvc_nmds_no_legend <-
  ggplot() + 
  geom_point(data = nmds_uvc3, aes(x = NMDS1, y = NMDS2, fill = Location), 
             colour = c("black"), shape = 21, alpha=0.8, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  geom_text(data = nmds_uvc3, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 0.5, y= -0.7, label = "stress = 0.0591",size = 5,colour ="black")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        #legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        #legend.text = element_text(size = 11, face ="bold", colour ="black"),
        legend.position = "none",
        #legend.background = element_blank(),
        #legend.box.background = element_rect(colour = "black"),
        #legend.key=element_blank(), #hides grey background in legend
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", colour = "Depth range (m)", y = "NMDS2", size = "Richness")  + 
  labs(x = "NMDS1", y = "NMDS2", fill="Location") +
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) +
  xlim(-0.7,0.7) +
  ylim(-0.7,0.7) 


## eDNA Bray-Curtis NMDS

data.scores.edna.nooutlier.b <- as.data.frame(scores(nmds.bray.ord.nooutlier, display ="sites"))  
data.scores.edna.nooutlier.b$id_seq <- rownames(data.scores.edna.nooutlier.b)
nmds_edna2 <- merge(data.scores.edna.nooutlier.b, mdata, by.x = "id_seq")


edna_nmds2 <-
  ggplot() + 
  geom_point(data = nmds_edna2, aes(x = NMDS1, y = NMDS2, fill = location), 
             colour = c("black"), shape = 21, alpha=0.8, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  geom_text(data = nmds_edna2, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 1.1, y= -0.8, label = "stress = 0.0693",size = 5,colour ="black")+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"),
        legend.direction = "vertical", 
        legend.box = "vertical",
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) +
  labs(x = "NMDS1", colour = "Depth range (m)", y = "NMDS2", size = "Richness")  + 
  labs(x = "NMDS1", y = "NMDS2", fill="Location") +
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) +
  xlim(-1.5,1.5) +
  ylim(-0.8,0.85) 




