#######################################
## Diversity indices for Greece data ##
## box plots for results             ##
#######################################

library(vegan)
library(BiodiversityR)
library(tidyverse)

#####
## Format data
#####

meta <- read.csv("greece_eDNA_sponge_metatdata.csv")
motus3.5 <- read.csv("motus_contam_removed_over100reads_edit_UVCadded.csv") 

# Remove Parablennius yatabei - only 93% to genus and species level...
motus3.5 <- motus3.5[!motus3.5$scientific_name_final== "Parablennius yatabei",]

## just samples

samp <- motus3.5[c(25:48,50:74)]

## Remove samples with only one species
# sample.DBE_A_BROWN1
# sample.DBE_C_BROWN3
# sample.DBM_B_BROWN2et
# sample.KE_B_PURPLE2
# sample.KE_C_2eDNA
# sample.KM_A_PURPLE1
# sample.KM_A_PURPLE3
# sample.KM_C_PURPLE3

samp <- subset(samp, select = -c(sample.DBE_A_BROWN1,
                                 sample.DBE_C_BROWN3,
                                 sample.DBM_B_BROWN2et,
                                 sample.KE_B_PURPLE2,
                                 sample.KE_C_2eDNA,
                                 sample.KM_A_PURPLE1,
                                 sample.KM_A_PURPLE3,
                                 sample.KM_C_PURPLE3))


## organise metadata in same order as dataframe 
sampheader <- as.data.frame(colnames(samp))
colnames(sampheader)[1] <- "id_seq"
mdata <- inner_join(sampheader, meta)

## Remove rows that sum to zero (get rid of snappers, no longer represented after sample removal)
samp <- samp[rowSums(samp[])>0,] 

## make species table, make sure it's numeric

spec <- as.data.frame(t(motus3.5[c(24:48,50:74)]))

## Remove samples with only one species
row_names_remove<-c("sample.DBE_A_BROWN1",
                    "sample.DBE_C_BROWN3",
                    "sample.DBM_B_BROWN2et",
                    "sample.KE_B_PURPLE2",
                    "sample.KE_C_2eDNA",
                    "sample.KM_A_PURPLE1",
                    "sample.KM_A_PURPLE3",
                    "sample.KM_C_PURPLE3")
spec <- spec[!(row.names(spec) %in% row_names_remove),]


colnames(spec) <- spec[1,]
spec <- spec[c(2:42),]
spec <- mutate_all(spec, function(x) as.numeric(as.character(x)))

## Remove Lutjanus sp. because it have no observations
spec <- spec[-c(27)]

#####
## Diversity index calculations with greece data
#####

## BiodiversityR package lets you compare indices by a grouping variable

mdata$sample_type <- as.factor(mdata$sample_type)
mdata$location <- as.factor(mdata$location)

diversitycomp(spec, y = mdata, factor1 = "sample_type", index = "Shannon")

#sample_type  n   Shannon
#eDNA   10 2.4776590
#sponge  7 1.3092507
#UVC    24 2.1574934

diversitycomp(spec, y = mdata, factor1 = "location", index = "Shannon")

#location       n   Shannon
#Dafni Beach 18 2.0310613
#Korakonissi 23 2.0319370

diversitycomp(spec, y = mdata, factor1 = "sample_type", factor2 = "location", index = "Shannon")

#, ,  = n

#location
#sample_type Dafni Beach Korakonissi
#eDNA             4           6
#sponge           2           5
#UVC             12          12

#, ,  = Shannon

#location
#sample_type Dafni Beach Korakonissi
#eDNA    2.09233456  2.04426514
#sponge  0.63290373  0.84383304
#UVC     1.94984346  1.86928225


## See results
# sites/samples seperate; individual data points for edge
#diversityresult(spec, y = mdata, factor ="sample_type",level = "UVC", index = "Shannon", method = "each site")

#diversityresult(spec, y = mdata, factor ="sample_type",level = "eDNA", index = "Shannon", method = "each site")

shan <- diversityresult(spec, index = "Shannon", method = "each site")

## calculate the J-evenness
even <- diversityresult(spec, index = "Shannon", method = "each site") / log(specnumber(spec))
colnames(even)[1] <- "Evenness"

## calculate the Richness
rich <- diversityresult(spec, index = "richness", method = "each site")
colnames(rich)[1] <- "Richness"

shanplot <- cbind(mdata, shan, even, rich)


#####
## Box plots
#####

library(ggplot2)
library(cowplot)

##Figure 5A
splot <-
  ggplot(shanplot, aes(x=sample_type, y=Shannon, fill = sample_type)) + 
  scale_fill_manual(values = c("#004D40", "#FFC107", "#1E88E5")) +
  geom_boxplot(shape=16, colour = "black") +
  facet_grid(. ~ location) +
  labs(x = "", y ="Shannon Index") +
  theme_light() +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000", face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"), 
        strip.text.y = element_text(colour = "#000000", angle = 360, face = "bold"),
        text = element_text(size = 14),
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.position = "none")

## richness plot, not used in paper
rplot <-
  ggplot(shanplot, aes(x=sample_type, y=Richness, fill = sample_type)) + 
  scale_fill_manual(values = c("#004D40", "#FFC107", "#1E88E5")) +
  geom_boxplot(shape=16, colour = "black") +
  facet_grid(. ~ location) +
  labs(x = "", y ="Richness") +
  theme_light() +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"), 
        axis.ticks = element_line(colour = "#000000"),
        strip.text.y = element_text(colour = "#000000", angle = 360),
        axis.text.x = element_text(colour = "#000000"),
        axis.text.y = element_text(colour = "#000000"),
        text = element_text(size = 14),
        legend.position = "none")

## Figure 5B
eplot <-
  ggplot(shanplot, aes(x=sample_type, y=Evenness, fill = sample_type)) + 
  scale_fill_manual(name = "Sample Type",
                    values = c("#004D40", "#FFC107", "#1E88E5")) +
  geom_boxplot(shape=16, colour = "black") +
  facet_grid(. ~ location) +
  labs(x = "", y ="J-Evenness") +
  theme_light() +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000", face = "bold"),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"), 
        strip.text.y = element_text(colour = "#000000", angle = 360, face = "bold"),
        text = element_text(size = 14),
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        legend.text = element_text(size = 11, face ="bold", colour ="black"),
        legend.direction = "vertical", 
        legend.box = "vertical",
        axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"))  


#####
## Test for significant differences
#####
library(FSA)
###
## Dafni Beach
###

shanplot_db <- subset(shanplot, location=="Dafni Beach")

shanplot_db$sample_type <- ordered(shanplot_db$sample_type,
                                   levels = c("eDNA", "sponge", "UVC"))

levels(shanplot_db$sample_type)

## Shannon index
group_by(shanplot_db, sample_type) %>%
  summarise(
    count = n(),
    mean = mean(Shannon, na.rm = TRUE),
    sd = sd(Shannon, na.rm = TRUE),
    median = median(Shannon, na.rm = TRUE),
    IQR = IQR(Shannon, na.rm = TRUE)
  )

## Kruskal-wallis
#kruskal.test(Shannon ~ sample_type, data = shanplot_db) #not significant, p value = 0.06
## Wilcoxin rank-sum
#pairwise.wilcox.test(shanplot_db$Shannon, shanplot_db$sample_type,
#                     p.adjust.method = "BH")
#
## Mann-Whitney U test eDNA and UVC
shanplot_db2 <- subset(shanplot_db, sample_type=="eDNA"| sample_type=="UVC")
wilcox.test(Shannon ~ sample_type, data=shanplot_db2, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE) #p = 0.4306


## J-evenness
group_by(shanplot_db, sample_type) %>%
  summarise(
    count = n(),
    mean = mean(Evenness, na.rm = TRUE),                           
    sd = sd(Evenness, na.rm = TRUE),
    median = median(Evenness, na.rm = TRUE),
    IQR = IQR(Evenness, na.rm = TRUE)
  )

## Mann-Whitney U test eDNA and UVC
wilcox.test(Evenness ~ sample_type, data=shanplot_db2, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE) #p = 0.5853


###
## Korakonissi
###
shanplot_k <- subset(shanplot, location=="Korakonissi")

levels(shanplot_k$sample_type)

## Shannon index
group_by(shanplot_k, sample_type) %>%
  summarise(
    count = n(),
    mean = mean(Shannon, na.rm = TRUE),
    sd = sd(Shannon, na.rm = TRUE),
    median = median(Shannon, na.rm = TRUE),
    IQR = IQR(Shannon, na.rm = TRUE)
  )

## Kruskal-wallis
kruskal.test(Shannon ~ sample_type, data = shanplot_k) #significant, p value = 0.0004456
## Wilcoxin rank-sum
pairwise.wilcox.test(shanplot_k$Shannon, shanplot_k$sample_type,
                     p.adjust.method = "BH", paired = FALSE)
# UVC and sponge; significantly different

## Mann-Whitney U test eDNA and UVC
shanplot_k2 <- subset(shanplot_k, sample_type=="eDNA"| sample_type=="UVC")
wilcox.test(Shannon ~ sample_type, data=shanplot_k2, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE) #p = 0.0832


## J-evenness
group_by(shanplot_k, sample_type) %>%
  summarise(
    count = n(),
    mean = mean(Evenness, na.rm = TRUE),
    sd = sd(Evenness, na.rm = TRUE),
    median = median(Evenness, na.rm = TRUE),
    IQR = IQR(Evenness, na.rm = TRUE)
  )

## Kruskal-wallis
kruskal.test(Evenness ~ sample_type, data = shanplot_k) #not significant, p value = 0.1432
## Wilcoxin rank-sum
pairwise.wilcox.test(shanplot_k$Evenness, shanplot_k$sample_type,
                     p.adjust.method = "BH")










