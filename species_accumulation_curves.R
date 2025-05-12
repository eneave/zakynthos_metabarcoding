########################################################
## Species accumulation curve with Greece data
########################################################

setwd("C:/your/directory/beta_diversity")

library(BiodiversityR)

#####
## Prepare data for specie accumulation curves
#####

## dataframe for x component of accumcomp
samp_edna_abund <- data.frame(t(samp_edna_2))
colnames(samp_edna_abund) <- samp_edna1$scientific_name_final

## dataframe for y component of accumcomp
samp_edna_abund2 <- samp_edna_abund
samp_edna_abund2$id_seq <- rownames(samp_edna_abund)
# Remove Parablennius yatabei
samp_edna_abund2 <- samp_edna_abund2[-c(23)]
meta_abund <- merge(samp_edna_abund2, mdata, by.x = "id_seq")

accum.ourstudy <- accumcomp(samp_edna_abund, y = meta_abund, factor = "location", plotit = F)

#accumresult(samp_edna_abund, y = meta_abund, factor = "location", method="rarefaction", permutations=100,
#            plotit = T)

#####
## Add Aglieri data to specie accumulation curves
#####

## https://stackoverflow.com/questions/48771601/how-to-plot-multiple-species-accumulation-curves-in-one-plot-using-r

## might be worth trying this method as an alternative to BiodiversityR

ag2020 <- read.csv("Aglieri_2020_zak_edna.csv")


## dataframe for x component of accumcomp
samp_ag_abund <- data.frame(t(ag2020[5:16]))
colnames(samp_ag_abund) <- ag2020$scientific_name

## dataframe for y component of accumcomp
samp_ag_abund2 <- samp_ag_abund
samp_ag_abund2$id_seq <- rownames(samp_ag_abund)

agmdata <- read.csv("Aglieri_2020_mdata.csv")
meta_abund_ag <- merge(samp_ag_abund2, agmdata, by.x = "id_seq")
as.character(meta_abund_ag$mpa)

accum.agstudy <- accumcomp(samp_ag_abund, y = meta_abund_ag, factor = "mpa", plotit = F)

#####
## Plot Aglieri data and our data together
#####

accum.ourstudy.long <- accumcomp.long(accum.ourstudy, ci=NA, label.freq=2)
accum.ourstudy.long$uperror <- accum.ourstudy.long$Richness + accum.ourstudy.long$SD
accum.ourstudy.long$lowerror <- accum.ourstudy.long$Richness - accum.ourstudy.long$SD


accum.agstudy.long <- accumcomp.long(accum.agstudy, ci=NA, label.freq=2)
accum.agstudy.long$uperror <- accum.agstudy.long$Richness + accum.agstudy.long$SD
accum.agstudy.long$lowerror <- accum.agstudy.long$Richness - accum.agstudy.long$SD

accum.data <- rbind(accum.ourstudy.long, accum.agstudy.long)


#####
## Code for plotting from North Atlantic paper
#####

library(ggplot2)

accumcurve <-
  ggplot(data=accum.data, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), linewidth=2) +
  geom_point( aes(colour=Grouping), size=4, alpha=0.9) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, alpha=0.2)), ymax = uperror, ymin = lowerror), 
              show.legend=FALSE) + 
  ylim(0,35) +
  xlim(0.9,8) +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = c(.97, .05),
        legend.justification = c("right", "bottom"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.key=element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2)) +
  scale_color_manual(values = c("MPA" = "#FFB000",
                                "Dafni Beach" = "#FE6100",
                                "outside_MPA" = "#648FFF",
                                "Korakonissi" = "#785EF0"),
                     labels = c("Dafni Beach - NMPZ",
                                "Korakonissi - Outside NMPZ",
                                "Aglieri et al. - NMPZ",
                                "Aglieri et al. - Outside NMPZ")) +
  labs(x = "No. of eDNA samples", y = "No. of Fish Species", colour = "Location")

ggsave(filename="C:/accumcurve.jpg", 
       plot = accumcurve, width = 10, height = 5, dpi = 300, units = "in")

