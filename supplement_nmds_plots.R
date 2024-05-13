##################################################
## Final NMDS and stress plots for Greece paper ##
##################################################

## Stessplot for jaccard disstance NMDS for UVC

stressplot(nmds.jaccard.ord.2)

## Stessplot for bray/jaccard disstance NMDS for eDNA before/after outlier removal

stressplot(nmds.jaccard.ord.edna)
stressplot(nmds.bray.ord.edna)

stressplot(nmds.jaccard.ord.edna.nooutlier)
stressplot(nmds.bray.ord.nooutlier)

## Need to show the eDNA plots before and after outlier removal,
## for both Jaccard and Bray-curtis distance

# eDNA All samples Bray-curtis
data.scores.edna.orig.b <- as.data.frame(scores(nmds.bray.ord.edna, display ="sites"))  
data.scores.edna.orig.b$id_seq <- rownames(data.scores.edna.orig.b)
nmds_edna_b <- merge(data.scores.edna.orig.b, mdata, by.x = "id_seq")

edna_nmds_b <-
  ggplot() + 
  geom_point(data = nmds_edna_b, aes(x = NMDS1, y = NMDS2, fill = location), 
             colour = c("black"), shape = 21, alpha=0.5, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  #geom_text(data = nmds_edna, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 4500, y= -2000, label = "stress = 4.66e-05",size = 5,colour ="black")+
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
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) 

# eDNA All samples Jaccard
data.scores.edna.orig.j <- as.data.frame(scores(nmds.jaccard.ord.edna, display ="sites"))  
data.scores.edna.orig.j$id_seq <- rownames(data.scores.edna.orig.j)
nmds_edna_j <- merge(data.scores.edna.orig.j, mdata, by.x = "id_seq")

edna_nmds_j <-
  ggplot() + 
  geom_point(data = nmds_edna_j, aes(x = NMDS1, y = NMDS2, fill = location), 
             colour = c("black"), shape = 21, alpha=0.5, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  #geom_text(data = nmds_edna, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 1750, y= -0.5, label = "stress = 9.02e-05",size = 5,colour ="black")+
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
  scale_fill_manual(values = c("#E66100", "#5D3A9B")) 

# eDNA No Outliers Bray-curtis

data.scores.edna.noout.b <- as.data.frame(scores(nmds.bray.ord.nooutlier, display ="sites"))  
data.scores.edna.noout.b$id_seq <- rownames(data.scores.edna.noout.b)
nmds_edna_noout_b <- merge(data.scores.edna.noout.b, mdata, by.x = "id_seq")

edna_nmds_noout_b <-
  ggplot() + 
  geom_point(data = nmds_edna_noout_b, aes(x = NMDS1, y = NMDS2, fill = location), 
             colour = c("black"), shape = 21, alpha=0.5, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  #geom_text(data = nmds_edna, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 0.5, y= -0.5, label = "stress = 0.0693",size = 5,colour ="black")+
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
  scale_fill_manual(values = c("#E66100", "#5D3A9B"))

# eDNA No Outliers Jaccard

edna_nmds_noout_j <-
  ggplot() + 
  geom_point(data = nmds_edna, aes(x = NMDS1, y = NMDS2, fill = location), 
             colour = c("black"), shape = 21, alpha=0.5, size = 10)+ 
  #geom_richtext(data = fig2a_legend, aes(x=x, y=y, label=label), colour = c("#009E73", "#CC79A7", "#D55E00"), size = 6.5) +
  #geom_text(data = nmds_edna, aes(x = NMDS1, y= NMDS2, label = Richness),size = 4.5,colour ="white")+
  annotate(geom = "text", x = 0.5, y= -1.5, label = "stress = 0.0443",size = 5,colour ="black")+
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
  scale_fill_manual(values = c("#E66100", "#5D3A9B"))


library(cowplot)

supp_nmds <- plot_grid(edna_nmds_j, edna_nmds_b, edna_nmds_noout_j, edna_nmds_noout_b,
                  rel_widths = c(1,1,1,1), 
                  labels = c('A','B','C','D'), 
                  label_size = 20,
                  nrow = 2)







