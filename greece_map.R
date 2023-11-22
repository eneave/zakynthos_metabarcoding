####################
## Map for Greece ##
####################

setwd("C:/YourDirectoryHere/map")

## load API key
# youtube tutortial on google API: https://www.youtube.com/watch?v=Of_M4kcE9yM&t=0s

library("rstudioapi")

register_google(key = "YourAPIKeyHere")

library("ggplot2")
library("ggmap")

# Korakonissi 37.72106, 20.72644
# Dafni Beach 37.720392, 20.954633

sites <- read.csv("map.csv")

# center of map 37.781411, 20.785463

## need to check if this is correct
bigzak <- get_map(location = c(lon= 20.785463, lat= 37.781411),
                  zoom = 6,
                  maptype = "satellite",
                  color = "color")


zak <- get_map(location = c(lon= 20.785463, lat= 37.781411),
        zoom = 11,
        maptype = "satellite",
        color = "color")

## Inset map

bigmap <-
ggmap(bigzak, extent = "panel", legend = "bottomright") +
  geom_rect(xmin = 20.57, xmax = 21,   ymin = 37.63, ymax = 37.95, fill = "red", alpha = 0.8) +
  geom_text(x = 18.5, y= 36, label = "Ionian Sea",
            size = 4.5,color ="white") +
  geom_text(x = 23, y= 39, label = "Greece",
            size = 4.5,color ="white") +
  geom_text(x = 20.5, y= 37.4, label = "Study area",
            size = 4,color ="white") +
  labs(x = "Longitude", y ="Latitude") +
  theme(axis.text.y = element_text(color = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(color = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave(filename="C:/YourDirectoryHere/map/bigmap.jpg", 
       plot = bigmap, width = 5, height = 5, dpi = 300, units = "in")



bigmap2 <-
    ggmap(bigzak, extent = "panel", legend = "bottomright") +
    geom_rect(xmin = 20.57, xmax = 21,   ymin = 37.63, ymax = 37.95, fill = "red", alpha = 0.8) +
    geom_text(x = 19.5, y= 36, label = "Ionian Sea",
              size = 4.5,color ="white") +
    geom_text(x = 23, y= 39, label = "Greece",
              size = 5,color ="white") +
    geom_text(x = 21.5, y= 37.2, label = "Study area",
              size = 4,color ="red") +
    labs(x = "", y ="") +
    theme(axis.text.y = element_text(color = "white", size = 16, face = "bold"), 
          axis.text.x = element_text(color = "white", face = "bold", size = 16), 
          axis.line = element_line(color = "white"),
          axis.line.x.top = element_line(color = "white"),
          axis.line.y.right = element_line(color = "white"),
          axis.ticks = element_line(color = "white"),
          panel.background = element_rect(fill = "transparent",
                                          colour = NA_character_), # necessary to avoid drawing panel outline
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          plot.background = element_rect(fill = "transparent",
                                         colour = NA_character_), # necessary to avoid drawing plot outline
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"),
          legend.key = element_rect(fill = "transparent"))
  
ggsave(filename="C:/YourDirectoryHere/map/bigmap_transparent.png", 
         plot = bigmap2, width = 3, height = 3, dpi = 300, units = "in")

  

## Main map

mainmap <-
ggmap(zak, extent = "panel", legend = "bottomright") +
  geom_point(aes(x = Longitude, y = Latitude), data = sites, size = 7, fill = c("#5D3A9B", "#E66100"),
             pch = 21, color = "#FFFFFF") +
  geom_text(data = sites, aes(x = Longitude-0.05, y= Latitude-0.01, label = Location),
            size = 4.5,color ="white",fontface = "bold") +
  scale_x_continuous(limits = c(20.57, 21.0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(37.63, 37.95), expand = c(0, 0)) +
  labs(x = "Longitude", y ="Latitude") +
  theme(axis.text.y = element_text(color = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(color = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave(filename="C:/YourDirectoryHere/map/mainmap.jpg", 
       plot = mainmap, width = 6.5, height = 6.5, dpi = 300, units = "in")








