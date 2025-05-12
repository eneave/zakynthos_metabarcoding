#######################
## Map for Zakynthos ##
#######################

setwd("C:/your/directory/map")

library("rstudioapi")

register_google(key = "your-API-key")

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

bigzak2 <- get_map(location = c(lon= 20.785463, lat= 37.781411),
                  zoom = 6,
                  maptype = "roadmap",
                  color = "bw",
                  labels = FALSE)

bigzak3 <- get_map(location = c(lon= 20.785463, lat= 37.781411),
                   zoom = 7,
                   maptype = "roadmap",
                   color = "bw",
                   labels = FALSE)

zak <- get_map(location = c(lon= 20.785463, lat= 37.781411),
        zoom = 11,
        maptype = "satellite",
        color = "color")

head(sites)

# Korakonissi
k <- get_map(location = c(lon= 20.72950, lat= 37.71750),
             zoom = 18,
             maptype = "satellite",
             color = "color")

k2 <- get_map(location = c(lon= 20.72950, lat= 37.71720),
             zoom = 19,
             maptype = "satellite",
             color = "color")

# Dafni Beach
b <- get_map(location = c(lon= 20.95363, lat= 37.72039),
             zoom = 18,
             maptype = "satellite",
             color = "color")

b2 <- get_map(location = c(lon= 20.95300, lat= 37.72010),
             zoom = 19,
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

ggsave(filename="C:/bigmap.jpg", 
       plot = bigmap, width = 5, height = 5, dpi = 300, units = "in")

bigmap2 <-
  ggmap(bigzak3, extent = "panel", legend = "bottomright") +
  geom_rect(xmin = 20.57, xmax = 21,   ymin = 37.63, ymax = 37.95, color = "black",
            fill = NA, size = 2) +
  geom_text(x = 20.5, y= 37.4, label = "Study area",
            size = 4,color ="black") +
  labs(x = "Longitude", y ="Latitude") +
  theme(axis.text.y = element_text(color = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(color = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
ggsave(filename="C:/bigmap2.jpg", 
       plot = bigmap2, width = 5, height = 5, dpi = 300, units = "in")


bigmap3 <-
  ggmap(bigzak3, extent = "panel", legend = "bottomright") +
  geom_rect(xmin = 20.57, xmax = 21,   ymin = 37.63, ymax = 37.95, color = "black",
            fill = NA, size = 2) +
  geom_text(x = 20.5, y= 37.4, label = "Study area",
            size = 10,color ="black") +
  labs(x = "", y ="") +
  theme(axis.text.y = element_text(color = "black", size = 30, face = "bold"), 
        axis.text.x = element_text(color = "black", face = "bold", size = 30), 
        axis.title.y = element_text(face = "bold", size = 16), 
        axis.title.x = element_text(face = "bold", size = 16, color = "black"), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))
ggsave(filename="C:/bigmap3.jpg", 
       plot = bigmap3, width = 5, height = 5, dpi = 300, units = "in")

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

ggsave(filename="C:/mainmap.jpg", 
       plot = mainmap, width = 6.5, height = 6.5, dpi = 300, units = "in")


mainmap2 <-
  ggmap(zak, extent = "panel", legend = "bottomright") +
  geom_point(aes(x = Longitude, y = Latitude), data = sites, size = 4, fill = c("#5D3A9B", "#E66100"),
             pch = 24, color = "#FFFFFF") +
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

ggsave(filename="C:/mainmap2.jpg", 
       plot = mainmap2, width = 5.2, height = 5.2, dpi = 300, units = "in")

## Site maps

ggmap(k)

kmap <-
ggmap(k2) +
  labs(x = "Longitude", y ="Latitude") +
  geom_point(x= 20.72990, y= 37.71730, size = 4, fill = "#5D3A9B",
             pch = 21, color = "#FFFFFF") +
  geom_point(x= 20.72950, y= 37.71725, size = 4, fill = "#5D3A9B",
             pch = 21, color = "#FFFFFF") +
  geom_point(x= 20.72930, y= 37.71750, size = 4, fill = "#5D3A9B",
             pch = 21, color = "#FFFFFF") +
  geom_text(x= 20.72982, y= 37.71730, size = 4.5, label = "4",
            color ="white", fontface = "bold") +
  geom_text(x= 20.72942, y= 37.71725, size = 4.5, label = "5",
            color ="white", fontface = "bold") +
  geom_text(x= 20.72922, y= 37.71750, size = 4.5, label = "6",
            color ="white", fontface = "bold") +
  theme(axis.text.y = element_text(color = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(color = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave(filename="C:/kmap.jpg", 
       plot = kmap, width = 5.2, height = 5.2, dpi = 300, units = "in")

ggmap(b)

bmap <-
ggmap(b2) +
  labs(x = "Longitude", y ="Latitude") +
  geom_point(x= 20.95270, y= 37.72010, size = 4, fill = "#E66100",
             pch = 21, color = "#FFFFFF") +
  geom_point(x= 20.95310, y= 37.719980, size = 4, fill = "#E66100",
             pch = 21, color = "#FFFFFF") +
  geom_point(x= 20.95340, y= 37.71965, size = 4, fill = "#E66100",
             pch = 21, color = "#FFFFFF") +
  geom_text(x= 20.95262, y= 37.72010, size = 4.5, label = "1",
            color ="white", fontface = "bold") +
  geom_text(x= 20.95302, y= 37.719980, size = 4.5, label = "2",
            color ="white", fontface = "bold") +
  geom_text(x= 20.95332, y= 37.71965, size = 4.5, label = "3",
           color ="white", fontface = "bold") +
  theme(axis.text.y = element_text(color = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(color = "black", face = "bold", size = 12), 
        axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, color = "black"), 
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"))

ggsave(filename="C:/bmap.jpg", 
       plot = bmap, width = 5.2, height = 5.2, dpi = 300, units = "in")

