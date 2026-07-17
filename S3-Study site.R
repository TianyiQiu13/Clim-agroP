library(ggplot2)
library(rworldmap)
library(scatterpie)
library(ggrepel)
library(ggspatial)
library(ggforce)
library(cowplot)
library(ggalluvial)
library(raster)
library(rmapshaper)
library(sf)

setwd("E:/Data/Clim-agroP")
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12),plot.subtitle = element_text(size=13))#11.6inches

a <- read.csv("Study site-new.csv")

worldMap <- fortify(map_data("world"), region = "subregion")
worldMap <- worldMap[worldMap$region != "Antarctica",]

cropland <- raster("E:/Data/Raster/Cropland area/Cropland2000_5m.tif")
cropland.df <- as.data.frame(cropland,xy=TRUE)
cropland.df$Cropland2000_5m[which(cropland.df$Cropland2000_5m=="0")] <- NA

koppen <- rast("E:/Data/Raster/Koppen/Beck_KG_V1_present_0p5.tif")
reclass_matrix <- matrix(c(
  1, 1,   2, 1,   3, 1, 
  4, 2,   5, 2,   6, 2,   7, 2,
  8, 3,   11, 3,  14, 3,     
  9, 3,   10, 3,  12, 3,  13, 3, 15, 3, 16, 3,
  17, 4,  18, 4,  21, 4,  22, 4,
  19, 4,  20, 4,  23, 4,  24, 4, 
  25, 4,  26, 4,  27, 4,  28, 4, 
  29, 4,  30, 4             
), ncol = 2, byrow = TRUE)
kg_reclass <- classify(koppen, reclass_matrix)
kg_vector <- as.polygons(kg_reclass)
kg_poly <- st_as_sf(kg_vector)


a$methods <- factor(a$methods,levels = c("FACE","OTCs","chamber","infrared heater"),
                    labels = c("FACE","OTC","Chamber","Infrared heater"))
a$type <- factor(a$type,levels = c("CO2","Warming","CO2 + Warming"))

ggplot()+
  geom_polygon(data=worldMap,aes(x=long,y=lat,group=group),
               fill="gray",color="white",size=0.2)+
  geom_sf(data = kg_poly,aes(fill = factor(Beck_KG_V1_present_0p5)), color = "NA",alpha=0.85) +
  scale_fill_manual(values = c("1" = "#1F78B4", "2" = "#EDBB6B", "3" = "#5B8F4C", "4" = "gray"),na.value = "transparent",
                    labels = c("Tropical","Arid","Temperate","Cold"),
                    name = "Climate zone") +
  geom_point(data=a,aes(x=longitude,y=latitude,shape=methods,color=type),size=3.4)+
  theme_cowplot()+
  scale_x_continuous(limits= c(-175,195))+
  scale_y_continuous(limits= c(-60,85))+
  scale_color_manual(values = c("#16a085","#c0392b", "#4198AC"),
                     labels = c(expression(CO[2]),"Warming",expression(CO[2]~+~warming)))+
  xlab(expression(paste("Longitude",sep="")))+
  ylab(expression(paste("Latitude",sep="")))+font+
  guides(fill = guide_legend("Climate classification",order = 1),
         color = guide_legend("Climate treatment",order = 2),
         shape = guide_legend("Method",order = 3))+
  theme(legend.direction = "vertical",legend.position = "bottom",legend.spacing = unit(2,"cm"))#11.6*6.2
