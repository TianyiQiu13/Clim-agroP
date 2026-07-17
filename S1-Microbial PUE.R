library(cowplot)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(terra)
library(tidyterra)
library(patchwork)
library(raster)

setwd("E:/Data/Clim-agroP")
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12),plot.subtitle = element_text(size=13))#11.6inches

pH <- rast("E:/Data/Raster/CC-raster/pH-0.1.tif")
maizePUE <- rast("E:/Data/Raster/Crop PUE/Maize_PUE.tif")
maizePUE <- resample(maizePUE,pH,method="bilinear")
microbialPUE <- rast("E:/Data/Raster/Microbial NUE and PUE/PUE.tif")
microbialPUE <- resample(microbialPUE,pH,method="bilinear")

find_low_pue <- function(pue_raster) {
  thresh <- global(pue_raster, fun=function(x) quantile(x, 0.25, na.rm=TRUE))[[1]]
  return(pue_raster < thresh)
}
low_maize <- find_low_pue(maizePUE)
low_microbial <- find_low_pue(microbialPUE)
pH_mask <- pH < 6.5
crisis_hotspots <- pH_mask & low_maize & low_microbial
crisis_map <- mask(maizePUE, crisis_hotspots, maskvalues = 0)

worldMap <- fortify(map_data("world"), region = "subregion")
worldMap <- worldMap[worldMap$region != "Antarctica",]
cropland <- raster("E:/Data/Raster/Cropland area/Cropland2000_5m.tif")
cropland.df <- as.data.frame(cropland,xy=TRUE)
cropland.df$Cropland2000_5m[which(cropland.df$Cropland2000_5m=="0")] <- NA
crisis_map.df <- as.data.frame(crisis_map,xy=TRUE)
crisis_map.df <- crisis_map.df[!is.na(crisis_map.df[[3]]), ]
crisis_map.df$Category <- " "

cropland1 <- rast("E:/Data/Raster/Cropland area/Cropland2000_5m.tif")
cropland1 <- resample(cropland1,pH,method="bilinear")
cell_size_km2 <- cellSize(cropland1, unit = "km") 
total_crop_area <- sum(values(cropland1 * cell_size_km2), na.rm = TRUE)
crisis_crop_area <- sum(values(crisis_hotspots * cropland1 * cell_size_km2), na.rm = TRUE)
ratio <- (crisis_crop_area / total_crop_area) * 100##2.37%

microbialPUE.df <- as.data.frame(microbialPUE,xy=TRUE)

p1 <- ggplot()+
  geom_polygon(data=worldMap,aes(x=long,y=lat,group=group),
               fill="darkgray",color="white",size=0.2)+
  geom_raster(data=microbialPUE.df,aes(x,y,fill=PUE)) +
  scale_fill_gradientn(colors = hcl.colors(9,"RdYlBu",rev = T)[1:9],na.value = NA,name = "Soil microbial PUE")+
  theme_cowplot()+
  scale_x_continuous(limits= c(-175,195))+
  scale_y_continuous(limits= c(-60,85))+font+
  theme(legend.direction = "horizontal",legend.position = c(0.5,0.1),
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())

p2 <- ggplot()+
  geom_polygon(data=worldMap,aes(x=long,y=lat,group=group),
               fill="darkgray",color="white",size=0.2)+
  geom_raster(data=cropland.df,aes(x,y,fill=Cropland2000_5m*100)) +
  scale_fill_gradientn(colors = hcl.colors(9,"YlGn",rev = T)[1:7],na.value = NA,name = "Cropland area (%)")+
  new_scale_fill() +
  geom_tile(data = crisis_map.df, aes(x, y, fill = Category)) +
  scale_fill_manual(values = c(" " = "#c0392b"), name = "Hotspots with P crisis\n(~2.37% of global croplands)") +
  theme_cowplot()+
  scale_x_continuous(limits= c(-175,195))+
  scale_y_continuous(limits= c(-60,85))+font+
  theme(legend.direction = "horizontal",legend.position = c(0.5,0.1),
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())

ggdraw()+ 
  draw_plot(p1, x=0, y=0.5, width = 1, height = 0.5)+
  draw_plot(p2, x=0, y=0, width = 1, height = 0.5)+ 
  draw_plot_label(label = c("a","b"), size = 15,
                  x=c(0.01,0.01),
                  y=c(1,0.5))##11.6*9.4
