library(brms)
library(tidyverse)
library(bayesplot)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(rstanarm)
library(coda)
library(rstan)
library(broom)
library(DHARMa)
library(tidybayes)
library(ggeffects)
library(broom.mixed)
library(dplyr)
library(HDInterval)
library(performance)
library(metafor)
library(ggpmisc)
library(ggpubr)
library(terra)
library(tidyterra)
library(ggnewscale)
library(ggdist)
library(multcompView)
library(patchwork)
library(raster)
library(ggtext)

setwd("E:/Data/Clim-agroP")
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12),plot.subtitle = element_text(size=13))#11.6inches

#####===Overall effects===#####
a <- read.csv("Effect sizes.csv")
a$study <- as.factor(a$study)
a$pH2 <- ifelse(a$pH<6.5,"acid",ifelse(a$pH<=7.5,"neutral","alkaline"))
a$duration2 <- ifelse(a$duration>=5,"long","short")
a$group <- interaction(a$pH2,a$duration2)

prior_interaction <- c(
  prior(normal(0, 1), class = "b",coef = "groupalkaline.long"),
  prior(cauchy(0, 0.5), class = "sd")
)

fit_AP_interaction <- brm(
  formula = rrAP | se(sqrt(varAP)) ~ -1 + group + (1|study),
  data = a,
  seed = 1234,
  prior = prior_interaction,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

models_interaction <- list(
  AP = fit_AP_interaction
)
posteriors_interactiongroupacid.long <- map2_dfr(models_interaction, names(models_interaction), ~ {
  spread_draws(.x, b_groupacid.long)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(interaction="groupacid.long")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_interactiongroupneutral.long <- map2_dfr(models_interaction, names(models_interaction), ~ {
  spread_draws(.x, b_groupneutral.long)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(interaction="groupneutral.long")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_interactiongroupalkaline.long <- map2_dfr(models_interaction, names(models_interaction), ~ {
  spread_draws(.x, b_groupalkaline.long)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(interaction="groupalkaline.long")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_interactiongroupacid.short <- map2_dfr(models_interaction, names(models_interaction), ~ {
  spread_draws(.x, b_groupacid.short)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(interaction="groupacid.short")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_interactiongroupneutral.short <- map2_dfr(models_interaction, names(models_interaction), ~ {
  spread_draws(.x, b_groupneutral.short)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(interaction="groupneutral.short")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_interactiongroupalkaline.short <- map2_dfr(models_interaction, names(models_interaction), ~ {
  spread_draws(.x, b_groupalkaline.short)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(interaction="groupalkaline.short")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)

posteriors_interaction <- rbind(posteriors_interactiongroupacid.long,posteriors_interactiongroupneutral.long,
                                posteriors_interactiongroupalkaline.long,posteriors_interactiongroupacid.short,
                                posteriors_interactiongroupneutral.short,posteriors_interactiongroupalkaline.short)

gg_sum_interaction <- group_by(as.data.frame(posteriors_interaction), interaction) %>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)%>%
  mean_qi(b_Intercept2, .width = .9)%>%
  mutate(pH2=c(rep("Acid",2),rep("Alkaline",2),rep("Neutral",2)),
         duration2=c("Long-term<br>(≥ 5 years)","Short-term<br>(< 5 years)","Long-term<br>(≥ 5 years)",
                     "Short-term<br>(< 5 years)","Long-term<br>(≥ 5 years)","Short-term<br>(< 5 years)"))
posteriors_interaction <- posteriors_interaction %>%
  complete(interaction, fill = list(b_Intercept = NA))%>%
  mutate(pH2=c(rep("Acid",24000),rep("Alkaline",24000),rep("Neutral",24000)),
         duration2=c(rep("Long-term<br>(≥ 5 years)",12000),rep("Short-term<br>(< 5 years)",12000),
                     rep("Long-term<br>(≥ 5 years)",12000),rep("Short-term<br>(< 5 years)",12000),
                     rep("Long-term<br>(≥ 5 years)",12000),rep("Short-term<br>(< 5 years)",12000)))
gg_sum_interaction <- gg_sum_interaction%>%mutate(n_sub=c(13,11,1,10,22,26))

p1 <- ggplot(posteriors_interaction, aes(x = b_Intercept2, 
                                 y = factor(duration2,levels = c("Short-term<br>(< 5 years)","Long-term<br>(≥ 5 years)")),
                                 fill =  factor(pH2,levels = c("Acid","Neutral","Alkaline")),
                                 color = factor(pH2,levels = c("Acid","Neutral","Alkaline")))) +
  stat_halfeye(
    .width = c(0.9),
    point_interval = median_qi,
    interval_size = 2.4,
    slab_alpha = 0.5,
    position = position_dodge(width = 0.9)
  ) +
  theme_cowplot()+
  geom_vline(xintercept = 0, linetype = "dotted") +
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  geom_text(
    data = mutate_if(gg_sum_interaction, is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept2} ({.lower},{.upper})"), x = Inf),
    position = position_dodge(width = 0.9),
    hjust = "inward", size=3.8,show.legend = FALSE)+
  geom_richtext(
    data = mutate_if(gg_sum_interaction, is.numeric, round, 2),
    aes(label = str_glue("*n*={n_sub}"), x = -Inf),
    position = position_dodge(width = 0.9),
    hjust = "inward", size=3.8,fill = NA, label.color = NA,show.legend = FALSE)+
  labs(
    x = expression(paste("Available P changes compared to ambient (%)",sep = "")),
    y = NULL
  ) +
  font+
  theme(legend.position = c(0.44,1),
        legend.title = element_blank(),legend.direction = "horizontal",axis.text.y = element_markdown(lineheight = 1.1))

#####===Relationship between pH and Fe-associated OC===#####
a$pH2 <- factor(a$pH2,levels = c("acid","neutral","alkaline"),labels = c("Acid","Neutral","Alkaline"))
model <- aov(OC_Fe ~ pH2, data = a)
hsd <- TukeyHSD(model)
tukey_letters <- multcompLetters4(model, hsd)$pH2$Letters
label_data <- data.frame(
  pH2 = names(tukey_letters),
  Letter = as.character(tukey_letters),
  y_pos = c(9.5,7,4.5),
  n_sub = c(94,131,78)
)

p2_main <-ggplot(a, aes(x = pH2, y = OC_Fe, fill = pH2)) +
  stat_halfeye(
    adjust = 0.5, 
    width = 0.6, 
    .width = 0,
    slab_alpha = 0.5,
    justification = -0.3, 
    point_colour = NA
  )+
  geom_boxplot(
    aes(color = pH2),
    width = 0.15, 
    outlier.shape = NA, 
    alpha = 0.5
  ) +
  geom_richtext(
    data = mutate_if(label_data, is.numeric, round, 1),
    aes(label = str_glue("{Letter}<br>*n*={n_sub}"), y = -Inf),lineheight=1.5,
    vjust = "inward", size=3.8,fill = NA, label.color = NA,show.legend = FALSE)+
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  scale_y_continuous(limits = c(-1,10))+
  theme_cowplot() +
  font+
  labs(x = "Soil pH", y = "Fe-bound OC concentration (g/kg)")+
  theme(legend.position = "none")

p2_inset <- ggplot(a, aes(x = pH, y = OC_Fe)) +
  geom_point(alpha = 0.3, size = 1.2, color = "gray40") +
  geom_smooth(method = "lm", color = "black", fill = "gray80", size = 0.8) +
  scale_y_continuous(limits = c(0,10.5))+
  theme_bw(base_size = 9) + 
  stat_poly_eq(aes(pH,OC_Fe,
                   label = paste(..p.value.label..,..n.label..,sep="*\", \"*")),
               formula = y ~ poly(x, 1, raw = TRUE), parse = T,label.x = "right",label.y = "top",size=3.4)+
  theme(panel.grid = element_blank(),
        background = element_rect(fill = "transparent"),axis.text = element_text(size = 10,colour = "black"),
        axis.title = element_text(size = 10)) +
  labs(x = NULL, y = NULL)

p2 <- ggdraw(p2_main) +
  draw_plot(p2_inset, x = 0.48, y = 0.58, width = 0.45, height = 0.4)

#####===Phosphorus crisis===#####
pH <- rast("E:/Data/Raster/CC-raster/pH-0.1.tif")
maizePUE <- rast("E:/Data/Raster/Crop PUE/Maize_PUE.tif")
maizePUE <- resample(maizePUE,pH,method="bilinear")

find_low_pue <- function(pue_raster) {
  thresh <- global(pue_raster, fun=function(x) quantile(x, 0.25, na.rm=TRUE))[[1]]
  return(pue_raster < thresh)
}
low_maize <- find_low_pue(maizePUE)
pH_mask <- pH < 6.5
crisis_hotspots <- pH_mask & low_maize
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
ratio <- (crisis_crop_area / total_crop_area) * 100##12.92%

wang_sites.df <- data.frame(lon=120.55,lat=31.5)

p3 <- ggplot()+
  geom_polygon(data=worldMap,aes(x=long,y=lat,group=group),
               fill="darkgray",color="white",size=0.2)+
  geom_raster(data=cropland.df,aes(x,y,fill=Cropland2000_5m*100)) +
  scale_fill_gradientn(colors = hcl.colors(9,"YlGn",rev = T)[1:7],na.value = NA,name = "Cropland area (%)")+
  new_scale_fill() +
  geom_tile(data = crisis_map.df, aes(x, y, fill = Category)) +
  scale_fill_manual(values = c(" " = "#c0392b"), name = "Hotspots with P crisis\n(~12.9% of global croplands)") +
  theme_cowplot()+
  scale_x_continuous(limits= c(-175,195))+
  scale_y_continuous(limits= c(-60,85))+font+
  theme(legend.direction = "horizontal",legend.position = c(0.5,0.1),
        axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())

ggdraw()+ 
  draw_plot(p1, x=0.01, y=0.52, width = 0.49, height = 0.475)+
  draw_plot(p2, x=0.55, y=0.52, width = 0.44, height = 0.475)+ 
  draw_plot(p3, x=0, y=0, width = 1, height = 0.52)+ 
  draw_plot_label(label = c("a","b","c"), size = 15,
                  x=c(0.01,0.51,0.01),
                  y=c(1,1,0.52))##11.6*9.4
