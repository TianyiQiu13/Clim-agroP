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

setwd("E:/Data/Clim-agroP")
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12),plot.subtitle = element_text(size=13))#11.6inches

a <- read.csv("Effect sizes.csv")
a$pH2 <- ifelse(a$pH<6.5,"acid",ifelse(a$pH<=7.5,"neutral","alkaline"))

prior_fit <- c(
  prior(normal(0, 0.5), class = "b")
)

###Overall###
fit_AP_OC_Fe <- brm(
  formula = rrAP | se(sqrt(varAP)) ~ OC_Fe:pH2 -1,
  data = a,
  seed = 1234,
  prior = prior_fit,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

###Long-term###
fit_AP_OC_Fe_long_term <- brm(
  formula = rrAP | se(sqrt(varAP)) ~ OC_Fe:pH2 -1,
  data = filter(a,duration>=5),
  seed = 1234,
  prior = prior_fit,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

models_fit <- list(
  Overall = fit_AP_OC_Fe,
  Long_term = fit_AP_OC_Fe_long_term
)

posteriors_fit <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, `b_OC_Fe:pH2acid`, `b_OC_Fe:pH2neutral`, `b_OC_Fe:pH2alkaline`) %>%
    pivot_longer(
      cols = starts_with("b_OC_Fe:"),
      names_to = "pH",
      values_to = "b_Intercept"
    ) %>%
    mutate(
      pH2 = str_to_title(str_remove(pH, "b_OC_Fe:pH2")),
      type = .y
    )
}) %>%
  mutate(
    pH3 = factor(pH2, levels = c("Acid", "Neutral", "Alkaline")),
    type = factor(type, levels = c("Overall", "Long_term"),labels = c("Overall","Long-term"))
  )

gg_sum_fit <- posteriors_fit %>%
  group_by(type, pH3) %>%
  median_qi(b_Intercept, .width = 0.90) %>%
  mutate_if(is.numeric, round, 2)

p1 <- ggplot(filter(posteriors_fit,type=="Overall"), aes(x = b_Intercept, y = pH3, fill = pH3, color = pH3)) +
  stat_halfeye(
    .width = c(0.9),
    point_interval = median_qi,
    interval_size = 1.2,
    slab_alpha = 0.5
  ) +
  facet_wrap(~type)+
  scale_x_continuous(limits = c(-0.05,0.035))+
  scale_y_discrete(labels=c("Acid<br>(*n*=25)","Neutral<br>(*n*=48)","Alkaline<br>(*n*=11)"))+
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_text(
    data = mutate_if(filter(gg_sum_fit,type=="Overall"), is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept} ({.lower},{.upper})"), x = Inf),
    hjust = "inward", size=3.8,show.legend = FALSE)+
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3") +
  paletteer::scale_color_paletteer_d("ggsci::category10_d3") +
  theme_cowplot() +
  labs(
    x = expression(paste("Effect sizes of Fe-bound OC on available P changes (regression coefficients)", sep = "")),
    y = NULL
  ) +font+
  theme(legend.position = "none",axis.text.y = element_markdown(lineheight = 1.2))

p2 <- ggplot(filter(posteriors_fit,type=="Long-term"), aes(x = b_Intercept, y = pH3, fill = pH3, color = pH3)) +
  stat_halfeye(
    .width = c(0.9),
    point_interval = median_qi,
    interval_size = 1.2,
    slab_alpha = 0.5
  ) +
  facet_wrap(~type)+
  scale_y_discrete(labels=c("Acid<br>(*n*=13)","Neutral<br>(*n*=22)","Alkaline<br>(*n*=1)"))+
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_text(
    data = mutate_if(filter(gg_sum_fit,type=="Long-term"), is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept} ({.lower},{.upper})"), x = Inf),
    hjust = "inward", size=3.8,show.legend = FALSE)+
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3") +
  paletteer::scale_color_paletteer_d("ggsci::category10_d3") +
  theme_cowplot() +
  labs(
    x = expression(paste("Effect sizes of Fe-bound OC on available P changes (regression coefficients)", sep = "")),
    y = NULL
  ) +font+
  theme(legend.position = "none",axis.text.y = element_markdown(lineheight = 1.2))

plot_grid( p1,p2,
           align = 'hv', 
           nrow = 2,
           ncol=1,labels = "auto",
           label_size = 15,label_x = 0.01
)##7.0*10.2
