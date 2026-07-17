library(brms)
library(tidyverse)
library(bayesplot)
library(cowplot)
library(ggplot2)
library(rworldmap)
library(ggrepel)
library(ggspatial)
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
library(ggdist)
library(HDInterval)
library(performance)
library(metafor)
library(glmulti)
library(metaforest)
library(pdp)

setwd("E:/Data/Clim-agroP")
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12),plot.subtitle = element_text(size=13))#11.6inches

#####===Redox validation===#####
a <- read.csv("Effect sizes.csv")
a$study <- as.factor(a$study)
a$system <- ifelse(a$crop=="rice","paddy","nonpaddy")
a$redox <- interaction(a$system,a$irrigation)

###paddy vs. nonpaddy system###
prior_system <- c(
  prior(normal(0, 1), class = "b",coef = "systempaddy"),
  prior(cauchy(0, 0.5), class = "sd")
)

fit_AP_system <- brm(
  formula = rrAP | se(sqrt(varAP)) ~ -1 + system + (1|study),
  data = a,
  seed = 1234,
  prior = prior_system,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

models_system <- list(
  system = fit_AP_system
)

posteriors_systemnonpaddy <- map2_dfr(models_system, names(models_system), ~ {
  spread_draws(.x, b_systemnonpaddy)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(system="systemnonpaddy")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_systempaddy <- map2_dfr(models_system, names(models_system), ~ {
  spread_draws(.x, b_systempaddy)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(system="systempaddy")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)

posteriors_system <- rbind(posteriors_systemnonpaddy,posteriors_systempaddy)
gg_sum_system <- group_by(as.data.frame(posteriors_system), system) %>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)%>%
  mean_qi(b_Intercept2, .width = .9)
posteriors_system <- posteriors_system %>%
  complete(system, fill = list(b_Intercept = NA))
gg_sum_system <- gg_sum_system%>%mutate(n_sub=c(71,12))

p1 <- ggplot(posteriors_system, aes(x = b_Intercept2, 
                                         y = factor(system,levels = c("systempaddy","systemnonpaddy"),
                                                    labels = c("Paddy<br>(*n*=12)","Non-paddy<br>(*n*=71)")),
                                         fill =  factor(system,levels = c("systempaddy","systemnonpaddy"),
                                                        labels = c("Paddy<br>(*n*=12)","Non-paddy<br>(*n*=71)")),
                                         color = factor(system,levels = c("systempaddy","systemnonpaddy"),
                                                        labels = c("Paddy<br>(*n*=12)","Non-paddy<br>(*n*=71)")))) +
  stat_halfeye(
    .width = c(0.9),
    point_interval = median_qi,
    interval_size = 2.4,
    slab_alpha = 0.5,
    position = position_dodge(width = 0.9)
  ) +
  theme_cowplot()+
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_manual(values=c("#16a085","#c0392b"))+
  scale_color_manual(values=c("#16a085","#c0392b"))+
  geom_text(
    data = mutate_if(gg_sum_system, is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept2} ({.lower},{.upper})"), x = Inf),
    position = position_dodge(width = 0.9),
    hjust = "inward", size=3.8,show.legend = FALSE)+
  labs(
    x = expression(paste("Available P changes compared to ambient (%)",sep = "")),
    y = NULL
  ) +
  font+
  theme(legend.position = "none",
        legend.title = element_blank(),legend.direction = "horizontal",axis.text.y = element_markdown(lineheight = 1.3))

###redox###
prior_redox <- c(
  prior(normal(0, 1), class = "b",coef = "redoxpaddy.conventional")
)

fit_AP_redox <- brm(
  formula = rrAP | se(sqrt(varAP)) ~ -1 + redox,
  data = a,
  seed = 1234,
  prior = prior_redox,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

models_redox <- list(
  redox = fit_AP_redox
)

posteriors_redoxnonpaddy.conventional <- map2_dfr(models_redox, names(models_redox), ~ {
  spread_draws(.x, b_redoxnonpaddy.conventional)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(redox="redoxnonpaddy.conventional")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_redoxpaddy.conventional <- map2_dfr(models_redox, names(models_redox), ~ {
  spread_draws(.x, b_redoxpaddy.conventional)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(redox="redoxpaddy.conventional")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_redoxnonpaddy.rainfed <- map2_dfr(models_redox, names(models_redox), ~ {
  spread_draws(.x, b_redoxnonpaddy.rainfed)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(redox="redoxnonpaddy.rainfed")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_redoxnonpaddy.saving <- map2_dfr(models_redox, names(models_redox), ~ {
  spread_draws(.x, b_redoxnonpaddy.saving)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(redox="redoxnonpaddy.saving")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_redoxpaddy.saving <- map2_dfr(models_redox, names(models_redox), ~ {
  spread_draws(.x, b_redoxpaddy.saving)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(redox="redoxpaddy.saving")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)

posteriors_redox <- rbind(posteriors_redoxnonpaddy.conventional,posteriors_redoxnonpaddy.rainfed,posteriors_redoxnonpaddy.saving,
                          posteriors_redoxpaddy.conventional,posteriors_redoxpaddy.saving)
gg_sum_redox <- group_by(as.data.frame(posteriors_redox), redox) %>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)%>%
  mean_qi(b_Intercept2, .width = .9)%>%
  mutate(system=c(rep("nonpaddy",3),rep("paddy",2)),
         irrigation=c("conventional","rainfed","saving","conventional","saving"))
posteriors_redox <- posteriors_redox %>%
  complete(redox, fill = list(b_Intercept = NA))%>%
  mutate(system=c(rep("nonpaddy",36000),rep("paddy",24000)),
         irrigation=c(rep("conventional",12000),rep("rainfed",12000),rep("saving",12000),
                      rep("conventional",12000),rep("saving",12000)))
  
gg_sum_redox <- gg_sum_redox%>%mutate(n_sub=c(23,43,5,10,2))

p2 <- ggplot(posteriors_redox, aes(x = b_Intercept2, 
                                    y = factor(system,levels = c("paddy","nonpaddy"),
                                               labels = c("Paddy","Non-paddy")),
                                    fill =  factor(irrigation,levels = c("rainfed","conventional","saving"),
                                                   labels = c("Rainfed","Conventional irrigation","Water-saving irrigation")),
                                    color = factor(irrigation,levels = c("rainfed","conventional","saving"),
                                                   labels = c("Rainfed","Conventional irrigation","Water-saving irrigation")))) +
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
    data = mutate_if(gg_sum_redox, is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept2} ({.lower},{.upper})"), x = Inf),
    position = position_dodge(width = 0.9),
    hjust = "inward", size=3.8,show.legend = FALSE)+
  geom_richtext(
    data = mutate_if(gg_sum_redox, is.numeric, round, 2),
    aes(label = str_glue("*n*={n_sub}"), x = -Inf),
    position = position_dodge(width = 0.9),
    hjust = "inward", size=3.8,fill = NA, label.color = NA,show.legend = FALSE)+
  labs(
    x = expression(paste("Available P changes compared to ambient (%)",sep = "")),
    y = NULL
  ) +
  font+
  theme(legend.position = c(0.05,1),
        legend.title = element_blank(),legend.direction = "horizontal")

ggdraw()+ 
  draw_plot(p1, x=0.01, y=0.5, width = 0.98, height = 0.5)+
  draw_plot(p2, x=0.01, y=0, width = 0.98, height = 0.5)+ 
  draw_plot_label(label = c("a","b"), size = 15,
                  x=c(0.01,0.01),
                  y=c(1,0.5))##7.0*10.2

#####===PH importance after accounting for redox===#####
AP2 <- a[,c("rrAP","varAP","study","MAT","MAP",
            "OC_Fe","TP","pH","redox","duration","CO2_magnitude","Temp_magnitude")]
AP2 <- na.omit(AP2)
colnames(AP2) <- c("yi","vi","study","MAT","MAP",
                   "OC_Fe","TP","pH","redox","duration","CO2_magnitude","Temp_magnitude")
AP2[,-c(1,2,3,9)] <- scale(AP2[,-c(1,2,3,9)],scale = T,center = T)

extractRVI <- function(x) {
  ww = exp(-(x@crits - x@crits[1])/2)
  ww = ww/sum(ww)
  clartou = function(x) {
    pieces <- sort(strsplit(x, ":")[[1]])
    if (length(pieces) > 1)
      paste(pieces[1], ":", pieces[2], sep = "")
    else x
  }
  tet = lapply(x@formulas, function(x) sapply(attr(delete.response(terms(x)), "term.labels"), clartou))
  allt <- unique(unlist(tet))
  imp <- sapply(allt, function(x) sum(ww[sapply(tet, function(t) x %in% t)]))
  return(sort(imp))
}

rma.glmulti <- function(formula, data, ...) {
  rma(formula, vi, data=data, method="ML", control=list(stepadj=0.5, optimizer="optimParallel",ncpus=20),...)
}

res <- glmulti(yi ~ MAT + MAP + OC_Fe + TP + pH + redox + duration + 
                     CO2_magnitude + Temp_magnitude, 
                   data=AP2,
                   level=1, fitfunction=rma.glmulti, crit="aicc", confsetsize=128
)

imp <- data.frame(predictor=names(extractRVI(res)), importance=extractRVI(res)) %>%
  mutate(type="Water-regime incorporated")%>%
  mutate(sig=c("Climate","Experiment","Soil","Climate","Experiment","Soil","Management","Soil","Experiment"))

p3 <- ggplot(imp, aes(x = reorder(predictor, importance), y = importance,
                            fill = factor(sig,levels = c("Soil","Climate","Experiment","Management")),
                            color=factor(sig,levels = c("Soil","Climate","Experiment","Management")))) + 
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) + 
  scale_x_discrete(labels=c("MAT",expression(eCO[2]~magnitude),"Total P stocks","MAP",
                            "Warming magnitude","Fe-bound OC","Water regime","pH","Duration"))+
  coord_flip() + 
  facet_wrap(~type)+
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  geom_hline(yintercept = 0.8,linetype="dashed",color="#c0392b")+
  labs(y = "AICc-based importance for available P changes")+
  theme_cowplot() + 
  font+
  theme(axis.title.y = element_blank(),legend.title = element_blank())

###control redox###
rma.glmulti_control <- function(formula, data, ...) {
  form_str <- paste(deparse(formula), collapse = " ")
  new_formula <- as.formula(paste(form_str, "+ redox"))
  rma(new_formula, vi, data=data, method="ML", control=list(stepadj=0.5, optimizer="optimParallel",ncpus=20),...)
}
res_control <- glmulti(yi ~ MAT + MAP + OC_Fe + TP + pH + duration + 
                 CO2_magnitude + Temp_magnitude, 
               data=AP2,
               level=1,fitfunction=rma.glmulti_control, crit="aicc", confsetsize=128
)
imp_control <- data.frame(predictor=names(extractRVI(res_control)), importance=extractRVI(res_control))%>%
  mutate(type="Water-regime controlled")%>%
  mutate(sig=c("Climate","Experiment","Soil","Climate","Experiment","Soil","Soil","Experiment"))

p4 <- ggplot(imp_control, aes(x = reorder(predictor, importance), y = importance,
                      fill = factor(sig,levels = c("Soil","Climate","Experiment")),
                      color=factor(sig,levels = c("Soil","Climate","Experiment")))) + 
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) + 
  scale_x_discrete(labels=c("MAT",expression(eCO[2]~magnitude),"Total P stocks","MAP",
                            "Warming magnitude","Fe-bound OC","pH","Duration"))+
  coord_flip() + 
  facet_wrap(~type)+
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  geom_hline(yintercept = 0.8,linetype="dashed",color="#c0392b")+
  labs(y = "")+
  theme_cowplot() + 
  font+
  theme(axis.title.y = element_blank(),legend.position = "none")

###redox residuals###
base_model <- rma(yi, vi, mods = ~ redox, data = AP2, method = "ML")
AP2$yi_residual <- residuals(base_model)

rma.glmulti_res <- function(formula, data, ...) {
  rma(formula, vi, data=data, method="ML", control=list(stepadj=0.5, optimizer="optimParallel",ncpus=20),...)
}

res_res <- glmulti(yi_residual ~ MAT + MAP + OC_Fe + TP + pH + duration + 
                     CO2_magnitude + Temp_magnitude, 
                   data=AP2,
                   level=1, fitfunction=rma.glmulti_res, crit="aicc", confsetsize=128
)

imp_res <- data.frame(predictor=names(extractRVI(res_res)), importance=extractRVI(res_res)) %>%
  mutate(type="Water-regime residuals")%>%
  mutate(sig=c("Experiment","Climate","Soil","Climate","Experiment","Soil","Soil","Experiment"))

p5 <- ggplot(imp_res, aes(x = reorder(predictor, importance), y = importance,
                              fill = factor(sig,levels = c("Soil","Climate","Experiment")),
                              color=factor(sig,levels = c("Soil","Climate","Experiment")))) + 
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) + 
  scale_x_discrete(labels=c(expression(eCO[2]~magnitude),"MAT","Total P stocks","MAP",
                            "Warming magnitude","Fe-bound OC","pH","Duration"))+
  coord_flip() + 
  facet_wrap(~type)+
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  geom_hline(yintercept = 0.8,linetype="dashed",color="#c0392b")+
  labs(y = "")+
  theme_cowplot() + 
  font+
  theme(axis.title.y = element_blank(),legend.position = "none")

ggdraw()+ 
  draw_plot(p3, x=0.199, y=0.5, width = 0.65, height = 0.5)+
  draw_plot(p4, x=0.01, y=0, width = 0.48, height = 0.5)+ 
  draw_plot(p5, x=0.51, y=0, width = 0.48, height = 0.5)+ 
  draw_text("              AICc-based importance for available P changes",x=0.5,y=0.02,size=13)+
  draw_plot_label(label = c("a","b","c"), size = 15,
                  x=c(0.2,0.01,0.51),
                  y=c(1,0.5,0.5))##11.6*11.0
