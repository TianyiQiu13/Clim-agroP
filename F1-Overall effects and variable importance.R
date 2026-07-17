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
library(ggtext)

setwd("E:/Data/Clim-agroP")
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12),plot.subtitle = element_text(size=13))#11.6inches

#####===Overall effects===#####
a <- read.csv("Effect sizes.csv")
a$study <- as.factor(a$study)

prior_conditions <- c(
  prior(normal(0, 1), class = "b",coef = "conditionsCO2"),
  prior(normal(0, 1), class = "b",coef = "conditionsTemp"),
  prior(cauchy(0, 0.5), class = "sd")
)

fit_AP_conditions <- brm(
  formula = rrAP | se(sqrt(varAP)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_TP_conditions <- brm(
  formula = rrTP | se(sqrt(varTP)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_phosphatase_conditions <- brm(
  formula = rrphosphatase | se(sqrt(varphosphatase)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_MBP_conditions <- brm(
  formula = rrMBP | se(sqrt(varMBP)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_AMF_conditions <- brm(
  formula = rrAMF | se(sqrt(varAMF)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_shoots_P_uptake_conditions <- brm(
  formula = rrshoots_P_uptake | se(sqrt(varshoots_P_uptake)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_roots_P_uptake_conditions <- brm(
  formula = rrroots_P_uptake | se(sqrt(varroots_P_uptake)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_biomass_conditions <- brm(
  formula = rrbiomass | se(sqrt(varbiomass)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

fit_yield_conditions <- brm(
  formula = rryield | se(sqrt(varyield)) ~ -1 + conditions + (1|study),
  data = a,
  seed = 1234,
  prior = prior_conditions,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

models_conditions <- list(
  AP = fit_AP_conditions,
  TP = fit_TP_conditions,
  phosphatase = fit_phosphatase_conditions,
  MBP = fit_MBP_conditions,
  AMF = fit_AMF_conditions,
  shoots_P_uptake = fit_shoots_P_uptake_conditions,
  roots_P_uptake = fit_roots_P_uptake_conditions,
  biomass = fit_biomass_conditions,
  yield = fit_yield_conditions
)

posteriors_conditionsCO2 <- map2_dfr(models_conditions, names(models_conditions), ~ {
  spread_draws(.x, b_conditionsCO2) %>%
    mutate(Response = .y)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept","Response"))%>%mutate(conditions="CO2")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_conditionsTemp <- map2_dfr(models_conditions, names(models_conditions), ~ {
  spread_draws(.x, b_conditionsTemp) %>%
    mutate(Response = .y)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept","Response"))%>%mutate(conditions="Temp")%>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)
posteriors_conditions <- rbind(posteriors_conditionsCO2,posteriors_conditionsTemp)

level_order <- c("AP","TP","phosphatase","MBP","AMF","shoots_P_uptake","roots_P_uptake","biomass","yield")
vector_labels <- c("Available P (*n*=83)","Total P (*n*=44)","Phosphatase (*n*=17)",
                   "Microbial biomass P (*n*=18)","AMF (*n*=35)","Shoot P uptake (*n*=169)","Root P uptake (*n*=27)",
                   "Shoot biomass (*n*=81)","Yield (*n*=55)")

gg_sum_conditions <- group_by(as.data.frame(posteriors_conditions), conditions,Response) %>%
  mutate(b_Intercept2 = (exp(b_Intercept)-1)*100)%>%
  mean_qi(b_Intercept2, .width = .9)
posteriors_conditions <- posteriors_conditions %>%
  complete(Response = level_order, conditions, fill = list(b_Intercept = NA))
posteriors_conditions$conditions <- factor(posteriors_conditions$conditions,levels = c("CO2","Temp"),
                                           labels = c("CO2","Warming"))
gg_sum_conditions$conditions <- factor(gg_sum_conditions$conditions,levels = c("CO2","Temp"),
                                       labels = c("CO2","Warming"))
gg_sum_conditions <- gg_sum_conditions%>%mutate(n_sub=c(20,63,65,14,11,21,126,39,45,
                                                        15,20,16,4,6,6,43,5,10))

p1 <- ggplot(posteriors_conditions, aes(x = b_Intercept2, 
                                  y = factor(Response,levels = rev(level_order),labels = rev(vector_labels)),
                                  fill = conditions, color = conditions)) +
  stat_halfeye(
    .width = c(0.9),
    point_interval = median_qi,
    interval_size = 1.2,
    slab_alpha = 0.5,
    position = position_dodge(width = 0.6)
  ) +
  theme_cowplot()+
  geom_vline(xintercept = 0, linetype = "dotted") +
  scale_fill_manual(values=c("#16a085","#c0392b"),labels=c(expression(eCO[2]),"Warming"))+
  scale_color_manual(values=c("#16a085","#c0392b"),labels=c(expression(eCO[2]),"Warming"))+
  geom_text(
    data = mutate_if(gg_sum_conditions, is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept2} ({.lower},{.upper})"), x = Inf, color = conditions),
    position = position_dodge(width = 0.6),
    hjust = "inward", size=3.8,show.legend = FALSE)+
  geom_richtext(
    data = mutate_if(gg_sum_conditions, is.numeric, round, 0),
    aes(label = str_glue("*n*={n_sub}"), x = -Inf, color = conditions),
    position = position_dodge(width = 0.6),
    hjust = "inward", size=3.8,fill = NA, label.color = NA,show.legend = FALSE)+
  labs(
    x = expression(paste("Relative changes compared to ambient (%)",sep = "")),
    y = NULL
  ) +
  font+
  guides(fill=guide_legend("conditions"),
         color=guide_legend("conditions"))+
  theme(legend.position = c(0.05,1),
        legend.title = element_blank(),legend.direction = "horizontal",axis.text.y = element_markdown())

#####===Variable importance===#####
AP2 <- a[,c("rrAP","varAP","study","MAT","MAP",
               "OC_Fe","TP","pH","duration","CO2_magnitude","Temp_magnitude")]
AP2 <- na.omit(AP2)
colnames(AP2) <- c("yi","vi","study","MAT","MAP",
                   "OC_Fe","TP","pH","duration","CO2_magnitude","Temp_magnitude")
AP2[,-c(1,2,3)] <- scale(AP2[,-c(1,2,3)],scale = T,center = T)

###Glmulti###
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
res <- glmulti(yi ~ MAT + MAP + OC_Fe + TP + pH + duration + 
                 CO2_magnitude + Temp_magnitude, 
               data=AP2,
               level=1,fitfunction=rma.glmulti, crit="aicc", confsetsize=128
)
imp <- data.frame(predictor=names(extractRVI(res)), importance=extractRVI(res))%>%
  mutate(type="Glmulti")%>%mutate(sig=c(rep("Climate",2),rep("Experiment",2),rep("Soil",3),"Experiment"))
imp$predictor <- factor(imp$predictor,levels = c("pH","TP","OC_Fe","MAT","MAP",
                                      "duration","CO2_magnitude","Temp_magnitude"),
                           labels = c("pH","Total P stocks","Fe-bound OC","MAT","MAP",
                                      "Duration","CO2 magnitude","Warming magnitude"))

p2 <- ggplot(imp, aes(x = reorder(predictor, importance), y = importance,
                   fill = factor(sig,levels = c("Soil","Climate","Experiment")),
                   color=factor(sig,levels = c("Soil","Climate","Experiment")))) + 
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) + 
  scale_x_discrete(labels=c("MAP","MAT","Warming magnitude",expression(eCO[2]~magnitude),
                            "Total P stocks","Fe-bound OC","pH","Duration"))+
  coord_flip() + 
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  geom_hline(yintercept = 0.8,linetype="dashed",color="#c0392b")+
  labs(y = "AICc-based importance for available P changes")+
  theme_cowplot() + 
  font+
  theme(axis.title.y = element_blank(),legend.position = "none")

###MetaForest###
set.seed(36326)
MF <- MetaForest(yi ~ MAT + MAP + OC_Fe + TP + pH + duration +
                   CO2_magnitude + Temp_magnitude,
                 data = AP2,vi = "vi",
                 whichweights = "random",
                 num.trees = 10000)
imp_MF <- data.frame(MF$forest$variable.importance)%>%mutate(type="MetaForest")%>%
  `colnames<-`(c("importance","type"))%>%mutate(sig=c(rep("Climate",2),rep("Soil",3),rep("Experiment",3)))
imp_MF$predictor <- rownames(imp_MF)%>%factor(levels = c("pH","TP","OC_Fe","MAT","MAP",
                                                          "duration","CO2_magnitude","Temp_magnitude"),
                                              labels = c("pH","Total P stocks","Fe-bound OC","MAT","MAP",
                                                         "Duration","CO2 magnitude","Warming magnitude"))

p3 <- ggplot(imp_MF, aes(x = reorder(predictor, importance), y = importance,
                   fill = factor(sig,levels = c("Soil","Climate","Experiment")),
                   color=factor(sig,levels = c("Soil","Climate","Experiment")))) + 
  geom_bar(stat = "identity", width = 0.7, alpha = 0.5) + 
  scale_x_discrete(labels=c("Warming magnitude","MAT","Fe-bound OC","MAP",
                            "Total P stocks","pH",expression(eCO[2]~magnitude),"Duration"))+
  coord_flip() + 
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  labs(y = "MetaForest-based importance for available P changes                ")+
  theme_cowplot() + 
  font+
  theme(axis.title.y = element_blank(),legend.position = c(0.65,0.15), 
        legend.title = element_blank())

###Bayesian regression###
prior_fit <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(cauchy(0, 0.5), class = "sd")
)

fit_AP_fit <- brm(
  formula = yi | se(sqrt(vi)) ~ MAT + MAP + OC_Fe + TP + pH + duration +
    CO2_magnitude + Temp_magnitude +(1|study),
  data = AP2,
  seed = 1234,
  prior = prior_fit,
  cores = 4,
  chains = 4,
  iter = 6000,
  control = list(stepsize=0.01, adapt_delta = 0.99, max_treedepth = 15)
)

models_fit <- list(
  AP = fit_AP_fit
)
posteriors_fitMAT <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_MAT)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="MAT")
posteriors_fitMAP <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_MAP)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="MAP")
posteriors_fitOC_Fe <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_OC_Fe)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="OC_Fe")
posteriors_fitTP <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_TP)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="TP")
posteriors_fitpH <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_pH)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="pH")
posteriors_fitduration <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_duration)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="duration")
posteriors_fitCO2_magnitude <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_CO2_magnitude)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="CO2_magnitude")
posteriors_fitTemp_magnitude <- map2_dfr(models_fit, names(models_fit), ~ {
  spread_draws(.x, b_Temp_magnitude)
})%>%`colnames<-`(c(".chain",".iteration","draw","b_Intercept"))%>%mutate(fit="Temp_magnitude")

posteriors_fit <- rbind(posteriors_fitMAT,posteriors_fitMAP,posteriors_fitOC_Fe,
                        posteriors_fitTP,posteriors_fitpH,
                        posteriors_fitduration,
                        posteriors_fitCO2_magnitude,posteriors_fitTemp_magnitude)

gg_sum_fit <- group_by(as.data.frame(posteriors_fit), fit) %>%
  mean_qi(b_Intercept, .width = .9)%>%
  mutate(sig=c(rep("Experiment",2),rep("Climate",2),rep("Soil",2),"Experiment","Soil"))
posteriors_fit <- posteriors_fit %>%
  complete(fit, fill = list(b_Intercept = NA))%>%
  mutate(sig=c(rep("Experiment",24000),rep("Climate",24000),rep("Soil",24000),
               rep("Experiment",12000),rep("Soil",12000)))
posteriors_fit$fit <- factor(posteriors_fit$fit,levels = c("pH","TP","OC_Fe","MAT","MAP",
                                                           "duration","CO2_magnitude","Temp_magnitude"),
                             labels = c("pH","Total P stocks","Fe-bound OC","MAT","MAP",
                                        "Duration","CO2 magnitude","Warming magnitude"))
gg_sum_fit$fit <- factor(gg_sum_fit$fit,levels = c("pH","TP","OC_Fe","MAT","MAP",
                                                   "duration","CO2_magnitude","Temp_magnitude"),
                         labels = c("pH","Total P stocks","Fe-bound OC","MAT","MAP",
                                    "Duration","CO2 magnitude","Warming magnitude"))

p4 <- ggplot(posteriors_fit, aes(x = b_Intercept, 
                           y = factor(fit,levels = rev(c("pH","Total P stocks","Fe-bound OC","MAT","MAP",
                                                         "Duration","CO2 magnitude","Warming magnitude"))),
                           fill = factor(sig,levels = c("Soil","Climate","Experiment")), 
                           color = factor(sig,levels = c("Soil","Climate","Experiment")))) +
  stat_halfeye(
    .width = c(0.9),
    point_interval = median_qi,
    interval_size = 2.4,
    slab_alpha = 0.5
  ) +
  scale_y_discrete(labels=rev(c("pH","Total P stocks","Fe-bound OC","MAT","MAP",
                            "Duration",expression(eCO[2]~magnitude),"Warming magnitude")))+
  theme_cowplot()+
  geom_vline(xintercept = 0, linetype = "dotted") +
  paletteer::scale_fill_paletteer_d("ggsci::category10_d3")+
  paletteer::scale_color_paletteer_d("ggsci::category10_d3")+
  geom_text(
    data = mutate_if(gg_sum_fit, is.numeric, round, 2),
    aes(label = str_glue("{b_Intercept} ({.lower},{.upper})"), x = Inf),
    hjust = "inward", size=3.8,show.legend = FALSE)+
  labs(
    x = expression(paste("Regression coefficients with available P changes        ",sep = "")),
    y = NULL
  ) +
  font+
  theme(legend.position = c(0.05,1),
        legend.title = element_blank(),legend.direction = "horizontal")

plot_grid( p1,p4,p2,p3,
           align = 'hv', 
           nrow = 2,
           ncol=2,labels = "auto",
           label_size = 15,label_x = 0.01
)###11.6*11.0
