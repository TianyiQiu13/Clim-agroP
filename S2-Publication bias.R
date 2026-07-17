library(metafor)


setwd("E:/Data/Clim-agroP")
font=theme(axis.title=element_text(size=13),axis.text = element_text(size=12,colour = 'black'),
           strip.text = element_text(size=12),legend.title = element_text(size = 12),
           legend.text = element_text(size = 12),plot.subtitle = element_text(size=13))#11.6inches

a <- read.csv("Effect sizes.csv")

## FSN
fsn(a$rrAP,a$varAP,alpha = 0.05)
fsn(a$rrTP,a$varTP,alpha = 0.05)
fsn(a$rrphosphatase,a$varphosphatase,alpha = 0.05)
fsn(a$rrMBP,a$varMBP,alpha = 0.05)
fsn(a$rrAMF,a$varAMF,alpha = 0.05)
fsn(a$rrshoots_P_uptake,a$varshoots_P_uptake,alpha = 0.05)
fsn(a$rrroots_P_uptake,a$varroots_P_uptake,alpha = 0.05)
fsn(a$rrbiomass,a$varbiomass,alpha = 0.05)
fsn(a$rryield,a$varyield,alpha = 0.05)
