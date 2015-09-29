library(ggplot2)
library(plyr)
library(vegan)
library(lme4)
library(grid)
library(nlme)


setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank())

###################################################################################################
###################################################################################################
###################################################################################################

#read in experiment info
info <- read.csv('exp_info072015.csv')

#read in community dataset (multivariate and basic diversity stats data)
community <- read.csv('dispersion_means_press_bayesian_09092015.csv')
names(community)[names(community)=='ANPP'] <- 'ANPP_sitelevel'

#read in anpp data
anpp <- read.csv('AllAnppData_02172015.csv')

###################################################################################################
#calculate mean anpp for each treatment and trt year at each site
anppMean <- ddply(anpp, c('site_code', 'project_name', 'community_type', 'treatment_year', 'treatment'), summarise,
                  ANPP_trtmean=mean(anpp),
                  ANPP_trtvar=var(anpp))
names(anppMean)[names(anppMean)=='treatment_year'] <- 'trt.year'

#merge anpp and community data
anppComm <- merge(anppMean, community, by=c('site_code', 'project_name', 'community_type', 'trt.year', 'treatment'))

#calculate response ratios
trt <- subset(anppComm, subset=(plot_mani>0))
ctl <- subset(anppComm, subset=(plot_mani==0))
keep <- c('expt', 'trt.year', 'dispersion', 'S_mean', 'J_mean', 'ANPP_trtmean', 'ANPP_trtvar')
ctl <- ctl[,colnames(ctl) %in% keep]
names(ctl)[names(ctl)=='dispersion'] <- 'dispersion_ctl'
names(ctl)[names(ctl)=='S_mean'] <- 'S_ctl'
names(ctl)[names(ctl)=='J_mean'] <- 'J_ctl'
names(ctl)[names(ctl)=='ANPP_trtmean'] <- 'ANPP_trtmean_ctl'
names(ctl)[names(ctl)=='ANPP_trtvar'] <- 'ANPP_trtvar_ctl'
RR <- merge(trt, ctl, by=c('expt', 'trt.year'))
#note that mean_change is already a response ratio of sorts
RR$dispersion_lnRR <- with(RR, log(dispersion/dispersion_ctl))
RR$S_lnRR <- with(RR, log(S_mean/S_ctl))
RR$J_lnRR <- with(RR, log(J_mean/J_ctl))
RR$ANPP_lnRR <- with(RR, log(ANPP_trtmean/ANPP_trtmean_ctl))
RR$varANPP_lnRR <- with(RR, log(ANPP_trtvar/ANPP_trtvar_ctl))

#get final year of data for each experiment
finalYear <- RR[with(RR, ave(trt.year, label, FUN=max)==trt.year),]

#mixed effects models for mean ANPP change
#these models are not correct, need more work to properly account for random effects of site and experiment
summary(richnessModel <- lmer(ANPP_lnRR ~ S_lnRR + (1|expt), data=finalYear))
summary(evennessModel <- lmer(ANPP_lnRR ~ J_lnRR + (1|expt), data=finalYear))
summary(meanchangeModel <- lmer(ANPP_lnRR ~ mean_change + (1|expt), data=finalYear))
summary(dispersionModel <- lmer(ANPP_lnRR ~ dispersion_lnRR + (1|expt), data=finalYear))

summary(richness <- lm(S_lnRR~ANPP_lnRR, finalYear))
richnessText <- grobTree(textGrob(expression(paste(R^2, " = 0.268, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
S_plot <- ggplot(finalYear, aes(x=S_lnRR, y=ANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='ln RR Richness') +
  scale_y_continuous(name='ln RR ANPP') +
  coord_cartesian(ylim=c(-1,1.5)) +
  annotation_custom(richnessText)
summary(evenness <- lm(J_lnRR~ANPP_lnRR, finalYear))
evennessText <- grobTree(textGrob(expression(paste(R^2, " = 0.254, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
J_plot <- ggplot(finalYear, aes(x=J_lnRR, y=ANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='ln RR Evenness') +
  scale_y_continuous(name='ln RR ANPP') +
  coord_cartesian(ylim=c(-1,1.5)) +
  annotation_custom(evennessText)
summary(meanchange <- lm(mean_change~ANPP_lnRR, finalYear))
meanchangeText <- grobTree(textGrob(expression(paste(R^2, " = 0.444, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
meanchange_plot <- ggplot(finalYear, aes(x=mean_change, y=ANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='ln RR ANPP') +
  coord_cartesian(ylim=c(-1,1.5)) +
  annotation_custom(meanchangeText)
summary(dispersion <- lm(dispersion_lnRR~ANPP_lnRR, finalYear))
dispersionText <- grobTree(textGrob(expression(paste(R^2, " = 0.176, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
dispersion_plot <- ggplot(finalYear, aes(x=dispersion_lnRR, y=ANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='ln RR Multivariate Dispersion') +
  scale_y_continuous(name='ln RR ANPP') +
  coord_cartesian(ylim=c(-1,1.5)) +
  annotation_custom(dispersionText)
pushViewport(viewport(layout=grid.layout(2,2)))
print(S_plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(J_plot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanchange_plot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(dispersion_plot, vp=viewport(layout.pos.row=2, layout.pos.col=2))



#mixed effects models for change in ANPP variance
#these models are not correct, need more work to properly account for random effects of site and experiment
summary(richnessModel <- lmer(varANPP_lnRR ~ S_lnRR + (1|expt), data=finalYear))
summary(evennessModel <- lmer(varANPP_lnRR ~ J_lnRR + (1|expt), data=finalYear))
summary(meanchangeModel <- lmer(varANPP_lnRR ~ mean_change + (1|expt), data=finalYear))
summary(dispersionModel <- lmer(varANPP_lnRR ~ dispersion_lnRR + (1|expt), data=finalYear))

summary(richness2 <- lm(S_lnRR~varANPP_lnRR, finalYear))
richnessText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.147, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
S_plot2 <- ggplot(finalYear, aes(x=S_lnRR, y=varANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='ln RR Richness') +
  scale_y_continuous(name='ln RR ANPP variance') +
  coord_cartesian(ylim=c(-6,4)) +
  annotation_custom(richnessText2)
summary(evenness2 <- lm(J_lnRR~varANPP_lnRR, finalYear))
evennessText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.113, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
J_plot2 <- ggplot(finalYear, aes(x=J_lnRR, y=varANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='ln RR Evenness') +
  scale_y_continuous(name='ln RR ANPP variance') +
  coord_cartesian(ylim=c(-6,4)) +
  annotation_custom(evennessText2)
summary(meanchange2 <- lm(mean_change~varANPP_lnRR, finalYear))
meanchangeText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.052, p = 0.002", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
meanchange_plot2 <- ggplot(finalYear, aes(x=mean_change, y=varANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='ln RR ANPP variance') +
  coord_cartesian(ylim=c(-6,4)) +
  annotation_custom(meanchangeText2)
summary(dispersion2 <- lm(dispersion_lnRR~varANPP_lnRR, finalYear))
dispersionText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.034, p = 0.015", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
dispersion_plot2 <- ggplot(finalYear, aes(x=dispersion_lnRR, y=varANPP_lnRR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='ln RR Multivariate Dispersion') +
  scale_y_continuous(name='ln RR ANPP variance') +
  coord_cartesian(ylim=c(-6,4)) +
  annotation_custom(dispersionText2)
pushViewport(viewport(layout=grid.layout(2,2)))
print(S_plot2, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(J_plot2, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanchange_plot2, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(dispersion_plot2, vp=viewport(layout.pos.row=2, layout.pos.col=2))









