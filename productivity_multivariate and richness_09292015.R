library(ggplot2)
library(plyr)
library(vegan)
library(lme4)
library(grid)
library(nlme)
library(reshape2)


setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank())

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

###################################################################################################
###################################################################################################
###################################################################################################

#read in experiment info
info <- read.csv('exp_info072015.csv')

#read in community dataset (multivariate and basic diversity stats data)
community <- read.csv('dispersion_means_press_bayesian_11172015.csv')
names(community)[names(community)=='ANPP'] <- 'ANPP_sitelevel'

#read in anpp data
anpp <- read.csv('AllAnppData_02172015.csv')

#treatment types
treatments <- read.csv('treatments.csv')

###################################################################################################
#calculate mean anpp for each treatment and trt year at each site
anppMean <- ddply(anpp, c('site_code', 'project_name', 'community_type', 'treatment_year', 'treatment'), summarise,
                  ANPP_trtmean=mean(anpp),
                  ANPP_trtvar=var(anpp))
names(anppMean)[names(anppMean)=='treatment_year'] <- 'trt.year'

#merge anpp and community data
anppComm <- merge(anppMean, community, by=c('site_code', 'project_name', 'community_type', 'trt.year', 'treatment'))

#calculate response ratios
anppComm$eH <- exp(anppComm$H_mean)
trt <- subset(anppComm, subset=(plot_mani>0))
ctl <- subset(anppComm, subset=(plot_mani==0))
keep <- c('expt', 'trt.year', 'dispersion', 'S_mean', 'evenness_mean', 'eH',  'ANPP_trtmean', 'ANPP_trtvar')
ctl <- ctl[,colnames(ctl) %in% keep]
names(ctl)[names(ctl)=='dispersion'] <- 'dispersion_ctl'
names(ctl)[names(ctl)=='S_mean'] <- 'S_ctl'
names(ctl)[names(ctl)=='evenness_mean'] <- 'evenness_ctl'
names(ctl)[names(ctl)=='ANPP_trtmean'] <- 'ANPP_trtmean_ctl'
names(ctl)[names(ctl)=='ANPP_trtvar'] <- 'ANPP_trtvar_ctl'
names(ctl)[names(ctl)=='eH'] <- 'H_ctl'
RR <- merge(trt, ctl, by=c('expt', 'trt.year'))
#note that mean_change is already a response ratio of sorts
RR$dispersion_lnRR <- with(RR, log(dispersion/dispersion_ctl))
RR$dispersion_percent <- with(RR, abs(100*(dispersion-dispersion_ctl)/dispersion_ctl))
RR$S_lnRR <- with(RR, log(S_mean/S_ctl))
RR$S_percent <- with(RR, abs(100*(S_mean-S_ctl)/S_ctl))
RR$evenness_lnRR <- with(RR, log(evenness_mean/evenness_ctl))
RR$evenness_percent <- with(RR, abs(100*(evenness_mean-evenness_ctl)/evenness_ctl))
RR$eH_percent <- with(RR, abs(100*(eH-H_ctl)/H_ctl))
RR$ANPP_lnRR <- with(RR, log(ANPP_trtmean/ANPP_trtmean_ctl))
RR$ANPP_percent <- with(RR, 100*(ANPP_trtmean-ANPP_trtmean_ctl)/ANPP_trtmean_ctl)
RR$varANPP_lnRR <- with(RR, log(ANPP_trtvar/ANPP_trtvar_ctl))
RR$varANPP_percent <- with(RR, 100*(ANPP_trtvar-ANPP_trtvar_ctl)/ANPP_trtvar_ctl)
RR$mean_change_percent <- with(RR, 100*mean_change)


RRtrt <- merge(RR, treatments, by=c('expt', 'label', 'treatment'))

#get final year of data for each experiment
finalYear <- RRtrt[with(RRtrt, ave(trt.year, label, FUN=max)==trt.year),]

#mixed effects models for mean ANPP change - ln response ratio
#these models are not correct, need more work to properly account for random effects of site and experiment
summary(richnessModel <- lmer(ANPP_lnRR ~ S_lnRR + (1|expt), data=finalYear))
summary(evennessModel <- lmer(ANPP_lnRR ~ evenness_lnRR + (1|expt), data=finalYear))
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


#mixed effects models for mean ANPP change - percent change (following Leuzinger's figures)
#these models are not correct, need more work to properly account for random effects of site and experiment
summary(richnessModel <- lmer(ANPP_percent ~ S_percent + (1|expt), data=finalYear))
summary(evennessModel <- lmer(ANPP_percent ~ J_percent + (1|expt), data=finalYear))
summary(meanchangeModel <- lmer(ANPP_percent ~ mean_change + (1|expt), data=finalYear))
summary(dispersionModel <- lmer(ANPP_percent ~ dispersion_lnRR + (1|expt), data=finalYear))

summary(richness <- lm(S_percent~ANPP_percent, finalYear))
richnessText <- grobTree(textGrob(expression(paste(R^2, " = 0.254, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
S_plot <- ggplot(finalYear, aes(x=S_percent, y=ANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Change in Richness (%)') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(richnessText)
summary(evenness <- lm(J_percent~ANPP_percent, finalYear))
evennessText <- grobTree(textGrob(expression(paste(R^2, " = 0.315, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
J_plot <- ggplot(finalYear, aes(x=J_percent, y=ANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Change in Evenness (%)') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(evennessText)
summary(meanchange <- lm(mean_change~ANPP_percent, finalYear))
meanchangeText <- grobTree(textGrob(expression(paste(R^2, " = 0.482, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
meanchange_plot <- ggplot(finalYear, aes(x=mean_change, y=ANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(meanchangeText)
summary(dispersion <- lm(dispersion_lnRR~ANPP_percent, finalYear))
dispersionText <- grobTree(textGrob(expression(paste(R^2, " = 0.210, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
dispersion_plot <- ggplot(finalYear, aes(x=dispersion_lnRR, y=ANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='ln RR Multivariate Dispersion') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(dispersionText)
pushViewport(viewport(layout=grid.layout(2,2)))
print(S_plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(J_plot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanchange_plot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(dispersion_plot, vp=viewport(layout.pos.row=2, layout.pos.col=2))

#mean change and anpp percent by factor manipulated
finalYear$order <- factor(finalYear$factors_mani, levels=c('CO2 (C)', 'water (W)', 'nutrients (N)', 'non-resource (O)', 'C+O', 'W+C', 'W+O', 'W+C+O', 'N+C', 'N+W', 'N+O', 'N+W+C', 'N+W+O', 'N+W+C+O'))
summary(meanchange <- lm(mean_change~ANPP_percent, finalYear))
meanchangeText <- grobTree(textGrob(expression(paste(R^2, " = 0.482, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(finalYear, aes(x=mean_change, y=ANPP_percent)) +
  geom_point(aes(colour=order), size=3) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(meanchangeText)

summary(meanchange <- lm(mean_change~ANPP_percent, finalYear))
meanchangeText <- grobTree(textGrob(expression(paste(R^2, " = 0.482, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(finalYear, aes(x=mean_change, y=ANPP_percent)) +
  geom_point(size=3) +
  geom_smooth(method=lm, formula=y~x^3, colour='black') +
  scale_x_continuous(name='Distance Between Treatment Centroids') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(meanchangeText)
summary(richness <- lm(S_percent~ANPP_percent, finalYear))
meanchangeText <- grobTree(textGrob(expression(paste(R^2, " = 0.254, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(finalYear, aes(x=S_percent, y=ANPP_percent)) +
  geom_point(size=3) +
  geom_smooth(method=lm, formula=y~x^3, colour='black') +
  scale_x_continuous(name='Change in Richness (%)') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(meanchangeText)

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




#productivity mean change and variance predicted by change in mean and variance in plant community
#final year only
summary(meanchangeModel <- lmer(ANPP_percent ~ mean_change + (1|expt), data=finalYear))
summary(dispersionModel <- lmer(ANPP_percent ~ dispersion_percent + (1|expt), data=finalYear))
summary(meanchangevarModel <- lmer(varANPP_percent ~ mean_change + (1|expt), data=finalYear))
summary(dispersionvarModel <- lmer(varANPP_percent ~ dispersion_percent + (1|expt), data=finalYear))

summary(meanchange <- lm(mean_change~ANPP_percent, finalYear))
meanchangeText <- grobTree(textGrob(expression(paste(R^2, " = 0.482, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
meanchange_plot <- ggplot(finalYear, aes(x=mean_change, y=ANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(meanchangeText)
summary(dispersion <- lm(dispersion_percent~ANPP_percent, finalYear))
dispersionText <- grobTree(textGrob(expression(paste(R^2, " = 0.187, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
dispersion_plot <- ggplot(finalYear, aes(x=dispersion_percent, y=ANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Change in Multivariate Dispersion (%)') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(dispersionText)
summary(meanchange2 <- lm(mean_change~varANPP_percent, finalYear))
meanchangeText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.020, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
meanchange_plot2 <- ggplot(finalYear, aes(x=mean_change, y=varANPP_percent)) +
  geom_point(shape=1) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='Change in ANPP Variance (%)') 
summary(dispersion2 <- lm(dispersion_percent~varANPP_percent, finalYear))
dispersionText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.210, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
dispersion_plot2 <- ggplot(finalYear, aes(x=dispersion_percent, y=varANPP_percent)) +
  geom_point(shape=1) +
  scale_x_continuous(name='Change in Multivariate Dispersion (%)') +
  scale_y_continuous(name='Change in ANPP Variance (%)') 
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanchange_plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersion_plot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanchange_plot2, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(dispersion_plot2, vp=viewport(layout.pos.row=2, layout.pos.col=2))



#productivity mean change and variance predicted by change in mean and variance in plant community
#all years
summary(meanchangeModel <- lmer(ANPP_percent ~ mean_change + (1|expt.year), data=RRtrt))
summary(dispersionModel <- lmer(ANPP_percent ~ dispersion_percent + (1|expt.year), data=RRtrt))
summary(meanchangevarModel <- lmer(varANPP_percent ~ mean_change + (1|expt.year), data=RRtrt))
summary(dispersionvarModel <- lmer(varANPP_percent ~ dispersion_percent + (1|expt.year), data=RRtrt))

summary(meanchange <- lm(mean_change~ANPP_percent, RRtrt))
meanchangeText <- grobTree(textGrob(expression(paste(R^2, " = 0.149, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
meanchange_plot <- ggplot(RRtrt, aes(x=mean_change, y=ANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='Change in ANPP (%)') +
  annotation_custom(meanchangeText)
summary(dispersion <- lm(dispersion_percent~ANPP_percent, RRtrt))
# dispersionText <- grobTree(textGrob(expression(paste(R^2, " = 0.187, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
dispersion_plot <- ggplot(RRtrt, aes(x=dispersion_percent, y=ANPP_percent)) +
  geom_point(shape=1) +
  scale_x_continuous(name='Change in Multivariate Dispersion (%)') +
  scale_y_continuous(name='Change in ANPP (%)')
summary(meanchange2 <- lm(mean_change~varANPP_percent, RRtrt))
# meanchangeText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.002, p < 0.001", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
meanchange_plot2 <- ggplot(RRtrt, aes(x=mean_change, y=varANPP_percent)) +
  geom_point(shape=1) +
  scale_x_continuous(name='Mean Multivariate Shift') +
  scale_y_continuous(name='Change in ANPP Variance (%)') 
summary(dispersion2 <- lm(dispersion_percent~varANPP_percent, RRtrt))
dispersionText2 <- grobTree(textGrob(expression(paste(R^2, " = 0.002, p = 0.045", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
dispersion_plot2 <- ggplot(RRtrt, aes(x=dispersion_percent, y=varANPP_percent)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Change in Multivariate Dispersion (%)') +
  scale_y_continuous(name='Change in ANPP Variance (%)') +
  annotation_custom(dispersionText2)
pushViewport(viewport(layout=grid.layout(2,2)))
print(meanchange_plot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispersion_plot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanchange_plot2, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(dispersion_plot2, vp=viewport(layout.pos.row=2, layout.pos.col=2))


summary(ANPPvar <- lm(dispersion~ANPP_trtvar, RRtrt))
# dispersionText3 <- grobTree(textGrob(expression(paste(R^2, " = 0.002, p = 0.045", sep="")), x=0.05, y=0.05, hjust=0, gp=gpar(col="black", fontsize=15)))
ggplot(RRtrt, aes(x=dispersion, y=ANPP_trtvar)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Multivariate Dispersion') +
  scale_y_continuous(name='Variance in ANPP')




ggplot(RRtrt, aes(x=S_percent, y=mean_change)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Richness Change') +
  scale_y_continuous(name='Mean Change')




ggplot(RRtrt, aes(x=J_percent, y=mean_change)) +
  geom_point(shape=1) +
  geom_smooth(method=lm) +
  scale_x_continuous(name='Evenness Change') +
  scale_y_continuous(name='Mean Change')



ggplot(RRtrt, aes(x=plot_mani, y=ANPP_percent)) +
  geom_point()
ggplot(RRtrt, aes(x=plot_mani, y=mean_change)) +
  geom_point()


keep <- c('label', 'treatment', 'S_percent', 'eH_percent', 'evenness_percent', 'mean_change_percent', 'dispersion_percent')
RRpercentChange <- RRtrt[,colnames(RRtrt) %in% keep]
RRpercentChangeLong <- melt(RRpercentChange, id.vars=c('label', 'treatment'))

ggplot(barGraphStats(data=RRpercentChangeLong, variable='value', byFactorNames=c('variable')), aes(x=variable, y=mean)) +
         geom_bar(stat='identity', fill='white', colour='black') +
         geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2))


