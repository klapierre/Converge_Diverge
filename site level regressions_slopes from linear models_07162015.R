library(ggplot2)
library(plyr)
library(lme4)
library(nlme)
library(colorspace)
library(colorRamps)
library(visreg)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=32, vjust=-0.35), axis.text.x=element_text(size=24),
             axis.title.y=element_text(size=32, angle=90, vjust=1), axis.text.y=element_text(size=24),
             plot.title = element_text(size=32, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.text=element_text(size=20))


setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

#read in experiment info
info <- read.csv('exp_info072015.csv')

#read in full change in mean and dispersion dataset
all <- read.csv('dispersion_and_means_press_experiments_with_exp_info_08062015.csv')

#subset out change in means
means <- subset(all, subset=(mean.disp=='mean'))
means$dist_log <- log(means$dist)
means$dist_sqrt <- sqrt(means$dist)
means$trt.year_log <- log(means$trt.year)

#run mixed effects model with full dataset for change in mean
summary(meansFullModelExpt <- lmer(dist_sqrt ~ trt.year*plot_mani.x + (1|label), data=subset(means, trt.year>0), family=binomial(logit)))
# summary(meansFullModelExpt <- lmer(dist ~ trt.year_log*plot_mani.x + (1|expt/(trt.year_log*plot_mani.x)+label), data=subset(means, trt.year>0)))
# summary(meansFullModelExpt <- lmer(dist ~ trt.year_log*plot_mani.x + (1|expt/(trt.year_log*plot_mani.x)+label), data=subset(means, trt.year>0)))
# summary(meansFullModelExpt <- lme(dist ~ trt.year_log*plot_mani.x, random=~1|expt/label, data=subset(means, trt.year>0)))
# summary(meansFullModelExpt <- lmer(dist_sqrt~plot_mani.x*trt.year, random=~1|expt, data=subset(means, trt.year>0)))

plot(ranef(meansFullModelExpt)$label[,2], ranef(meansFullModelExpt)$label[,3])

#get slopes and intercepts of fit model
meansIntercepts <- fixef(meansFullModelExpt)[1] + unlist(ranef(meansFullModelExpt)$label)
meansSlopes <- fixef(meansFullModelExpt)[2] + unlist(ranef(meansFullModelExpt)$label)
meansModelData <- data.frame(intercepts=meansIntercepts, slopes=meansSlopes, label=unique(means$label))
labelPlotMani <- unique(means[c('label', 'plot_mani.x')])
meansModelDataLabel <- merge(meansModelData, labelPlotMani, by='label')


#figures for experiment-level data for change in mean
colorMani <- c('#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#253494', '#081d58')
colorManiGreens <- c('#c7e9c0', '#a1d99b', '#74c476', '#41ab5d', '#238b45', '#006d2c', '#00441b')
colorManiReds <- c('#fed976', '#feb24c', '#fd8d3c', '#fc4e2a', '#e31a1c', '#bd0026', '#800026')
colorManiDiverge <- c('#313695', '#4575b4', '#74add1', '#abd9e9', '#f46d43', '#d73027', '#a50026')

means$order <- factor(as.character(means$plot_mani.x), levels=as.character(c(0:7)))

#means mixed effect model visualization
ggplot(data=subset(means, trt.year>0), aes(x=trt.year, y=dist_sqrt)) +
  geom_point(pch=1) +
  facet_wrap(~expt) +
  geom_abline(aes(intercept=intercepts, slope=slopes), alpha=0.5, data=meansModelDataLabel) +
  scale_x_continuous('Treatment Year') +
  scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_y_continuous('Distance Between Centroids', limits=c(0,1)) +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24)) +
  guides(col=guide_legend(nrow=1, override.aes=list(size=1)))

#raw data figure
ggplot(data=subset(means, trt.year>0), aes(x=trt.year, y=dist_sqrt)) +
  #geom_smooth(aes(y=dist, colour=order, group=interaction(expt, order)), method=lm, formula=y~poly(x,2), se=F, size=0.25) +
  #geom_line(aes(y=dist, colour=order, group=interaction(expt, order))) +
  geom_smooth(aes(y=dist_sqrt, colour=order, group=interaction(expt, order)), method=lm, formula=y~log(x), se=F, size=0.25) +
  geom_smooth(aes(y=dist_sqrt, colour=order, fill=order, group=order), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
  scale_x_continuous('Treatment Year') +
  scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_y_continuous('Distance Between Centroids', limits=c(0,1)) +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24)) +
  guides(col=guide_legend(nrow=1, override.aes=list(size=1)))

#look at each site's response
ggplot(data=subset(means, trt.year>0), aes(x=trt.year, y=dist_sqrt)) +
  geom_smooth(method=loess, aes(y=dist, colour=order, group=interaction(expt, order)), se=F) +
  scale_x_continuous('Treatment Year') +
  scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_y_continuous('Distance Between Centroids', limits=c(0,1)) +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24)) +
  guides(col=guide_legend(override.aes=list(size=1))) +
  facet_wrap(~site_code, ncol=5)

#means with trt intensity
means$n_intensity <- with(means, n+p+k-soil_carbon)

#run mixed effects model with nutrient addition dataset for change in mean
summary(meansIntensityExpt <- lmer(dist_sqrt ~ trt.year*n_intensity + (1|label), data=subset(means, trt.year>0 & nutrients!=0), family=binomial(logit)))

ggplot(data=subset(means, trt.year>0 & nutrients!=0), aes(x=trt.year, y=dist_sqrt)) +
  #geom_smooth(aes(y=dist, colour=order, group=interaction(expt, order)), method=lm, formula=y~poly(x,2), se=F, size=0.25) +
  #geom_line(aes(y=dist, colour=order, group=interaction(expt, order))) +
  geom_smooth(aes(y=dist_sqrt, colour=n_intensity, group=interaction(expt, n_intensity)), method=lm, formula=y~log(x), se=F, size=0.25) +
  geom_smooth(aes(y=dist_sqrt, colour=n_intensity, fill=n_intensity, group=n_intensity), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
  scale_x_continuous('Treatment Year') +
#   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_y_continuous('Distance Between Centroids', limits=c(0,1)) +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24))

#run mixed effects model with precip trt dataset for change in mean
summary(meansIntensityExpt <- lmer(dist_sqrt ~ trt.year*precip + (1|label), data=subset(means, trt.year>0 & nutrients!=0), family=binomial(logit)))

ggplot(data=subset(means, trt.year>0 & water!=0), aes(x=trt.year, y=dist_sqrt)) +
  #geom_smooth(aes(y=dist, colour=order, group=interaction(expt, order)), method=lm, formula=y~poly(x,2), se=F, size=0.25) +
  #geom_line(aes(y=dist, colour=order, group=interaction(expt, order))) +
  geom_smooth(aes(y=dist_sqrt, colour=precip, group=interaction(expt, precip)), method=lm, formula=y~log(x), se=F, size=0.25) +
  geom_smooth(aes(y=dist_sqrt, colour=precip, fill=precip, group=precip), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
  scale_x_continuous('Treatment Year') +
#   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_colour_gradientn(colour=rainbow(3), name='Nutrient\nIntensity') +
  scale_y_continuous('Distance Between Centroids', limits=c(0,1)) +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24)) +
  guides(col=guide_legend(nrow=1, override.aes=list(size=1)))

###########################################
###########################################
###########################################


#subset out dispersion data
dispersion <- subset(all, subset=(mean.disp=='disp'))

#average dispersion across plots
dispersionAverage <- aggregate(data=dispersion, dist ~ site_code + project_name + community_type + label + expt.year + plot_mani.x + dist.1 + dist.2 + mean.disp + site + cal.year + trt.year + expt + dg + treatment + nutrients + n + p + k + soil_carbon + precip + water, mean)

#subset out controls and treatments and merge back with treatment data
dispersionTrt <- subset(dispersionAverage, subset=(plot_mani.x>0))
dispersionCtl <- subset(dispersionAverage, subset=(plot_mani.x==0))
keep <- c('site_code', 'project_name', 'community_type', 'cal.year', 'dist')
dispersionCtlCondensed <- dispersionCtl[,colnames(dispersionCtl) %in% keep]
names(dispersionCtlCondensed)[names(dispersionCtlCondensed)=='dist'] <- 'ctl_dist'
dispersionDiff <- merge(dispersionCtlCondensed, dispersionTrt, by=c("site_code", "project_name", "community_type", "cal.year"))

#calculate response ratio
dispersionDiff$dist_RR <- with(dispersionDiff, (dist-ctl_dist)/ctl_dist)

#run mixed effects model with full dataset for change in mean
summary(meansFullModelExpt <- lme(dist_RR~plot_mani.x*trt.year, random=~1|expt, data=subset(dispersionDiff, trt.year>0)))

#figure at experiment-level of disperison
dispersionDiff$order <- factor(as.character(dispersionDiff$plot_mani.x), levels=as.character(c(0:7)))

ggplot(data=subset(dispersionDiff, trt.year>0 & dist_RR<3), aes(x=trt.year, y=dist_RR)) +
  #geom_point(aes(y=dist_RR, colour=order)) +
  geom_smooth(aes(y=dist_RR, colour=order, group=interaction(expt, order)), method=lm, formula=y~log(x), se=F, size=0.25) +
  geom_smooth(aes(y=dist_RR, colour=order, fill=order, group=order), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
  scale_x_continuous('Treatment Year') +
  scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_y_continuous('Relative Difference in Dispersion\nbetween Treatment and Control') +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24)) +
  guides(col=guide_legend(nrow=1, override.aes=list(size=1)))

#no effect of number of factors
ggplot(data=subset(dispersionDiff, trt.year>0 & dist_RR<3), aes(x=trt.year, y=dist_RR)) +
  #geom_point(aes(y=dist_RR, colour=order)) +
  geom_smooth(aes(y=dist_RR, group=interaction(expt, order)), colour='gray40', method=lm, formula=y~log(x), se=F, size=0.25) +
  geom_smooth(aes(y=dist_RR), colour='black', method=lm, formula=y~log(x), size=3, se=T, alpha=0.9) +
  scale_x_continuous('Treatment Year') +
  scale_y_continuous('Relative Difference in Dispersion\nbetween Treatment and Control')

#look at each site's response
ggplot(data=subset(dispersionDiff, trt.year>0 & dist_RR<3 & site_code=='CDR'), aes(x=trt.year, y=dist_RR)) +
  geom_smooth(method=lm, aes(y=dist, colour=order, group=interaction(expt, order)), se=F) +
  scale_x_continuous('Treatment Year') +
  scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_y_continuous('Relative Difference in Dispersion\nbetween Treatment and Control') +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24)) +
  guides(col=guide_legend(override.aes=list(size=1))) +
  facet_wrap(~site_code, ncol=5) +
  facet_wrap(~expt)

#dispersion with trt intensity
dispersionDiff$n_intensity <- with(dispersionDiff, n+p+k-soil_carbon)

#run mixed effects model with nutrient addition dataset for change in mean
summary(dispersionIntensityExpt <- lmer(dist_RR ~ trt.year*n_intensity + (1|label), data=subset(dispersionDiff, trt.year>0 & nutrients!=0)))

ggplot(data=subset(dispersionDiff, trt.year>0 & nutrients!=0), aes(x=trt.year, y=dist_RR)) +
  #geom_smooth(aes(y=dist, colour=order, group=interaction(expt, order)), method=lm, formula=y~poly(x,2), se=F, size=0.25) +
  #geom_line(aes(y=dist, colour=order, group=interaction(expt, order))) +
  geom_smooth(aes(y=dist_RR, colour=n_intensity, group=interaction(n_intensity)), method=lm, formula=y~log(x), se=F, size=0.25) +
  scale_x_continuous('Treatment Year') +
  #   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  #   scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  scale_y_continuous('dist_RR') +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24))

#run mixed effects model with precip trt dataset for change in mean
summary(dispersionIntensityExpt <- lmer(dist_RR ~ trt.year*precip + (1|label), data=subset(dispersionDiff, trt.year>0 & water!=0)))

ggplot(data=subset(dispersionDiff, trt.year>0 & water!=0), aes(x=trt.year, y=dist_RR)) +
  #geom_smooth(aes(y=dist, colour=order, group=interaction(expt, order)), method=lm, formula=y~poly(x,2), se=F, size=0.25) +
  #geom_line(aes(y=dist, colour=order, group=interaction(expt, order))) +
  geom_smooth(aes(y=dist_RR, colour=precip, group=interaction(precip)), method=lm, formula=y~log(x), se=F, size=0.25) +
  scale_x_continuous('Treatment Year') +
  #   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
  #   scale_colour_gradientn(colour=rainbow(3), name='Nutrient\nIntensity') +
  scale_y_continuous('dist_RR') +
  theme(legend.position='right', legend.direction='vertical',
        legend.title=element_text(size=24)) +
  guides(col=guide_legend(nrow=1, override.aes=list(size=1)))

# #find experiments where diserpsion increased or decreased in last year
# finalYear <- aggregate(trt.year~expt, max, data=dispersionDiff)
# dispersionDiffFinal <- merge(finalYear, dispersionDiff, by=c('expt', 'trt.year'))
# dispersionDiffFinal$con_div <- with(dispersionDiffFinal, ifelse(dist_RR<=0, 'con', 'div'))
# keep <- c('label', 'con_div')
# dispersionDiffFinal <- dispersionDiffFinal[,colnames(dispersionDiffFinal) %in% keep]
# dispersionDiffCategorized <- merge(dispersionDiffFinal, dispersionDiff, by='label')
# 
# #only datasets that converged based on final year
# ggplot(data=subset(dispersionDiffCategorized, trt.year>0 & dist_RR<3 & con_div=='con'), aes(x=trt.year, y=dist_RR)) +
#   #geom_point(aes(y=dist_RR, colour=order)) +
#   geom_smooth(aes(y=dist_RR, colour=order, group=interaction(expt, order)), method=lm, formula=y~log(x), se=F, size=0.25) +
#   geom_smooth(aes(y=dist_RR, colour=order, fill=order, group=order), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
#   scale_x_continuous('Treatment Year') +
#   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_y_continuous('Relative Difference in Dispersion\nbetween Treatment and Control') +
#   theme(legend.position='right', legend.direction='vertical',
#         legend.title=element_text(size=24)) +
#   guides(col=guide_legend(nrow=1, override.aes=list(size=1)))
# 
# #only datasets that diverged based on final year
# ggplot(data=subset(dispersionDiffCategorized, trt.year>0 & dist_RR<3 & con_div=='div'), aes(x=trt.year, y=dist_RR)) +
#   #geom_point(aes(y=dist_RR, colour=order)) +
#   geom_smooth(aes(y=dist_RR, colour=order, group=interaction(expt, order)), method=lm, formula=y~log(x), se=F, size=0.25) +
#   geom_smooth(aes(y=dist_RR, colour=order, fill=order, group=order), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
#   scale_x_continuous('Treatment Year') +
#   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_y_continuous('Relative Difference in Dispersion\nbetween Treatment and Control') +
#   theme(legend.position='right', legend.direction='vertical',
#         legend.title=element_text(size=24)) +
#   guides(col=guide_legend(nrow=1, override.aes=list(size=1)))
# 
# #get datasets that converge or diverge based on slope
# slopesDisp <- ddply(dispersionDiff, c('expt', 'label'), function(x) {
#   model <- lm(dist_RR~log(trt.year+1), data=x)
#   coef(model)
# })
# names(slopesDisp)[names(slopesDisp)=="log(trt.year + 1)"] <- "slope"
# slopesDisp$con_div <- with(slopesDisp, ifelse(slope<=0, 'con', 'div'))
# keep <- c('label', 'con_div')
# slopesDisp <- slopesDisp[,colnames(slopesDisp) %in% keep]
# dispersionDiffSlopeCategorized <- merge(slopesDisp, dispersionDiff, by='label')

# #only datasets that converged based on slope
# ggplot(data=subset(dispersionDiffSlopeCategorized, trt.year>0 & dist_RR<3 & con_div=='con'), aes(x=trt.year, y=dist_RR)) +
#   #geom_point(aes(y=dist_RR, colour=order)) +
#   geom_smooth(aes(y=dist_RR, colour=order, group=interaction(expt, order)), method=lm, formula=y~log(x), se=F, size=0.25) +
#   geom_smooth(aes(y=dist_RR, colour=order, fill=order, group=order), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
#   scale_x_continuous('Treatment Year') +
#   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_y_continuous('Relative Difference in Dispersion\nbetween Treatment and Control') +
#   theme(legend.position='right', legend.direction='vertical',
#         legend.title=element_text(size=24)) +
#   guides(col=guide_legend(nrow=1, override.aes=list(size=1)))
# 
# #only datasets that diverged based on slope
# ggplot(data=subset(dispersionDiffSlopeCategorized, trt.year>0 & dist_RR<3 & con_div=='div'), aes(x=trt.year, y=dist_RR)) +
#   #geom_point(aes(y=dist_RR, colour=order)) +
#   geom_smooth(aes(y=dist_RR, colour=order, group=interaction(expt, order)), method=lm, formula=y~log(x), se=F, size=0.25) +
#   geom_smooth(aes(y=dist_RR, colour=order, fill=order, group=order), method=lm, formula=y~log(x), size=3, se=T, alpha=0.5) +
#   scale_x_continuous('Treatment Year') +
#   scale_fill_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_colour_manual(values=colorManiDiverge, name='Factors\nManipulated') +
#   scale_y_continuous('Relative Difference in Dispersion\nbetween Treatment and Control') +
#   theme(legend.position='right', legend.direction='vertical',
#         legend.title=element_text(size=24)) +
#   guides(col=guide_legend(nrow=1, override.aes=list(size=1)))




































# meanTrue <- ddply(maniType, c('mean.disp', 'plot_mani.x', 'trt.year', 'expt'), summarise,
#                   mean=mean(dist),
#                   N=length(dist),
#                   sd=sd(dist),
#                   se=sd/sqrt(N))
# 
# meanTrue$order <- factor(as.character(meanTrue$plot_mani.x), levels=as.character(c(1:7)))
# 
# meanPredicted <- ddply(maniType, c('mean.disp', 'plot_mani.x', 'trt.year'), summarise,
#                        mean=mean(dist),
#                        N=length(dist),
#                        sd=sd(dist),
#                        se=sd/sqrt(N))
# meanPredictedComplete <- meanPredicted[complete.cases(meanPredicted[,4]),]
# 
# meanPredictedComplete$order <- factor(as.character(meanPredictedComplete$plot_mani.x), levels=as.character(c(1:7)))

#using full dataset, get slopes of means through time for each treatment independantly
slopesMean <- ddply(subset(full, mean.disp=='mean'), c('expt', 'label'), function(x) {
  model <- lm(dist~trt.year, data=x)
  coef(model)
})
names(slopesMean)[names(slopesMean)=="trt.year"] <- "slope"
  
#merge slopes data with controls and extra info
slopesMeanControls <- merge(slopesMean, controlsInfo)
slopesMeanInfo <- merge(slopesMeanControls, info)

# summary(slopeMeansModel <- lm(slope ~ plot_mani.x, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=plot_mani.x, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ expt_length, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=expt_length, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ MAP, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=MAP, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ ANPP, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=ANPP, y=slope)) +
#   geom_point()
# 
# summary(slopeMeansModel <- lm(slope ~ species_num, data=slopesMeanInfo))
# ggplot(data=slopesMeanInfo, aes(x=species_num, y=slope)) +
#   geom_point()

#using full dataset, get slopes of dispersion through time for each treatment independantly
slopesDispersion <- ddply(subset(full, mean.disp=='disp'), c('expt', 'label'), function(x) {
  model <- lm(dist~trt.year, data=x)
  coef(model)
})
names(slopesDispersion)[names(slopesDispersion)=="trt.year"] <- "slope"

#merge slopes data with controls and extra info
slopesDispersionInfo <- merge(slopesDispersion, controlsAggregate)

#get treatment and control effects as columns
slopesDispersionInfoTrt <- subset(slopesDispersionInfo, subset=(plot_mani.x!=0))
slopesDispersionInfoCtl <- subset(slopesDispersionInfo, subset=(plot_mani.x==0))
slopesDispersionInfoCtlSlim <- slopesDispersionInfoCtl[,-c(1, 3, 5)]
names(slopesDispersionInfoCtlSlim)[names(slopesDispersionInfoCtlSlim)=='slope'] <- 'ctl_slope'
slopesDispersionInfoRR <- merge(slopesDispersionInfoTrt, slopesDispersionInfoCtlSlim)

#separate out convergence vs divergence datasets
slopesDispersionRRconverge <- subset(slopesDispersionInfoRR, subset=(slope<0))
slopesDispersionRRconverge$converge_diverge <- 'converge'
slopesDispersionRRdiverge <- subset(slopesDispersionInfoRR, subset=(slope>0))
slopesDispersionRRdiverge$converge_diverge <- 'diverge'

# #get response ratios
# slopesDispersionRRconverge$RR <- with(slopesDispersionRRconverge, ifelse(ctl_slope>0, -1*(slope - ctl_slope) / ctl_slope,                                                                      (slope - ctl_slope) / ctl_slope))
# slopesDispersionRRdiverge$RR <- with(slopesDispersionRRdiverge, ifelse(ctl_slope<0, -1*(slope - ctl_slope) / ctl_slope,
#                                                            (slope - ctl_slope) / ctl_slope))
# 
# #merge experiment info with the response ratio dataframes
# slopesDispersionInfoRRconverge <- merge(slopesDispersionRRconverge, info)
# slopesDispersionInfoRRdiverge <- merge(slopesDispersionRRdiverge, info)
# 
# #number of factors model
# summary(lm(RR ~ plot_mani.x, data=slopesDispersionInfoRRconverge))
# ggplot(data=slopesDispersionInfoRRconverge, aes(x=plot_mani.x, y=RR)) +
#   geom_point()
# 
# summary(lm(RR ~ plot_mani.x, data=slopesDispersionInfoRRdiverge))
# ggplot(data=slopesDispersionInfoRRdiverge, aes(x=plot_mani.x, y=RR)) +
#   geom_point()
# 
# #gamma diversity model
# summary(lm(RR ~ species_num, data=slopesDispersionInfoRRconverge))
# ggplot(data=slopesDispersionInfoRRconverge, aes(x=species_num, y=RR)) +
#   geom_point()
# 
# summary(lm(RR ~ species_num, data=slopesDispersionInfoRRdiverge))
# ggplot(data=slopesDispersionInfoRRdiverge, aes(x=species_num, y=RR)) +
#   geom_point()
# 
# #MAP model
# summary(lm(RR ~ MAP, data=slopesDispersionInfoRRconverge))
# ggplot(data=slopesDispersionInfoRRconverge, aes(x=MAP, y=RR)) +
#   geom_point()
# 
# summary(lm(RR ~ MAP, data=slopesDispersionInfoRRdiverge))
# ggplot(data=slopesDispersionInfoRRdiverge, aes(x=MAP, y=RR)) +
#   geom_point()
# 
# #ANPP model
# summary(lm(RR ~ ANPP, data=slopesDispersionInfoRRconverge))
# ggplot(data=slopesDispersionInfoRRconverge, aes(x=ANPP, y=RR)) +
#   geom_point()
# 
# summary(lm(RR ~ ANPP, data=slopesDispersionInfoRRdiverge))
# ggplot(data=slopesDispersionInfoRRdiverge, aes(x=ANPP, y=RR)) +
#   geom_point()






#make dataset indicating what control treatments are, plus all other relevant information
controls <- all[,-c(1:5, 7:8, 10:74)]
controlsAggregate <- aggregate(plot_mani.x ~ label, mean, data=all)
exptLength <- aggregate(trt.year ~ label, max, data=all)
names(exptLength)[names(exptLength)=='trt.year'] <- 'expt_length'
controlsInfo <- merge(controlsAggregate, exptLength)

#make column of manipulation type
full$num <- with(full, nut2+pp+car+heat+other)
maniType <- full
maniType$mani_type <- with(maniType, ifelse(nut2==1, 'nutrients',
                                            ifelse(pp==1, 'water',
                                                   ifelse(car==1, 'carbon',
                                                          ifelse(num==0, 'control', 'other')))))
maniType$mani_type <- with(maniType, ifelse(num>1, 'multiple', mani_type))








