setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

library(vegan)
library(reshape2)
library(ggplot2)
library(gtools)
library(plyr)
library(grid)
library(nlme)

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"),
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

#read in the dispersion/mean dataset
allData<-read.csv('dispersion_and_means_press_experiments_with_exp_info_08062015.csv')

#read in experiment information dataset
expInfo <- read.csv('exp_info072015.csv')

############################################################################################
#subset out only mean comparison
changeInMean <- subset(allData, mean.disp=="mean")

#subset out only dispersion comparison for controls and treatments
dispersion <- subset(allData, mean.disp=="disp")

#means of dispersion across all plots within a year and treatment
dispersionMeans <- ddply(dispersion, 
                         c("expt.year", "plot_mani.x", "dist.1", "mean.disp", "site",
                           "cal.year", "trt.year", "expt", "dg", "label", "water", "carbon", "nutrients", 
                           "light", 'precip', 'n', 'p', 'k'), 
                         summarise,
                         disp=mean(dist))

#create columns for dispersion in control plots and dispersion in treatment plots
dispersionControls <- subset(dispersionMeans, dg=="c_c")
dispersionTreatments <- subset(dispersionMeans, dg=="t_t")

dispersionControls <- rename(dispersionControls, c("disp"="dispersionControl", "dist.1"="ctl"))
dispersionControls <- dispersionControls[,-c(2:4, 9:10)]
dispersionTreatments <- rename(dispersionTreatments, c("disp"="dispersionTreatment", "dist.1"="trt"))
dispersionTreatments <- dispersionTreatments[,-9]

#merge together the dispersion treatment and controls
changeInDispersion <- merge(dispersionControls, dispersionTreatments, by=c("expt.year", "site", "cal.year", "trt.year", "expt"))

#get the response ratio between trt and ctl dispersion
changeInDispersion$dist_RR <- with(changeInDispersion, (dispersionTreatment - dispersionControl)/dispersionControl)

######################################################################################################
###just picking one year, the final year
finalYear <- aggregate(changeInMean["trt.year"], by=changeInMean[c("label", "plot_mani.x", "site")], FUN=max)

finalMean <- merge(changeInMean, finalYear, by=c("label", "trt.year", "plot_mani.x", "site"))

finalDispersion <- merge(changeInDispersion, finalYear, by=c("label", "trt.year", "plot_mani.x", "site"))

#get positive and negative dispersion responses only
finalDispersionPositive <- subset(finalDispersion, dist_RR>0)
finalDispersionNegative <- subset(finalDispersion, dist_RR<=0)

#################################################################################################
#merge basic exp info with mean and variance data
finalMeanPressInfo <- merge(finalMean, expInfo, by="expt")
# finalDivergePressInfo <- merge(finalDispersionNegative, expInfo, by="expt")
# finalConvergePressInfo <- merge(finalDispersionPositive, expInfo, by="expt")
finalDispersionInfo <- merge(finalDispersion, expInfo, by="expt")
finalDispersionInfo$con_div <- ifelse(finalDispersionInfo$dist_RR<=0, 'div', 'con')

#by manipulation type
waterMean <- subset(finalMeanPressInfo, subset=(water==1))
waterMean$mani_type <- "water"
waterDiv <- subset(finalDispersionInfo, subset=(water.x==1 & plot_mani.x!=0))
waterDiv$mani_type <- "water"
nutsMean <- subset(finalMeanPressInfo, subset=(nutrients==1))
nutsMean$mani_type <- "nutrients"
nutsDiv <- subset(finalDispersionInfo, subset=(nutrients.x==1 & plot_mani.x!=0))
nutsDiv$mani_type <- "nutrients"
carbonMean <- subset(finalMeanPressInfo, subset=(carbon==1))
carbonMean$mani_type <- "carbon"
carbonDiv <- subset(finalDispersionInfo, subset=(carbon.x==1 & plot_mani.x!=0))
carbonDiv$mani_type <- "carbon"
meanManiType <- rbind(waterMean, nutsMean, carbonMean)
dispersionManiType <- rbind(waterDiv, nutsDiv, carbonDiv)
dispersionManiType$broad.ecosystem.type <- as.character(dispersionManiType$broad.ecosystem.type)


shapiro.test(meanManiType$dist)
qqnorm(meanManiType$dist)
meanManiType$dist_log <- log(meanManiType$dist)
shapiro.test(meanManiType$dist_log)
qqnorm(meanManiType$dist_log)
meanManiType$dist_sqrt <- sqrt(meanManiType$dist)
shapiro.test(meanManiType$dist_sqrt)
qqnorm(meanManiType$dist_sqrt)

#################################################################################################
#linear model

summary(meansMultipleRegression <- lm(dist_sqrt ~ mani_type + MAP.x + ANPP.x + species_num.x + broad.ecosystem.type.y, data=meanManiType))
anova(meansMultipleRegression)

summary(dispersionMultipleRegression <- lm(dist_RR ~ mani_type + MAP + ANPP + species_num + broad.ecosystem.type, data=subset(dispersionManiType, dist_RR<3)))
anova(dispersionMultipleRegression)

summary(divergeMultipleRegression <- lm(dist_RR ~ mani_type + MAP + ANPP + species_num + broad.ecosystem.type, data=subset(dispersionManiType, con_div=='div' & dist_RR<3)))
anova(divergeMultipleRegression)

summary(convergeMultipleRegression <- lm(dist_RR ~ mani_type + MAP + ANPP + species_num + broad.ecosystem.type, data=subset(dispersionManiType, con_div=='con' & dist_RR<3)))
anova(convergeMultipleRegression)

#################################################################################################
#means figures
summary(meansMani <- aov(dist_sqrt~mani_type, meanManiType))
TukeyHSD(meansMani)
meanManiPlot <- ggplot(barGraphStats(data=meanManiType, variable="dist_sqrt", byFactorNames=c("mani_type")), aes(x=mani_type, y=mean)) +
  geom_bar(stat="identity", fill='white', colour='black') +
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se), width=0.2)) +
  scale_y_continuous(name="Change in Mean", breaks=seq(0,0.6,0.2)) +
  coord_cartesian(ylim=c(0,0.6)) +
  xlab("Manipulation Type") +
  annotate('text', x=1 ,y=0.57, label='ab', size=7) +
  annotate('text', x=2 ,y=0.55, label='a', size=7) +
  annotate('text', x=3 ,y=0.58, label='b', size=7)
summary(meansMAP <- lm(MAP.x~dist_RR, subset(dispersionManiType, dist_RR<3)))
meanMAPplot <- ggplot(subset(dispersionManiType, dist_RR<3), aes(x=MAP.x, y=dist_RR)) +
  geom_point(shape=1) +
  scale_x_continuous(name="Mean Annual Precipitation") +
  scale_y_continuous(breaks=seq(0,1,0.1), name=element_blank()) +
  coord_cartesian(ylim=c(0,1))
summary(meansANPP <- lm(ANPP.x~dist_RR, subset(dispersionManiType, dist_RR<3)))
meansANPPText <- grobTree(textGrob(expression(paste(R^2, " = 0.092, p < 0.001", sep="")), x=0.03, y=0.95, hjust=0, gp=gpar(col="black", fontsize=15)))
meanANPPplot <- ggplot(subset(dispersionManiType, dist_RR<3), aes(x=ANPP.x, y=dist_RR)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  scale_x_continuous(breaks=seq(0,1000,200), name="ANPP") +
  scale_y_continuous(breaks=seq(0,1,0.1), name=element_blank()) +
  coord_cartesian(ylim=c(0,1)) +
  annotation_custom(meansANPPText)
summary(meansEcosystem <- aov(dist_RR~broad.ecosystem.type.y, subset(dispersionManiType, dist_RR<3)))
TukeyHSD(meansEcosystem)
meanEcosystemPlot <- ggplot(barGraphStats(data=subset(dispersionManiType, dist_RR<3), variable="dist_RR", byFactorNames=c("broad.ecosystem.type.y")), aes(x=broad.ecosystem.type.y, y=mean)) +
  geom_bar(stat="identity", fill='white', colour='black') +
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se), width=0.2)) +
  scale_y_continuous(name="Change in Mean", breaks=seq(0,0.6,0.2)) +
  coord_cartesian(ylim=c(0,0.6)) +
  xlab("Ecosystem Type") +
  annotate('text', x=1 ,y=0.52, label='ab', size=7) +
  annotate('text', x=2 ,y=0.59, label='a', size=7) +
  annotate('text', x=3 ,y=0.40, label='b', size=7) +
  annotate('text', x=4 ,y=0.53, label='a', size=7) +
  annotate('text', x=5 ,y=0.55, label='ab', size=7)
summary(meansSpp <- lm(species_num.x~dist_RR, subset(dispersionManiType, dist_RR<3)))
meanSppPlot <- ggplot(subset(dispersionManiType, dist_RR<3), aes(x=species_num.x, y=dist_RR)) +
  geom_point(shape=1) +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,1,0.1), name=element_blank()) +
  coord_cartesian(ylim=c(0,1))
pushViewport(viewport(layout=grid.layout(2,3))) 
print(meanManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(meanMAPplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(meanANPPplot, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(meanEcosystemPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(meanSppPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))


#################################################################################################
#dispersion figures
color <- c('#339900', '#CC0000')

summary(dispMani <- aov(dist_RR~mani_type, subset(dispersionManiType, dist_RR<3)))
TukeyHSD(dispMani)
dispManiPlot <- ggplot(barGraphStats(data=subset(dispersionManiType, dist_RR<3), variable="dist_RR", byFactorNames=c("mani_type", "con_div")), aes(x=mani_type, y=mean, fill=con_div)) +
  geom_bar(stat="identity", position=position_dodge(0)) +
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se)), width=0.2, position=position_dodge(0)) +
  scale_y_continuous(name="Relative Difference in Dispersion\nbetween Treatment and Control") +
  #coord_cartesian(ylim=c(0,0.6)) +
  xlab("Manipulation Type") +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=color) +
  annotate('text', x=1 ,y=0.18, label='a', size=7) +
  annotate('text', x=2 ,y=0.31, label='b', size=7) +
  annotate('text', x=3 ,y=0.32, label='b', size=7) +
  theme(legend.position='none')
summary(conMAP <- lm(MAP~dist_RR, subset(dispersionManiType, dist_RR<3 & con_div=='con')))
summary(divMAP <- lm(MAP~dist_RR, subset(dispersionManiType, dist_RR<3 & con_div=='div')))
dispMAPplot <- ggplot(subset(dispersionManiType, dist_RR<3), aes(x=MAP, y=dist_RR, colour=con_div)) +
  geom_point(shape=1) +
  scale_x_continuous(name="Mean Annual Precipitation") +
  scale_y_continuous(name=element_blank()) +
  #coord_cartesian(ylim=c(0,1)) +
  scale_colour_manual(values=color) +
  theme(legend.position='none')
summary(conANPP <- lm(ANPP~dist_RR, subset(dispersionManiType, dist_RR<3 & con_div=='con')))
summary(divANPP <- lm(ANPP~dist_RR, subset(dispersionManiType, dist_RR<3 & con_div=='div')))
conANPPText <- grobTree(textGrob(expression(paste(R^2, " = 0.54, p = 0.012", sep="")), x=700, y=0.95, hjust=0, gp=gpar(col="black", fontsize=15)))
dispANPPplot <- ggplot(subset(dispersionManiType, dist_RR<3), aes(x=ANPP, y=dist_RR, colour=con_div)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="#339900") +
  scale_x_continuous(breaks=seq(0,1000,200), name="ANPP") +
  scale_y_continuous(name=element_blank()) +
  #coord_cartesian(ylim=c(0,1)) +
  scale_colour_manual(values=color) +
  theme(legend.position='none') +
  annotation_custom(meansANPPText)
summary(conEcosystem <- aov(dist_RR~broad.ecosystem.type, subset(dispersionManiType, dist_RR<3 & con_div=='con')))
summary(divEcosystem <- aov(dist_RR~broad.ecosystem.type, subset(dispersionManiType, dist_RR<3 & con_div=='div')))
dispEcosystemPlot <- ggplot(barGraphStats(data=subset(dispersionManiType, dist_RR<3), variable="dist_RR", byFactorNames=c("broad.ecosystem.type", "con_div")), aes(x=broad.ecosystem.type, y=mean, fill=con_div)) +
  geom_bar(stat="identity", position=position_dodge(0)) +
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se)), width=0.2, position=position_dodge(0)) +
  scale_y_continuous(name="Relative Difference in Dispersion\nbetween Treatment and Control") +
  #coord_cartesian(ylim=c(0,0.6)) +
  xlab("Ecosystem Type") +
  geom_hline(aes(yintercept=0)) +
  scale_fill_manual(values=color) +
  theme(legend.position='none')
summary(conSpp <- lm(species_num~dist_RR, subset(dispersionManiType, dist_RR<3 & con_div=='con')))
summary(divSpp <- lm(species_num~dist_RR, subset(dispersionManiType, dist_RR<3 & con_div=='div')))
dispSppPlot <- ggplot(subset(dispersionManiType, dist_RR<3), aes(x=species_num, y=dist_RR, colour=con_div)) +
  geom_point(shape=1) +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(name=element_blank()) +
  #coord_cartesian(ylim=c(0,1)) +
  scale_colour_manual(values=color) +
  theme(legend.position='none')
pushViewport(viewport(layout=grid.layout(2,3))) 
print(dispManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(dispMAPplot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(dispANPPplot, vp=viewport(layout.pos.row=1, layout.pos.col=3))
print(dispEcosystemPlot, vp=viewport(layout.pos.row=2, layout.pos.col=1))
print(dispSppPlot, vp=viewport(layout.pos.row=2, layout.pos.col=2))



























#plots with number of manipulations

# finalDispersionPressWithControlNegativeLength <- ggplot(finalDivergePressInfo, aes(x=species_num, y=dist_RR)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T, color="black") +
#   ggtitle("Divergence") +
#   scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
#   scale_y_continuous(breaks=seq(-1.0,0.0,0.2), name="Dispersion") +
#   #coord_cartesian(ylim=c(0,0.6)) +
#   annotation_custom(divergeRegText)
finalDispersionPressSpp <- ggplot(finalDispersionInfo, aes(x=species_num, y=dist_RR, colour=con_div)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, aes(y=dist_RR, colour=con_div)) +
  ggtitle("Change in Dispersion") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(-1.0,1.5,0.25), name="Relative Difference in Dispersion\nbetween Treatment and Control") +
  coord_cartesian(ylim=c(-1.0,1.5)) +
  annotation_custom(divergeRegText) +
  annotation_custom(convergeRegText)
# finalDispersionPressWithControlPositiveLength <- ggplot(finalConvergePressInfo, aes(x=species_num, y=dist_RR)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T, color="black") +
#   ggtitle("Convergence") +
#   scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
#   scale_y_continuous(breaks=seq(0,2,0.5), name="Dispersion") +
#   #coord_cartesian(ylim=c(0,0.6)) +
#   annotation_custom(convergeRegText)
pushViewport(viewport(layout=grid.layout(1,2))) 
print(finalMeanPressSpp, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(finalDispersionPressSpp, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(finalDispersionPressWithControlPositiveLength, vp=viewport(layout.pos.row=1, layout.pos.col=3))

#regessions with ANPP
summary(meansFit <- lm(ANPP.x~dist, finalMeanPressInfo))
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.002, p = 0.476", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
summary(divergeFit <- lm(ANPP~dist_RR, finalDivergePressInfo))
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.005, p = 0.380", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
summary(convergeFit <- lm(ANPP~dist_RR, finalConvergePressInfo))
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.002, p = 0.524", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))

#plots with ANPP
finalMeanPressANPP <- ggplot(finalMeanPressInfo, aes(x=ANPP.y, y=dist)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Change in Mean") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,1,0.1), name="Change in Mean") +
  coord_cartesian(ylim=c(0,1)) +
  annotation_custom(meansRegText)
finalDispersionPressWithControlPositiveANPP <- ggplot(finalDivergePressInfo, aes(x=ANPP, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Divergence") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(divergeRegText)
finalDispersionPressWithControlNegativeANPP <- ggplot(finalConvergePressInfo, aes(x=ANPP, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Convergence") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(convergeRegText)
pushViewport(viewport(layout=grid.layout(1,3))) 
print(finalMeanPressANPP, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(finalDispersionPressWithControlPositiveANPP, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(finalDispersionPressWithControlNegativeANPP, vp=viewport(layout.pos.row=1, layout.pos.col=3))

#regessions with MAP
meansFit <- lm(MAP.x~dist, finalMeanPressInfo)
summary(meansFit)
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.007, p = 0.103", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
summary(divergeFit <- lm(MAP~dist_RR, finalDivergePressInfo))
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.001, p = 0.600", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
summary(convergeFit <- lm(MAP~dist_RR, finalConvergePressInfo))
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " < 0.001, p = 0.91", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))

# #plots with number of manipulations
# finalMeanPressMAP <- ggplot(finalMeanPressInfo, aes(x=MAP.x, y=dist)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T, color="black") +
#   ggtitle("Change in Mean") +
#   scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
#   scale_y_continuous(breaks=seq(0,1,0.1), name="Change in Mean") +
#   coord_cartesian(ylim=c(0,1)) +
#   annotation_custom(meansRegText)
# finalDispersionPressWithControlPositiveMAP <- ggplot(finalDivergePressInfo, aes(x=MAP, y=disp)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T, color="black") +
#   ggtitle("Divergence") +
#   scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
#   scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
#   coord_cartesian(ylim=c(0,0.6)) +
#   annotation_custom(divergeRegText)
# finalDispersionPressWithControlNegativeMAP <- ggplot(finalConvergePressInfo, aes(x=MAP, y=disp)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T, color="black") +
#   ggtitle("Convergence") +
#   scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
#   scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
#   coord_cartesian(ylim=c(0,0.6)) +
#   annotation_custom(convergeRegText)
# pushViewport(viewport(layout=grid.layout(1,3))) 
# print(finalMeanPressMAP, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionPressWithControlPositiveMAP, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(finalDispersionPressWithControlNegativeMAP, vp=viewport(layout.pos.row=1, layout.pos.col=3))


#################################################################################################
#plotting by manipulation type

divergeManiPlot <- ggplot(barGraphStats(data=divergeManiType, variable="disp", byFactorNames=c("mani_type")), aes(x=mani_type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.5, 0.1)), name="Dispersion") +
  xlab("Manipulation Type") +
  coord_cartesian(ylim=c(0,0.4)) +
  ggtitle("Divergence")
convergeManiPlot <- ggplot(barGraphStats(data=convergeManiType, variable="disp", byFactorNames=c("mani_type")), aes(x=mani_type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.5, 0.1)), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.4)) +
  xlab("Manipulation Type") +
  ggtitle("Convergence")
pushViewport(viewport(layout=grid.layout(1,3)))
print(meanManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(divergeManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(convergeManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))

#################################################################################################
#by biome type
#plotting by manipulation type
meanManiPlot <- ggplot(barGraphStats(data=finalMeanPressInfo, variable="dist", byFactorNames=c("broad.ecosystem.type.y")), aes(x=broad.ecosystem.type.y, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.4, 0.1)), name="Change in Mean") +
  coord_cartesian(ylim=c(0,0.4)) +
  xlab("Manipulation Type") +
  ggtitle("Change in Mean")
divergeManiPlot <- ggplot(barGraphStats(data=finalDivergePressInfo, variable="disp", byFactorNames=c("broad.ecosystem.type")), aes(x=broad.ecosystem.type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.5, 0.1)), name="Dispersion") +
  xlab("Manipulation Type") +
  coord_cartesian(ylim=c(0,0.4)) +
  ggtitle("Divergence")
convergeManiPlot <- ggplot(barGraphStats(data=finalConvergePressInfo, variable="disp", byFactorNames=c("broad.ecosystem.type")), aes(x=broad.ecosystem.type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.5, 0.1)), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.4)) +
  xlab("Manipulation Type") +
  ggtitle("Convergence")
pushViewport(viewport(layout=grid.layout(1,3)))
print(meanManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(divergeManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(convergeManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))


#################################################################################################









