setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\07_16_2015')

library(vegan)
library(reshape2)
library(ggplot2)
library(gtools)
library(plyr)
library(grid)
library(lme4)

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
allData<-read.csv('dispersion_and_means_press_experiments_with_exp_info_03232015.csv')

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
                           "light"), 
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

#get the difference between trt and ctl dispersion
changeInDispersion$dispersionDifference <- changeInDispersion$dispersionTreatment - changeInDispersion$dispersionControl

######################################################################################################
###just picking one year, the final year
finalYear <- aggregate(changeInMean["trt.year"], by=changeInMean[c("label", "plot_mani.x", "site")], FUN=max)

finalMean <- merge(changeInMean, finalYear, by=c("label", "trt.year", "plot_mani.x", "site"))

finalDispersion <- merge(changeInDispersion, finalYear, by=c("label", "trt.year", "plot_mani.x", "site"))

# ###are the max years also the final years?
# yearsMean <- merge(maxMeanYear, finalMean, by=c("label", "site", "plot_mani"))
# yearsMean <- yearsMean[,-c(5:9, 11:18, 20, 22:23)]
# yearsDispersion <- merge(maxDispersion, finalDispersion, by=c("label", "site", "plot_mani"))
# yearsDispersion <- yearsDispersion[,-c(5, 7:12, 14:15, 17:21, 23)]
# 
# #remove treatments that are pulses, keeping only the press treatments (for which we want to look at the last year only) #be sure to add back in CDR e002 and SGS ESA
# yearsMeanPress <- subset(yearsMean, label!="CUL::Culardoch::0::burn" & label!="CUL::Culardoch::0::burnclip" & label!="CUL::Culardoch::0::N10burn" & label!="CUL::Culardoch::0::N10burnclip" & label!="CUL::Culardoch::0::N20burn" & label!="CUL::Culardoch::0::N20burnclip" & label!="CUL::Culardoch::0::N50burn" & label!="CUL::Culardoch::0::N50burnclip" & label!="KBS::T7::0::T1F1" & label!="KBS::T7::0::T1F0" & site!="dcgs" & site!="KAEFS")
# 
# yearsMeanPress$Equal <- ifelse(yearsMeanPress$cal.year.x==yearsMeanPress$cal.year.y, "yes", "no")
# yearsMeanPress$difference <- yearsMeanPress$dist.x - yearsMeanPress$dist.y
# summary(yearsMeanPress$difference)
# yearsMeanPress <- yearsMeanPress[order(yearsMeanPress$difference),]
# hist(yearsMeanPress$difference)
# 
# #remove treatments that are pulses, keeping only the press treatments (for which we want to look at the last year only) #be sure to add back in CDR e002 and SGS ESA
# yearsDispersionPress <- subset(yearsDispersion, label!="CUL::Culardoch::0::burn" & label!="CUL::Culardoch::0::burnclip" & label!="CUL::Culardoch::0::N10burn" & label!="CUL::Culardoch::0::N10burnclip" & label!="CUL::Culardoch::0::N20burn" & label!="CUL::Culardoch::0::N20burnclip" & label!="CUL::Culardoch::0::N50burn" & label!="CUL::Culardoch::0::N50burnclip" & label!="KBS::T7::0::T1F1" & label!="KBS::T7::0::T1F0" & site!="dcgs" & site!="KAEFS")
# yearsDispersionPress$Equal <- ifelse(yearsDispersionPress$cal.year.x==yearsDispersionPress$cal.year.y, "yes", "no")
# yearsDispersionPress$difference <- yearsDispersionPress$dispersionDifference.x - yearsDispersionPress$dispersionDifference.y
# summary(yearsDispersionPress$difference)
# yearsDispersionPress$absDifference <- abs(yearsDispersionPress$dispersionDifference.x) - abs(yearsDispersionPress$dispersionDifference.y)
# summary(yearsDispersionPress$absDifference)
# yearsDispersion <- yearsDispersionPress[order(yearsDispersionPress$absDifference),]
# hist(yearsDispersionPress$absDifference)
# hist(yearsDispersionPress$difference)
# 
# #plotting the mean and dispersion by plot_mani, with ci
# finalMeanPlotCI <- ggplot(barGraphStats(data=finalMean, variable="dist", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(breaks=seq(0,1,0.2), name="Change in Mean") +
#   coord_cartesian(ylim=c(0,1))
# finalDispersionPlotCI <- ggplot(barGraphStats(data=finalDispersion, variable="dispersionDifference", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_y_continuous(name="Dispersion") +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(finalMeanPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# 
# #plotting the mean and dispersion by plot_mani, with se
# finalMeanPlotSE <- ggplot(barGraphStats(data=finalMean, variable="dist", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(breaks=seq(0,0.6,0.2), name="Change in Mean") +
#   coord_cartesian(ylim=c(0,0.6))
# finalDispersionPlotSE <- ggplot(barGraphStats(data=finalDispersion, variable="dispersionDifference", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(name="Dispersion") +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(finalMeanPlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionPlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=2))

#same but with abs of dispersion
finalDispersion$absDispersion <- abs(finalDispersion$dispersionDifference)
# 
# finalDispersionAbsValuePlotCI <- ggplot(barGraphStats(data=finalDispersion, variable="absDispersion", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_y_continuous(name="Dispersion") +
#   scale_x_continuous(breaks=seq(1,7,1), "Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(finalMeanPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionAbsValuePlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# 
# finalDispersionAbsValuePlotSE <- ggplot(barGraphStats(data=finalDispersion, variable="absDispersion", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(name="Dispersion") +
#   scale_x_continuous(breaks=seq(1,7,1), "Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(finalMeanPlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionAbsValuePlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=2))

###plot regressions with dispersion separated into converge and diverge panels
finalDispersionPressPositive <- subset(finalDispersion, dispersionDifference>0)
finalDispersionPressNegative <- subset(finalDispersion, dispersionDifference<=0)

# finalMeanPressScatter <- ggplot(finalMean, aes(x=plot_mani.x, y=dist)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Change in Mean") +
#   ggtitle("Change in Mean")
# finalDispersionPressPositiveScatter <- ggplot(finalDispersionPressPositive, aes(x=plot_mani.x, y=dispersionDifference)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Difference in Dispersion") +
#   ggtitle("Divergence")
# finalDispersionPressNegativeScatter <- ggplot(finalDispersionPressNegative, aes(x=plot_mani.x, y=dispersionDifference)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Difference in Dispersion") +
#   ggtitle("Convergence")
# pushViewport(viewport(layout=grid.layout(1,3)))
# print(finalMeanPressScatter, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionPressPositiveScatter, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(finalDispersionPressNegativeScatter, vp=viewport(layout.pos.row=1, layout.pos.col=3))

############################################################################
#making dispersion graphs with controls as plot_mani=0
finalYearDispersionMeansPress <- aggregate(dispersionMeans["trt.year"], by=dispersionMeans[c("label", "plot_mani.x", "site")], FUN=max)

finalDispersionMeansPress <- merge(finalYearDispersionMeansPress, dispersionMeans, by=c("label", "trt.year", "plot_mani.x", "site"))

# #plot regression
# finalMeanPressScatter <- ggplot(finalMeanPress, aes(x=plot_mani, y=dist)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T) +
#   scale_y_continuous(name="Change In Mean") +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# finalDispersionPressScatter <- ggplot(finalDispersionMeansPress, aes(x=plot_mani, y=disp)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=T) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Dispersion")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(finalMeanPressScatter, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionPressScatter, vp=viewport(layout.pos.row=1, layout.pos.col=2))

#merge with dataframe that has difference between control and treatments to pick out positives vs negatives
finalDispersionPressPosNeg <- merge(finalDispersionMeansPress, finalDispersion, by=c("label", "trt.year", "cal.year", "expt", "site", "plot_mani.x"), all=TRUE)

#get positive and negative dispersion responses only
finalDispersionPressWithControlPositive <- subset(finalDispersionPressPosNeg, dispersionDifference>0 | plot_mani.x==0)
finalDispersionPressWithControlNegative <- subset(finalDispersionPressPosNeg, dispersionDifference<=0 | plot_mani.x==0)

#################################################################################################
###number of manipulations

#regessions with just number of manipulations
meansFit <- lm(plot_mani.x~dist, finalMean)
summary(meansFit)
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.236, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(plot_mani.x~disp, finalDispersionPressWithControlNegative)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.112, p < 0.001", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))
divergeFit <- lm(plot_mani.x~disp, finalDispersionPressWithControlPositive)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.128, p < 0.001", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))

#plots with number of manipulations
finalMeanPressScatter <- ggplot(finalMean, aes(x=plot_mani.x, y=dist)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Change in Mean") +
  scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
  scale_y_continuous(breaks=seq(0,1,0.1), name="Change in Mean") +
  coord_cartesian(ylim=c(0,1)) +
  annotation_custom(meansRegText)
finalDispersionPressWithControlPositiveScatter <- ggplot(finalDispersionPressWithControlPositive, aes(x=plot_mani.x, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Divergence") +
  scale_x_continuous(breaks=seq(0,7,1), name="Number of Manipulations") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(divergeRegText)
finalDispersionPressWithControlNegativeScatter <- ggplot(finalDispersionPressWithControlNegative, aes(x=plot_mani.x, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Convergence") +
  scale_x_continuous(breaks=seq(0,7,1), name="Number of Manipulations") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(convergeRegText)
pushViewport(viewport(layout=grid.layout(1,3))) 
print(finalMeanPressScatter, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(finalDispersionPressWithControlPositiveScatter, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(finalDispersionPressWithControlNegativeScatter, vp=viewport(layout.pos.row=1, layout.pos.col=3))

#################################################################################################
###length of experiment

#regessions with just number of manipulations
meansFit <- lm(trt.year~dist, finalMean)
summary(meansFit)
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.183, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(trt.year~disp, finalDispersionPressWithControlNegative)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.005, p = 0.253", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))
divergeFit <- lm(trt.year~disp, finalDispersionPressWithControlPositive)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.076, p < 0.001", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))

#plots with number of manipulations
finalMeanPressLength <- ggplot(finalMean, aes(x=trt.year, y=dist)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Change in Mean") +
  scale_x_continuous(breaks=seq(0,30,5), name="Experiment Length") +
  scale_y_continuous(breaks=seq(0,1,0.1), name="Change in Mean") +
  coord_cartesian(ylim=c(0,1)) +
  annotation_custom(meansRegText)
finalDispersionPressWithControlPositiveLength <- ggplot(finalDispersionPressWithControlPositive, aes(x=trt.year, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Divergence") +
  scale_x_continuous(breaks=seq(0,30,5), name="Experiment Length") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(divergeRegText)
finalDispersionPressWithControlNegativeLength <- ggplot(finalDispersionPressWithControlNegative, aes(x=trt.year, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Convergence") +
  scale_x_continuous(breaks=seq(0,30,5), name="Experiment Length") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(convergeRegText)
pushViewport(viewport(layout=grid.layout(1,3))) 
print(finalMeanPressLength, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(finalDispersionPressWithControlPositiveLength, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(finalDispersionPressWithControlNegativeLength, vp=viewport(layout.pos.row=1, layout.pos.col=3))

#################################################################################################
#merge basic exp info with mean and variance data
finalMeanPressInfo <- merge(finalMean, expInfo, by="expt")
finalDivergePressInfo <- merge(finalDispersionPressWithControlPositive, expInfo, by="expt")
finalConvergePressInfo <- merge(finalDispersionPressWithControlNegative, expInfo, by="expt")
#################################################################################################

###community richness

#regessions with just number of manipulations
meansFit <- lm(species_num.x~dist, finalMeanPressInfo)
summary(meansFit)
AIC(meansFit)
# meansFitExp <- lm(species_num~(-dist^2), finalMeanPressInfo)
# summary(meansFitExp)
# AIC(meansFitExp)
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.057, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
divergeFit <- lm(species_num~disp, finalDivergePressInfo)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.217, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(species_num~disp, finalConvergePressInfo)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.162, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))

#plots with number of manipulations
finalMeanPressLength <- ggplot(finalMeanPressInfo, aes(x=species_num.x, y=dist)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Change in Mean") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,1,0.1), name="Change in Mean") +
  coord_cartesian(ylim=c(0,1)) +
  annotation_custom(meansRegText)
finalDispersionPressWithControlPositiveLength <- ggplot(finalDivergePressInfo, aes(x=species_num, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Divergence") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(divergeRegText)
finalDispersionPressWithControlNegativeLength <- ggplot(finalConvergePressInfo, aes(x=species_num, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Convergence") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(convergeRegText)
pushViewport(viewport(layout=grid.layout(1,3))) 
print(finalMeanPressLength, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(finalDispersionPressWithControlPositiveLength, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(finalDispersionPressWithControlNegativeLength, vp=viewport(layout.pos.row=1, layout.pos.col=3))

#regessions with ANPP
meansFit <- lm(ANPP.x~dist, finalMeanPressInfo)
summary(meansFit)
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.002, p = 0.476", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
divergeFit <- lm(ANPP~disp, finalDivergePressInfo)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.012, p = 0.084", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(ANPP~disp, finalConvergePressInfo)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.013, p = 0.082", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))

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
divergeFit <- lm(MAP~disp, finalDivergePressInfo)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.001, p = 0.600", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(MAP~disp, finalConvergePressInfo)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " < 0.001, p = 0.91", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))

#plots with number of manipulations
finalMeanPressMAP <- ggplot(finalMeanPressInfo, aes(x=MAP.x, y=dist)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Change in Mean") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,1,0.1), name="Change in Mean") +
  coord_cartesian(ylim=c(0,1)) +
  annotation_custom(meansRegText)
finalDispersionPressWithControlPositiveMAP <- ggplot(finalDivergePressInfo, aes(x=MAP, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Divergence") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(divergeRegText)
finalDispersionPressWithControlNegativeMAP <- ggplot(finalConvergePressInfo, aes(x=MAP, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Convergence") +
  scale_x_continuous(breaks=seq(0,185,50), name="Gamma Diversity") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(convergeRegText)
pushViewport(viewport(layout=grid.layout(1,3))) 
print(finalMeanPressMAP, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(finalDispersionPressWithControlPositiveMAP, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(finalDispersionPressWithControlNegativeMAP, vp=viewport(layout.pos.row=1, layout.pos.col=3))


#################################################################################################
#by manipulation type
waterMean <- subset(finalMeanPressInfo, subset=(water==1))
waterMean$mani_type <- "water"
waterDiv <- subset(finalDivergePressInfo, subset=(water==1 & plot_mani.x!=0))
waterDiv$mani_type <- "water"
waterCon <- subset(finalConvergePressInfo, subset=(water==1 & plot_mani.x!=0))
waterCon$mani_type <- "water"
nutsMean <- subset(finalMeanPressInfo, subset=(nutrients==1))
nutsMean$mani_type <- "nutrients"
nutsDiv <- subset(finalDivergePressInfo, subset=(nutrients==1 & plot_mani.x!=0))
nutsDiv$mani_type <- "nutrients"
nutsCon <- subset(finalConvergePressInfo, subset=(nutrients==1 & plot_mani.x!=0))
nutsCon$mani_type <- "nutrients"
carbonMean <- subset(finalMeanPressInfo, subset=(carbon==1))
carbonMean$mani_type <- "carbon"
carbonDiv <- subset(finalDivergePressInfo, subset=(carbon==1 & plot_mani.x!=0))
carbonDiv$mani_type <- "carbon"
carbonCon <- subset(finalConvergePressInfo, subset=(carbon==1 & plot_mani.x!=0))
carbonCon$mani_type <- "carbon"
meanManiType <- rbind(waterMean, nutsMean, carbonMean)
divergeManiType <- rbind(waterDiv, nutsDiv, carbonDiv)
divergeManiType$broad.ecosystem.type <- as.character(divergeManiType$broad.ecosystem.type)
convergeManiType <- rbind(waterCon, nutsCon, carbonCon)
convergeManiType$broad.ecosystem.type <- as.character(convergeManiType$broad.ecosystem.type)

#plotting by manipulation type
meanManiPlot <- ggplot(barGraphStats(data=meanManiType, variable="dist", byFactorNames=c("mani_type")), aes(x=mani_type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(1.96*se), ymax=mean+(1.96*se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.4, 0.1)), name="Change in Mean") +
  coord_cartesian(ylim=c(0,0.4)) +
  xlab("Manipulation Type") +
  ggtitle("Change in Mean")
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
#mixed effects model
###need to add in manipulation type

summary(meansMultipleRegression <- lm(dist ~ plot_mani.x + mani_type + trt.year + MAP.x + ANPP.x + species_num.x + broad.ecosystem.type.y, data=meanManiType))
anova(meansMultipleRegression)

summary(divergeMultipleRegression <- lm(disp ~ plot_mani.x + mani_type + trt.year + MAP + ANPP + species_num + broad.ecosystem.type, data=divergeManiType))
anova(divergeMultipleRegression)

summary(convergeMultipleRegression <- lm(disp ~ plot_mani.x + mani_type + trt.year + MAP + ANPP + species_num + broad.ecosystem.type, data=convergeManiType))
anova(convergeMultipleRegression)











