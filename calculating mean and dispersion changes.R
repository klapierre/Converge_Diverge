library(vegan)
library(reshape2)
library(gtools)
library(plyr)
library(grid)
library(tidyr)
library(dplyr)

#kim
setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

# #meghan
# setwd("~/Dropbox/converge_diverge/datasets/LongForm")

#read in the merged dataset
alldata<-read.csv("SpeciesRelativeAbundance_11232015.csv")%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))

<<<<<<< HEAD
expinfo<-read.csv("ExperimentInformation_11202015.csv")%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))%>%
  select(exp_year, plot_mani, treatment)

alldata2a<-merge(alldata, expinfo, by=c("exp_year","treatment"), all=F)
=======
expinfo<-read.csv("ExperimentInformation_Nov2015.csv")%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))%>%
  select(exp_year, plot_mani, treatment)

alldata2<-merge(alldata, expinfo, by=c("exp_year","treatment"),all=T)%>%
  filter(project_name!="e001"&project_name!="e002")#there are known problems with cdr_e001
>>>>>>> bed5ea4202e6bcd396c0a4be89ca8ca456c71d9d

#ESA has dozens of duplicate species
alldata2<-subset(alldata2a, project_name!='ESA')


#make a new dataframe with just the label;
exp_year=alldata2%>%
  select(exp_year)%>%
  unique()

#makes an empty dataframe
for.analysis=data.frame(row.names=1) 


###first, get bray curtis dissimilarity values for each year within each experiment between all combinations of plots
###second, get distance of each plot within a trt to the trt centroid 
###third: mean_change is the distance between trt and control centriods
####fourth: dispersion is the average dispersion of plots within a treatment to treatment centriod
for(i in 1:length(exp_year$exp_year)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=alldata2[alldata2$exp_year==as.character(exp_year$exp_year[i]),]%>%
    select(exp_year, treatment, plot_mani, genus_species, relcov, plot_id)
  
  #need this to keep track of plot mani
  labels=subset%>%
    select(plot_mani, treatment)%>%
    unique()
  
  #transpose data
  species=subset%>%
    spread(genus_species, relcov, fill=0)

  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,5:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances among treatment centroids; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #####there is an error here that it isn't giving us the distances, just the treatment names
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(exp_year=exp_year$exp_year[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="mean_change" 
  
  #merging back with labels to get back plot_mani
  centroid=merge(cent_C_T, labels, by="treatment")
  
  #collecting and labeling distances to centroid from betadisper
  trt_disp=data.frame(data.frame(exp_year=exp_year$exp_year[i], 
                                              plot_id=species$plot_id, 
                                              treatment=species$treatment,
                                              dist=disp$distances))#%>%
#     tbl_df%>%
#     group_by(exp_year, treatment)%>%
#     summarize(dispersion=mean(dist))
  
  distances<-merge(centroid, trt_disp, by=c("exp_year","treatment"))
  
#   #getting diversity indixes
#     H<-diversity(species[,5:ncol(species)])
#     S<-specnumber(species[,5:ncol(species)])
#     InvD<-diversity(species[,5:ncol(species)],"inv")
#     SimpEven<-InvD/S
#     out1<-cbind(H, S)
#     output<-cbind(out1, SimpEven)
#     divmeasure<-cbind(species, output)%>%
#       select(exp_year, treatment, H, S, SimpEven)
#     
#     ##merging all measures of diversity
#     alldiv<-merge(distances, divmeasure, by=c("exp_year","treatment"))

      #pasting dispersions into the dataframe made for this analysis
#     for.analysis=rbind(alldiv, for.analysis)  
}

#####################################################

#making column to label whether the comparison is for dispersions or means/centroids
for(i in 1:length(for.analysis$dist.1)) if(for.analysis$dist.1[i]==for.analysis$dist.2[i]) for.analysis$mean.disp[i]="disp" else for.analysis$mean.disp[i]="mean" 
expt.year.m=matrix(unlist(strsplit(levels(expt.year.list$expt.year), "::")), ncol=5, byrow=T)
expt.labels=data.frame(labels=expt.year.list, 
                       site=expt.year.m[,1], 
                       cal.year=expt.year.m[,4], 
                       trt.year=as.numeric(expt.year.m[,5]), 
                       expt=as.factor(paste(expt.year.m[,1], expt.year.m[,2], expt.year.m[,3], sep="::"))) #adding labels
for.analysis2=merge(for.analysis, expt.labels, by="expt.year")

# #write the dissimilarity data to a csv file
# write.csv(for.analysis2, "C:\\Users\\Kim\\Dropbox\\converge_diverge\\datasets\\FINAL_SEPT2014\\R files\\analysis files\\test_results.csv", row.names=FALSE)

#sets for.analysis2 (the bray curtis results) as a new dataframe named df
df <- for.analysis2

#creates column dg, which indicates what treatments are being compared (ctl vs ctl=c_c, ctl vs trt=c_t, trt vs trt=t_t)
df$dg <- NA
df$dg[df$mean.disp=="mean"] <- "c_t"
df$dg[df$mean.disp=="disp" & df$plot_mani==0] <- "c_c"
df$dg[df$mean.disp=="disp" & df$plot_mani>0] <- "t_t"
df$dg <- as.factor(df$dg)

#create a label with the site, exp, community, and trt 
df$label=as.factor(paste(df$expt, df$dist.2, sep="::"))

###picking out only the press experiments, for now; discuss how to include these pulse treatments with the group!
dfPress<-subset(df, subset=(label!="CUL::Culardoch::0::burn" & label!="CUL::Culardoch::0::burnclip" & label!="CUL::Culardoch::0::N10burn" & label!="CUL::Culardoch::0::N10burnclip" & label!="CUL::Culardoch::0::N20burn" & label!="CUL::Culardoch::0::N20burnclip" & label!="CUL::Culardoch::0::N50burn" & label!="CUL::Culardoch::0::N50burnclip" & label!="KBS::T7::0::T1F1" & label!="KBS::T7::0::T1F0" & site!="dcgs" & label!="ASGA::Exp1::a::2_0_PA" & label!="ASGA::Exp1::a::2_0_UN" & label!="ASGA::Exp1::a::2_0_CO" & label!="ASGA::Exp1::a::2_1_PA" & label!="ASGA::Exp1::a::2_1_UN" & label!="ASGA::Exp1::a::2_1_CO"))

# write.csv(dfPress, 'dfpresstest.csv')

#combine with exp info x2
dfPressInfo1<-merge(dfPress, expInfo, by="label")
dfPressInfoFinal<-merge(dfPressInfo1, expSiteInfo, by=c("site_code", "project_name", "community_type"))

###EXPORT dispersion data
write.csv(dfPressInfoFinal, 'dispersion_and_means_press_experiments_with_exp_info_08062015.csv')


#final year only
finalYear <- aggregate(dfPressInfoFinal["trt.year"], by=dfPressInfoFinal[c("expt")], FUN=max)
finalYearPressInfoFinal <- merge(finalYear, dfPressInfoFinal, by=c("expt", "trt.year"))
write.csv(finalYearPressInfoFinal, "dispersion and means_FINAL YEAR_press experiments_with exp info_08062015.csv")

############################################################################################
#subset out only mean comparison
changeInMean <- subset(dfPress, mean.disp=="mean")

#subset out only dispersion comparison for controls and treatments
dispersion <- subset(dfPress, mean.disp=="disp")

#means of dispersion across all plots within a year
dispersionMeans <- ddply(dispersion, 
                         c("expt.year", "plot_mani", "dist.1", "mean.disp", "site",
                           "cal.year", "trt.year", "expt", "dg", "label"), 
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

############################################################################################
# #just picking one year, the maximum for mean
# maxMean <- aggregate(changeInMean["dist"], by=changeInMean[c("label", "plot_mani", "site")], FUN=max)
# #merge maxMean back with original dataset to get what year the max mean change occured in
# maxMeanYear <- merge(maxMean, changeInMean, by=c("dist", "label", "plot_mani", "site"))
# 
# #just picking one year, the maximum for dispersion
# #getting the absolute values of the dispersion
# changeInDispersion$absDispersion <- abs(changeInDispersion$dispersionDifference)
# #finding the max of the absolute values of dispersion
# maxAbsDispersion <- aggregate(changeInDispersion["absDispersion"], by=changeInDispersion[c("label", "plot_mani", "site")], FUN=max)
# #merging back with the main diseprsion data
# maxDispersion <- merge(changeInDispersion, maxAbsDispersion, by=c("label", "plot_mani", "site", "absDispersion"))
# 
# #plotting the mean and dispersion by plot_mani, with CI
# maxMeanPlotCI <- ggplot(barGraphStats(data=maxMean, variable="dist", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_y_continuous(breaks=(seq(0, 1.0, 0.2)), name="Change in Mean") +
#   coord_cartesian(ylim=c(0,1.1)) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") 
# maxDispersionPlotCI <- ggplot(barGraphStats(data=maxDispersion, variable="dispersionDifference", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_y_continuous(breaks=(seq(-0.5, 0.5, 0.25)), name="Dispersion") +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(maxMeanPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(maxDispersionPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# 
# #plotting the mean and dispersion by plot_mani, with se
# maxMeanPlotSE <- ggplot(barGraphStats(data=maxMean, variable="dist", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=(seq(0, 1.0, 0.2)), name="Change in Mean") +
#   coord_cartesian(ylim=c(0,0.8)) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# maxDispersionPlotSE <- ggplot(barGraphStats(data=maxDispersion, variable="dispersionDifference", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=(seq(-0.1, 0.2, 0.05)), name="Dispersion") +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(maxMeanPlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(maxDispersionPlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# 
# #same but with abs of dispersion
# maxDispersion$absDispersion <- abs(maxDispersion$dispersionDifference)
# 
# maxDispersionAbsValuePlotCI <- ggplot(barGraphStats(data=maxDispersion, variable="absDispersion", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_y_continuous(breaks=(seq(-0.1, 0.35, 0.05)), name="Dispersion") +
#   coord_cartesian(ylim=c(0,0.35)) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(maxMeanPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(maxDispersionAbsValuePlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# 
# maxDispersionAbsValuePlotSE <- ggplot(barGraphStats(data=maxDispersion, variable="absDispersion", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
#   scale_y_continuous(breaks=(seq(-0.1, 0.25, 0.05)), name="Dispersion") +
#   coord_cartesian(ylim=c(0,0.25)) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations")
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(maxMeanPlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(maxDispersionAbsValuePlotSE, vp=viewport(layout.pos.row=1, layout.pos.col=2))

######################################################################################################
###just picking one year, the final year
finalYear <- aggregate(changeInMean["trt.year"], by=changeInMean[c("label", "plot_mani", "site")], FUN=max)

finalMean <- merge(changeInMean, finalYear, by=c("label", "trt.year", "plot_mani", "site"))

finalDispersion <- merge(changeInDispersion, finalYear, by=c("label", "trt.year", "plot_mani", "site"))


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

#same but with the pulse experiments gone
# finalMeanPress <- subset(finalMean, label!="CUL::Culardoch::0::burn" & label!="CUL::Culardoch::0::burnclip" & label!="CUL::Culardoch::0::N10burn" & label!="CUL::Culardoch::0::N10burnclip" & label!="CUL::Culardoch::0::N20burn" & label!="CUL::Culardoch::0::N20burnclip" & label!="CUL::Culardoch::0::N50burn" & label!="CUL::Culardoch::0::N50burnclip" & label!="KBS::T7::0::T1F1" & label!="KBS::T7::0::T1F0" & site!="dcgs" & site!="KAEFS")
# 
# finalDispersionPress <- subset(finalDispersion, label!="CUL::Culardoch::0::burn" & label!="CUL::Culardoch::0::burnclip" & label!="CUL::Culardoch::0::N10burn" & label!="CUL::Culardoch::0::N10burnclip" & label!="CUL::Culardoch::0::N20burn" & label!="CUL::Culardoch::0::N20burnclip" & label!="CUL::Culardoch::0::N50burn" & label!="CUL::Culardoch::0::N50burnclip" & label!="KBS::T7::0::T1F1" & label!="KBS::T7::0::T1F0" & site!="dcgs" & site!="KAEFS")
# 
# finalMeanPressPlotCI <- ggplot(barGraphStats(data=finalMeanPress, variable="dist", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Change in Mean") +
#   coord_cartesian(ylim=c(0,0.9))
# finalDispersionPressPlotCI <- ggplot(barGraphStats(data=finalDispersionPress, variable="absDispersion", byFactorNames=c("plot_mani")), aes(x=plot_mani, y=mean)) +
#   geom_bar(stat="identity") +
#   geom_errorbar(aes(ymin=mean-(1.96*sd), ymax=mean+(1.96*sd), width=0.2)) +
#   scale_y_continuous(name="Dispersion") +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   coord_cartesian(ylim=c(0,0.22))
# pushViewport(viewport(layout=grid.layout(1,2)))
# print(finalMeanPressPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionPressPlotCI, vp=viewport(layout.pos.row=1, layout.pos.col=2))

###plot regressions with dispersion separated into converge and diverge panels
finalDispersionPressPositive <- subset(finalDispersion, dispersionDifference>0)
finalDispersionPressNegative <- subset(finalDispersion, dispersionDifference<=0)

# finalMeanPressScatter <- ggplot(finalMeanPress, aes(x=plot_mani, y=dist)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=F) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Change in Mean") +
#   ggtitle("Change in Mean")
# finalDispersionPressPositiveScatter <- ggplot(finalDispersionPressPositive, aes(x=plot_mani, y=dispersionDifference)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=F) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Difference in Dispersion") +
#   ggtitle("Divergence")
# finalDispersionPressNegativeScatter <- ggplot(finalDispersionPressNegative, aes(x=plot_mani, y=dispersionDifference)) +
#   geom_point(shape=1) +
#   geom_smooth(method=lm, se=F) +
#   scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
#   scale_y_continuous(name="Difference in Dispersion") +
#   ggtitle("Convergence")
# pushViewport(viewport(layout=grid.layout(1,3)))
# print(finalMeanPressScatter, vp=viewport(layout.pos.row=1, layout.pos.col=1))
# print(finalDispersionPressPositiveScatter, vp=viewport(layout.pos.row=1, layout.pos.col=2))
# print(finalDispersionPressNegativeScatter, vp=viewport(layout.pos.row=1, layout.pos.col=3))

############################################################################
#making dispersion graphs with controls as plot_mani=0
# dispersionMeansPress <- subset(dispersionMeans, label!="CUL::Culardoch::0::burn" & label!="CUL::Culardoch::0::burnclip" & label!="CUL::Culardoch::0::N10burn" & label!="CUL::Culardoch::0::N10burnclip" & label!="CUL::Culardoch::0::N20burn" & label!="CUL::Culardoch::0::N20burnclip" & label!="CUL::Culardoch::0::N50burn" & label!="CUL::Culardoch::0::N50burnclip" & label!="KBS::T7::0::T1F1" & label!="KBS::T7::0::T1F0" & site!="dcgs" & label!="ASGA::Exp1::a::2_0_PA" & label!="ASGA::Exp1::a::2_0_UN" & label!="ASGA::Exp1::a::2_0_CO" & label!="ASGA::Exp1::a::2_1_PA" & label!="ASGA::Exp1::a::2_1_UN" & label!="ASGA::Exp1::a::2_1_CO")

finalYearDispersionMeansPress <- aggregate(dispersionMeans["trt.year"], by=dispersionMeans[c("label", "plot_mani", "site")], FUN=max)

finalDispersionMeansPress <- merge(finalYearDispersionMeansPress, dispersionMeans, by=c("label", "trt.year", "plot_mani", "site"))

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
finalDispersionPressPosNeg <- merge(finalDispersionMeansPress, finalDispersion, by=c("label", "trt.year", "cal.year", "expt", "site", "plot_mani"), all=TRUE)

#get positive and negative dispersion responses only
finalDispersionPressWithControlPositive <- subset(finalDispersionPressPosNeg, dispersionDifference>0 | plot_mani==0)
finalDispersionPressWithControlNegative <- subset(finalDispersionPressPosNeg, dispersionDifference<=0 | plot_mani==0)

#################################################################################################
###number of manipulations

#regessions with just number of manipulations
meansFit <- lm(plot_mani~dist, finalMean)
summary(meansFit)
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.224, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(plot_mani~disp, finalDispersionPressWithControlNegative)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.110, p < 0.001", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))
divergeFit <- lm(plot_mani~disp, finalDispersionPressWithControlPositive)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.127, p < 0.001", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))

#plots with number of manipulations
finalMeanPressScatter <- ggplot(finalMean, aes(x=plot_mani, y=dist)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Change in Mean") +
  scale_x_continuous(breaks=seq(1,7,1), name="Number of Manipulations") +
  scale_y_continuous(breaks=seq(0,1,0.1), name="Change in Mean") +
  coord_cartesian(ylim=c(0,1)) +
  annotation_custom(meansRegText)
finalDispersionPressWithControlPositiveScatter <- ggplot(finalDispersionPressWithControlPositive, aes(x=plot_mani, y=disp)) +
  geom_point(shape=1) +
  geom_smooth(method=lm, se=T, color="black") +
  ggtitle("Divergence") +
  scale_x_continuous(breaks=seq(0,7,1), name="Number of Manipulations") +
  scale_y_continuous(breaks=seq(0,0.6,0.1), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.6)) +
  annotation_custom(divergeRegText)
finalDispersionPressWithControlNegativeScatter <- ggplot(finalDispersionPressWithControlNegative, aes(x=plot_mani, y=disp)) +
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
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.177, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(trt.year~disp, finalDispersionPressWithControlNegative)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.005, p = 0.276", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))
divergeFit <- lm(trt.year~disp, finalDispersionPressWithControlPositive)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.075, p < 0.001", sep="")), x=0.97, y=0.97, hjust=1, gp=gpar(col="black", fontsize=15)))

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
finalMeanPressInfo <- merge(finalMean, expInfo, by="label")
finalDivergePressInfo <- merge(finalDispersionPressWithControlPositive, expInfo, by="label")
finalConvergePressInfo <- merge(finalDispersionPressWithControlNegative, expInfo, by="label")
#################################################################################################

###community richness

#regessions with just number of manipulations
meansFit <- lm(species_num~dist, finalMeanPressInfo)
summary(meansFit)
AIC(meansFit)
# meansFitExp <- lm(species_num~(-dist^2), finalMeanPressInfo)
# summary(meansFitExp)
# AIC(meansFitExp)
meansRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.070, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
divergeFit <- lm(species_num~disp, finalDivergePressInfo)
summary(divergeFit)
divergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.177, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))
convergeFit <- lm(species_num~disp, finalConvergePressInfo)
summary(convergeFit)
convergeRegText <- grobTree(textGrob(expression(paste(R^2, " = 0.097, p < 0.001", sep="")), x=0.03, y=0.97, hjust=0, gp=gpar(col="black", fontsize=15)))

#plots with number of manipulations
finalMeanPressLength <- ggplot(finalMeanPressInfo, aes(x=species_num, y=dist)) +
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
convergeManiType <- rbind(waterCon, nutsCon, carbonCon)

#plotting by manipulation type
meanManiPlot <- ggplot(barGraphStats(data=meanManiType, variable="dist", byFactorNames=c("mani_type")), aes(x=mani_type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.4, 0.1)), name="Change in Mean") +
  coord_cartesian(ylim=c(0,0.4)) +
  xlab("Manipulation Type") +
  ggtitle("Change in Mean")
divergeManiPlot <- ggplot(barGraphStats(data=divergeManiType, variable="disp", byFactorNames=c("mani_type")), aes(x=mani_type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.5, 0.1)), name="Dispersion") +
  xlab("Manipulation Type") +
  coord_cartesian(ylim=c(0,0.4)) +
  ggtitle("Divergence")
convergeManiPlot <- ggplot(barGraphStats(data=convergeManiType, variable="disp", byFactorNames=c("mani_type")), aes(x=mani_type, y=mean)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se), width=0.2)) +
  scale_y_continuous(breaks=(seq(0, 0.5, 0.1)), name="Dispersion") +
  coord_cartesian(ylim=c(0,0.4)) +
  xlab("Manipulation Type") +
  ggtitle("Convergence")
pushViewport(viewport(layout=grid.layout(1,3)))
print(meanManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(divergeManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=2))
print(convergeManiPlot, vp=viewport(layout.pos.row=1, layout.pos.col=3))













