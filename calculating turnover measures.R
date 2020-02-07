library(vegan)
library(gtools)
library(grid)
library(codyn)
library(tidyverse)

#kim
setwd("C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm")

# #meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=40, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=34),
             axis.title.y=element_text(size=40, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=34),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


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




#read in the merged dataset
alldata<-read.csv("SpeciesRelativeAbundance_March2019.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_code=paste(site_code, project_name, community_type, treatment, sep="::"))

expinfo<-read.csv("ExperimentInformation_March2019.csv")%>%
  mutate(exp_code=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
  select(exp_code, plot_mani, calendar_year, treatment_year)

alldata2<-merge(alldata, expinfo, by=c("exp_code","calendar_year"), all=F)%>%
  filter(treatment_year>0)%>%
  select(-treatment_year)%>%
  filter(exp_code!='GVN::FACE::0::A'&exp_code!='GVN::FACE::0::E')


# codyn function modification ---------------------------------------------
#modifying codyn functions to output integer numbers of species appearing and disappearing, plus total spp number over two year periods, rather than ratios
turnover_allyears <- function(df, 
                              time.var, 
                              species.var, 
                              abundance.var, 
                              metric=c("total", "disappearance","appearance")) {
  
  # allows partial argument matching
  metric = match.arg(metric) 
  
  # sort and remove 0s
  df <- df[order(df[[time.var]]),]
  df <- df[which(df[[abundance.var]]>0),]
  
  ## split data by year
  templist <- split(df, df[[time.var]])
  
  ## create two time points (first year and each other year)
  t1 <- templist[1]
  t2 <- templist[-1]
  
  ## calculate turnover for across all time points
  out <- Map(turnover_twoyears, t1, t2, species.var, metric)
  output <- as.data.frame(unlist(out))
  names(output)[1] = metric
  
  ## add time variable column
  alltemp <- unique(df[[time.var]])
  output[time.var] =  alltemp[2:length(alltemp)]
  
  # results
  return(output)
}

turnover_twoyears <- function(d1, d2, 
                              species.var, 
                              metric=c("total", "disappearance","appearance")){
  
  # allows partial argument matching
  metric = match.arg(metric)
  
  # create character vectors of unique species from each df
  d1spp <- as.character(unique(d1[[species.var]]))
  d2spp <- as.character(unique(d2[[species.var]]))
  
  # ID shared species
  commspp <- intersect(d1spp, d2spp)
  
  # count number not present in d2
  disappear <- length(d1spp)-length(commspp)
  
  # count number that appear in d2
  appear <- length(d2spp)-length(commspp)
  
  # calculate total richness
  totrich <- sum(disappear, appear, length(commspp))
  
  # output based on metric 
  if(metric == "total"){
    output <- totrich
  } else {
    if(metric == "appearance"){
      output <- appear
    } else {
      if(metric == "disappearance"){
        output <- disappear
      }
    }
  }
  
  # results
  return(output)
}







# generating appearances and disappearances for each experiment ---------------------------------------------
#make a new dataframe with just the label;
exp_code=alldata2%>%
  select(exp_code)%>%
  unique()

#makes an empty dataframe
for.analysis=data.frame(row.names=1)

for(i in 1:length(exp_code$exp_code)) {
  
  #creates a dataset for each unique trt, exp combo
  subset=alldata2[alldata2$exp_code==as.character(exp_code$exp_code[i]),]%>%
    select(exp_code, calendar_year, plot_mani, genus_species, relcov, plot_id)%>%
    group_by(exp_code, calendar_year, plot_mani, plot_id, genus_species)%>%
    summarise(relcov=max(relcov))%>%
    ungroup()
  
  #need this to keep track of experiment labels
  labels=subset%>%
    select(calendar_year, exp_code)%>%
    unique()
  
  #calculating appearances and disappearances
  appear<-turnover_allyears(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', metric='appearance')
  disappear<-turnover_allyears(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', metric='disappearance')
  total<-turnover_allyears(df=subset, time.var='calendar_year', species.var='genus_species', abundance.var='relcov', metric='total')

  #merging back with labels to get back experiment labels
  turnover<-merge(appear, disappear, by=c('calendar_year'))
  turnoverAll<-merge(turnover, total, by=c('calendar_year'))
  turnoverLabel<-merge(turnoverAll, labels, by=c('calendar_year'), all=T)

  #pasting into the dataframe made for this analysis
  for.analysis=rbind(turnoverLabel, for.analysis)  
}

# write.csv(for.analysis, 'appear_disappear_Mar2019.csv')

turnoverData <- for.analysis%>%
  mutate(turnover=appearance+disappearance)%>%
  filter(!is.na(turnover))%>%
  #relativize by species richness
  mutate(turnover_relative=turnover/total)%>%
  #filter down to just final year of data for each treatment
  group_by(exp_code)%>%
  filter(calendar_year==max(calendar_year))%>%
  separate(exp_code, c('site_code','project_name','community_type','treatment'), sep='::')%>%
  left_join(read.csv('stdtimebytrt_shape_classification_04072019.csv'))%>%
  mutate(variable=ifelse(is.na(variable), 'control', as.character(variable)))%>%
  filter(variable!='richness')%>%
  mutate(effect=ifelse(variable=='control', 'control', ifelse(N01_shape==0, 'non_sig', 'sig')))
  
  

ggplot(barGraphStats(data=turnoverData, variable="turnover_relative", byFactorNames=c("effect")), aes(x=effect, y=mean)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2)
