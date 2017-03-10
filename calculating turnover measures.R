library(vegan)
library(reshape2)
library(gtools)
library(plyr)
library(grid)
library(tidyr)
library(dplyr)
library(codyn)

#kim
setwd("C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm")

# #meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

#read in the merged dataset
alldata<-read.csv("SpeciesRelativeAbundance_Dec2016.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_code=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
  # #get rid of duplicate species within a plot and year for IMGERS_Yu 2006 plot 107 and 2008 plot 401 duplicate species is salsola collina pall.
  tbl_df()%>%
  group_by(exp_code, site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species)%>%
  summarise(relcov=mean(relcov))%>%
  tbl_df()

expinfo<-read.csv("ExperimentInformation_Dec2016.csv")%>%
  mutate(exp_code=paste(site_code, project_name, community_type, treatment, sep="::"))%>%
  select(exp_code, plot_mani, calendar_year)

alldata2<-merge(alldata, expinfo, by=c("exp_code","calendar_year"), all=F)


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
  
  ## create consecutive pairs of time points
  t1 <- templist[-length(templist)]
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


#write csv

write.csv(for.analysis, 'DiversityMetrics_Dec2016.csv')
