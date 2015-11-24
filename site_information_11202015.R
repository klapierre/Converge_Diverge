setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)
library(vegan)

#import the list of all experiments site information
ExpInfo <- read.csv("SpeciesRelativeAbundance_Nov2015.csv")%>%
  select(-X)

ExpList<-ExpInfo%>%
  select(site_code, project_name)%>%
  unique()
#write.csv(ExpList, "Experiment_List.csv")

#Getting control ANPP
# ANPP<-read.csv("ANPP_Nov2015.csv")
# 
# Experiment_Info<-read.csv("ExperimentInformation_Nov2015.csv")%>%
#   select(site_code, project_name, community_type, treatment, plot_mani)%>%
#   unique()
# 
# ExpList<-Experiment_Info%>%
#   select(site_code, project_name, community_type)%>%
#   unique()
# 
# controlANPP<-merge(ANPP, Experiment_Info, by=c("site_code","project_name","community_type","treatment"))%>%
#   filter(plot_mani==0)%>%
#   na.omit%>%
#   tbl_df()%>%
#   group_by(site_code, project_name, community_type, treatment_year)%>%
#   summarize(anpp=mean(anpp))%>%
#   tbl_df()%>%
#   group_by(site_code, project_name, community_type)%>%
#   summarize(anpp=mean(anpp))
# 
# ExpANPP<-merge(controlANPP, ExpList, by=c("site_code","project_name","community_type"), all=T)
# 
# write.csv(ExpANPP, "ExperimentANPP_toadd.csv")

# siteList<-ExpInfo%>%
#   select(site_code)%>%
#   unique()
#write.csv(siteList, "SiteList_LatLong.csv")

SiteClimate<-read.csv("siteList_climate.csv")

ExpANPP<-read.csv("ExperimentANPP.csv")

ExpLength<-ExpInfo%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(experiment_length=max(treatment_year))

##need to rarefy species. not sure how to do this with cover data.


##calculate chao richness and rarefied richness for each site
#import species data
species <- read.csv("SpeciesRawAbundance_Nov2015.csv")%>%
  select(site_code, project_name, community_type, plot_id, calendar_year, genus_species, abundance)%>%
  mutate(exp=paste(site_code, project_name, community_type, sep='::'))%>%
  #get rid of duplicate species within a plot and year in the dataset; once we contact the dataowners, this step will no longer be needed
  tbl_df()%>%
  group_by(exp,site_code, project_name, community_type, calendar_year, plot_id, genus_species)%>%
  summarise(abundance=mean(abundance))

SampleIntensity<-species%>%
  tbl_df()%>%
  group_by(exp, plot_id, calendar_year)%>%
  summarize(SampleIntensity=length(abundance))%>%
  tbl_df()%>%
  group_by(exp)%>%
  summarize(SampleIntensity=length(SampleIntensity))

exp=species%>%
  select(exp)%>%
  unique()

#create empty dataframe for loop
estimatedRichness=data.frame(row.names=1) 

for(i in 1:length(exp$exp)) {

  #creates a dataset for each unique experiment
  subset <- species%>%
    filter(exp==exp$exp[i])%>%
    select(exp, plot_id, calendar_year, genus_species, abundance)
  
  #transpose data into wide form
  speciesData <- subset%>%
    spread(genus_species, abundance, fill=0)
  
  #calculate species accumulation curves
  pool <- poolaccum(speciesData[,4:ncol(speciesData)], permutations=100)
  chao <- as.data.frame(as.matrix(pool$chao))#this gives us estimated richness from 1-X samples
  chao$aveChao<-rowMeans(chao)
  chao$n<-row.names(chao)
  chao$exp<-exp$exp[i]
  chao2<-chao%>%
    select(exp,n, aveChao)
  
  #rbind back
  estimatedRichness<-rbind(chao2, estimatedRichness)
}

  #calculate rarefied value
  #we should get all of the chao richness estimates for samples 1-X for each experiment, cbind them all together, find the min sample number across all experiments, and then use that to get rarefied richness for each of the experiments
  
  #in theory we could also get estimated richness from our pool/species accumulation curves by extrapolating out beyond the number of samples collected

pplots<-species%>%
  filter(project_name=="pplots")%>%
  spread(genus_species, abundance, fill=0)

specpool(pplots[,7:ncol(pplots)])
pool<-poolaccum(pplots[,7:ncol(pplots)])
chao <- as.data.frame(as.matrix(pool$chao))
chao$average<-rowMeans(chao)
chao$n<-row.names(chao)
chao2<-chao%>%
  select(n, average)
specpool<-as.data.frame(as.matrix(specpool$chao))




















