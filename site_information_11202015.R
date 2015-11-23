setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)
library(vegan)

#import the site information
ExpInfo <- read.csv("SpeciesRelativeAbundance_11232015.csv")%>%
  select(-X)

#Getting control ANPP
ANPP<-read.csv("ANPP_11202015.csv")

Experiment_Info<-read.csv

controlANPP<-merge(ANPP, Experiment_Info, by=c("site_code","project_name","community_type","treatment_year","calendar_year","treatment"))%>%
  select(plot_mani==0)%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type, treatment_year)%>%
  summarize(anpp=mean(anpp))%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(anpp=mean(anpp))

ExpANPP<-merge(controlANPP, Experiment_Info, by=c("site_code","project_name","community_type", all=T))

# siteList<-ExpInfo%>%
#   select(site_code)%>%
#   unique()
#write.csv(siteList, "SiteList_LatLong.csv")

ExpList<-ExpInfo%>%
   select(site_code, project_name)%>%
   unique()
write.csv(ExpList, "Experiment_List.csv")

experiment_length<-ExpInfo%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(experiment_length=max(treatment_year))

siteinfo2<-merge(siteinfo, experiment_length, by=c("site_code","project_name","community_type"), all=T)

species_num<-species%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type, genus_species)%>%
  summarize(rich=length(abundance))%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(rich=length(rich))

site_info3<-merge(siteinfo2, species_num, by=c("site_code","project_name","community_type"), all=T)

##need to rarefy species. not sure how to do this with cover data.


##calculate chao richness and rarefied richness for each site
#import species data
species <- read.csv("SpeciesRawAbundance_11202015.csv")%>%
  select(site_code, project_name, community_type, plot_id, calendar_year, genus_species, abundance)%>%
  mutate(exp=paste(site_code, project_name, community_type, sep='::'))

exp=species%>%
  select(exp)%>%
  unique()

#create empty dataframe for loop
estimatedRichness=data.frame(row.names=1) 

for(i in 1:length(exp$exp)) {

  #creates a dataset for each unique experiment
  subset <- species[species$exp==as.character(site$exp[1]),]%>%
    select(exp, plot_id, calendar_year, genus_species, abundance)
  
  #transpose data into wide form
  speciesData <- subset%>%
    spread(genus_species, abundance, fill=0)
  
  #calculate species accumulation curves
  pool <- poolaccum(speciesData[4:ncol(speciesData)], permutations=100)
  chao <- summary(pool, display='chao') #this gives us estimated richness from 1-X samples, but then how do we get that into a dataframe?
  chao <- as.data.frame(chao) #this doesn't remotely work; I also tried making a matrix first
  
  #calculate rarefied value
  #we should get all of the chao richness estimates for samples 1-X for each experiment, cbind them all together, find the min sample number across all experiments, and then use that to get rarefied richness for each of the experiments
  
  #in theory we could also get estimated richness from our pool/species accumulation curves by extrapolating out beyond the number of samples collected

}




#sample from vegan documentation - not remotely helpful
data(dune)
data(dune.env)
attach(dune.env)
pool <- specpool(dune, Management)
pool
op <- par(mfrow=c(1,2))
boxplot(specnumber(dune) ~ Management, col="hotpink", border="cyan3",
        notch=TRUE)
boxplot(specnumber(dune)/specpool2vect(pool) ~ Management, col="hotpink",
        border="cyan3", notch=TRUE)
par(op)
data(BCI)
## Accumulation model
pool <- poolaccum(BCI)
summary(pool, display = "chao")
plot(pool)
## Quantitative model
estimateR(BCI[1:5,])





















