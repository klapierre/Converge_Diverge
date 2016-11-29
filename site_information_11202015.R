#kim's
setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#meghan's
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)
library(vegan)

#import the list of all experiments site information
ExpInfo <- read.csv("SpeciesRelativeAbundance_March2016.csv")%>%
  select(-X)

ExpList<-ExpInfo%>%
  select(site_code, project_name, community_type)%>%
  unique()

#Getting ANPP
ANPP<-read.csv("ANPP_March2016.csv")

Experiment_Info<-read.csv("ExperimentInformation_Mar2016.csv")%>%
  select(site_code, project_name, community_type, treatment, plot_mani, public)%>%
  unique()

ExpList2<-Experiment_Info%>%
  select(site_code, project_name, community_type, public)%>%
  unique()
# write.csv(ExpList2, "Experiment_List_Feb2016a.csv")

controlANPP<-merge(ANPP, Experiment_Info, by=c("site_code","project_name","community_type","treatment"))%>%
  filter(plot_mani==0)%>%
  na.omit%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type, treatment_year)%>%
  summarize(anpp=mean(anpp))%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(anpp=mean(anpp))

ANPP_nocont<-read.csv("ANPP_noControls.csv")

AllANPP<-rbind(ANPP_nocont, controlANPP)

ExpANPP<-merge(AllANPP, ExpList2, by=c("site_code","project_name","community_type"), all=T)

# siteList<-ExpInfo%>%
#   select(site_code)%>%
#   unique()
#write.csv(siteList, "SiteList_LatLong.csv")

SiteClimate<-read.csv("siteList_climate_Feb2016.csv")%>%
  mutate(MAP=ifelse(site_code=="Finse", 1030, MAP))%>%
  select(site_code, MAP, MAT)
#for Finse_WarmNut there is a big differnce between this and what they published, and their coordinates were VERY vauge. I am replacing with their value: 1030 mm.

ExpLength<-ExpInfo%>%
  tbl_df()%>%
  group_by(site_code, project_name, community_type)%>%
  summarize(experiment_length=max(treatment_year))


##calculate chao richness and rarefied richness for each site

species <- read.csv("SpeciesRawAbundance_March2016.csv")%>%
  select(site_code, project_name, community_type, plot_id, calendar_year, genus_species, abundance)%>%
  mutate(exp=paste(site_code, project_name, community_type, sep='::'))%>%
  #get rid of duplicate species within a plot and year in the dataset; once we contact the dataowners, this step will no longer be needed
  tbl_df()%>%
  group_by(exp,site_code, project_name, community_type, calendar_year, plot_id, genus_species)%>%
  summarise(abundance=mean(abundance))%>%
  filter(genus_species!="")


SampleIntensity<-species%>%
  tbl_df()%>%
  group_by(exp, plot_id, calendar_year)%>%
  summarize(SampleIntensity=length(abundance))%>%
  tbl_df()%>%
  group_by(exp)%>%
  summarize(SampleIntensity=length(SampleIntensity))#how many plots were sampled over the course of the experiment

exp<-SampleIntensity%>%
  select(exp)

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

ExpRichness<-estimatedRichness%>%
  filter(n==34)%>%#the lowest sampling intensity, not including GVN_FACE
  separate(exp, c("site_code", "project_name", "community_type"), sep="::")%>%
  mutate(rrich=aveChao)%>%
  select(-n, -aveChao)

gface<-data.frame(site_code="GVN", project_name="FACE", community_type=0, rrich=30.85)

ExpRichness<-rbind(ExpRichness, gface)

ExpDetails1<-merge(ExpRichness, ExpLength, by=c("site_code","project_name","community_type"))
ExpDetails<-merge(ExpDetails1, ExpANPP, by=c("site_code","project_name","community_type"))

SiteExpDetails<-merge(SiteClimate, ExpDetails, by="site_code")

write.csv(SiteExpDetails, "SiteExperimentDetails_March2016.csv")

pairs(SiteExpDetails[,c(2,3,6:8)])
with(SiteExpDetails,cor.test(MAP, anpp))
with(SiteExpDetails,plot(MAP, anpp))
with(SiteExpDetails,cor.test(rrich, experiment_length))


#####subsetting the data

Experiment_Info<-read.csv("ExperimentInformation_Mar2016.csv")
###Look at experiments that do N only
pplots<-Experiment_Info%>%
  filter(project_name=="pplots")%>%
  filter(p==0)%>%
  select(site_code, project_name, community_type, treatment, n)%>%
  unique()

Nitrogenonly<-Experiment_Info%>%
  filter(nutrients==1)%>%
  filter(p==0&
        k==0&
        CO2==0&
        precip==0&
        temp==0&
        mow_clip==0&
        burn==0&
        herb_removal==0&
        trt_details==0&
        other_trt==0&
        successional==0&
        plant_mani==0&
        pulse==0)%>%
  select(site_code, project_name, community_type, treatment, n)%>%
  unique()

# Nitrogenonly2<-rbind(Nitrogenonly, pplots)%>%
#   select(site_code,project_name, community_type)%>%
#   unique()

Nitrogenonly2<-rbind(Nitrogenonly, pplots)
write.csv(Nitrogenonly2, "nitrogen.csv")



###Look at experiments that do N&P only
pplots_np<-Experiment_Info%>%
  filter(project_name=="pplots")%>%
  select(site_code, project_name, community_type, treatment, n, p)%>%
  unique()

NitPhos<-Experiment_Info%>%
  filter(nutrients==1)%>%
  filter(p>0&
           k==0&
           CO2==0&
           precip==0&
           temp==0&
           mow_clip==0&
           burn==0&
           herb_removal==0&
           trt_details==0&
           other_trt==0&
           successional==0&
           plant_mani==0&
           pulse==0)%>%
  select(site_code, project_name, community_type, treatment, n, p)%>%
  unique()

Nitphos2<-rbind(NitPhos, pplots_np)

write.csv(Nitphos2, "nitrogen_phosphorus.csv")

###Drought

drought<-Experiment_Info%>%
  filter(water==1)%>%
  filter(n==0&
        p==0&
           k==0&
           CO2==0&
           precip<0&
           temp==0&
           mow_clip==0&
           burn==0&
           herb_removal==0&
           trt_details==0&
           other_trt==0&
           successional==0&
           plant_mani==0&
           pulse==0)%>%
  select(site_code, project_name, community_type, treatment, precip)%>%
  unique()

##Irrigation
irrigation<-Experiment_Info%>%
  filter(water==1)%>%
  filter(n==0&
           p==0&
           k==0&
           CO2==0&
           precip>0&
           temp==0&
           mow_clip==0&
           burn==0&
           herb_removal==0&
           trt_details==0&
           other_trt==0&
           successional==0&
           plant_mani==0&
           pulse==0)%>%
  select(site_code, project_name, community_type, treatment, precip)%>%
  unique()

