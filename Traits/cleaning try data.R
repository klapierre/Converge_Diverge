library(tidyverse)
library(data.table)

#meghan's
setwd("C://Users/mavolio2/Dropbox/converge_diverge/datasets/Traits/Try Data Nov 2019")

#kim's desktop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\Traits\\Try Data Nov 2019')

dat<-fread("7764.txt",sep = "\t",data.table = FALSE,stringsAsFactors = FALSE,strip.white = TRUE)

#removing trait outliers
dat2<-dat%>%
  select(DatasetID, SpeciesName, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, ValueKindName, OrigUncertaintyStr, UncertaintyName, Replicates, StdValue, UnitName, OrigObsDataID, ErrorRisk)%>%
  mutate(ErrorRisk2=ifelse(is.na(ErrorRisk), 0, ErrorRisk))%>%
  filter(ErrorRisk2<8)%>%
  filter(!is.na(TraitID))

#mering corre with try
key<-read.csv("corre2trykey.csv")%>%
  select(species_matched, AccSpeciesID, AccSpeciesName)%>%
  unique()

length(unique(key$species_matched))

dat3<-dat2%>%
  left_join(key)

test<-dat3%>%
  filter(TraitID==597)

table(test$OrigValueStr)

length(unique(test$DatasetID))

##how many traits for sp?
sdivtrt<-read.csv("TRY_traits_type_11252019.csv")%>%
  filter(sdiv_trait==1)

traitnum<-dat3%>%
  select(species_matched, TraitID)%>%
  unique()%>%
  group_by(TraitID)%>%
  summarise(nusp=length(species_matched))%>%
  right_join(sdivtrt)

write.csv(traitnum, "try_traits_touse_nov2019.csv", row.names=F)
##categorical traits dataset
cattraits<-dat3%>%
  filter(t)

#removing outliers by species and genus
dat4<-dat3%>%
  group_by(species_matched, TraitID)%>%
  summarize(outlier_sp=(StdValue))
