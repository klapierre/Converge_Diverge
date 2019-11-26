library(tidyverse)
library(data.table)

#meghan's
setwd("C://Users/mavolio2/Dropbox/converge_diverge/datasets/Traits/Try Data Nov 2019")

#kim's desktop
setwd('C:\\Users\\komatsuk\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\Traits\\Try Data Nov 2019')

dat<-fread("7764.txt",sep = "\t",data.table = FALSE,stringsAsFactors = FALSE,strip.white = TRUE)

#removing trait outliers
dat2<-dat%>%
  select(DatasetID,ObsDataID, AccSpeciesID, AccSpeciesName, TraitID, OriglName, TraitName, OrigValueStr, OrigUnitStr, StdValue, UnitName, ErrorRisk)%>%
  mutate(ErrorRisk2=ifelse(is.na(ErrorRisk), 0, ErrorRisk))%>%
  filter(ErrorRisk2<8)%>%
  filter(!is.na(TraitID))

#mering corre with try
key<-read.csv("corre2trykey.csv")%>%
  select(species_matched, AccSpeciesID, AccSpeciesName)%>%
  unique()

length(unique(key$species_matched))

dat3<-dat2%>%
  right_join(key)%>%
  select(-ErrorRisk, -ErrorRisk2)


##how many traits for sp?
sdivtrt<-read.csv("TRY_traits_type_11252019.csv")%>%
  filter(sdiv_trait==1)

traitnum<-dat3%>%
  select(species_matched, TraitID)%>%
  unique()%>%
  group_by(TraitID)%>%
  summarise(nusp=length(species_matched))%>%
  right_join(sdivtrt)

# write.csv(traitnum, "try_traits_export_nov2019.csv", row.names=F)



###cleaning life history traits
trait59<-dat3%>%
  filter(TraitID==59&OrigValueStr!="")%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=="1"|OrigValueStr=="always annual"|OrigValueStr=="ann"|OrigValueStr=="annual"|OrigValueStr=="Annual"|OrigValueStr=="annual-winter annual"|OrigValueStr=="annuals"|OrigValueStr=="winter annual"|OrigValueStr=="winter annuals"|OriglName=="Plant phenology: Annual"&OrigValueStr=="yes"|OriglName=="Plant phenology: Perennial"&OrigValueStr=="no"|OrigValueStr=="summer annuals", "Annual", 
        ifelse(OrigValueStr=="2"|OrigValueStr=="1, 2"|OrigValueStr=="1,2"|OrigValueStr=="1-2"|OriglName=="Plant phenology: Biennial"&OrigValueStr=="yes"|OrigValueStr=="always annual, always biennial"|OrigValueStr=="always biennial"|OrigValueStr=="annual-winter annual, biennial"|OrigValueStr=="annual, Biennial"|OrigValueStr=="Annual, Biennial"|OrigValueStr=="annual/bieenial"|OrigValueStr=="annual/biennial"|OrigValueStr=='annual/bisannual'|OrigValueStr=="biannual"|OrigValueStr=="biasannual"|OrigValueStr=="biennial"|OrigValueStr=="sometimes annual, always biennial"|OrigValueStr=="winter annual-biennial"|OrigValueStr=="always annual, always biennial, always pluriennial-hapaxanthic", "Biennial", 
        ifelse(OriglName=="Plant phenology: Biennial"&OrigValueStr=="no", NA, OrigValueStr))))

table(trait59$CleanTraitValue)

trait59_test<-dat3%>%
  filter(TraitID==59)%>%
  filter(OrigValueStr=="always biennial, always pluriennial-hapaxanthic")%>%
  select(OriglName, OrigValueStr, species_matched)%>%
  unique()


#leaf area - merging different traits that all correspond to leaf area
#do everything with the StdValue, which is converted to mm2
trait3108_3109_3110_3111_3112_3113_3114<-dat3%>%
  filter(TraitID %in% c(3114, 3108, 3110, 3112, 3109, 3111, 3113)) #all data related to leaf areas

#getting averages within each species for each trait type, are they comparable?
traitLeafAreaTest <- trait3108_3109_3110_3111_3112_3113_3114%>%
  group_by(DatasetID, species_matched, TraitID)%>%
  summarise(DatasetValue=mean(StdValue))%>% #averaging by trait and species within each dataset
  ungroup()%>%
  group_by(species_matched, TraitID)%>%
  summarise(SppValue=mean(DatasetValue))%>% #averaging by trait and species across datasets
  ungroup()





#removing outliers by species and genus
dat4<-dat3%>%
  group_by(species_matched, TraitID)%>%
  summarize(outlier_sp=(StdValue))
