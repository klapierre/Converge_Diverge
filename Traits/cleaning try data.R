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
        ifelse(OrigValueStr=="2"|OrigValueStr=="1, 2"|OrigValueStr=="1,2"|OrigValueStr=="1-2"|OriglName=="Plant phenology: Biennial"&OrigValueStr=="yes"|OrigValueStr=="always annual, always biennial"|OrigValueStr=="always biennial"|OrigValueStr=="annual-winter annual, biennial"|OrigValueStr=="annual, Biennial"|OrigValueStr=="Annual, Biennial"|OrigValueStr=="annual/bieenial"|OrigValueStr=="annual/biennial"|OrigValueStr=='annual/bisannual'|OrigValueStr=="biannual"|OrigValueStr=="biasannual"|OrigValueStr=="biennial"|OrigValueStr=="sometimes annual, always biennial"|OrigValueStr=="winter annual-biennial"|OrigValueStr=="always annual, always biennial, always pluriennial-hapaxanthic"|OrigValueStr=="strict monocarpic bi-annuals and poly-annuals"|OrigValueStr=="Biennial", "Biennial", 
        ifelse(OriglName=="Plant phenology: Biennial"&OrigValueStr=="no", NA, "Perennial"))))%>%
  filter(!is.na(CleanTraitValue))

table(trait59$CleanTraitValue)

trait59_test<-trait59%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  select(Annual, Biennial, Perennial)%>%
  unique

trait59_clean<-trait59%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(is.na(Biennial)&is.na(Perennial), "Annual", ifelse(is.na(Annual)&is.na(Perennial)|Annual=="Annual"&Biennial=="Biennial"&is.na(Perennial), "Biennial", "Perennial")))%>%
  select(species_matched, CleanTraitValue)%>%
  mutate(CleanTraitName="Lifespan", CleanTraitUnit=NA)


trait42<-dat3%>%
  filter(TraitID==42&OrigValueStr!=""&OrigValueStr!="?"&OriglName!="CONSENSUS")%>%
  mutate(CleanTraitValue=ifelse(OriglName=="aquatic"|OriglName=="carnivory"|OriglName=="Crop"|OriglName=="Ecological type"|OriglName=="final growth form 4 categories"&OrigValueStr=="herb"|OrigValueStr=="amphibiousubmerged"|OrigValueStr=="angiosperm"|OrigValueStr=="annual"|OrigValueStr=="Annual"|OrigValueStr=="aquatic"|OrigValueStr=="Aquatic"|OrigValueStr=="aquatic fresh water"|OrigValueStr=="aquatic, fresh water, floating"|OrigValueStr=="Aquatic_Epiphyte_Shrub_Tree_Vine_Herb"|OrigValueStr=="Aquatic_Herb"|OrigValueStr=="Aquatic_Shrub_Herb"|OrigValueStr=="Aquatic_Tree_Herb"|OriglName=="GGF"&OrigValueStr=="H"|OriglName=="SGF"&OrigValueStr=="H"|OrigValueStr=="carnivore"|OrigValueStr=="CARNIVORE"|OrigValueStr=="carnivorous"|OrigValueStr=="chasmophyte"|OrigValueStr=="crop"|OrigValueStr=="crops"|OrigValueStr=="Colonizing"|OrigValueStr=="Columnar"|OrigValueStr=="Conical"|OrigValueStr=="Decumbent"|OrigValueStr=="Emergent attached to the substrate"|OrigValueStr=="Epiphyte"|OrigValueStr=="Erect"|OriglName=="shoot growth form"|OriglName=="Shape and Orientation"|OrigValueStr=="free"|OrigValueStr=="geophyte"|OrigValueStr=="gymnosperm"|OriglName=="parasite"|OriglName=="Parasitic"|OrigValueStr=="hydrohalophyte"|OrigValueStr=="Hydrophytes"|OrigValueStr=="Liana_Herb"|OrigValueStr=="HEMI-PARASITE"|OriglName=="Growth form"&OrigValueStr=="herb"|OriglName=="Growth Form (herb,shrub,tree,herbaceous vine,liana/woody vine)"&OrigValueStr=="herb", NA,
    ifelse(OrigValueStr=="b H"|OrigValueStr=="a T"|OrigValueStr=="forb"|OrigValueStr=="herbaceous legume"|OrigValueStr=="annual forb"|OriglName=="GF"&OrigValueStr=="H"|OrigValueStr=="Forb"|OrigValueStr=="forb (herbaceous, with or without woody base)"|OrigValueStr=="h"|OrigValueStr=="D"|OrigValueStr=="HS"|OrigValueStr=="HSA"|OrigValueStr=="HSL"|OrigValueStr=="HSLT"|OrigValueStr=="HST"|OrigValueStr=="club moss"|OrigValueStr=="CLUBMOSS"|OrigValueStr=="Club moss"|OrigValueStr=="Epiphyte_Herb"|OrigValueStr=="Epiphyte_Liana_Tree_Shrub_Herb"|OrigValueStr=="Epiphyte_Tree_Vine_Shrub_Herb"|OrigValueStr=="Forb/herb"|OrigValueStr=="Forbs"|OrigValueStr=='Forb/herb, Subshrub'|OrigValueStr=="Forb/herb, Vine"|OrigValueStr=="Forb/herb, Shrub, Subshrub"|OrigValueStr=="Forb/herb, Shrub, Subshrub, Vine"|OrigValueStr=="H"|OrigValueStr=="herbs"|OrigValueStr=="perennial leguminous herb", "Forb", 
     ifelse(OrigValueStr=="Tree"|OrigValueStr=="tree"|OrigValueStr=="Absence"|OrigValueStr=="shrub"|OrigValueStr=="woody plant"|OrigValueStr=="f P"|OrigValueStr=="S"|OrigValueStr=="SH"|OrigValueStr=="ST"|OrigValueStr=="T"|OrigValueStr=="t"|OrigValueStr=="Woody"|OrigValueStr=="W"|OrigValueStr=="Shrub"|OrigValueStr=="subshrub (woody <1m)"|OrigValueStr=="sh"|OrigValueStr=="t"|OrigValueStr=="P"|OrigValueStr=="c C"|OrigValueStr=="Chaemaephyte"|OrigValueStr=="conifer"|OrigValueStr=="Conifers"|OrigValueStr=="d z"|OrigValueStr=="d Z"|OrigValueStr=="Deciduous shrub or tree"|OrigValueStr=="Dwarf shrub"|OrigValueStr=="e N"|OrigValueStr=="Epiphyte_Shrub_Herb"|OrigValueStr=="Epiphyte_Vine_Tree_Shrub_Herb"|OrigValueStr=="erect dwarf shrub"|OrigValueStr=="evergreen shrub or tree"|OrigValueStr=="large shrub"|OrigValueStr=="low to high shrub", "Woody", 
     ifelse(OriglName=="Succulence"|OriglName=="Succulence index"|OriglName=="Succulent"|OrigValueStr=="C"|OrigValueStr=="succulent"|OrigValueStr=="Succulent"|OrigValueStr=="succulent/non-woody"|OrigValueStr=="C"|OrigValueStr=="cactus", "Succulent", 
     ifelse(OrigValueStr=="grass"|OrigValueStr=="sedge"|OrigValueStr=="G"|OrigValueStr=="annual grass"|OrigValueStr=="Bunch"|OrigValueStr=="grass (Poaceae only)"|OrigValueStr=="g"|OrigValueStr=="se"|OrigValueStr=="rus"|OrigValueStr=="C3 grass"|OrigValueStr=="C4 grass"|OrigValueStr=="Caesp"|OrigValueStr=="caesp"|OrigValueStr=="cereal"|OrigValueStr=="forage grass"|OrigValueStr=="graminoid"|OrigValueStr=="GRAMINOID"|OrigValueStr=="graminoid/aquatic"|OrigValueStr=="graminoid/aquatic/non-woody"|OrigValueStr=="graminoid/non-woody"|OrigValueStr=="Graminoids"|OrigValueStr=="Graminoids Tussock"|OrigValueStr=="Graminoid"|OrigValueStr=="Grass"|OrigValueStr=="grass (clonal)"|OrigValueStr=="grasslike", "Graminoid",
     ifelse(OriglName=="climber"|OriglName=="ClimbingMode"|OrigValueStr=="climber"|OrigValueStr=="Vine"|OrigValueStr=="V"|OrigValueStr=="twiner/climber."|OrigValueStr=="L"|OrigValueStr=="C+Sc"|OrigValueStr=="climber"|OrigValueStr=="climber or creeper"|OrigValueStr=="climber/non-woody"|OrigValueStr=="climber/parasitic"|OrigValueStr=="climber/vine"|OrigValueStr=="climber/woody"|OrigValueStr=="Climbing"|OrigValueStr=="Climber"|OrigValueStr=="Lianas and climbers"|OrigValueStr=="g L"|OrigValueStr=="liana"|OrigValueStr=="Liana"|OrigValueStr=="lianas"|OrigValueStr=="Lianas (wody climbers)"|OrigValueStr=="lianas/Woody Liana"|OrigValueStr=="Lianna"|OrigValueStr=="Liana_Shrub"|OrigValueStr=="Liana_Vine"|OrigValueStr=="Liana_Vine_Herb", "Vine",
     ifelse(OrigValueStr=="F"|OrigValueStr=="fern"|OrigValueStr=="Fern"|OrigValueStr=="FERN"|OrigValueStr=="FERN ALLY"|OrigValueStr=="Fern or fern ally"|OrigValueStr=="fern/non-woody"|OrigValueStr=="FERNALLY"|OrigValueStr=="Ferns"|OrigValueStr=="Ferns and allies (Lycophytes)"|OrigValueStr=="M", "Fern",
            OrigValueStr))))))))

table(trait42$CleanTraitValue)

trait42_test<-trait42%>%
  select(OriglName, CleanTraitValue)%>%
  unique()

trait42_test_messy<-trait42%>%
  filter(TraitID==42&OrigValueStr!=""&OrigValueStr!="?"&OriglName!="CONSENSUS")%>%
  select(OriglName, OrigValueStr, species_matched, CleanTraitValue)%>%
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
