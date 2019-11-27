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
  mutate(CleanTraitValue=ifelse(OriglName=="aquatic"|OriglName=="carnivory"|OriglName=="Crop"|OriglName=="Ecological type"|OriglName=="final growth form 4 categories"&OrigValueStr=="herb"|OrigValueStr=="amphibiousubmerged"|OrigValueStr=="angiosperm"|OrigValueStr=="annual"|OrigValueStr=="Annual"|OrigValueStr=="aquatic"|OrigValueStr=="Aquatic"|OrigValueStr=="aquatic fresh water"|OrigValueStr=="aquatic, fresh water, floating"|OrigValueStr=="Aquatic_Epiphyte_Shrub_Tree_Vine_Herb"|OrigValueStr=="Aquatic_Herb"|OrigValueStr=="Aquatic_Shrub_Herb"|OrigValueStr=="Aquatic_Tree_Herb"|OriglName=="GGF"&OrigValueStr=="H"|OriglName=="SGF"&OrigValueStr=="H"|OrigValueStr=="carnivore"|OrigValueStr=="CARNIVORE"|OrigValueStr=="carnivorous"|OrigValueStr=="chasmophyte"|OrigValueStr=="crop"|OrigValueStr=="crops"|OrigValueStr=="Colonizing"|OrigValueStr=="Columnar"|OrigValueStr=="Conical"|OrigValueStr=="Decumbent"|OrigValueStr=="Emergent attached to the substrate"|OrigValueStr=="Epiphyte"|OrigValueStr=="Erect"|OriglName=="shoot growth form"|OriglName=="Shape and Orientation"|OrigValueStr=="free"|OrigValueStr=="geophyte"|OrigValueStr=="gymnosperm"|OriglName=="parasite"|OriglName=="Parasitic"|OrigValueStr=="hydrohalophyte"|OrigValueStr=="Hydrophytes"|OrigValueStr=="Liana_Herb"|OrigValueStr=="HEMI-PARASITE"|OriglName=="Growth form"&OrigValueStr=="herb"|OriglName=="Growth Form (herb,shrub,tree,herbaceous vine,liana/woody vine)"&OrigValueStr=="herb"|OriglName=="Plant growth form"&OrigValueStr=="herb"|OriglName=="Plant growth form"&OrigValueStr=="herbaceous"|OriglName=="Plant growth form"&OrigValueStr=="herb/sub-shrub"|OriglName=="Plant growth form"&OrigValueStr=="herb/shrub"|OriglName=="GrowthForm"&OrigValueStr=="herb"|OriglName=="GROWTHFORM_STD"&OrigValueStr=="Herb"|OriglName=="GROWTHFORM_DIV"&OrigValueStr=="Herb"|OriglName=="GROWTHFORM_DIV"&OrigValueStr=="Shrub_Herb"|OriglName=="Life form"&OrigValueStr=="herbaceous"|OriglName=="Life form"&OrigValueStr=="Herbaceous perennial"|OrigValueStr=="herbaceous monocotyl"|OriglName=="plant growth form"&OrigValueStr=="herb"|OriglName=="GROWTHFORM_ORG"&OrigValueStr=="Herb"|OrigValueStr=="no"|OrigValueStr=="No"|OrigValueStr=="Nano-chamaephyte"|OrigValueStr=="non-woody"|OrigValueStr=="Nontree"|OrigValueStr=="Nonvascular"|OrigValueStr=="Parasite"|OrigValueStr=="Parasite_Herb"|OrigValueStr=="parasitic"|OrigValueStr=="perennial"|OrigValueStr=="Perennial"|OriglName=="Life form: geophyte"|OriglName=="Life form: epiphyte/parasite"|OriglName=="Plant form: Non-distinctive"|OriglName=="Plant form: prostrate"|OriglName=="Plant form: cushion"|OriglName=="Plant form: open"|OrigValueStr=="xerohalophyte"|OrigValueStr=="xerophyte"|OrigValueStr=="weed"|OrigValueStr=="weedy"|OrigValueStr=="Shrub_Vine_Herb"|OrigValueStr=="Tree_Vine_Shrub_Herb"|OrigValueStr=="Vine_Herb"|OrigValueStr=="Vine_Shrub_Tree_Herb"|OrigValueStr=="Vine_Shrub_Tree"|OrigValueStr=="Vine_Shrub_Herb"|OrigValueStr=="terrestrial"|OrigValueStr=="Terrestrial Herb"|OrigValueStr=="therophyte"|OrigValueStr=="Thicket Forming"|OrigValueStr=="moss"|OrigValueStr=="Multiple Stem"|OrigValueStr=="psammophile"|OrigValueStr=="pteridophyte"|OrigValueStr=="Single Crown"|OrigValueStr=="Single Stem"|OrigValueStr=="PS"|OrigValueStr=="rhiz"|OrigValueStr=="rhizomatous"|OrigValueStr=="Rhizomatous"|OrigValueStr=="SC"|OrigValueStr=="Scap"|OrigValueStr=="scrub"|OrigValueStr=="Rhiz"|OrigValueStr=="Herb/Aquatic"|OrigValueStr=="herb/non-woody"|OrigValueStr=="herb/palmoid/non-woody"|OrigValueStr=="Herb/Shrub"|OrigValueStr=="herb/shrub/climber/non-woody/woody"|OrigValueStr=="herb/shrub/non-woody/woody"|OrigValueStr=="herb/shrub/palmoid/non-woody/woody"|OrigValueStr=="Herb/Shrub/Subshrub"|OrigValueStr=="Herb/Shrub/Tree"|OrigValueStr=="Herb/Shrub/Vine"|OrigValueStr=="non-succulent", NA,
    ifelse(OrigValueStr=="b H"|OrigValueStr=="a T"|OrigValueStr=="forb"|OrigValueStr=="herbaceous legume"|OrigValueStr=="annual forb"|OriglName=="GF"&OrigValueStr=="H"|OrigValueStr=="Forb"|OrigValueStr=="forb (herbaceous, with or without woody base)"|OrigValueStr=="h"|OrigValueStr=="D"|OrigValueStr=="HS"|OrigValueStr=="HSA"|OrigValueStr=="HSL"|OrigValueStr=="HSLT"|OrigValueStr=="HST"|OrigValueStr=="club moss"|OrigValueStr=="CLUBMOSS"|OrigValueStr=="Club moss"|OrigValueStr=="Epiphyte_Herb"|OrigValueStr=="Epiphyte_Liana_Tree_Shrub_Herb"|OrigValueStr=="Epiphyte_Tree_Vine_Shrub_Herb"|OrigValueStr=="Forb/herb"|OrigValueStr=="Forbs"|OrigValueStr=='Forb/herb, Subshrub'|OrigValueStr=="Forb/herb, Vine"|OrigValueStr=="Forb/herb, Shrub, Subshrub"|OrigValueStr=="Forb/herb, Shrub, Subshrub, Vine"|OrigValueStr=="H"|OrigValueStr=="herbs"|OrigValueStr=="perennial leguminous herb"|OriglName=="PlantGrowthFormConsolidated"&OrigValueStr=="herb"|OriglName=="PlantGrowthFormConsolidated"&OrigValueStr=="herb/shrub"|OriglName=="GrowthFormCleartext"&OrigValueStr=="herb"|OriglName=="GrowthformCleartext"&OrigValueStr=="herb"|OriglName=="Life form"&OrigValueStr=="herb"|OriglName=='Life Form'&OrigValueStr=="herb"|OrigValueStr=="herbaceous dicotyl"|OrigValueStr=="Herbaceous Dicot"|OriglName=="type"&OrigValueStr=="herb"|OriglName=="growth_form   TRY"&OrigValueStr=="herb"|OriglName=="Plant Growth Form"&OrigValueStr=="Herb"|OriglName=="Plant Growth Form"&OrigValueStr=="herbaceous"|OriglName=="Plant growth forb"&OrigValueStr=="Herbaceous"|OrigValueStr=="perennial forb"|OriglName=="Life form: forb"|OrigValueStr=="Tree_Vine_Aquatic_Shrub_Herb"|OrigValueStr=="Tree_Vine_Herb"|OrigValueStr=="variable forb"|OrigValueStr=="n hyd"|OrigValueStr=="n Hyd"|OrigValueStr=="m Hel"|OrigValueStr=="rept"|OrigValueStr=="Rept"|OrigValueStr=="herb/aquatic/non-woody"|OrigValueStr=="herb/hemiparasitic/non-woody", "Forb", 
     ifelse(OrigValueStr=="Tree"|OrigValueStr=="tree"|OrigValueStr=="Absence"|OrigValueStr=="shrub"|OrigValueStr=="woody plant"|OrigValueStr=="f P"|OrigValueStr=="S"|OrigValueStr=="SH"|OrigValueStr=="ST"|OrigValueStr=="T"|OrigValueStr=="t"|OrigValueStr=="Woody"|OrigValueStr=="W"|OrigValueStr=="Shrub"|OrigValueStr=="subshrub (woody <1m)"|OrigValueStr=="sh"|OrigValueStr=="t"|OrigValueStr=="P"|OrigValueStr=="c C"|OrigValueStr=="Chaemaephyte"|OrigValueStr=="conifer"|OrigValueStr=="Conifers"|OrigValueStr=="d z"|OrigValueStr=="d Z"|OrigValueStr=="Deciduous shrub or tree"|OrigValueStr=="Dwarf shrub"|OrigValueStr=="e N"|OrigValueStr=="Epiphyte_Shrub_Herb"|OrigValueStr=="Epiphyte_Vine_Tree_Shrub_Herb"|OrigValueStr=="erect dwarf shrub"|OrigValueStr=="evergreen shrub or tree"|OrigValueStr=="large shrub"|OrigValueStr=="low to high shrub"|OrigValueStr=="palmoid"|OrigValueStr=="prostrate dwarf shrub"|OriglName=="Life form: erect dwarf shrub"|OriglName=="Life form: prostrate dwarf shrub"|OriglName=="Life form: shrub"|OriglName=="Life form: tree"|OrigValueStr=="Woody Liana"|OrigValueStr=="Woody"|OrigValueStr=="Woody evergreen"|OrigValueStr=="Woody deciduous"|OrigValueStr=="woody at base"|OrigValueStr=="woody"|OrigValueStr=="TREE"|OrigValueStr=="Tree (deciduous)"|OrigValueStr=="Tree (evergreen)"|OrigValueStr=="tree / shrub"|OrigValueStr=="Tree shrub intermediate"|OrigValueStr=="Tree, Shrub"|OrigValueStr=="Tree, Subshrub, Shrub"|OrigValueStr=="tree/palmoid/woody"|OrigValueStr=="Tree/Treelet"|OrigValueStr=="tree/woody"|OrigValueStr=="Tree_Shrub"|OrigValueStr=="trees"|OrigValueStr=="trees/tree"|OrigValueStr=="trees/tree"|OrigValueStr=="Shrub, Subshrub"|OrigValueStr=="Shrub, Subshrub, Tree"|OrigValueStr=="Shrub, Tree"|OrigValueStr=="Shrub,Subshrub"|OrigValueStr=="shrub/palmoid/woody"|OrigValueStr=="Shrub/Subshrub"|OrigValueStr=="shrub/tree"|OrigValueStr=="Shrub/Tree"|OrigValueStr=="Shrub/Tree intermediate"|OrigValueStr=="shrub/tree/palmoid/woody"|OrigValueStr=="Shrub/Tree/Subshrub"|OrigValueStr=="	Shrub/Tree/Treelet"|OrigValueStr=="shrub/tree/woody"|OrigValueStr=="shrub/woody"|OrigValueStr=="Shrub_Tree"|OrigValueStr=="shrub|tree"|OrigValueStr=="shrubs"|OrigValueStr=="small tree"|OrigValueStr=="Small_Tree"|OrigValueStr=="sub-shrub"|OrigValueStr=="Sub-Shrub (Chamaephyte)"|OrigValueStr=="subshrub"|OrigValueStr=="Subshrub"|OrigValueStr=="Subshrub, Shrub"|OrigValueStr=="Subshrub, Shrub, Tree", "Woody", 
     ifelse(OriglName=="Succulence"|OriglName=="Succulence index"|OriglName=="Succulent"|OrigValueStr=="C"|OrigValueStr=="succulent"|OrigValueStr=="Succulent"|OrigValueStr=="succulent/non-woody"|OrigValueStr=="C"|OrigValueStr=="cactus"|OriglName=="Stem succulent"|OriglName=="Leaf succulence"|OrigValueStr=="stem-succulent", "Succulent", 
     ifelse(OrigValueStr=="grass"|OrigValueStr=="sedge"|OrigValueStr=="G"|OrigValueStr=="annual grass"|OrigValueStr=="Bunch"|OrigValueStr=="grass (Poaceae only)"|OrigValueStr=="g"|OrigValueStr=="se"|OrigValueStr=="rus"|OrigValueStr=="C3 grass"|OrigValueStr=="C4 grass"|OrigValueStr=="Caesp"|OrigValueStr=="caesp"|OrigValueStr=="cereal"|OrigValueStr=="forage grass"|OrigValueStr=="graminoid"|OrigValueStr=="GRAMINOID"|OrigValueStr=="graminoid/aquatic"|OrigValueStr=="graminoid/aquatic/non-woody"|OrigValueStr=="graminoid/non-woody"|OrigValueStr=="Graminoids"|OrigValueStr=="Graminoids Tussock"|OrigValueStr=="Graminoid"|OrigValueStr=="Grass"|OrigValueStr=="grass (clonal)"|OrigValueStr=="grasslike"|OrigValueStr=="Herbaceous Monocot"|OrigValueStr=="pasture grass"|OrigValueStr=="perennial graminoid"|OrigValueStr=="perennial grass"|OrigValueStr=="Perennial grass"|OrigValueStr=="prairie grass"|OriglName=="Life form: graminoid"|OriglName=="Low Growing Grass"|OrigValueStr=="Sedge"|OrigValueStr=="SEDGE","Graminoid",
     ifelse(OriglName=="climber"|OriglName=="ClimbingMode"|OrigValueStr=="climber"|OrigValueStr=="Vine"|OrigValueStr=="V"|OrigValueStr=="twiner/climber."|OrigValueStr=="L"|OrigValueStr=="C+Sc"|OrigValueStr=="climber"|OrigValueStr=="climber or creeper"|OrigValueStr=="climber/non-woody"|OrigValueStr=="climber/parasitic"|OrigValueStr=="climber/vine"|OrigValueStr=="climber/woody"|OrigValueStr=="Climbing"|OrigValueStr=="Climber"|OrigValueStr=="Lianas and climbers"|OrigValueStr=="g L"|OrigValueStr=="liana"|OrigValueStr=="Liana"|OrigValueStr=="lianas"|OrigValueStr=="Lianas (wody climbers)"|OrigValueStr=="lianas/Woody Liana"|OrigValueStr=="Lianna"|OrigValueStr=="Liana_Shrub"|OrigValueStr=="Liana_Vine"|OrigValueStr=="Liana_Vine_Herb"|OrigValueStr=="vine"|OriglName=="Life form: climber"|OriglName=="Life form: liana"|OriglName=="Plant form: climbing"|OrigValueStr=="Vines (non-woody climbers)"|OrigValueStr=="Herb_Liana_Vine"|OrigValueStr=="Shrub_Liana_Vine"|OrigValueStr=="Shrub_Vine"|OrigValueStr=="Tree_Shrub_Liana_Herb_Vine"|OrigValueStr=="Vine"|OrigValueStr=="herb/climber/non-woody"|OrigValueStr=="herb/climber/parasitic/non-woody"|OrigValueStr=="herb/climber/woody"|OrigValueStr=="Herb/Liana", "Vine",
     ifelse(OrigValueStr=="F"|OrigValueStr=="fern"|OrigValueStr=="Fern"|OrigValueStr=="FERN"|OrigValueStr=="FERN ALLY"|OrigValueStr=="Fern or fern ally"|OrigValueStr=="fern/non-woody"|OrigValueStr=="FERNALLY"|OrigValueStr=="Ferns"|OrigValueStr=="Ferns and allies (Lycophytes)"|OrigValueStr=="M"|OriglName=="Life form: fern/fern ally", "Fern",
            NA))))))))%>%
  filter(!is.na(CleanTraitValue))

table(trait42$CleanTraitValue)

trait42_test<-trait42%>%
  select(CleanTraitValue)%>%
  unique()

trait42_clean<-trait42%>%
  select(species_matched, CleanTraitValue)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue=ifelse(Fern=="Fern", "Fern", ifelse(Forb=="Forb"&is.na(Graminoid)&is.na(Succulent)&is.na(Vine)&is.na(Woody), 'Forb', 999)))

ifelse(Graminoid=="Graminoid"&is.na(Forb)&is.na(Succulent)&is.na(Woody)&is.na(Vine),"Graminiod", ifelse(Woody=="Woody"&is.na(Forb)&is.na(Vine)&is.na(Graminoid)&is.na(Succulent),"Woody", ifelse(Vine=="Vine"&is.na(Forb)&is.na(Graminoid)&is.na(Woody)&is.na(Succulent),"Vine", ifelse(Succulent=="Succulent"&is.na(Forb)&is.na(Graminoid)&is.na(Vine)&is.na(Woody),"Succulent", "Help")))))))

trait42_test_messy<-trait42%>%
  filter(TraitID==42&OrigValueStr!=""&OrigValueStr!="?"&OriglName!="CONSENSUS")%>%
  select(OriglName, OrigValueStr, species_matched, CleanTraitValue)%>%
  unique()

##leaf area - merging different traits that all correspond to leaf area
#do everything with the StdValue, which is converted to mm2
#DECISION: combine all leaf area data into one clean variable (see below); some of the regressions are very poor, however most are good (>0.7)

#filter outliers
traitLeafAreaGenus <- dat3%>%
  filter(TraitID %in% c(3114, 3108, 3110, 3112, 3109, 3111, 3113))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
  group_by(genus, TraitID)%>%
  summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
  ungroup()

traitLeafAreaSpp <- dat3%>%
  filter(TraitID %in% c(3114, 3108, 3110, 3112, 3109, 3111, 3113))%>% #all data related to leaf areas
  group_by(species_matched, TraitID)%>%
  summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
  ungroup()

trait3108_3109_3110_3111_3112_3113_3114_clean <- dat3%>%
  filter(TraitID %in% c(3114, 3108, 3110, 3112, 3109, 3111, 3113))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
  left_join(traitLeafAreaGenus)%>%
  left_join(traitLeafAreaSpp)%>%
  mutate(genus_zscore=abs((StdValue-genus_mean)/genus_sd), spp_zscore=abs((StdValue-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
  filter(genus_zscore<4, spp_zscore<2)%>%
  group_by(DatasetID, species_matched, TraitID)%>%
  summarise(DatasetValue=mean(StdValue))%>% #averaging by ways to measure leaf area and species within each dataset
  ungroup()%>%
  group_by(species_matched, TraitID)%>%
  summarise(SppValue=mean(DatasetValue))%>% #averaging by ways to measure leaf area and species across datasets
  ungroup()%>%
  group_by(species_matched)%>%
  summarise(CleanTraitValue=mean(SppValue))%>% #averaging by species across datasets and ways to measure leaf area
  ungroup()%>%
  mutate(CleanTraitName='leaf_area', CleanTraitUnit='mm2')
  
  # #getting averages within each species for each trait type, are they comparable?
  # traitLeafAreaTest <- trait3108_3109_3110_3111_3112_3113_3114%>% #note that there are a small number of datasets where there is the same number for leaf vs "if leaflet" (so obv it was not a species with compound leaves, but now the data is in there twice); this needs fixing
  #   group_by(DatasetID, species_matched, TraitID2)%>%
  #   summarise(StdValue=mean(StdValue))%>%
  #   spread(key=TraitID2, value=StdValue, fill=NA)%>%
  #   group_by(DatasetID, species_matched, TraitID)%>%
  #   summarise(DatasetValue=mean(StdValue))%>% #averaging by trait and species within each dataset
  #   ungroup()%>%
  #   group_by(species_matched, TraitID)%>%
  #   summarise(SppValue=mean(DatasetValue))%>% #averaging by trait and species across datasets
  #   ungroup()%>%
  #   mutate(TraitID2=paste("Trait",TraitID, sep=''))%>%
  #   select(DatasetID, species_matched, TraitID2, StdValue)
  # 
  # #regressions to check for comparability of different leaf area variables
  # #for whole leaves
  # ggplot(data=subset(traitLeafAreaTest, Trait3108!=Trait3110), aes(x=Trait3108, y=Trait3110)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,2500) + ylim(0,2500) #with vs without petiole for whole leaves
  # cor(traitLeafAreaTest$Trait3108,traitLeafAreaTest$Trait3110, use = "complete.obs") #r=0.86
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3108, y=Trait3112)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #without vs undefined petiole for whole leaves
  # cor(traitLeafAreaTest$Trait3108,traitLeafAreaTest$Trait3112, use = "complete.obs") #r=0.64
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3108, y=Trait3114)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #without petiole and whole leaf vs undefined petiole and undefined leaf
  # cor(traitLeafAreaTest$Trait3108,traitLeafAreaTest$Trait3114, use = "complete.obs") #r=0.26
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3110, y=Trait3112)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #with vs undefined petiole for whole leaves
  # cor(traitLeafAreaTest$Trait3110,traitLeafAreaTest$Trait3112, use = "complete.obs") #r=0.78
  # 
  # #for leaflets
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3109, y=Trait3111)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,50000) + ylim(0,50000) #with vs without petiole for leaflets
  # with(subset(traitLeafAreaTest, Trait3108<50000),cor(Trait3108,Trait3111, use = "complete.obs")) #r=0.49
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3109, y=Trait3113)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) #+ ylim(0,2500) #without vs undefined petiole for leaflets
  # cor(traitLeafAreaTest$Trait3109,traitLeafAreaTest$Trait3113, use = "complete.obs") #r=0.11
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3109, y=Trait3114)) + geom_point() + geom_abline(intercept=0, slope=1) #+ xlim(0,100000) #+ ylim(0,2500) #without petiole and leaflet vs undefined petiole and undefined leaf
  # cor(traitLeafAreaTest$Trait3109,traitLeafAreaTest$Trait3114, use = "complete.obs") #r=0.47
  # 
  # #whole leaves vs leaflets
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3108, y=Trait3109)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,100000) + ylim(0,50000) #whole leaves vs leaflets without petiole
  # with(subset(traitLeafAreaTest, Trait3108<50000),cor(Trait3108,Trait3109, use = "complete.obs")) #r=0.29
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3110, y=Trait3111)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,20000) + ylim(0,20000) #whole leaves vs leaflets with petiole
  # cor(traitLeafAreaTest$Trait3110,traitLeafAreaTest$Trait3111, use = "complete.obs") #r=0.71
  # 
  # ggplot(data=traitLeafAreaTest, aes(x=Trait3112, y=Trait3113)) + geom_point() + geom_abline(intercept=0, slope=1) #+ xlim(0,100000) #+ ylim(0,2500) ##whole leaves vs leaflets undefined petiole
  # cor(traitLeafAreaTest$Trait3112,traitLeafAreaTest$Trait3113, use = "complete.obs") #r=0.80

  
  
  
##specific leaf area (SLA) - merging different traits that all correspond to SLA
  #DECISION: combine all SLA data into one clean variable (see below)
  #do everything with the StdValue, which is converted to mm2/mg

#filter outliers
traitSLAGenus <- dat3%>%
  filter(TraitID %in% c(3115, 3116, 3117))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
  group_by(genus, TraitID)%>%
  summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
  ungroup()

traitSLASpp <- dat3%>%
  filter(TraitID %in% c(3115, 3116, 3117))%>% #all data related to leaf areas
  group_by(species_matched, TraitID)%>%
  summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
  ungroup()

  trait3115_3116_3117_clean<-dat3%>%
    filter(TraitID %in% c(3115, 3116, 3117))%>% #all data related to SLA
    separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
    left_join(traitSLAGenus)%>%
    left_join(traitSLASpp)%>%
    mutate(genus_zscore=abs((StdValue-genus_mean)/genus_sd), spp_zscore=abs((StdValue-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
    filter(genus_zscore<4, spp_zscore<2)%>%
    group_by(DatasetID, species_matched, TraitID)%>%
    summarise(DatasetValue=mean(StdValue))%>% #averaging by ways to measure SLA and species within each dataset
    ungroup()%>%
    group_by(species_matched, TraitID)%>%
    summarise(SppValue=mean(DatasetValue))%>% #averaging by ways to measure SLA and species across datasets
    ungroup()%>%
    group_by(species_matched)%>%
    summarise(CleanTraitValue=mean(SppValue))%>% #averaging by species across datasets and ways to measure SLA
    ungroup()%>%
    mutate(CleanTraitName='SLA', CleanTraitUnit='mm2 mg-1')
  #hist(trait3115_3116_3117$CleanTraitValue)
  
  # #getting averages within each species for each trait type, are they comparable?
  # traitSLAtest <- trait3115_3116_3117%>%
  #   filter(TraitID %in% c(3115, 3116, 3117))%>%
  #   group_by(DatasetID, species_matched, TraitID)%>%
  #   summarise(DatasetValue=mean(StdValue))%>% #averaging by trait and species within each dataset
  #   ungroup()%>%
  #   group_by(species_matched, TraitID)%>%
  #   summarise(SppValue=mean(DatasetValue))%>% #averaging by trait and species across datasets
  #   ungroup()%>%
  #   mutate(TraitID2=paste("Trait",TraitID, sep=''))%>%
  #   select(-TraitID)%>%
  #   spread(key=TraitID2, value=SppValue, fill=NA)
  # 
  # #regressions to check for comparability of different leaf area variables
  # ggplot(data=traitSLAtest, aes(x=Trait3115, y=Trait3116)) + geom_point() + geom_abline(intercept=0, slope=1) #with vs without petiole; close enough, can be combined
  # lm(data=traitSLAtest$Trait3115,traitSLAtest$Trait3116, use = "complete.obs") #r=0.46
  # 
  # ggplot(data=traitSLAtest, aes(x=Trait3115, y=Trait3117)) + geom_point() + geom_abline(intercept=0, slope=1) #with vs undefined petiole; close enough, can be combined
  # cor(traitSLAtest$Trait3115,traitSLAtest$Trait3117, use = "complete.obs") #r=0.73
  # 
  # ggplot(data=traitSLAtest, aes(x=Trait3116, y=Trait3117)) + geom_point() + geom_abline(intercept=0, slope=1) #without vs undefined petiole; close enough, can be combined
  # cor(traitSLAtest$Trait3116,traitSLAtest$Trait3117, use = "complete.obs") #r=0.57

  

##leaf water content - merging different traits that all correspond to leaf water content
#DECISION: combine all leaf water content data into one clean variable (see below)
#do everything with the StdValue, which is converted to mm2
  
  #filter outliers
  traitLeafWaterGenus <- dat3%>%
    filter(TraitID %in% c(3120, 3121, 3122))%>% #all data related to leaf areas
    separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
    group_by(genus, TraitID)%>%
    summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
    ungroup()
  
  traitLeafWaterSpp <- dat3%>%
    filter(TraitID %in% c(3120, 3121, 3122))%>% #all data related to leaf areas
    group_by(species_matched, TraitID)%>%
    summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
    ungroup()
  
  trait3120_3121_3122<-dat3%>%
    filter(TraitID %in% c(3120, 3121, 3122))%>% #all data related to SLA
    separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
    left_join(traitLeafWaterGenus)%>%
    left_join(traitLeafWaterSpp)%>%
    mutate(genus_zscore=abs((StdValue-genus_mean)/genus_sd), spp_zscore=abs((StdValue-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
    filter(genus_zscore<4, spp_zscore<2)
  
  trait3120_3121_3122_clean <- trait3120_3121_3122%>%
    group_by(DatasetID, species_matched, TraitID)%>%
    summarise(DatasetValue=mean(StdValue))%>% #averaging by ways to measure SLA and species within each dataset
    ungroup()%>%
    group_by(species_matched, TraitID)%>%
    summarise(SppValue=mean(DatasetValue))%>% #averaging by ways to measure SLA and species across datasets
    ungroup()%>%
    mutate(AltValue=ifelse(TraitID==3122, (-0.3513 + 1.1673*SppValue), SppValue))%>% #convert saturated to unsaturated (highly correlated)  
    filter(TraitID!=3121)%>% #drop undefined saturation state (not correlated with other two values and least abundant measurement)
    group_by(species_matched)%>%
    summarise(CleanTraitValue=mean(SppValue))%>% #averaging by species across datasets and ways to measure SLA
    ungroup()%>%
    mutate(CleanTraitName='leaf_water_content', CleanTraitUnit='g(W)/g(DM)')
  #hist(trait3115_3116_3117$CleanTraitValue)
  
  # #getting averages within each species for each trait type, are they comparable?
  # traitLeafWatertest <- trait3120_3121_3122%>%
  #   group_by(DatasetID, species_matched, TraitID)%>%
  #   summarise(DatasetValue=mean(StdValue))%>% #averaging by trait and species within each dataset
  #   ungroup()%>%
  #   group_by(species_matched, TraitID)%>%
  #   summarise(SppValue=mean(DatasetValue))%>% #averaging by trait and species across datasets
  #   ungroup()%>%
  #   mutate(TraitID2=paste("Trait",TraitID, sep=''))%>%
  #   select(-TraitID)%>%
  #   spread(key=TraitID2, value=SppValue, fill=NA)
  # 
  # #regressions to check for comparability of different leaf area variables
  # ggplot(data=traitLeafWatertest, aes(x=Trait3120, y=Trait3121)) + geom_point() + geom_abline(intercept=0, slope=1) + xlim(0,5) #undefined vs not saturated
  # summary(lm(data=traitLeafWatertest, Trait3120~Trait3121))
  # 
  # ggplot(data=traitLeafWatertest, aes(x=Trait3120, y=Trait3122)) + geom_point() + geom_abline(intercept=0, slope=1) #saturated vs not saturated
  # summary(lm(data=traitLeafWatertest, Trait3120~Trait3122))
  # 
  # ggplot(data=traitLeafWatertest, aes(x=Trait3121, y=Trait3122)) + geom_point() + geom_abline(intercept=0, slope=1) #undefinted vs saturated
  # summary(lm(data=traitLeafWatertest, Trait3121~Trait3122))
  
  
  
##filtering continuous traits that TRY has already standardized
traitStandardGenus <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>% #all data related to leaf areas
  separate(species_matched, into=c('genus', 'species'), sep=' ')%>%
  group_by(genus, TraitID)%>%
  summarise(genus_mean=mean(StdValue), genus_sd=sd(StdValue))%>%
  ungroup()

traitStandardSpp <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>% #all data related to leaf areas
  group_by(species_matched, TraitID)%>%
  summarise(spp_mean=mean(StdValue), spp_sd=sd(StdValue))%>%
  ungroup()

#get list of trait units and types
traitStandardContinuousList <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>%
  select(TraitID, UnitName)%>%
  unique()%>%
  mutate(remove=ifelse(TraitID==48&UnitName=='', 1, ifelse(TraitID==3107&UnitName=='cm', 1, 0)))%>% #remove two trait unit names that have been fixed (below) or are redundant
  filter(remove==0)%>%
  select(-remove)

traitStandardContinuous_clean <- dat3%>%
  filter(TraitID %in% c(6,9,12,26,45,46,47,48,55,56,77,80,82,83,84,95,131,145,146,200,363,403,683,1111,2809,3106,3107))%>%
  mutate(StdValue2=ifelse(TraitID==3107&UnitName=='cm', StdValue/100, StdValue))%>%  #fix 195 cases where height was cm, but should be standardized to m
  separate(species_matched, into=c('genus', 'species'), sep=' ', remove=F)%>%
  left_join(traitStandardGenus)%>%
  left_join(traitStandardSpp)%>%
  mutate(genus_zscore=abs((StdValue2-genus_mean)/genus_sd), spp_zscore=abs((StdValue2-spp_mean)/spp_sd))%>% #calculate z-scores for genera and species
  filter(genus_zscore<4, spp_zscore<2)%>%
  group_by(DatasetID, species_matched, TraitID)%>%
  summarise(DatasetValue=mean(StdValue2))%>% #averaging by species within each dataset
  ungroup()%>%
  group_by(species_matched, TraitID)%>%
  summarise(CleanTraitValue=mean(DatasetValue))%>% #averaging by across datasets
  ungroup()%>%
  left_join(traitStandardContinuousList)%>%
  mutate(CleanTraitName=ifelse(TraitID==6, 'rooting_depth', ifelse(TraitID==9, 'root:shoot', ifelse(TraitID==12, 'leaf_longevity', ifelse(TraitID==26, 'seed_dry_mass', ifelse(TraitID==45, 'stomata_conductance', ifelse(TraitID==46, 'leaf_thickness', ifelse(TraitID==47, 'LDMC', ifelse(TraitID==48, 'leaf_density', ifelse(TraitID==55, 'leaf_dry_mass', ifelse(TraitID==56, 'leaf_N:P', ifelse(TraitID==77, 'RGR', ifelse(TraitID==80, 'root_N', ifelse(TraitID==82, 'root_density', ifelse(TraitID==83, 'root_diameter', ifelse(TraitID==84, 'root_C', ifelse(TraitID==95, 'germination_efficiency', ifelse(TraitID==131, 'seed_number', ifelse(TraitID==145, 'leaf_width', ifelse(TraitID==146, 'leaf_C:N', ifelse(TraitID==200, 'number_floristic_zones', ifelse(TraitID==363, 'root_dry_mass', ifelse(TraitID==403, 'shoot_dry_mass', ifelse(TraitID==683, 'root_P', ifelse(TraitID==1111, 'seedbank_density', ifelse(TraitID==2809, 'seedbank_duration', ifelse(TraitID==3106, 'plant_height_vegetative', ifelse(TraitID==3107, 'plant_height_regenerative', TraitID))))))))))))))))))))))))))))%>%
  rename(CleanTraitUnit=UnitName)%>%
  #dropping some traits
  filter(CleanTraitName!='germination_efficiency')

# ggplot(data=traitStandardContinuous_clean, aes(x=CleanTraitValue)) + geom_histogram() + facet_wrap(~CleanTraitName, scales='free')
# #some traits have very skewed distributions, but looking at the data they seem ok (leaf_C:N, leaf_dry_mass, root_dry_mass, seed_dry_mass, seed_number, seedbank_density, seedbank_duration)


##clean up categorical traits that need little cleaning

#heterotrophy
trait201_clean <- dat3%>%
  filter(TraitID==201)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=='always carnivorous', 'carnivorous', ifelse(OrigValueStr %in% c('always hemiparasitic', 'hemi-parasitic'), 'hemiparasitic', ifelse(OrigValueStr=='mycotrophic', 'mycotrophic', 'autotrophic'))))%>%
  mutate(remove=ifelse(species_matched=='Drosera rotundifolia'&CleanTraitValue=='autotrophic', 1, ifelse(species_matched=='Bartsia alpina'&CleanTraitValue=='autotrophic', 1, ifelse(species_matched=='Botrychium lunaria'&CleanTraitValue=='autotrophic', 1, 0))))%>% #fix three overlapping species (if anything in addition to autotroph, listed as such)
  filter(remove==0)%>%
  mutate(CleanTraitName='heterotrophy', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()

#chemical plant defense
trait346_clean <- dat3%>%
  filter(TraitID==346)%>%
  filter(!is.na(OrigValueStr))%>%
  unique()
  
trait357_clean <- dat3%>%
  filter(TraitID==357)%>%
  filter(OrigValueStr!='')%>%
  rename(CleanTraitValue=OrigValueStr)%>%
  mutate(CleanTraitName='clonal_organ_role', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue, CleanTraitValue)%>%
  mutate(CleanTraitValue2=paste(necessary, additive, regenerative, none, sep=','))%>%
  mutate(CleanTraitValue=ifelse(CleanTraitValue2 %in% c('NA,additive,NA,NA', 'NA,additive,NA,none'), 'additive', ifelse(CleanTraitValue2 %in% c('NA,additive,regenerative,NA', 'NA,additive,regenerative,none'), 'additive,regenerative', ifelse(CleanTraitValue2=='NA,NA,NA,none', 'none', ifelse(CleanTraitValue2 %in% c('NA,NA,regenerative,NA', 'NA,NA,regenerative,none'), 'regenerative', ifelse(CleanTraitValue2 %in% c('necessary,additive,NA,NA', 'necessary,additive,NA,none'), 'necessary,additive', ifelse(CleanTraitValue2 %in% c('necessary,additive,regenerative,NA', 'necessary,additive,regenerative,none'), 'necessary,additive,regenerative', ifelse(CleanTraitValue2 %in% c('necessary,NA,NA,NA', 'necessary,NA,NA,none'), 'necessary', 'necessary,regenerative'))))))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait597_clean <- dat3%>%
  filter(TraitID==597)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue2=OrigValueStr)%>%
  mutate(CleanTraitName='flowering_requirement', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(is.na(High)&is.na(Medium), 'Low', ifelse(is.na(High), 'Medium', 'High')))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait613_clean <- dat3%>%
  filter(TraitID==613)%>%
  filter(!is.na(OrigValueStr))%>%
  rename(CleanTraitValue2=OrigValueStr)%>%
  mutate(CleanTraitName='vegetative_spread_rate', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue2, CleanTraitUnit)%>%
  unique()%>%
  spread(CleanTraitValue2, CleanTraitValue2)%>%
  mutate(CleanTraitValue=ifelse(is.na(Rapid)&is.na(Moderate)&is.na(Slow), 'None', ifelse(is.na(Rapid)&is.na(Moderate), 'Slow', ifelse(is.na(Rapid), 'Moderate', 'Rapid'))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

trait1187_clean <- dat3%>%
  filter(TraitID==1187)%>%
  rename(CleanTraitValue=OrigValueStr)%>%
  mutate(CleanTraitName='stem_longevity', CleanTraitUnit='year')%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()
    
trait1188_clean <- dat3%>%
  filter(TraitID==1188)%>%
  mutate(CleanTraitValue=ifelse(OrigValueStr=='self', 'self-supporting', OrigValueStr))%>%
  mutate(CleanTraitName='stem_support', CleanTraitUnit=NA)%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)%>%
  unique()%>%
  mutate(CleanTraitValue2=ifelse(CleanTraitValue=='self-supporting', 'self_supporting', CleanTraitValue))%>%
  select(-CleanTraitValue)%>%
  spread(CleanTraitValue2, CleanTraitValue2, fill='fix')%>%
  mutate(CleanTraitValue2=paste(creeping, decumbent, procumbent, scrambling, self_supporting, tendrils, twining, sep=','))%>%
  mutate(CleanTraitValue=ifelse(CleanTraitValue2 %in% c('creeping,fix,fix,fix,fix,fix,fix','creeping,decumbent,fix,fix,fix,fix,fix'), 'creeping', ifelse(CleanTraitValue2 %in% c('fix,decumbent,fix,fix,fix,fix,fix', 'fix,decumbent,fix,scrambling,fix,fix,fix'), 'decumbent', ifelse(CleanTraitValue2 %in% c('fix,decumbent,procumbent,fix,fix,fix,fix', 'fix,fix,procumbent,fix,fix,fix,fix'), 'procumbent', ifelse(CleanTraitValue2=='fix,fix,fix,fix,fix,fix,twining', 'twining', ifelse(CleanTraitValue2=='fix,fix,fix,fix,fix,tendrils,fix', 'tendrils', ifelse(CleanTraitValue2=='fix,fix,fix,scrambling,fix,fix,fix', 'scrambling', 'self-supporting')))))))%>%
  select(species_matched, CleanTraitName, CleanTraitValue, CleanTraitUnit)

    
#combining traits
traits <- rbind(trait59_clean, trait3115_3116_3117_clean, trait3108_3109_3110_3111_3112_3113_3114_clean, traitStandardContinuous_clean, trait201_clean, trait346_clean, trait357_clean, trait597_clean, trait613_clean, trait1187_clean, trait1188_clean)
