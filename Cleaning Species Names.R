#Meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm/")
dir()

library(gtools)
library(reshape2)
library(tidyr)
library(dplyr)
library(taxize)

splist_raw<-read.csv("SpeciesList_Sept2016.csv")

splist<-splist_raw%>%
  select(genus_species)%>%
  unique()

Vsplist<-as.vector(splist$genus_species)
V1<-Vsplist[1:200]
V2<-Vsplist[201:400]
V3<-Vsplist[401:600]
V4<-Vsplist[601:800]
V5<-Vsplist[801:1000]
V6.1<-Vsplist[1001:1050]
V6.2<-Vsplist[1051:1100]
V6.3<-Vsplist[1101:1150]
V6.4<-Vsplist[1151:1200]
V7<-Vsplist[1201:1400]
V8.1<-Vsplist[1401:1500]
V8.2<-Vsplist[1501:1600]
V9<-Vsplist[1601:1800]
V10<-Vsplist[1801:2000]
V11<-Vsplist[2001:2200]
V12<-Vsplist[2201:2400]
V13<-Vsplist[2401:2600]
V14<-Vsplist[2601:2671]

#corrected_names1<-tnrs(V1, source="iPlant_TNRS", sleep = 0)
#corrected_names2<-tnrs(V2, source="iPlant_TNRS", sleep = 0)
#corrected_names3<-tnrs(V3, source="iPlant_TNRS", sleep = 0)
#corrected_names4<-tnrs(V4, source="iPlant_TNRS", sleep = 0)
#corrected_names5<-tnrs(V5, source="iPlant_TNRS", sleep = 0)
#corrected_names6.1<-tnrs(V6.1, source="iPlant_TNRS", sleep = 0)
#corrected_names6.2<-tnrs(V6.2, source="iPlant_TNRS", sleep = 0)
#corrected_names6.3<-tnrs(V6.3, source="iPlant_TNRS", sleep = 0)
#corrected_names6.4<-tnrs(V6.4, source="iPlant_TNRS", sleep = 0)
#corrected_names7<-tnrs(V7, source="iPlant_TNRS", sleep = 0)
#corrected_names8.1<-tnrs(V8.1, source="iPlant_TNRS", sleep = 0)
#corrected_names8.2<-tnrs(V8.2, source="iPlant_TNRS", sleep = 0)
#corrected_names9<-tnrs(V9, source="iPlant_TNRS", sleep = 0)
#corrected_names10<-tnrs(V10, source="iPlant_TNRS", sleep = 0)
#corrected_names11<-tnrs(V11, source="iPlant_TNRS", sleep = 0)
#corrected_names12<-tnrs(V12, source="iPlant_TNRS", sleep = 0)
#corrected_names13<-tnrs(V13, source="iPlant_TNRS", sleep = 0)
#corrected_names14<-tnrs(V14, source="iPlant_TNRS", sleep = 0)

correct_names<-rbind(corrected_names1, corrected_names2, corrected_names3, corrected_names4, corrected_names5, corrected_names6.1,corrected_names6.2,corrected_names6.3,corrected_names6.4, corrected_names7, corrected_names8.1,corrected_names8.2, corrected_names9, corrected_names10, corrected_names11, corrected_names12, corrected_names13, corrected_names14)

correct_names$match<-ifelse(correct_names$acceptedname==correct_names$matchedname, 1, 0)

cn<-merge(correct_names, splist, by.x="submittedname", by.y="genus_species", all=T)
write.csv(cn, "correct_names_toFIX.csv")


clean<-read.csv("correct_names_Corrected_April2016.csv")%>%
  mutate(genus_species=submittedname)%>%
  select(genus_species, acceptedname, type)%>%
  filter(genus_species!="caribou feces")%>%#this drops 1 line of code in MAT2
  filter(genus_species!="standing dead")%>%
  filter(genus_species!="nostoc sp.")%>%
  unique()

write.csv(clean, "cleanspecieslist_April2016_submitted_accepted.csv")
