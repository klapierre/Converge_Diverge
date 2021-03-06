library(rnoaa)
library(ggmap)
library(tidyverse)
library(maps)
library(mapdata)
library(mapproj)
library(geosphere)
library(lubridate)

#get list of global historical climate network daily stations
#stations<-ghcnd_stations()
#write.csv(stations, "~/Dropbox/converge_diverge/datasets/LongForm/climate/ghcn_stations.csv")

#read in data
stations<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/ghcn_stations.csv")
stations<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\ghcn_stations.csv")

stations_subset<-stations%>%
  select(id, latitude, longitude, name, element, first_year, state, last_year)%>%
  filter(element=="PRCP")%>%
  unique()%>%
  filter(first_year<1983&last_year>2015)

station_list<-stations_subset%>%
  select(id, latitude, longitude, name)

correloc<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\siteList_LatLong.csv")

##ploting stations
map <- borders("world", colour="gray80", fill="gray80") # create a layer of borders

ggplot() +
  map+
  geom_point(data=stations_subset, aes(x=longitude, y=latitude))+
  geom_point(data=correloc, aes(x=longitude,y=latitude),col="red")

#pairwise ditances
tocompare<-rbind(station_list, correloc)

long<-tocompare[,3]
lat<-tocompare[,2]
longlat<-cbind(long,lat)
dist1<-as.data.frame(distm(longlat,fun=distCosine))/1000

colnames(dist1)<-tocompare$id
dist2<-cbind(tocompare, dist1)

dist3<-dist2[c(13968:14018),]

dist4<-dist3%>%
  gather(station, distance, 5:14022)

closest_stations<-dist4%>%
  filter(distance<200&distance!=0)%>%
  group_by(name)%>%
  mutate(rank=rank(distance))%>%
  filter(rank<4)%>%
  filter(station!="corre5"&station!='corre6'&station!="corre12")

#write and then check data quality, pick 1 station for each site.
#write.csv(closest_stations, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\close_stations.csv")

touse<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\close_stations.csv")%>%
  filter(use!=0)

touse<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/close_stations.csv")%>%
  filter(use!=0)

token <-'qyXXjlaezYAyeQyiMoKPbIYEQYgCfrch'
options(noaakey=token)

ppt.all<-data.frame()
stations_touse<-touse$station

for (i in 1:length(stations_touse)){
  station<-stations_touse[i]
  ppt.station<-ghcnd_search(station, var = "PRCP")
  ppt.station<-as.data.frame(ppt.station)
  ppt.all<-rbind(ppt.all, ppt.station)
  
}

#create columns for year, month, day and sort
ppt.all$year<-year(ppt.all$prcp.date)
ppt.all$month<-month(ppt.all$prcp.date)
ppt.all$day<-day(ppt.all$prcp.date)

write.csv(ppt.all, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\precip_data.csv")
ppt.all<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\precip_data.csv")
ppt.all<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/precip_data.csv")

##drop failed quality control data
ppt.all2<-ppt.all%>%
  filter(prcp.qflag==" ")%>%
  na.omit()


#drop years with less than 90% of data
check<-ppt.all2%>%
  group_by(prcp.id, year)%>%
  summarize(num=length(prcp.prcp))%>%
  filter(num>330)%>%# each year needs >90% of data
  filter(year>1980)%>%
  summarise(num2=length(num))

check1<-ppt.all2%>%
  group_by(prcp.id, year)%>%
  summarize(num=length(prcp.prcp))%>%
  filter(num>330)

anppsites<-touse%>%
  filter(name=="ANG"|name=="CDR"|name=="KNZ"|name=="maerc"|name=="NWT")%>%
  mutate(site_code=name, prcp.id=station)%>%
  select(site_code, prcp.id)

anpp_sites_data1<-merge(anppsites, ppt.all2, by="prcp.id")#get data for only anpp datasets
anpp_site_data<-merge(anpp_sites_data1, check1, by=c("prcp.id","year"))#filter only years with 90% of data

anpp_yearly<-anpp_site_data%>%
  group_by(prcp.id, year, site_code)%>%
  summarize(precip=sum(prcp.prcp)/10)

write.csv(anpp_yearly, "C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\yearly_precip.csv")

###i went through and flagged sites for which there is okay data (btwn 1-3 years of missing data or probably when less than 90% of the data were not there, and bad data, most of the data is missing and good data. This is called close_stations_R2.

####for the ANPP paper
anpp.precip<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\yearly_precip.csv")%>%
  select(-prcp.id, -X)

sev<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\sev_black grama_00-12.csv")%>%
  group_by(year)%>%
  summarize(precip=sum(precp))%>%
  mutate(site_code='SEV')

dl<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\DL_NSFC.csv")%>%
  mutate(year=ï..year)%>%
  select(year, precip, site_code)

imgers<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\IMGERS.csv")%>%
  mutate(year=ï..year)%>%
  select(year, precip, site_code)

kbs<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\KBS_weather.csv")%>%
  group_by(year)%>%
  summarize(precip=sum(precipitation_mm))%>%
  mutate(site_code='KBS')

klu1<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\KLU-monthly.csv")%>%
  group_by(Year)%>%
  summarize(precip=sum(Total.Precip..mm.))%>%
  mutate(site_code='KLU', year=Year)%>%
  select(-Year)%>%
  filter(year!=2007)

klu2<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\KLU-2007.csv")%>%
  group_by(Year)%>%
  summarize(precip=sum(Total.Precip..mm.))%>%
  mutate(site_code='KLU', year=Year)%>%
  select(-Year)

klu3<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\KLU-2009.csv")%>%
  group_by(Year)%>%
  summarize(precip=sum(Total.Precip..mm.))%>%
  mutate(site_code='KLU', year=Year)%>%
  select(-Year)
klu<-rbind(klu1, klu2, klu3)%>%
  ungroup


serc1<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\SERC\\SERC_Precip_1976-2000.csv")%>%
  mutate(year=ï..year,site_code="SERC")%>%
  group_by(year, site_code)%>%
  summarize(precip=sum(precip_mm, na.rm=T))
serc2<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\SERC\\SERC_Precip_2002_2016.csv")%>%
  mutate(site_code="SERC")%>%
  group_by(year, site_code)%>%
  summarize(precip=sum(precip_mm, na.rm=T))%>%
  filter(year!=2017)
serc3<-read.csv("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\site specific\\SERC\\SERC_PRECIP2000-2009.csv")%>%
  mutate(year=ï..year,site_code="SERC")%>%
  group_by(year, site_code)%>%
  summarize(precip=sum(precip_mm, na.rm=T))%>%
  filter(year!=2000)

serc_test1<-rbind(serc1, serc2)
serc_test2<-rbind(serc1, serc3)
serc_test<-merge(serc_test1, serc3, by="year", all=T)
serc_test3<-merge(serc_test1, serc3, by="year")
##use serc1 and serc3 data as is. 
## for serc2 model the precip data to make it higher

lm(precip.y~precip.x, data=serc_test3)
serc2.2<-serc2%>%
  mutate(precip_corr=precip*1.258+82.33)%>%
  select(-precip)%>%
  mutate(precip=precip_corr)%>%
  select(-precip_corr)%>%
  filter(year>2009)

serc<-rbind(serc_test2, serc2.2)%>%
  ungroup%>%
  mutate(site_code=as.factor(site_code))%>%
  mutate(year=as.integer(year))%>%
  mutate(precip=as.numeric(precip))


precip.all<-rbind(anpp.precip, serc, klu, sev, kbs, imgers, dl)
write.csv(precip.all,"C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\real_precip_anppSites.csv" )
