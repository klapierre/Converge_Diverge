library(tidyverse)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\datasets\\LongForm')

expInfo <- read.csv('ExperimentInformation_March2019.csv')%>%
  select(-X)%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-site_code, -project_name, -community_type)

siteData <- read.csv('SiteExperimentDetails_March2019.csv')%>%
  select(-X)%>%
  rename(site_productivity=anpp)%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::'))

coords <- read.csv('siteList_LatLong.csv')%>%
  rename(site_code=name)%>%
  select(site_code, latitude, longitude)

ecosystem <- read.csv('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\CoRRE\\converge_diverge\\Leuzinger working group\\ecosystem_type.csv')

anppPublic <- read.csv('ANPP_Oct2017.csv')%>%
  mutate(site_proj_comm=paste(site_code, project_name, community_type, sep='::'))%>%
  select(-X, -site_code, -project_name, -community_type, -calendar_year, -block)%>%
  group_by(site_proj_comm, treatment_year, treatment)%>%
  summarise(mean=mean(anpp), sd=sd(anpp), len=length(anpp))%>%
  ungroup()

anppPublicCtl <- anppPublic%>%
  left_join(expInfo)%>%
  filter(plot_mani==0)%>%
  select(site_proj_comm, treatment_year, mean, sd, len)%>%
  rename(x_c=mean, sd_c=sd, rep_c=len)

anppPublicFormat <- anppPublic%>%
  rename(x_t=mean, sd_t=sd, rep_t=len)%>%
  left_join(anppPublicCtl)%>%
  left_join(siteData)%>%
  left_join(expInfo)%>%
  left_join(coords)%>%
  left_join(ecosystem)%>%
  filter(public==1)%>%
  rename(trt_corre=treatment, sampling_year=treatment_year, citation=site_proj_comm)%>%
  mutate(response='ANPP', elevation='', other_clim='', soil_texture='', sampling_depth='', response_level='community', plant_organ='', s_t=ifelse(project_name=='BioCON', 16, ''), experiment_type='field', community_type=ifelse(project_name=='BioCON', 'planted', 'natural'), start_year='', dominant_species='', comments='', d_t=ifelse(precip<0, precip, ''), n_t=ifelse(n>0, n, ''), p_t=ifelse(p>0, p, ''), k_t=ifelse(k>0, k, ''), c_added=ifelse(CO2>0, CO2, ''), i_added=ifelse(precip>0, precip, ''), w_t1=ifelse(temp>0, temp, ''), n_c=ifelse(n>0, 0, ''), p_c=ifelse(p>0, 0, ''), k_c=ifelse(k>0, 0, ''), c_c='', c_t='', i_c='', i_t='', r_t1='', r_t2='', w_t2='', w_t3='', x_units='gm-1')%>%
  rename(lat=latitude, lon=longitude, mat=MAT, map=MAP, duration=experiment_length)%>%
  mutate(c=ifelse(CO2>0, 'c', ''), d=ifelse(precip<0, 'd', ''), f=ifelse(n>0|p>0|k>0, 'f', ''), i=ifelse(precip>0, 'i', ''), l= ifelse(trt_type=='light', 'l', ''), r=ifelse(burn>0|mow_clip>0|herb_removal>0, 'r', ''), s=ifelse(plant_mani>0, 's', ''), w=ifelse(temp>0, 'w', ''), treatment=paste(c,d,f,i,l,r,s,w, sep=''))%>%
  filter(treatment!='')%>%
  mutate(n_trt=ifelse(n>0, 1, 0), p_trt=ifelse(p>0, 1, 0), k_trt=ifelse(k>0, 1, 0), npk=paste(n_trt, p_trt, k_trt, sep=''))%>%
  mutate(mow=ifelse(mow_clip>0, 1, 0), till=ifelse(trt_type=='till', 1, 0), defol=0, graz=ifelse(trt_type %in% c('P*burn*graze', 'N*P*burn*graze', 'N*burn*graze', 'burn*graze')|herb_removal>0, 1, 0), disturbance_type=paste(mow, till, defol, graz, burn, sep=''))%>%
  select(citation,response,lat,lon,elevation,mat,map,other_clim,soil_texture,ecosystem,sampling_depth,response_level,plant_organ,treatment,npk,disturbance_type,c_c,c_t,c_added,d_t,n_c,n_t,p_c,p_t,k_c,k_t,i_c,i_t,i_added,r_t1,r_t2,s_t,w_t1,w_t2,w_t3,experiment_type,community_type,sampling_year,start_year,duration,dominant_species,x_c,x_t,x_units,sd_c,sd_t,rep_c,rep_t,comments)
  

# write.csv(anppPublicFormat, 'CoRRE_anpp_public_formatted_02182020.csv', row.names=F)