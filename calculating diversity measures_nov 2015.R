library(vegan)
library(reshape2)
library(gtools)
library(plyr)
library(grid)
library(tidyr)
library(dplyr)

#kim's laptop
setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#kim's desktop
setwd('C:\\Users\\la pierrek\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

# #meghan
setwd("~/Dropbox/converge_diverge/datasets/LongForm")

#read in the merged dataset
alldata<-read.csv("SpeciesRelativeAbundance_March2019.csv")%>%
  select(site_code, project_name, community_type, calendar_year, treatment, block, plot_id, genus_species, relcov)%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))

expinfo<-read.csv("ExperimentInformation_March2019.csv")%>%
  mutate(exp_year=paste(site_code, project_name, community_type, calendar_year, sep="::"))%>%
  select(exp_year, plot_mani, treatment)

alldata2<-merge(alldata, expinfo, by=c("exp_year","treatment"), all=F)

#make a new dataframe with just the label;
exp_year=alldata2%>%
  select(exp_year)%>%
  unique()

#makes an empty dataframe
for.analysis=data.frame(row.names=1) 


###first, get bray curtis dissimilarity values for each year within each experiment between all combinations of plots
###second, get distance of each plot within a trt to the trt centroid 
###third: mean_change is the distance between trt and control centriods
####fourth: dispersion is the average dispersion of plots within a treatment to treatment centriod
for(i in 1:length(exp_year$exp_year)) {
  
  #creates a dataset for each unique year, trt, exp combo
  subset=alldata2[alldata2$exp_year==as.character(exp_year$exp_year[i]),]%>%
    select(exp_year, treatment, plot_mani, genus_species, relcov, plot_id)
  
  #need this to keep track of plot mani
  labels=subset%>%
    select(plot_mani, treatment)%>%
    unique()
  
  #transpose data
  species=subset%>%
   spread(genus_species, relcov, fill=0)

  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,5:ncol(species)], method="bray")
  
  #calculate distances of each plot to treatment centroid (i.e., dispersion)
  disp=betadisper(bc, species$treatment, type="centroid")
  
  #getting distances among treatment centroids; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.data.frame(as.matrix(vegdist(disp$centroids, method="euclidean"))) 
  
  #extracting only the distances we need and adding labels for the comparisons;
  cent_C_T=data.frame(exp_year=exp_year$exp_year[i],
                      treatment=row.names(cent_dist),
                      mean_change=t(cent_dist[names(cent_dist)==labels$treatment[labels$plot_mani==0],]))
  
  #not sure why the name didn't work in the previous line of code, so fixing it here
  names(cent_C_T)[3]="mean_change" 
  
  #merging back with labels to get back plot_mani
  centroid=merge(cent_C_T, labels, by="treatment")
  
  #collecting and labeling distances to centroid from betadisper
  trt_disp=data.frame(data.frame(exp_year=exp_year$exp_year[i], 
                                 plot_id=species$plot_id,
                                 treatment=species$treatment,
                                 dist=disp$distances))%>%
    tbl_df()%>%
    group_by(exp_year, treatment)%>%
    summarize(dispersion=mean(dist))
  
  #merge together change in mean and dispersion data
  distances<-merge(centroid, trt_disp, by=c("exp_year","treatment"))
  
  #getting diversity indixes
    H<-diversity(species[,5:ncol(species)])
    S<-specnumber(species[,5:ncol(species)])
    InvD<-diversity(species[,5:ncol(species)],"inv")
    SimpEven<-InvD/S
    out1<-cbind(H, S)
    output<-cbind(out1, SimpEven)
    divmeasure<-cbind(species, output)%>%
      select(exp_year, treatment, H, S, SimpEven)%>%
      tbl_df()%>%
      group_by(exp_year, treatment)%>%
      summarise(H=mean(H), S=mean(S), SimpEven=mean(SimpEven))
    
    ##merging all measures of diversity
    alldiv<-merge(distances, divmeasure, by=c("exp_year","treatment"))

    #pasting dispersions into the dataframe made for this analysis
    for.analysis=rbind(alldiv, for.analysis)  
}


#write csv

# write.csv(for.analysis, 'DiversityMetrics_Nov2017.csv')
