library(data.table)
library(ggplot2)
library(plyr)


#### load heteroplasmies - conservative (one mother and two children families only) #####

hq.counts<-fread("../../../files/Analysis/Fixing_heteroplasmy_table/hqcounts_cleaned_conservative_09272018.txt",header=T,sep="\t")

fam<-fread("../../../files/Analysis/Fixing_heteroplasmy_table/famfile_cleared_conservative_09272018.txt",header=T,sep="\t")

hq.counts<-base:::merge(hq.counts,fam[which(duplicated(fam$individual_id)=="FALSE"),c("individual_id","twin_info")],by="individual_id",all.x=T,sort=F)

############ Adjust allele frequency for each sample such that it is consistent (for an allele) across all samples in a family for the same heteroplasmy ###########

#create a table of 'reference' alleles for each site in a family based on the heteroplasmy table
#this reference table will be used to adjust the allele frequency for each sample accordingly

#dhets is a table which contains the 'reference' major and minor alleles
dhets<-ddply(hq.counts,.(FID,position,fam_het_id),function(x){
  major.list<-x$major[which(x$heteroplasmy=="yes")]
  major<-sample(major.list,1)
  minor.list<-x$minor[which(x$heteroplasmy=="yes" & x$major==major)]
  minor<-sample(minor.list,1)
  return(data.frame(major=major,minor=minor))
})


#define function to adjust allele frequency of all individuals in a family based on this reference
adj.af<-function(x){
  fid<-as.character(x$FID[1])
  position<-x$position[1]
  major<-as.character(dhets$major[which(dhets$FID==fid &dhets$position==position)])
  minor<-as.character(dhets$minor[which(dhets$FID==fid &dhets$position==position)])
  same<-paste(major,minor,sep="")
  flip1<-paste(minor,major,sep="")
  flip2<-paste(minor,".",sep="")
  lost<-paste(major,".",sep="")
  x$adj.f<-NA
  for(i in 1:nrow(x)){
    if(paste(x$major[i],x$minor[i],sep="")%in%c(same,flip1,flip2,lost)){
      if(paste(x$major[i],x$minor[i],sep="")==same){x$adj.f[i]<-x$maf[i]}
      if(paste(x$major[i],x$minor[i],sep="")%in%c(flip1,flip2)){x$adj.f[i]<-(1-x$maf[i])}
      if(paste(x$major[i],x$minor[i],sep="")==lost){x$adj.f[i]<-0}  
    }
    if(!(paste(x$major[i],x$minor[i],sep="")%in%c(same,flip1,flip2,lost))){
      x$adj.f[i]<-0
      #if(x$major[i]==major & x$maf[i]< 0.01){x$adj.f[i]<-0}
      #else{x$adj.f[i]<-NA}
    }
  }
  return(x)}

#go through each family and 'fix' the allele frequency according to the reference alleles
hq.counts<-ddply(hq.counts,.(FID,position),adj.af)

#check whether this has been done correctly - the plot should look like half of an X - no off-diagonals should be observed
ggplot(hq.counts,aes(maf,adj.f))+geom_point()+theme_bw()

fwrite(hq.counts,"../../../files/Database/Heteroplasmy_tables/hq_counts_adj_frequency_09272018.txt",sep="\t",col.names=T,row.names=F,quote=F)

#write major/minor allele codes for each heteroplasmy
fwrite(dhets,"../../../files/Database/hq_ref_alt_code_11102018.txt")
