library(data.table)
library(dplyr)
library(ggplot2)
library(here)

#setwd("/Users/Azaidi/Documents/M2_new/files/Analysis/Fixing_heteroplasmy_table")

#load fam file
fam<-fread(here("Data_files/famfile_09272018.txt"),header=T,sep="\t")
#fam<-fam[-which(fam$FID==c("Zclonal")),]
fam<-fam[-which(fam$run%in%c("TR13n","TR29n","TR13","TR29")),]

#load heteroplasmy table
hets<-fread(here("Data_files/Allruns_2_33_08222018.header.trim.rd.hets"),header=T,sep="\t")

#starting number of heteroplasmies: 1473

#merge hets with fam file
hets2<-merge(hets,fam,by="fqid",all.x=T)

#some rows/heteroplasmies (N = 9) were not merged with fam file - Fam file has missing info. following are the fastq IDs with explanations for each
#hets2[which(is.na(hets2$FID))]
#M117C1-ch: This family was sequenced for the PNAS paper but was removed for the M2 project becasue we did not have a second child for this family. recoverable if needed
#M211-ch/M211C1-bl/M211C2-ch: This family was sequenced for the PNAS paper but was removed for the M2 project because we did not have a second child for this family. also recoverable if needed
#m308G-bl-0330_S8 and M475c2_S17: These samples were not sequenced for a second tissue so were removed from the fam file

#remove these rows
hets2<-hets2[which(is.na(hets2$FID)=="FALSE"),]
hets2<-hets2[-which(hets2$individual_id%in%c("F254g1","F308m1c2")),]

#remaining hets: 1473

#check if all families are represented in the heteroplasmy file - i.e. does every family have at least one heteroplasmy?

unique(fam$FID)[-which(unique(fam$FID)%in%unique(hets2$FID))]
#F326 does not have detectable heteroplasmies

#remove duplicated samples

#load duplicates pre-TR33
dups<-fread(here("Data_files/duplicates_removed_08092018.txt"),header=F)

#load duplicates sequenced in TR33
dups.33<-fread(here("Data_files/TR33_duplicates.txt"))

#remove duplicated sequences from heteroplasmy file
hets3<-hets2[-which(hets2$fqid%in%dups$V1),]

#881 heteroplasmies left after removing duplicates

#remove older samples sequenced in TR33
hets3<-hets3[-which(hets3$fqid%in%dups.33$fqid_old),]

#797 heteroplasmies remaining

#remove duplicated samples from fam file as well 
fam2<-fam[-which(fam$fqid%in%dups$V1|fam$fqid%in%dups.33$fqid_old),]

# check how many tissues sequenced per person
dfam.2t<-fam2%>%
  group_by(individual_id)%>%
  summarize(ntissue=length(unique(tissue_id)))

# > dfam.2t[which(dfam.2t$ntissue<2),]
# individual_id ntissue
# 70         F171g1       1
# 161        F239g1       1
# 185        F254g1       1
# 188    F254g1m1c2       1
# 202      F261m1c1       1
# 291      F308m1c2       1
# 329      F330m1c1       1
# 332      F332m1c1       1

#remove these individuals from heteroplasmy file
hets4<-hets3[-which(hets3$individual_id%in%c("F171g1",
                                            "F239g1",
                                            "F254g1",
                                            "F254g1m1c2",
                                            "F261m1c1",
                                            "F308m1c2",
                                            "F330m1c2",
                                            "F332m1c1")),]

#789 heteroplasmies remaining

fam3<-fam2[-which(fam2$individual_id%in%c("F171g1",
                                          "F239g1",
                                          "F254g1",
                                          "F254g1m1c2",
                                          "F261m1c1",
                                          "F308m1c2",
                                          "F330m1c2",
                                          "F332m1c1")),]


#how many mothers have 2 children
dfam.2c<-fam3%>%
  group_by(mother_id,fam_str)%>%
  summarize(nchildren=length(unique(individual_id)))
#dfam.2c<-ddply(fam3,.(mother_id,fam_str),summarize,nchildren=length(which(duplicated(individual_id)=="FALSE")))

#dfam.2c[which(dfam.2c$nchildren<2),]
#mothers with less than two children
# mother_id    fam_str nchildren
# 1     F105g1    0-1-1-2         1
# 2     F113g1    0-0-1-2         1
# 3     F116g1    0-1-1-2         1
# 4     F140g1    0-1-1-2         1
# 5     F171g1    0-1-1-2         1
# 6    F172tg1    1-1-1-3         1
# 7     F239g1 0-1-1-3-1t         1
# 8    F249tg1    1-2-1-2         1
# 9     F254g1    0-1-1-2         1
# 10  F254g1m1    0-1-1-2         1
# 11    F261m1    0-0-1-2         1
# 12    F289g1    0-1-1-2         1
# 13    F291g1    0-1-1-3         1
# 14    F297g1    0-1-1-2         1
# 15    F330m1    0-0-1-2         1
# 16    F332m1    0-0-1-2         1

#remove these families
fam4<-fam3[-which(fam3$FID%in%c("F254","F261","F330","F332")),]
hets5<-hets4[-which(hets4$FID%in%c("F254","F261","F330","F332"))]

#775 heteroplasmies remaining

#change fam structure and mother category for these families
#F171 # both tissues for F171g1 not available 
#F239 # both tissues for F239g1 not available

fam4$fam_cat[which(fam4$FID=="F171")]<-"m1c2"
fam4$fam_str[which(fam4$FID=="F171")]<-"0-0-1-2"

fam4$fam_cat[which(fam4$FID=="F239")]<-"m1c3"
fam4$fam_str[which(fam4$FID=="F239")]<-"0-0-1-3-1t"

fam4$fam_cat[which(fam4$FID=="F308")]<-"m1c2"
fam4$fam_str[which(fam4$FID=="F239")]<-"0-0-1-2"

#F308 initially had three children. Since we removed the 2nd child, we have to relabel the 3rd child as c2 under "levels" column
fam4$level[which(fam4$individual_id=="F308m1c3")]<-"c2"

write.table(fam4,here("Data_files/famfile_cleared_conservative_09272018.txt"),sep="\t",col.names=T,row.names=F,quote=F)

#create fq_het_id in hets file
hets5$fq_het_id<-paste(hets5$fqid,hets5$position,sep="_")

#create fam_het_id in hets file
hets5$fam_het_id<-paste(hets5$FID,hets5$position,sep="_")



#### load count file for all family frequencies ####
counts<-fread("https://www.dropbox.com/s/gyk54cgek4kooo4/Allruns_08222018.trim.rd.counts?dl=1",header=F,sep="\t")
colnames(counts)<-c('fqid','ref','position','A','C','G','T','a','c','g','t','cvrg','nalleles','major','minor','maf')

#merge with fam file
counts2<-merge(counts,fam4,by="fqid",all.x=T)

#create fq_het_id in count file
counts2$fq_het_id<-paste(counts2$fqid,counts2$position,sep="_")

#remove duplicated sequences from count file
counts2<-counts2[-which(counts2$fqid%in%dups$V1 | counts2$fqid%in%dups.33$fqid_old),]

#create fam_het_id in count file to extract all frequencies from families
counts2$fam_het_id<-paste(counts2$FID,counts2$position,sep="_")

# #check which rows are missing family info after the merge
# unique(counts3$fqid[(which(is.na(counts3$FID)=="TRUE"))])
# [1] "M117-bl"            "M117-ch"            "M117C1-bl"          "M117C1-ch"          "M211-bl"            "M211-ch"           
# [7] "M211C1-bl"          "M211C1-ch"          "TR1329M_S24"        "TR23S24_S24"        "TR24S15_S15"        "TR24S24_S24"       
# [13] "TR25S12_S12"        "TR25S21_S21"        "TR26S9_S9"          "TR27S3_S3"          "TR27S6_S6"          "TR28S12_S12"       
# [19] "TR28S13_S13"        "TR28S15_S15"        "TR28S16_S16"        "TR28S17_S17"        "TR28S18_S18"        "TR28S19_S19"       
# [25] "TR28S1_S1"          "TR28S24_S24"        "TR28S3_S3"          "TR28S4_S4"          "TR28S6_S6"          "TR28S7_S7"         
# [31] "TR30S17_S17"        "TR30S1_S1"          "TR31S24_S24"        "TR33S10_S9"         "TR33S11_S10"        "TR33S12_S11"       
# [37] "TR33S13_S12"        "TR33S14_S13"        "TR33S15_S14"        "TR33S16_S15"        "TR33S17_S16"        "TR33S18_S17"       
# [43] "TR33S19_S18"        "TR33S1_S1"          "TR33S20_S19"        "TR33S21_S20"        "TR33S22_S21"        "TR33S23_S22"       
# [49] "TR33S24_S23"        "TR33S2_S2"          "TR33S4_S3"          "TR33S5_S4"          "TR33S6_S5"          "TR33S7_S6"         
# [55] "TR33S8_S7"          "TR33S9_S8"          "m205G-ch_S6"        "m287G-bl-0421_S20"  "m307-bl-0330_S6"    "m307-ch-0311_S1"   
# [61] "m307c1-bl-0330_S10" "m307c1-ch-0315_S22" "m307c2-bl-0330_S23" "m475c2_S17"         "m499-bl-0330_S19"   "m499-ch-0626_S22"  
# [67] "m499c1-bl-0330_S14" "m499c2-bl-0407_S7"  "m499c2-ch-0125_S6" 
#these are either clonses or samples that were removed b/c they did not have both tissues sequenced or two children
#remove these rows
#counts2<-counts2[which(is.na(counts2$FID)=="FALSE"),]

#remove families which are not present in hets3 file i.e. have no heteroplasmies after cleaning
counts3<-counts2[which(counts2$fam_het_id%in%hets5$fam_het_id),]

#add column to indicate whether heteroplasmy was discovered in that sample
counts3$heteroplasmy<-NA

counts3$heteroplasmy[which(counts3$fq_het_id%in%hets5$fq_het_id)]<-"yes"
counts3$heteroplasmy[which(is.na(counts3$heteroplasmy)==TRUE)]<-"no"
counts3$tissue_het_id<-paste(counts3$tissue_id,counts3$position,sep="_")

#tabulate for each heteroplasmy how many times in the family it was discovered as a heteroplasmy and how many times it had MAF>0.01
dcounts3<-counts3%>%
  group_by(position)%>%
  summarize(nmaf=length(which(maf>=0.01)),
            nhets=length(which(heteroplasmy=="yes")))
#dcounts3<-ddply(counts3,.(position),summarize,nmaf=length(which(maf>=0.01)),nhets=length(which(heteroplasmy=="yes")))

#add column listing the difference between these two numbers - this could highlight 'problematic' sites
dcounts3$diff<-with(dcounts3,nmaf-nhets)
dcounts3<-dcounts3[order(-dcounts3$diff),]

# #couple of observations from this:
# 1. nmaf is always equal to or higher than nhets:
#   - this means that often heteroplasmies are discovered at MAF>0.01 which do not always pass strand bias cutoffs
# 2. Sites whcih show a large difference between nhets and nmaf also tend to be previously known problematic sites (e.g. 185, 316, 204, 4233)

#tabulate number of families heteroplasic for the same site - to identify other problematic sites
dcounts.shared<-counts3%>%
  group_by(position)%>%
  summarize(nfams.hets=length(unique(FID[which(heteroplasmy=="yes")])),
            nfams.maf=length(unique(FID[which(maf>=0.01)])))

# dcounts.shared<-ddply(counts3,.(position),function(x){
#   nfams.hets=length(unique(x$FID[which(x$heteroplasmy=="yes")]))
#   nfams.maf=length(unique(x$FID[which(x$maf>=0.01)]))
#   return(data.frame(nfams.maf,nfams.hets))
# })


#dcounts.shared<-ddply(dcounts3,.(position),summarize,nmaf=sum(nmaf),nhets=sum(nhets),diff=sum(diff))
#dcounts.shared<-dcounts.shared[order(-dcounts.shared$nmaf),]

#identify sites with high tv/ts ration
hets5$tstv<-NA
for(i in 1:nrow(hets5)){
  if(paste(hets5$major[i],hets5$minor[i])%in%c("A G","G A","C T","T C")){
    hets5$tstv[i]<-"ts"
  }else{
    hets5$tstv[i]<-"tv"
  }
}

dhets5.tstv<-hets5%>%
  group_by(position)%>%
  summarize(nts=length(which(tstv=="ts")),
            ntv=length(which(tstv=="tv")))%>%
  arrange(-ntv)

# dhets5.tstv<-ddply(hets5,.(position),summarize,nts=length(which(tstv=="ts")),ntv=length(which(tstv=="tv")))
# dhets5.tstv<-dhets5.tstv[order(-dhets5.tstv$ntv),]


#sites that are "problematic" based on test1:
#dcounts3[which(dcounts3$diff>5),]
# position nmaf nhets diff
# 92     4233   43    12   31
# 30      567    8     1    7
# 10      185   12     6    6


#sites are that are "problematic" based on test2 (more than 9 fams share hets):
# dcounts.shared[dcounts.shared$nfams.maf>5,]
# position nfams.maf nfams.hets
# 7        152         7          7
# 10       185         6          6
# 11       189         6          6
# 15       204         8          8
# 20       215         6          6
# 25       316         6          6
# 92      4233        16         16
# 298    16093         6          6


#sites that have more transversions than transitions
#dhets5.tstv[which(dhets5.tstv$ntv-dhets5.tstv$nts>5),]
# position nts ntv
# 92      4233   0  19
# 42      1291   0   8
# 25       316   0   6
# 306    16182   0   6

#remove all these sites
#make list of unique sites from all three list
x<-sort(unique(c(dcounts3$position[which(dcounts3$diff>5)],
               dcounts.shared$position[dcounts.shared$nfams.maf>5],
               dhets5.tstv$position[which(dhets5.tstv$ntv-dhets5.tstv$nts>5)])))

# [1]  4233   567   185   152   189   204   215   316 16093 16182  1291

#REMOVE THESE HETEROPLASMIES
hets6<-hets5[-which(hets5$position%in%x),]

#REMOVE THESE POSITIONS FROM COUNT FILE
counts4<-counts3[-which(counts3$position%in%x),]

#make sure there are no families with no heteroplasmies
dcounts4<-counts4%>%
  group_by(fam_het_id)%>%
  summarize(nhets=length(which(heteroplasmy=="yes")))

#dcounts4<-ddply(counts4,.(fam_het_id),summarize,nhets=length(which(heteroplasmy=="yes")))
#which(dcounts4$nhets<1)
#no empty families. each famiy has at least one heteroplasmy in this table - good

#find out which people have less than 2 tissues
dcounts4.2t<-counts4%>%
  group_by(individual_id)%>%
  summarize(ntissues=length(unique(tissue_id)))

dcounts4.2c<-counts4%>%
  group_by(mother_id)%>%
  summarize(nchildren=length(unique(individual_id)))

#dcounts4.2t<-ddply(counts4,.(individual_id),summarize,ntissues=length(unique(tissue_id)))
#dcounts4.2c<-ddply(counts4,.(mother_id),summarize,nchildren=length(unique(individual_id)))

#write cleaned table to file
fwrite(counts4,
       here("Data_files/hqcounts_cleaned_conservative_09272018.txt"),
       sep="\t",
       col.names=T,
       row.names=F,
       quote=F)

#write cleaned heteroplasmy table to file
fwrite(hets6,
       here("Data_files/Fixing_heteroplasmy_table/hq_cleaned_conservative_09272018.txt"),
       sep="\t",
       col.names=T,
       row.names=F,
       quote=F)


############ Adjust allele frequency for each sample such that it is consistent (for an allele) across all samples in a family for the same heteroplasmy ###########

#create a table of 'reference' alleles for each site in a family based on the heteroplasmy table
#this reference table will be used to adjust the allele frequency for each sample accordingly

#dhets is a table which contains the 'reference' major and minor alleles
#write function to generate major and minor alleles for each FID, position combo
f_refalleles<-function(x){
  major.list<-x$major[which(x$heteroplasmy=="yes")]
  major<-sample(major.list,1)
  minor.list<-x$minor[which(x$heteroplasmy=="yes" & x$major==major)]
  minor<-sample(minor.list,1)
  return(data.table(major=major,minor=minor))
}

dhets<-counts4%>%
  group_by(FID,position,fam_het_id)%>%
  do(f_refalleles(.))

# dhets<-ddply(counts4,.(FID,position,fam_het_id),function(x){
#   major.list<-x$major[which(x$heteroplasmy=="yes")]
#   major<-sample(major.list,1)
#   minor.list<-x$minor[which(x$heteroplasmy=="yes" & x$major==major)]
#   minor<-sample(minor.list,1)
#   return(data.frame(major=major,minor=minor))
# })

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
hq.counts<-counts4%>%
  group_by(FID,position)%>%
  do(adj.af(.))

# hq.counts<-ddply(counts4,.(FID,position),adj.af)

#check whether this has been done correctly - the plot should look like half of an X - no off-diagonals should be observed
ggplot(hq.counts,aes(maf,adj.f))+geom_point()+theme_bw()

fwrite(hq.counts,here("Data_files/hq_counts_adj_frequency_09272018.txt")
       ,sep="\t",
       col.names=T,
       row.names=F,
       quote=F)

