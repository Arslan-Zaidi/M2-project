#Author: AA.Zaidi

library(ggplot2)
library(dplyr)

ddpcr<-read.table("../../../files/Analysis/Denovo_mutations/denovo_ddPCR_results_10122018.txt",header=T,sep="\t")

#calculate average heteroplasmy frequency
avg.ddpcr<-ddpcr%>%
  group_by(FID,individual,tissue,position)%>%
  summarize(ddpcr_f=mean(ddpcr_f),miseq_f=mean(miseq_f))

#calculate correlation b/w Miseq and ddPCR frequency
rcoef<-round(with(avg.ddpcr,cor(miseq_f,ddpcr_f,use="complete.obs")),2)

#calculate average deviation
avg.dev<-round(with(avg.ddpcr,mean(ddpcr_f-miseq_f,na.rm=T)),2)

cor.plt<-ggplot(avg.ddpcr,aes(miseq_f,ddpcr_f))+
  geom_point(size=0.6)+
  stat_smooth(method="lm")+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_bw()+
  geom_segment(aes(x=miseq_f,xend=miseq_f,y=miseq_f,yend=ddpcr_f),alpha=0.6)+
  annotate(geom="text",x=0.01,y=0.075,label=paste("r = ",rcoef,sep=""),size=6)+
  annotate(geom="text",x=0.01,y=0.07,label=paste("Avg. ddPCR - Miseq = ",avg.dev,sep=""),size=6)+
  theme(axis.title=element_text(size=15),axis.text=element_text(size=12))+
  labs(x="Miseq - MAF",y="ddPCR - frequency")

ggsave("../../../files/Analysis/validations/maf1/miseq_v_ddpcr_validations.pdf",cor.plt,height=7,width=7)

#calculate avg. per individual
avg.ddpcr.ind<-avg.ddpcr%>%
  group_by(FID,individual,position)%>%
  summarize(ddpcr_f=mean(ddpcr_f),miseq_f=mean(miseq_f,na.rm=T))

#in all four cases, only the individual in whom the heteroplasmy was discovered was found to be heteroplasmic (MAF>1%) using ddPCR
#no other individual in the family was heteroplasmic

#hair analysis

