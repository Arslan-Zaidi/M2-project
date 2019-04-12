library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)


hq.counts<-fread("files/Database/Heteroplasmy_tables/hq_counts_adj_frequency_09272018.txt",header=T)
hq.counts$individual_het_id<-paste(hq.counts$individual_id,hq.counts$position,sep="_")
hq.counts$mother_het_id<-paste(hq.counts$mother_id,hq.counts$position,sep="_")
#keep only heteroplasmies
hq<-hq.counts[,c('FID','mother_id','mot_cat','mother_het_id','individual_id','individual_het_id','tissue_id','position','level','tissue','major','minor','adj.f','age_collection','age_birth','mot_cat')]


#plot tissues against each other
dhq.2t<-dcast(hq,FID+mother_id+mot_cat+mother_het_id+individual_id+level+individual_het_id+age_collection+position~tissue,value.var="adj.f",fun.aggregate = mean,na.rm=T)
dhq.2t$diff.f<-abs(dhq.2t$bl-dhq.2t$ch)
dhq.2t$avg.f<-apply(dhq.2t[,c('bl','ch')],1,mean)
plt.2t<-ggplot(dhq.2t,aes(bl,ch))+
  geom_point(alpha=0.6)+
  theme_bw()+
  geom_abline(intercept=0,slope=1,color="red")+
  stat_smooth(method="lm",se=F)

#calculate correlation between bl and ch frequencies
#filter only individuals who are heteroplasmic in at least one tissue
dhq.2t.poly<-dhq.2t%>%
  filter((bl>=0.01 & bl<=0.99)|(ch>=0.01 & ch<=0.99))

cor1<-round(with(dhq.2t.poly,cor(bl,ch,method="pearson")),2)

#plot mother against child
dhq.mothers<-dhq.2t[which(dhq.2t$mot_cat!=""),]
dhq.children<-dhq.2t[which(dhq.2t$mot_cat==""),]
dhq.mc<-merge(dhq.mothers,dhq.children[,c('mother_id',"mother_het_id","individual_id","level","bl","ch","diff.f","avg.f")],by.x=c("individual_het_id"),by.y=c("mother_het_id"),all.y=T,suffixes=c(".m",".c"))

#correlation between heteroplasmy frequency b/w child's tissues
dhq.c.poly<-dhq.mc%>%
  distinct(FID,individual_id.c,position,.keep_all = T)%>%
  filter((bl.c>=0.01 & bl.c<=0.99) | (ch.c>=0.01 & ch.c<=0.99))
cor2<-round(with(dhq.c.poly,cor(bl.c,ch.c,method="pearson")),2)

#correlation between heteroplasmy frequency b/w mother's tissues
dhq.m.poly<-dhq.mc%>%
  distinct(FID,individual_id.m,position,.keep_all=T)%>%
  filter((bl.m>=0.01 & bl.m<=0.99) | (ch.m>=0.01 & ch.m<=0.99))
cor3<-round(with(dhq.m.poly,cor(bl.m,ch.m,method="pearson")),2)

#correlation between mother and child blood
dhq.mc.bl.poly<-dhq.mc%>%
  filter((bl.m>=0.01 & bl.m<=0.99)|(bl.c>=0.01 & bl.m<=0.99))
cor4<-round(with(dhq.mc.bl.poly,cor(bl.m,bl.c,method="pearson")),2)

#correlation between mother and child cheek
dhq.mc.ch.poly<-dhq.mc%>%
  filter((ch.m>=0.01 & ch.m<=0.99)|(ch.c>=0.01 & ch.m<=0.99))
cor5<-round(with(dhq.mc.ch.poly,cor(ch.m,ch.c,method="pearson")),2)

mother.plt<-
  ggplot(data=dhq.m.poly,aes(bl.m,ch.m))+
  geom_point(alpha=0.7)+
  theme_bw()+
  stat_smooth(method="lm",se=F)+
  geom_abline(intercept=0,slope=1,color="red")+
  labs(x="Mother - blood",y="Mother - cheek")+
  annotate("text",x=0.25,y=0.90,label=paste("r = ",cor2,sep=""),size=7)+
  theme(axis.title=element_text(size=14))

child.plt<-
  ggplot(dhq.c.poly,aes(x=bl.c,y=ch.c))+
  geom_point(alpha=0.7)+
  theme_bw()+
  stat_smooth(method="lm",se=F)+
  geom_abline(intercept=0,slope=1,color="red")+
  labs(x="Child - blood",y="Child - cheek")+
  annotate("text",x=0.25,y=0.90,label=paste("r = ",cor3,sep=""),size=7)+
  theme(axis.title=element_text(size=14))

blood.plt<-
  ggplot(data=dhq.mc.bl.poly,aes(bl.m,bl.c))+
  geom_point(alpha=0.7)+
  theme_bw()+
  stat_smooth(method="lm",se=F)+
  geom_abline(intercept=0,slope=1,color="red")+
  labs(x="Mother - blood",y="Child - blood")+
  annotate("text",x=0.25,y=0.90,label=paste("r = ",cor4,sep=""),size=7)+
  theme(axis.title=element_text(size=14))

cheek.plt<-
  ggplot(data=dhq.mc.ch.poly,aes(ch.m,ch.c))+
  geom_point(alpha=0.7)+
  theme_bw()+
  stat_smooth(method="lm",se=F)+
  geom_abline(intercept=0,slope=1,color="red")+
  labs(x="Mother - cheek",y="Child - cheek")+
  annotate("text",x=0.25,y=0.90,label=paste("r = ",cor5,sep=""),size=7)+
  theme(axis.title=element_text(size=14))

quartet.plt<-grid.arrange(mother.plt,child.plt,blood.plt,cheek.plt)

ggsave("files/Analysis/PNAS_styled_analyses/af_correlation_plots/quartet_plots_09142018.pdf",quartet.plt,height=10,width=10)

