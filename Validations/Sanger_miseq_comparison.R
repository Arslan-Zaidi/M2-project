library(ggplot2)
sanger<-read.table("../../../files/Analysis/validations/maf10/sanger_miseq_comparison.txt",header=T,sep="\t")

#calculate correlation b/w sanger and miseq frequency
f.cor<-with(sanger,cor(maf,adjusted_sanger_f))
#calculate average error between sanger and miseq frequency
f.error<-with(sanger,median(maf-adjusted_sanger_f))

cor.plt<-ggplot(sanger,aes(maf,adjusted_sanger_f))+
  geom_point(size=0.6)+
  stat_smooth(method="lm",se=T)+
  geom_abline(intercept=0,slope=1,color="red")+
  theme_bw()+
  labs(x="Miseq - MAF",y="Sanger - Frequency")+
  annotate(geom="text",x=0.2,y=0.5,label="r = 0.95",size=5)+
  annotate(geom="text",x=0.2,y=0.45,label="Avg. Sanger - Miseq = 0.01",size=5)+
  theme(axis.title=element_text(size=14),axis.text=element_text(size=14))+
  geom_segment(aes(x=maf,xend=maf,y=maf,yend=adjusted_sanger_f),alpha=0.5)

ggsave("../../../files/Analysis/validations/maf10/sanger_miseq_comparison_10122018.pdf",cor.plt,height=10,width=10)

