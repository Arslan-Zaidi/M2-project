

q.age_collection<-c(0,20,40,60,80,100)

#mothers
hq.mothers<-hq%>%
  filter(mot_cat!="")%>%
  mutate(mother_group=individual_het_id,
         age_birth_cat=NA,
         age_collect_cat=cut(age_collection/365,breaks=q.age_collection,ordered_result=T,include.lowest=T))

#their children
hq.children<-hq%>%
  filter(mother_id%in%hq.mothers$individual_id)

#get age of mother at childbirth for these children
q.age_birth<-c(0,20,30,40,50,60)

hq.children<-hq.children%>%
  mutate(mother_group=mother_het_id,
         age_birth_cat=cut(age_birth/365,breaks=q.age_birth,ordered_result=T,include.lowest=T),
         age_collect_cat=cut(age_collection/365,q.age_collection,ordered_result=T,include.lowest=T))

#rbind the two dataframes
hq.concat<-rbind(hq.mothers,hq.children)

#split hq data by mother_id
shq<-split(hq.concat,f=hq.concat$mother_group)
