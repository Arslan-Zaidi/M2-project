hq.syn<-fread(here("Data_files/hq_adjf_ancestral_syn_11282018.txt"),header=T,sep="\t")

hq.syn2<-hq.syn%>%
  filter(gene%in%c("ATP6","ATP8","COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND5","ND6"))%>%
  group_by(gene,syn)%>%
  distinct(position,major,minor)%>%
  summarize(n=length(position))%>%
  dcast(gene~syn,value.var="n")%>%
  replace_na(list(syn=0,nsyn=0))%>%
  rbind(data.frame(gene="total",nsyn=sum(.$nsyn),syn=sum(.$syn)))


fwrite(hq.syn2,here("Mutation_spectrum/syn_summary_table.txt"),
       sep="\t",
       col.names = T,
       row.names=F,
       quote=F)
