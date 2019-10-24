dhq.mothers<-dcast(hq.mothers,fam_het_id+FID+fam_cat+mother_id+mother_het_id+individual_id+individual_het_id+position~tissue,value.var="adj.f")

#write function to sample a single individual from within each family
#and output all heteroplasmies that exist in both tissues of these individuals

samplegerm<-function(){
  #sampled 1 individual at random from each family
  sampled_inds=allfam%>%
    distinct(FID,individual_id)%>%
    group_by(FID)%>%
    sample_n(size=1)
  
  #select these individuals and select heteroplasmies that are present in both tissues of these individuals
  allfam2<-allfam%>%
    filter(individual_id%in%sampled_inds$individual_id)%>%
    group_by(individual_id,position)%>%
    summarize(ngerm=length(which(adj.f>=0.01)))%>%
    filter(ngerm==2)
  
  #return the number of hetereoplasmies across all individuals present in both tissues
  return(nrow(allfam2))
  
}

num_germs<-replicate(1000,samplegerm())

write.table(num_germs,
            here("num_germline_mutations_random_08202019.txt"),
            col.names=F,
            row.names=F,
            quote=F)
  