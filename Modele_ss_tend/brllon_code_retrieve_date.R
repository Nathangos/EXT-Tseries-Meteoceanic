Indices_chosen<-read.csv("residus_clust/inds_extremes_donnees_S.csv")[,2]
Dates_extremes<-read.csv(file="ss_tend/HIVER/dates_prises.csv")[Indices_chosen,2]

TEMPS<-lapply(dataFiles,f_nom_variable,"Date")
intermed<-t(do.call(cbind.data.frame,TEMPS))
rownames(intermed)<-1:nrow(intermed)

# Date pic de la marée.  --------------------------------------------------
date_compar<-(intermed[,19])
df<-data.frame(date=date_compar)
Nber_obs<-as.data.frame(df %>% group_by(date) %>% summarise(n()))

resultat<-sapply(Dates_extremes,function(x,y){return(which(y==substr(x,1,10)))},y=date_compar)
resultat<-sapply(c(1:length(resultat)),FUN=function(x){
  cle<-names(resultat)
  valeur<-resultat[[x]]
  if(length(valeur)==1){
    return(valeur)
  }
  else{
    L<-nchar(cle[x])-1
    partie_j<-substr(cle[x],start = 13,stop = L)
    if(partie_j=="matin"){
      return(valeur[1])
    }
    else{
      return(valeur[2])
    }
  }
})
Indices_obtenus<-as.numeric(resultat)