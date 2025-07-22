#devtools::install_github("psyteachr/introdataviz")
require(introdataviz)
fonction_effet_lifting_periode_temps<-function(origine,sorties_lifting,periode_temps,nom){
  Copie<-origine[,periode_temps]
  N<-nrow(Copie)
  WHOLE<-melt(t(rbind.data.frame(sorties_lifting[,periode_temps],Copie)))
  colnames(WHOLE)<-c("Temps","Individu","Mesure")
  WHOLE$Temps<-as.factor(WHOLE$Temps)
  WHOLE$type<-ifelse(WHOLE$Individu<=N,"Post_l","Ante_l")
  TitreV<-paste0("Graphes Violin des données pour ",nom," avant et après lifting(",periode_temps[1],",",periode_temps[-1],")")
  PV<-ggplot2::ggplot(data = WHOLE,aes(Temps,Mesure,fill=type))+
    geom_split_violin()+
    ggtitle(TitreV)
  return(PV)
}
fonction_effet_lifting_norme<-function(origine,sorties,coeurs,nom){
  N_origine<-unlist(parApply(cl = coeurs,X = origine,MARGIN = 1,FUN = calcul_norme_L2))
  N_sortie<-unlist(parApply(cl = coeurs,X = sorties,MARGIN = 1,FUN = calcul_norme_L2))
  Changement_opere<-N_sortie/N_origine
  plot(density(Changement_opere),main=paste0("Modification de la norme post hissage pour ",nom))
  abline(v = 1)
}
