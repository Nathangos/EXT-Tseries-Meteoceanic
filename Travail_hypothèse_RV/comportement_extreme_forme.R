rm(list=ls())
set.seed(133)

require(ggplot2)
Taille_p<-50
Seuil<-4.22
Scale<-0.96
  
p_normale<-sort(seq.int(from = 0.008,to = 0.20,length.out = Taille_p)^(-1))

fonction_GPD_exp_quantile<-function(scale,shape,seuil,p){
  if(abs(shape)<10^{-4}){
    part1<-log(1/p)
    return(part1*scale+seuil)
  }
  part1<-p^(-shape)-1
  return((part1/shape)*scale+seuil)
}

l_gamma<-seq.int(from = -0.5,to = 0.5,by = 0.25)
valeurs_associees_theorique<-c()
val_gamma<-c()
val_1_p<-c()

for(gamma in l_gamma){
  valeurs_associees_theorique<-c(valeurs_associees_theorique,sapply(p_normale^(-1),FUN = fonction_GPD_exp_quantile,scale=Scale,shape=gamma,seuil=Seuil))
  val_gamma<-c(val_gamma,rep(gamma,Taille_p))
  val_1_p<-c(val_1_p,p_normale)
}

df_gamma_val<-cbind.data.frame(valeurs_associees_theorique,val_gamma,val_1_p)
colnames(df_gamma_val)<-c("niveau_retour","gamma","un_sur_proba_depassement")
df_gamma_val$gamma<-as.character(df_gamma_val$gamma)

ggplot(data = df_gamma_val,aes(x=un_sur_proba_depassement,y=niveau_retour,group=interaction(gamma),col=gamma))+
  geom_line()+
  geom_point()+
  scale_x_continuous(trans = "log10")+
  ggtitle(paste("Evolution du niveau de retour (log(x)) selon", "\u03B3", "pour des simulations"))+
  ylab("valeur")+
  xlab("1/p")
  