rm(list=ls())
set.seed(133)
require(reshape2)
require(ggplot2)
require(gridExtra)
require(grid)
require(dplyr)
require(threshr)
require(extRemes)
require(ggside)

M<-2000
vecteur_u<-runif(M)
qtile_ext<-0.90
mu<-3
sigma<-1
seuil<-qnorm(p = qtile_ext,mean = mu,sd = sigma)
print(seuil)
ech_ev<-1
gam_ev<--0.05

Generateur_beta<-function(u,quant,sigma,mu,seuil_ev,ech_ev,gam_ev){
  if(u<quant){
    return(qnorm(u,mean = mu,sd = sigma))
  }
  else{
    qprime<-(u-quant)/(1-quant)
    return(extRemes::qevd(p=qprime,
                          threshold =seuil_ev,shape = gam_ev,
                          scale=ech_ev,
                          type="GP"))
  }
}
vect_beta<-sapply(vecteur_u,Generateur_beta,quant=qtile_ext,
       gam_ev=gam_ev,ech_ev=ech_ev,sigma=sigma,
       mu=mu,seuil_ev=seuil)
plot(density(vect_beta))
abline(v=seuil,col="red")

# Generation courbes ------------------------------------------------------
c<-1
L_courbe<-37
begg<-1
end<-2

# change the peak position ------------------------------------------------

DPHASE<-TRUE
Generation_courbe<-function(beta,c,L_courbe,begg,end,dephase){
  alpha_b<--beta/3
  vect_time<-seq.int(from = begg,to = end,length.out = L_courbe)
  if(dephase==TRUE){
    delta<-rnorm(1,sd=2)/L_courbe
    vect_time<-vect_time+delta
  }
  return(alpha_b*vect_time^2+beta*vect_time+c)
}

vect_courbes<-t(sapply(vect_beta,FUN = Generation_courbe,c=c,
       L_courbe=L_courbe,begg=begg,end=end,dephase=DPHASE))

# Illustration ------------------------------------------------------------
Nshown<-50
percent_studied<-0.05
Q1<-as.numeric(apply(X = vect_courbes,MARGIN = 2,
          FUN = function(x){return(quantile(x,percent_studied/2))}))
Q2<-as.numeric(apply(X = vect_courbes,MARGIN = 2,
          FUN = function(x){return(quantile(x,1-percent_studied/2))}))

Mean<-colMeans(vect_courbes)
inds_taken<-sample(c(1:M),size = Nshown)
Curves_seen<-vect_courbes[inds_taken,]
Melteur<-melt(t(Curves_seen))
colnames(Melteur)<-c("Time","number","value")
Melteur$number<-as.character(Melteur$number)
Melteur$Q1<-Q1[Melteur$Time]
Melteur$Q2<-Q2[Melteur$Time]
Melteur$mean_curve<-Mean[Melteur$Time]


ggplot(data=Melteur,aes(x=Time,y=value,col=number,
                        group=interaction(number),
                        ))+
  geom_ribbon(mapping=aes(ymin=Q1,ymax=Q2,linetype="confidence_band"),
              alpha=0.02,
              fill="grey",col="darkblue")+
  ylab("value (m)")+
  geom_line(aes(y=mean_curve,linetype="mean"),col="red")+
  geom_line(alpha=0.7)+
  guides(col="none")+
  scale_linetype_manual("Legend",values=c("confidence_band"=2, 
                                          "mean"=5))+
  ggtitle(paste0(Nshown," curves shown among the ",M," simulated"))


