Inds_not_null<-which(resultat_Simul_Obs[,2]>quantile(resultat_Simul_Obs[,2],0.10))
Notnull<-resultat_Simul_Obs[Inds_not_null,2]
Q1<-0.10
Q2<-0.98
NTHS<-15
Min_u<-quantile(Notnull,Q1)
Max_u<-quantile(Notnull,Q2)
Result_threshr<-Outils_Threshr_choix(Q1 =Q1 ,
                                     Q2 =Q2,
                                     NT_ths = NTHS,
                                     variable = Notnull,
                                     N_v=2)
Result<-plot(Result_threshr)

Mean_found<-apply(Result$y,MARGIN = 1,FUN = mean)
Mean_found
Whole_q<-seq(Q1, Q2,length.out=NTHS)
qtile_ext<-0.78
NB_TH_for_graphs<-150
# vect_k<-c(20:500)
# Graphics_estimators_gamma(series = Notnull,vect_k = vect_k,
#                           Title_graphic  = "Extreme of L2")

POT::tcplot(Notnull,ask = FALSE,u.range = c(Min_u,Max_u),
            nt =NB_TH_for_graphs)
POT::mrlplot(data = Notnull,ask = FALSE,u.range = c(Min_u,Max_u),
             nt = NB_TH_for_graphs)
abline(v=quantile(Notnull,qtile_ext),col="red")
n.dens<-2^(14)
kernel_dens<-density(x=Notnull,n=n.dens)
threshold<-as.numeric(quantile(Notnull,qtile_ext))
threshold
extremes<-subset(Notnull,Notnull>threshold)
index<-which.max(threshold-kernel_dens$x<0)-1
value_fonc_u<-as.numeric(kernel_dens$y[index])
scale_fonc<-(p_u/value_fonc_u)
k<-which.max(Notnull-threshold<0)-1
gamma<-get.tail.index.estimate(Notnull,1-p_u,scale_fonc)
gap<-diff(kernel_dens$x)[1]

F_dens<-cumsum(kernel_dens$y)*gap
non_ext<-subset(Notnull,Notnull<=threshold)
intervals<-c(-Inf,kernel_dens$x+gap/2)
interval_vecteur<-findInterval(non_ext,vec = intervals)

indices_Excedent<-which(Notnull>threshold)
vector_unif<-rep(NA,length(Notnull))
Excesses<-Notnull[indices_Excedent]-threshold

Fcraft<-function(x,xi,sigma){
  return(pgp_craft(x=x,xi=xi,sigma=sigma))
}
Unif_excesses<-(qtile_ext)+(1-qtile_ext)*sapply(Excesses,Fcraft,
                                                sigma = scale_fonc,xi = gamma)
Unif_not_null<-Notnull
Unif_not_null[-indices_Excedent]<-F_dens[interval_vecteur]
Unif_not_null[indices_Excedent]<-Unif_excesses
plot(density(Unif_not_null))
summary(Unif_not_null)
goftest::ad.test(Notnull,null = "punif")