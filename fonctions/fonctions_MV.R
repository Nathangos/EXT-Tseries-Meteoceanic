fonction_Stheta_MV<-function(L,indices_select,Tau,Params_Matern,liste_noms,obs_P){
  S_theta<-0
  d<-length(Tau)
  for (i in 1:d){
    for(j in 1:d){
      Beta<-Params_Matern[1]
      Sigma_zero<-Params_Matern[2]
      Kappa<-Params_Matern[3]
      f_vario_nul<-function(h){
        numerateur<-Sigma_zero^(2)*(abs(h)/Kappa)^(2)
        denominateur<-((abs(h)/Kappa)^(2)+1)^(Beta)
        return(numerateur/denominateur)
      }
      vecteur_temps<-(1:L)/L
      a<-Params_Matern[4]
      f_vario_Bessel_ij<-function(h,i,j){
        idx<-h==0
        sigma_ij<-(Params_Matern[4+j]*Params_Matern[4+i])
        nu_ij<-(1/2)*(Params_Matern[6+j]+Params_Matern[6+i])
        kappa_ij<-nu_ij
        lambda_ij<-a/(2*sqrt(nu_ij))
        ans<-sigma_ij
        ans[!idx]<-sigma_ij * 2^(1 - kappa_ij)/gamma(kappa_ij) * (h[!idx]/lambda_ij)^kappa_ij * 
                 besselK(h[!idx]/lambda_ij, kappa_ij)
        return(ans)
      }
      
      f_vario_wind_gusts<-function(t,s){
        h<-abs(t-s)
        valeur_WG<-f_vario_nul(h)+(1/2)*(f_vario_Bessel_ij(0,i=i,j=i)+f_vario_Bessel_ij(0,i=j,j=j))-f_vario_Bessel_ij(h,i=i,j=j)
        return(valeur_WG)
      }
      Matrice_couples<-list()
      z<-1
      for(a in 1:L){
        for (b in 1:L){
          Matrice_couples[[z]]<-c(a,b)
          z<-z+1
        }
      }
      nom1<-liste_noms[i]
      nom2<-liste_noms[j]
      PIST_empirique<-extremogram_empirique_MV(Matrice_couples=Matrice_couples,indices = indices_select,liste_tau_MV = Tau,nom1=nom1,nom2=nom2,donnees=obs_P)
      PIST_theorique<-sapply(Matrice_couples,FUN=pi_s_t_theorique,vecteur_temps=vecteur_temps,f_variogram=f_vario_wind_gusts)
      S_theta<-sum((PIST_theorique-PIST_empirique)^2)+S_theta
    }
  }
  return(S_theta)
}
extremogram_empirique_MV<-function(Matrice_couples,indices,donnees,liste_tau_MV,nom1,nom2){
  PIST_empirique<-sapply(Matrice_couples,FUN =pi_s_t_MV,indices=indices,donnees=donnees,liste_tau_MV=liste_tau_MV,nom1=nom1,nom2=nom2)
  return(PIST_empirique)
}

pi_s_t_MV<-function(indices,donnees,liste_tau_MV,nom1,nom2,couples_indices){
  t<-couples_indices[1]
  s<-couples_indices[2]
  vecteur_u1<-liste_tau_MV[[nom1]]
  vecteur_u2<-liste_tau_MV[[nom2]]
  premier_processus_s<-donnees[[nom1]][indices,s]
  second_processus_t<-donnees[[nom2]][indices,t]
  numerateur<-sum(as.numeric((premier_processus_s>vecteur_u1[s])&(second_processus_t>vecteur_u2[t])))
  denominateur<-sum(as.numeric(premier_processus_s>vecteur_u1[s]))
  return(numerateur/denominateur)
}


Inner_product_MV<-function(liste_MV_simul,l_name){
  N_rows<-length(liste_MV_simul)
  sapply(c(1:N_rows),)
}