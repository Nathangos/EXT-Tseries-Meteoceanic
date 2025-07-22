# Modèles ML --------------------------------------------------------------

# # -----------------------------------------------------------------------

fnct_ML_realite_simulations_all_kernels<-function(ech_test,ech_train,variable){
  LISTE<-list()
  Matrice_perfomances<-matrix(NA,nrow=4,ncol=4)
  j<-1
  for(type_noyau in c("linear","sigmoid","polynomial","radial")){
    resultat<-fnct_ML_realite_simulations(ech_train = ech_train,ech_test=ech_test,
                                          type_noyau = type_noyau)
    LISTE[[type_noyau]]<-resultat
    Mat_confTest<-resultat$resultat_conf_test
    df<-matrix(0,nrow=2,ncol=2)
    if(ncol(Mat_confTest)!=2){
      indice<-as.numeric(colnames(Mat_confTest))+1
      df[,indice]<-Mat_confTest
    }
    else{
      df<-matrix(c(Mat_confTest),nrow = 2,ncol=2)
    }
    Total<-colSums(df)
    Proportion_simulations_predite<-sum(df[,2])/sum(df)
    Proportion_simulations_realite<-sum(df[2,])/sum(df)
    # ligne =vérité. ----------------------------------------------------------
    Sensilite<-df[2,2]/sum(df[2,])
    Specificite<-df[1,1]/sum(df[1,])
    
    # observer parmi les gens prédits positifs ---------------------------------
    Accuracy<-df[2,2]/sum(df[,2])
    df<-rbind(df,Total)
    df<-cbind.data.frame(df,rep(variable,nrow(df)),rep(type_noyau,nrow(df)))
    rownames(df)<-c("obs realite","simul realite","Total")
    colnames(df)<-c("obs pred","simul pred","variable","noyau")
    LISTE[[type_noyau]]$resultat_conf_test<-df
    
    F1_score<-2*(Accuracy*Sensilite)/(Accuracy+Sensilite)
    Matrice_perfomances[j,]<-c(type_noyau,round(F1_score,2),round(Proportion_simulations_predite,2),
                               round(Proportion_simulations_realite,2))
    j<-j+1
  }
  colnames(Matrice_perfomances)<-c("Noyau utilise","F1 score",
                                   "Prop_simulations_pred","Vraie_proportion")
  LISTE$performances_globales<- Matrice_perfomances
  return(LISTE)
}

fnct_logistigue_glm_mult_links<-function(ech_test,ech_train,variable,liste_fncts){
  LISTE_resultat<-list()
  Matrice_perfomances<-matrix(NA,nrow=length(liste_fncts),ncol=4)
  j<-1
  Matrice_perfomances[,1]<-liste_fncts
  for(fnct_link_name in liste_fncts){
    print(fnct_link_name)
    tryCatch(try(
      {resultat<-fnct_logistic_glm(fnct_lien = fnct_link_name,
                                     ech_test=ech_test,ech_train=ech_train)
      LISTE_resultat[[fnct_link_name]]<-resultat
      Mat_confTest<-resultat$resultat_conf_test
      df<-matrix(0,nrow=2,ncol=2)
      if(ncol(Mat_confTest)!=2){
        indice<-as.numeric(colnames(Mat_confTest))+1
        df[,indice]<-Mat_confTest
      }
      else{
        df<-matrix(c(Mat_confTest),nrow = 2,ncol=2)
      }
      Total<-colSums(df)
      Proportion_simulations_predite<-sum(df[,2])/sum(df)
      Proportion_simulations_realite<-sum(df[2,])/sum(df)
      # ligne =vérité. ----------------------------------------------------------
      Sensilite<-df[2,2]/sum(df[2,])
      Specificite<-df[1,1]/sum(df[1,])
      
      # observer parmi les gens prédits positifs ---------------------------------
      Accuracy<-df[2,2]/sum(df[,2])
      df<-rbind(df,Total)
      df<-cbind.data.frame(df,rep(variable,nrow(df)),
                           rep(fnct_link_name,nrow(df)))
      rownames(df)<-c("obs realite","simul realite","Total")
      colnames(df)<-c("obs pred","simul pred","variable","lien")
      F1_score<-2*(Accuracy*Sensilite)/(Accuracy+Sensilite)
      Matrice_perfomances[j,]<-c(fnct_link_name,round(F1_score,2),round(Proportion_simulations_predite,2),
                                 round(Proportion_simulations_realite,2))
      LISTE_resultat[[fnct_link_name]]$resultat_conf_test<-df
      }), error=function(e){
        print(e)
      }
    )
    j<-j+1
  }
  colnames(Matrice_perfomances)<-c("Noyau utilise","F1 score",
                                   "p pred","vrai p")
  LISTE_resultat$performances_globales<- Matrice_perfomances
  return(LISTE_resultat)
}

#' fnct_ML_realite_simulations
#'

#' @param ech_train : list. Training elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param ech_test : list. Testing elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param type_noyau : str. Chosen kernel. 
#' @return
#' @export
#'
#' @examples
fnct_ML_realite_simulations<-function(ech_train,ech_test,type_noyau){
  
  Modif<-FALSE
  simul_train<-ech_train$simul
  simul_test<-ech_test$simul
  real_train<-ech_train$real
  real_test<-ech_test$real
  
  if(is.null(nrow(simul_train))){
    simul_train<-matrix(simul_train,byrow = FALSE,ncol = 1,
                        nrow = length(simul_train))
    simul_test<-matrix(simul_test,byrow = FALSE,ncol = 1,
                       nrow = length(simul_test))
    real_train<-matrix(real_train,byrow = FALSE,ncol = 1,
                       nrow = length(real_train))
    real_test<-matrix(real_test,byrow = FALSE,ncol = 1,
                      nrow = length(real_test))
  }
  realisations<-rbind(real_train,simul_train)
  mean_realisations<-apply(X = realisations,MARGIN = 2, 
                           FUN = mean)
  sd_realisations<-apply(X = realisations,MARGIN = 2, 
                         FUN = sd)
  realisations<-as.data.frame(scale(as.matrix(realisations)))
  # Variable cible : indicatrice simulation ---------------------------------
  
  # Echantillon d'apprentissage ---------------------------------------------
  # y=1 for simulations et y=0 for observations
  realisations$y<-c(rep(0,nrow(real_train)),rep(1,nrow(simul_train)))
  realisations$y<-as.factor(realisations$y)
  realisations<-realisations[sample(nrow(realisations)),]
  
  modele_svm<-e1071::svm(y~.,data=realisations,kernel=type_noyau, 
                         probability=TRUE)
  
  predictions<-modele_svm$fitted
  DF_test<-rbind(real_test,simul_test)
  DF_test<-as.data.frame(scale(DF_test,center=mean_realisations, 
                 scale=sd_realisations))
  rownames(DF_test)<-c(1:nrow(DF_test))
  y<-c(rep(0,nrow(real_test)),rep(1,nrow(simul_test)))
  # predictions--SVM ---------------------------------------------------------
  predictions_modele<-predict(modele_svm,DF_test,
                              decision.values = TRUE)
  y_pred<-as.numeric(as.character(predictions_modele))
  scores<-attr(predictions_modele,"decision.values")
  probs_fonction_chosen<-as.numeric((1+exp((-1)*scores))^(-1))
  
  # Take into account the order ---------------------------------------------
  # sign(score) gives yi.  --------------------------------------------------
  if(y_pred[1]!=round(probs_fonction_chosen[1])){
    probs_fonction_chosen<-(1-probs_fonction_chosen)
  }
  return(list("Modèle"=modele_svm,"resultat_conf_train"=table(realisations$y,predictions),
              "resultat_conf_test"=table(y,predictions_modele),
              "type_noyau"=type_noyau, 
              "predictions_modele"=predictions_modele,
              "probs_y"=probs_fonction_chosen, 
              "y"=y_pred, 
              "y_test"=y))
}


# RandomForest ------------------------------------------------------------

#' fnct_RF
#'
#' @param ech_train : list. Training elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param ech_test : list. Testing elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param NB_trees : int. Number of trees used in the forest.
#'
#' @return
#' @export
#'
#' @examples
fnct_RF<-function(ech_train,ech_test,NB_trees){
  
  Modif<-FALSE
  simul_train<-ech_train$simul
  simul_test<-ech_test$simul
  real_train<-ech_train$real
  real_test<-ech_test$real
  
  if(is.null(nrow(simul_train))){
    simul_train<-matrix(simul_train,byrow = FALSE,ncol = 1,
                        nrow = length(simul_train))
    simul_test<-matrix(simul_test,byrow = FALSE,ncol = 1,
                       nrow = length(simul_test))
    real_train<-matrix(real_train,byrow = FALSE,ncol = 1,
                       nrow = length(real_train))
    real_test<-matrix(real_test,byrow = FALSE,ncol = 1,
                      nrow = length(real_test))
  }
  realisations<-rbind(real_train,simul_train)
  mean_realisations<-apply(X = realisations,MARGIN = 2, 
                           FUN = mean)
  sd_realisations<-apply(X = realisations,MARGIN = 2, 
                         FUN = sd)
  realisations<-as.data.frame(scale(as.matrix(realisations)))
  # Variable cible : indicatrice simulation ---------------------------------
  
  # Echantillon d'apprentissage ---------------------------------------------
  # y=1 for simulations et y=0 for observations
  L<-ncol(realisations)
  realisations$y<-c(rep("0",nrow(real_train)),rep("1",nrow(simul_train)))
  realisations$y<-as.factor(realisations$y)
  realisations<-realisations[sample(nrow(realisations)),]
  Input_rf<-realisations[,1:L]
  if(is.null(nrow(Input_rf))){
    Input_rf<-matrix(Input_rf,byrow = FALSE,ncol = 1,
                     nrow = length(Input_rf))
  }
  modele_rf<-randomForest::randomForest(Input_rf,y=realisations$y,
                                        data=realisations,
                                        ntree=NB_trees)
  predictions<-modele_rf$predicted
  DF_test<-rbind(real_test,simul_test)
  if(ncol(DF_test)==1){
    DF_test<-scale(DF_test,center=mean_realisations, 
          scale=sd_realisations)
  }
  else{
    DF_test<-as.data.frame(scale(DF_test,center=mean_realisations, 
                                 scale=sd_realisations))
  }
  y<-c(rep(0,nrow(real_test)),rep(1,nrow(simul_test)))
  # predictions--RF ---------------------------------------------------------
  probs_fonction_chosen<-predict(modele_rf,DF_test,
                              type = "prob")[,2]
  y_pred<-round(probs_fonction_chosen)
  return(list("Modèle"=modele_rf,"resultat_conf_train"=table(realisations$y,predictions),
              "resultat_conf_test"=table(y,y_pred),
              "probs_y"=probs_fonction_chosen, 
              "y"=y_pred, 
              "y_test"=y))
}

#' fnct_logistic_glm
#'
#' @param ech_train : list. Training elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param ech_test : list. Testing elements with simulated time series (or values)
#'                    and recorded time series (or values).
#' @param fnct_lien : str. Link function used. 
#'
#' @return List. Fitted model with main results. 
#' @export
#'
#' @examples
fnct_logistic_glm<-function(ech_train,ech_test,fnct_lien){
  
  Modif<-FALSE
  simul_train<-ech_train$simul
  simul_test<-ech_test$simul
  real_train<-ech_train$real
  real_test<-ech_test$real
  
  if(is.null(nrow(simul_train))){
    simul_train<-matrix(simul_train,byrow = FALSE,ncol = 1,
                        nrow = length(simul_train))
    simul_test<-matrix(simul_test,byrow = FALSE,ncol = 1,
                       nrow = length(simul_test))
    real_train<-matrix(real_train,byrow = FALSE,ncol = 1,
                       nrow = length(real_train))
    real_test<-matrix(real_test,byrow = FALSE,ncol = 1,
                      nrow = length(real_test))
  }
  realisations<-rbind(real_train,simul_train)
  mean_realisations<-apply(X = realisations,MARGIN = 2, 
                           FUN = mean)
  sd_realisations<-apply(X = realisations,MARGIN = 2, 
                         FUN = sd)
  realisations<-as.data.frame(scale(as.matrix(realisations)))
  # Variable cible : indicatrice simulation ---------------------------------
  
  # Echantillon d'apprentissage ---------------------------------------------
  
  # y=1 pour les simulations et y=0 pour les observations
  realisations$y<-c(rep(0,nrow(real_train)),rep(1,nrow(simul_train)))
  realisations$y<-as.factor(realisations$y)
  
  # Shuffle -----------------------------------------------------------------
  realisations<-realisations[sample(nrow(realisations)),]
  modele_logit<-glm(y~.,data=realisations,family = binomial(link = fnct_lien))
  predictions<-round(modele_logit$fitted)
  
  # Echantillon_test --------------------------------------------------------
  y<-c(rep(0,nrow(real_test)),rep(1,nrow(simul_test)))
  vect<-c(real_test,simul_test)
  DF_test<-rbind(real_test,simul_test)
  DF_test<-as.data.frame(scale(DF_test,center=mean_realisations, 
                                 scale=sd_realisations))
  rownames(DF_test)<-c(1:nrow(DF_test))
  # predictions--logit ---------------------------------------------------------
  probs_model<-as.numeric(predict(modele_logit,DF_test,
                                    type="response"))
  predictions_modele<-round(probs_model)
  return(list("Modèle"=modele_logit,"resultat_conf_train"=table(realisations$y,predictions),
              "resultat_conf_test"=table(y,predictions_modele),
              "type_fnct_lien"=fnct_lien,"y_test"=y,"probs_y"=probs_model))
}

#' Function_extract_elt
#'
#' @param modele_stat : model. Return object of craft GLM, SVM and RF functions.
#'
#' @return Key performance values.
#' @export
#'
#' @examples
Function_extract_elt<-function(modele_stat){
  
  # ROC ---------------------------------------------------------------------
  reality_<-modele_stat$y_test
  probs_algorithm_<-modele_stat$probs_y
  curve_roc_<-pROC::roc(reality_,
                         probs_algorithm_,plot=FALSE, 
                         quiet=TRUE)
  index_mid_<-which.min(abs(curve_roc_$thresholds-0.5))
  Liste_ROC_<-list("curve"=curve_roc_, 
                   "index_mid"=index_mid_)
  Mat_confTest_<-modele_stat$resultat_conf_test
  df_<-matrix(0,nrow=2,ncol=2)
  if(ncol(Mat_confTest_)!=2){
    indice<-as.numeric(colnames(Mat_confTest_))+1
    df_[,indice]<-Mat_confTest_
  }
  else{
    df_<-matrix(c(Mat_confTest_),nrow = 2,ncol=2)
  }
  Num<-2*(det(df_))
  Sum_lignes<-prod(apply(df_,MARGIN = 1,FUN = function(x){return(sum(x))}))
  
  Sum_cols<-prod(apply(df_,MARGIN = 2,FUN = function(x){return(sum(x))}))
  Denom<-Sum_lignes+Sum_cols

  # HSS near 0= no skill ----------------------------------------------------
  
  HSS_score<-Num/Denom
  Prop_<-t(apply(df_,MARGIN = 1,
                 function(x){return(x/sum(x))}))
  R_<-sum(diag(df_))/sum(df_)
 
  
  return(list("1"=Prop_,"2"=R_,"3"=Liste_ROC_,"4"=HSS_score))
}

# Uncertainty due to sampling of simulated +sampling for train ------------

#' fnct_calcul_correct_proportion
#'
#' @param Base_simul : dataframe of simulated time series (or L2 norm).
#' @param Base_data : dataframe of observed time series (or L2 norm).
#' @param hyp_param : vector(str). Link function and kernel used. 
#' @param K : int. K fold parameter. 
#'
#' @return
#' @export
#'
#' @examples
fnct_calcul_correct_proportion<-function(Base_simul, Base_data,hyp_param,
                                          type_sampling,K){
  if(!is.null(nrow(Base_data))){
    Nb_ech<-nrow(Base_data)
    Nb_simul<-nrow(Base_simul)
  }
  else{
    Nb_ech<-length(Base_data)
    Nb_simul<-length(Base_simul)
  }
  Base_ind_s<-c(1:Nb_simul)
  Base_ind_r<-c(1:Nb_ech)
  
  # sampling of real observations ------------------------------------------------------
  Indices_train_r<-sample(Base_ind_r,size = round(Nb_ech*0.7),replace = FALSE)
  Indices_test_r<-which(!(Base_ind_r%in% Indices_train_r))
  
  # Sampling of simulations, size  ------------------------------------------
  # because size simul >> size data. 
  Indices_s<-sample(Base_ind_s,size =Nb_ech,replace = FALSE)
  if(!is.null(nrow(Base_simul))){
    Sub_Simul<-Base_simul[Indices_s,]
  }
  else{
    Sub_Simul<-Base_simul[Indices_s]
  }
  if(type_sampling=="simple"){
    Inds_s_test<-sample(x = c(1:Nb_ech),size = round(Nb_ech*0.3))
    Inds_r_test<-sample(x = c(1:Nb_ech),size = round(Nb_ech*0.3))
    # Construction of the df --------------------------------------------------
    if(!is.null(nrow(Base_data))){
      ech_train<-list("simul"=Sub_Simul[-Inds_s_test,],
                      "real"=Base_data[-Inds_r_test,])
      ech_test<-list("simul"=Sub_Simul[Inds_s_test,],
                     "real"= Base_data[Inds_r_test,])
    }
    else{
      ech_train<-list("simul"=Sub_Simul[-Inds_s_test],
                      "real"=Base_data[-Inds_r_test])
      ech_test<-list("simul"=Sub_Simul[Inds_s_test],
                     "real"= Base_data[Inds_r_test])
    }
    
    modele_SVM<-fnct_ML_realite_simulations(ech_train = ech_train,ech_test=ech_test,
                                            type_noyau = hyp_param[1])
    SVM_results<-Function_extract_elt(modele_stat = modele_SVM)
    Prop_SVM<-SVM_results[[1]]
    R_SVM<-SVM_results[[2]]
    Liste_ROC_SVM<-SVM_results[[3]]
    HSS_SVM<-SVM_results[[4]]
    
    modele_GLM<-fnct_logistic_glm(fnct_lien = hyp_param[2],
                                    ech_test=ech_test,ech_train=ech_train)
    GLM_results<-Function_extract_elt(modele_stat = modele_GLM)
    Prop_GLM<-GLM_results[[1]]
    R_GLM<-GLM_results[[2]]
    Liste_ROC_GLM<-GLM_results[[3]]
    HSS_<-GLM_results[[4]]

    
    # RF ----------------------------------------------------------------------
    modele_RF<-fnct_RF(ech_train = ech_train,ech_test=ech_test,
                       NB_trees = as.numeric(hyp_param[3]))
    RF_results<-Function_extract_elt(modele_stat = modele_RF)
    Prop_RF<-RF_results[[1]]
    R_RF<-RF_results[[2]]
    Liste_ROC_RF<-RF_results[[3]]
    HSS_rf<-RF_results[[4]]
    
  }
  
  return(list("GLM"=list("accuracy"=R_GLM,
                         "modeles"=modele_GLM,
                         "ROC"=Liste_ROC_GLM,
                         "HSS"=HSS_),
              "SVM"=list("accuracy"=R_SVM,
                         "modeles"=modele_SVM,
                         "ROC"=Liste_ROC_SVM,
                         "HSS"=HSS_SVM),
              "RF"=list("accuracy"=R_RF,
                        "modeles"=modele_RF,
                        "ROC"=Liste_ROC_RF,
                        "HSS"=HSS_rf)))
}

#' Running_perfs_ML
#'
#' @param Base_simul : dataframe (or vector). Simulated time series (or values).
#' @param Base_data : dataframe (or vector). Observed time series (or values).
#' @param hyp_param : list. Hyper parameters such as the link function (GLM), the kernel (SVM)
#'                    and the number of trees (RF).
#' @param NB_times : int. Number of repetitions of the experiment (=number of models used).
#' @param alpha_prop : float. Confidence level used. 
#' @param K : int. Not used. Number of blocks used in cross validation.
#' @param type_sampling : str. Type of experiment used (simple or CV).
#' @param title_ROC : str. Titre used in ROC graphics. 
#' @param NB_shown_ROC : str. Number of ROC curves shown.
#'
#' @return
#' @export
#'
#' @examples
Running_perfs_ML<-function(Base_simul, Base_data,hyp_param,NB_times,alpha_prop,K,
                          type_sampling, title_ROC,NB_shown_ROC){
  
  
  Resultat<-replicate(n = NB_times,fnct_calcul_correct_proportion(Base_simul=Base_simul,
                                                                   Base_data = Base_data,
                                                                   hyp_param =hyp_param,type_sampling=type_sampling,
                                                                   K=K))
  Echs<-sample(c(1:NB_times),size = NB_shown_ROC)
  Whole_curve<-apply(Resultat,MARGIN = 2,
                             function(x){return(x$SVM$ROC$curve)})                                                                
  
  GGROC_several<-pROC::ggroc(Whole_curve[Echs],alpha=0.4)+
    guides(col="none")+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black")+
    ggtitle(paste0("SVM (", hyp_param[1],") ROC curves for ",title_ROC," (shown=",NB_shown_ROC,", calculated=",NB_times,")"))
  print(GGROC_several)
  Whole_curve_GLM<-apply(Resultat,MARGIN = 2,
                                 function(x){return(x$GLM$ROC$curve)})
  GGROC_several_2<-pROC::ggroc(Whole_curve_GLM[Echs],alpha=0.4)+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black")+
    guides(col="none")+
    ggtitle(paste0("GLM (", hyp_param[2],") ROC curves for ",title_ROC," (shown=",NB_shown_ROC,", calculated=",NB_times,")"))
  print(GGROC_several_2)
  
  Whole_curve_RF<-apply(Resultat,MARGIN = 2,function(x){return(x$RF$ROC$curve)})
  GGROC_several_RF<-pROC::ggroc(Whole_curve_RF[Echs],alpha=0.4)+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black")+
    guides(col="none")+
    ggtitle(paste0("Random Forest ROC curves for ",title_ROC," (shown=",NB_shown_ROC,", calculated=",NB_times,")"))
  print(GGROC_several_RF)

  # Accuracy ----------------------------------------------------------------

  GLM_ACCURACY<-apply(Resultat,MARGIN = 2,function(x){return(x$GLM$accuracy)})
  SVM_ACCURACY<-apply(Resultat,MARGIN = 2,function(x){return(x$SVM$accuracy)})
  RF_ACCURACY<-apply(Resultat,MARGIN = 2,function(x){return(x$RF$accuracy)})
  Accuracy<-cbind(SVM_ACCURACY,GLM_ACCURACY,RF_ACCURACY)
  Q1<-apply(X = Accuracy,MARGIN = 2,FUN = function(x){return(round(quantile(x,(alpha_prop/2)),2)*100)})
  Q2<-apply(X = Accuracy,MARGIN = 2,FUN = function(x){return(round(quantile(x,1-(alpha_prop/2)),2)*100)})
  df<-cbind.data.frame(Q1,Q2)
  resultat<-cbind.data.frame(c("SVM","GLM","RF"),paste0(df$Q1,"-",df$Q2))
  colnames(resultat)<-c("Modele","qu_accuracy")
  
  # HSS ----------------------------------------------------------------
  GLM_HSS<-apply(Resultat,MARGIN = 2,function(x){return(x$GLM$HSS)})
  SVM_HSS<-apply(Resultat,MARGIN = 2,function(x){return(x$SVM$HSS)})
  RF_HSS<-apply(Resultat,MARGIN = 2,function(x){return(x$RF$HSS)})
  HSS<-cbind(SVM_HSS,GLM_HSS,RF_HSS)
  Q1_h<-apply(X = HSS,MARGIN = 2,FUN = function(x){return(round(quantile(x,(alpha_prop/2)),2)*100)})
  Q2_h<-apply(X = HSS,MARGIN = 2,FUN = function(x){return(round(quantile(x,1-(alpha_prop/2)),2)*100)})
  df<-cbind.data.frame(Q1_h,Q2_h)
  resultat_h<-cbind.data.frame(c("SVM","GLM","RF"),paste0("(",df$Q1,")-(",df$Q2,")"))
  colnames(resultat_h)<-c("Modele","qu_hss")
  
  return(list("Accuracy"=resultat,
              "HSS"=resultat_h,
              "Ensemble_r"=Resultat))
}


Plot_ROC_curves_mult_links<-function(obj_link_GLM,titre_graphique,liste_cles,racine_repertoire,variable){
  liste_ROC<-list()
  liste_th_mid<-list()
  Mat_AUC<-matrix(NA,ncol=2,nrow=length(liste_cles))
  j<-1
  for(link in liste_cles){
    reality<-obj_link_GLM[[link]]$y_test
    probs_algorithm<-obj_link_GLM[[link]]$probs_y
    print(mean(probs_algorithm))
    curve_roc<-pROC::roc(reality,
                          probs_algorithm,plot=FALSE, 
                          quiet=TRUE)
    liste_ROC[[link]]<-curve_roc
    liste_th_mid[[link]]<-which.min(abs(curve_roc$thresholds-0.5))
    auc<-round(as.numeric(pROC::auc(curve_roc)),2)
    Mat_AUC[j,]<-c(link,auc)
    j<-j+1
  }
  GG<-pROC::ggroc(liste_ROC)+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black")+
    labs(color="Fonction")+
    ggtitle(titre_graphique)
  print(GG)
  for(link_fonction in names(liste_th_mid)){
    c_roc<-liste_ROC[[link_fonction]]
    indice_th_mid<-liste_th_mid[[link_fonction]]
    x_mid<-c_roc$specificities[indice_th_mid]
    y_mid<-c_roc$sensitivities[indice_th_mid]
    GG<-GG+
      annotate("point",x=x_mid,y=y_mid,pch=4,cex=3,col="red")
  }
  Mat_AUC<-as.data.frame(Mat_AUC)
  colnames(Mat_AUC)<-c("lien","AUC")
  Mat_AUC$variable<-rep(variable,nrow(Mat_AUC))
  # Exporter l'AUC --------------------------------------------------------------
  return(list("AUC"=Mat_AUC,"grap"=GG))
  
  # Fin code  ---------------------------------------------------------------
  
}
