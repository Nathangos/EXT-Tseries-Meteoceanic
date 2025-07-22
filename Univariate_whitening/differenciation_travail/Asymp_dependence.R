rm(list=ls())
nom_variable<-"Surcote"
repertoire<-"ss_tend/"
type_donnees<-"HIVER"
if(type_donnees=="HIVER"){
  repertoire<-paste0(repertoire,"HIVER/")
}
if(type_donnees=="HORS_HIVER"){
  repertoire<-paste0(repertoire,"HORS_HIVER/")
}
nom_recherche<-paste0(repertoire,nom_variable,"_ss_tend.csv")
nom_recherche2<-paste0("residus_clust/",nom_variable,"_residus.csv")
obs<-read.csv(nom_recherche)[,2:38]
Resid<-read.csv(nom_recherche2)[,2:38]
colnames(obs)<-c(1:ncol(obs))
colnames(Resid)<-colnames(obs)
t<-1
s<-37
nomG<-ifelse(nom_variable=="Surcote","S",nom_variable)
par(mfrow = c(1, 2),  # 1 row, 2 plots
    cex.lab = 1,    # axis label size
    cex.axis = 1,   # axis tick label size
    cex.main = 1,   # main title size
    cex = 1)
Two<-expression(bar(chi))
POT::chimeas(data = cbind.data.frame(obs[,t],obs[,s]),
             ask = FALSE,which=1,cex.lab=1.5)
POT::chimeas(data = cbind.data.frame(obs[,t],obs[,s]),
             ask = FALSE,
             which=2,
             ylabs = rep(NA,2),
             cex.lab=1.5)
mtext(expression(bar(chi)),cex=1.5,las=2,
      outer=TRUE,line=-10.5)
Time_s<-(t-19)*(1/6)
# Basenot<-expression("Asymptotic dependence between "~tilde(X)[M]^{t}~"and "~tilde(X)[M]^{s}~" for "~nom_variable)
#obj<-do.call("substitute", list(Basenot[[1]], 
# list(t = t,s=s,
#      nom_variable=nomG)))
# mtext(obj,outer = TRUE,line=-4)

# Dependences -------------------------------------------------------------
##############
Time_reference<-1
M<-20
h<-seq.int(1,to = M,by=1)
Corr_coeffs_resid<-rep(NA,M)
Corr_coeffs_data<-rep(NA,M)
Function_lagcorr<-function(series,lag){
  fin<-length(series)-lag
  x_<-series[c(1:fin)]
  debut<-lag+1
  x_lag1<-series[c(debut:length(series))]
  return(
    list("Pearson"=as.numeric(cor.test(x_,x_lag1,method="pearson")$estimate),
         "Kendall"=as.numeric(cor.test(x_,x_lag1,method="kendall")$estimate),
         "Spearman"=as.numeric(cor.test(x_,x_lag1,method="spearman")$estimate)))
}
Corr_coeffs_resid<-rbind(sapply(h,FUN = Function_lagcorr,
       series=Resid[,Time_reference]))
type<-c(rep("Pearson",M),rep("Kendall",M),
        rep("Spearman",M))
Corr_coeffs_data<-rbind(sapply(h,FUN = Function_lagcorr,
                          series=obs[,Time_reference])) 
df<-as.data.frame(t(Corr_coeffs_data))
DATA_CORR<-melt(t(sapply(df,unlist)))
df_resid<-as.data.frame(t(Corr_coeffs_resid))
RESID_CORR<-melt(t(sapply(df_resid,unlist)))

RESID_CORR$Legend<-rep("simulations",nrow(DATA_CORR))
DATA_CORR$Legend<-rep("data",nrow(DATA_CORR))

colnames(DATA_CORR)<-c("type","Lag","coefficient","Legend")
colnames(RESID_CORR)<-colnames(DATA_CORR)

WHOLE_corr<-rbind.data.frame(DATA_CORR,RESID_CORR)
WHOLE_corr$coefficient<-sapply(X = WHOLE_corr$coefficient,unlist)
WHOLE_corr
ggplot(data=WHOLE_corr,aes(x=Lag,y=coefficient,col=Legend,
                           shape=Legend))+
  facet_wrap(~type)+
  geom_point()+
  geom_hline(yintercept=0)+
  scale_shape_manual(values=c("data"=19,
                              "simulations"=17))+
  scale_color_manual(values=c("simulations"="orange",
                              "data"="blue"))+
  theme(axis.title=element_text(size=15),
        legend.text=element_text(size=10))

