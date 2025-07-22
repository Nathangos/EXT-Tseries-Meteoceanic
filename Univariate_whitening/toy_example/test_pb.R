rm(list=ls())
set.seed(133)
########## Understand the origin of the logistic behavior.
n<-250
# data
X1<-rnorm(n)
d<-10
Mat_data<-matrix(NA,nrow = n,ncol = d)
Mat_data[,1]<-X1
Mat_simul<-matrix(NA,nrow = n,ncol = d)

Mat_simul[,1]<-X1
Xnew<-X1
for(j in c(2:d)){
  X2<-1*Xnew
  Mat_data[,j]<-X2
  Xnew<-X2
}
qr(Mat_data)$rank
#simul
epsi<-rnorm(n = n)
for(j in c(2:7)){
  X2<-1*X1+epsi
  Mat_simul[,j]<-X2
  X1<-X2
}
Mat_simul[,8]<-rnorm(n)
Mat_simul[,9]<-rnorm(n)
Mat_simul[,10]<-rnorm(n)
plot(Mat_simul[,1],Mat_simul[,2])
Mat_simul<-as.data.frame(Mat_simul)
Mat_data<-as.data.frame(Mat_data)
# linearly dependent. 
Mat_data$y<-rep(0,nrow(Mat_data))
Mat_simul$y<-rep(1,nrow(Mat_simul))
colnames(Mat_simul)<-colnames(Mat_data)
Mat_xy<-rbind.data.frame(Mat_data,
                 Mat_simul)
qr(Mat_xy[,c(1:d)])$rank
qr(Mat_simul[,c(1:d)])$rank
qr(Mat_data[,c(1:d)])$rank
Model<-glm(formula = y~.,data = Mat_xy,
    family =binomial())
summary(Model)
table(round(Model$fitted.values),Mat_xy$y)
