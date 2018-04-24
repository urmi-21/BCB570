setwd("D:/courses/sem4/bcb570/hw/hw7")
source('general.R')
library(kernlab)
library(ROCR)
library(lattice)
library(readr)
rm(list=ls())
# Load the Protein Prediction dataset as text file
fdata <- read.csv('C_elegans_acc_gc.csv', header =FALSE, col.names = c("class","GCIntron(after)", "GC Exon"))
fdata_balanced<-fdata
#count num classes
num_pos<-length(which(fdata$class==1))
num_neg<-length(which(fdata$class==-1))
#over sample the positive classes
d<-0
while(d<num_neg-num_pos){
for(i in (which((fdata$class==1)))){
  if(d>=num_neg-num_pos){
    break
  }
  rand=runif(1, min = 0, max = 1)
  if(rand>0.4){
    fdata_balanced[length(fdata_balanced$class)+1,]<-fdata[i,]
    d=d+1
  }
}
}

##now fdata_balanced is our balanced class. In this example we are not doing cv so we oversampled before using a cv dataset
# Predict type of disease from expression data
y <- fdata_balanced$class
x <- fdata_balanced[c(2:3)]
x <- as.matrix(sapply(x, as.numeric))
n<-length(y)

tindex <- sample(n,round(n*0.7)) # indices of training samples to plot
xtrain <- x[tindex,]
xtest <- x[-tindex,]
ytrain <- y[tindex]
ytest <- y[-tindex]
istrain=rep(0,n)
istrain[tindex]=1

# Train a linear SVM ------------------------------------------------------

# train the SVM
svp <- ksvm(xtrain,ytrain,type="C-svc",kernel="rbf",C=10,scaled=c())

# General summary
svp

# Use the built-in function to pretty-plot the classifier
plot(svp,data=xtrain)

# Visualize decision boundary with home-made function (available in general.R)
#pdf('linearsvm.pdf',width=8)
plotlinearsvm2D(svp,xtrain)
#dev.off()

# Predict labels on test
ypred = predict(svp,xtest)
# Confusion table
cm<-table(ytest,ypred)
acc<-sum(diag(cm))/length(ytest)
cat("Accuraccy:",acc*100,"%")
pred <- prediction(ypred,ytest)


# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf,title="TPR vs. FPR")
auc<- performance(pred, measure = "auc")
cat("AUC ROC:",unlist(auc@y.values)[1])
# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf,title="Precision-Recall")





##Read data
X1625Data_encoded <- read_csv("1625Data_encoded.txt")
#firsc col is class
y <- X1625Data_encoded$`data11$V2`
x <- X1625Data_encoded[c(2:length(X1625Data_encoded))]
x <- as.matrix(sapply(x, as.numeric))
n<-length(y)
###use al data as train set
tindex <- sample(n,round(n*1)) # indices of training samples to plot
xtrain <- x[tindex,]
xtest <- x[-tindex,]
ytrain <- y[tindex]
ytest <- y[-tindex]
istrain=rep(0,n)
istrain[tindex]=1
#initialize c list
clist <- 2^seq(-1,10)
clist <- seq(1,1000,2)
errrbf <- numeric(length(clist))
errlin <- numeric(length(clist))

for (i in seq(length(clist))) {
  cat('.')
  svp <- ksvm(xtrain,ytrain,type="C-svc",kernel='vanilladot',C=clist[i],cross=5)
  errlin[i] <- cro
  ss(svp)
}

plot(clist,errlin,type='l',ylim=c(0,0.1),xlab="C",ylab="Error rate",col=1,lwd=2,main='CrossValidaton error')
grid()
legend("topright",c('Linear'),lwd=2,col=c(1))

bestC<--2
bestC<-clist[which.min(errlin)]
cat("Best C:",bestC,"from linear, Err:", min(errlin))
svp<- ksvm(xtrain,ytrain,type="C-svc",kernel='vanilladot',C=bestC)
#open test set
X746Data_encoded <- read_csv("746Data_encoded.txt")

testy<-X746Data_encoded$`data11$V2`
testx<-X746Data_encoded[c(2:length(X746Data_encoded))]
test_pred<- predict(svp,testx)
cm<-table(testy,test_pred)
acc<-sum(diag(cm))/length(testy)
cat("Accuraccy on 764 data:",acc*100,"%")

pred <- prediction(test_pred,testy)
# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf,title="TPR vs. FPR")
auc<- performance(pred, measure = "auc")
cat("AUC ROC:",unlist(auc@y.values)[1])
# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf,title="Precision-Recall")



##############solution 3

#read data
schillingData_encoded <- read_csv("schillingData_encoded.txt")
#balance class by under sampling
num_pos<-length(which(schillingData_encoded$`data11$V2`==1))
num_neg<-length(which(schillingData_encoded$`data11$V2`==-1))
schillingData_encoded_bal<-data.frame(matrix(ncol = 161, nrow = 0))
colnames(schillingData_encoded_bal)<-colnames(schillingData_encoded)
flag=-1
i=1
while(flag<0){
  if(length(which(schillingData_encoded_bal$`data11$V2`==-1)) == num_pos && length(which(schillingData_encoded_bal$`data11$V2`==1)) == num_pos){
    break
  }
  if(as.numeric(schillingData_encoded[i,]$`data11$V2`)==-1  ){
  if((length(which(schillingData_encoded_bal$`data11$V2`==-1)) < num_pos)){
    schillingData_encoded_bal[i,]<-schillingData_encoded[i,]  
    
  }
  }
  
  if(as.numeric(schillingData_encoded[i,]$`data11$V2`)==1  ){
    if((length(which(schillingData_encoded_bal$`data11$V2`==1)) < num_pos)){
      schillingData_encoded_bal[i,]<-as.data.frame(schillingData_encoded[i,]) 
      print(schillingData_encoded[i,]$`data11$V2` )
    }
  }
  
  i=i+1
}
#remove NAs
schillingData_encoded_bal<-na.omit(schillingData_encoded_bal)



#firsc col is class
y <- schillingData_encoded_bal$`data11$V2`
x <- schillingData_encoded_bal[c(2:length(schillingData_encoded_bal))]
x <- as.matrix(sapply(x, as.numeric))
n<-length(y)
###use al data as train set
tindex <- sample(n,round(n*1)) # indices of training samples to plot
xtrain <- x[tindex,]
xtest <- x[-tindex,]
ytrain <- y[tindex]
ytest <- y[-tindex]
istrain=rep(0,n)
istrain[tindex]=1
#initialize c list
clist <- seq(1,1000,2)
errrbf <- numeric(length(clist))
for (i in seq(length(clist))) {
  cat('.')
  svp <- ksvm(xtrain,ytrain,type="C-svc",kernel='rbf',C=clist[i],cross=5)
  errrbf[i] <- cross(svp)
}

plot(clist,errrbf,type='l',ylim=c(0,0.1),xlab="C",ylab="Error rate",col=1,lwd=2,main='CrossValidaton error')
grid()
legend("topright",c('RBF'),lwd=2,col=c(1))
bestC<--2
bestC<-clist[which.min(errrbf)]
cat("Best C:",bestC,"from RBF using schillingData, Err:", min(errrbf))
svp<- ksvm(xtrain,ytrain,type="C-svc",kernel='rbf',C=bestC)

#read test set
X1625Data_encoded <- read_csv("1625Data_encoded.txt")
testy<-X1625Data_encoded$`data11$V2`
testx<-X1625Data_encoded[c(2:length(X1625Data_encoded))]
test_pred<- predict(svp,testx)
cm<-table(testy,test_pred)
acc<-sum(diag(cm))/length(testy)
cat("Accuraccy on 764 data:",acc*100,"%")
pred <- prediction(test_pred,testy)
# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf,title="TPR vs. FPR")
auc<- performance(pred, measure = "auc")
cat("AUC ROC:",unlist(auc@y.values)[1])
# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf,title="Precision-Recall")



#read data
schillingData_encoded <- read_csv("schillingData_encoded.txt")
#balance class by under sampling
num_pos<-length(which(schillingData_encoded$`data11$V2`==1))
num_neg<-length(which(schillingData_encoded$`data11$V2`==-1))
#balance class by oversampling
schillingData_neg<-schillingData_encoded[which(schillingData_encoded$`data11$V2`==-1),]
schillingData_pos<-schillingData_encoded[which(schillingData_encoded$`data11$V2`==1),]
N<-5
#initialize c list
clist <- 2^seq(-1,10)
errrbf<-c()
bestMod<-c()
bestAcc<--1
for(c in clist){
cvErrs<-c()
for(i in(1:N)){
#create a balanced cross val dataset
cv<-data.frame(matrix(ncol = 161, nrow = 0))
colnames(cv)<-colnames(schillingData_encoded)
#sample 10% of data from pos class and equal amount from neg class to make cv set
n_cv<-as.integer(num_pos*0.1)
## set the seed to make your partition reproductible
set.seed(41*i+i**2)
cv_ind <- sample(seq_len(nrow(schillingData_pos)), size = n_cv)
cv_pos <- schillingData_pos[cv_ind, ]
train_pos <- schillingData_pos[-cv_ind, ]
cv_ind <- sample(seq_len(nrow(schillingData_neg)), size = n_cv)
cv_neg <- schillingData_neg[cv_ind, ]
train_neg <- schillingData_neg[-cv_ind, ]

#final CV
cv<-rbind(cv_pos,cv_neg)

#over sample the remaining pos class by replicating randomly
osize=length(train_neg$`data11$V2`)-length(train_pos$`data11$V2`)
oversamp_ind<-sample(seq_len(nrow(schillingData_pos)), size = osize,replace = T)
tr_pos_os<-train_pos[oversamp_ind,]
train_pos<-rbind(train_pos,tr_pos_os)
train_set<-rbind(train_pos,train_neg)

#####start training#######
y <- train_set$`data11$V2`
x <- train_set[c(2:length(train_set))]
x <- as.matrix(sapply(x, as.numeric))
n<-length(y)
#train SVM model
svp <- ksvm(x,y,type="C-svc",kernel="rbf",C=c,scaled=c())
# General summary
svp
xtest<-cv[c(2:length(train_set))]
ytest<-cv$`data11$V2`
ypred = predict(svp,xtest)
# Confusion table
cm<-table(ytest,ypred)
acc<-sum(diag(cm))/length(ytest)
#cat("Accuraccy:",acc*100,"%")
thisErr<-1-acc
#print(thisErr)
cvErrs<-c(cvErrs,thisErr)
}
errrbf<-c(errrbf, mean(cvErrs))
if(bestAcc<mean(cvErrs)){
  bestAcc<-mean(cvErrs)
  bestMod<-svp
}

}

plot(clist,errrbf,type='l',ylim=c(0,0.2),xlab="C",ylab="Error rate",col=1,lwd=2,main='CrossValidaton error')
grid()
legend("topright",c('RBF'),lwd=2,col=c(1))
bestC<--2
bestC<-clist[which.min(errrbf)]
cat("Best C:",bestC,"from RBF using balanced schillingData, Err:", min(errrbf))
bestMod

##predict on test data 
#read test set
X1625Data_encoded <- read_csv("1625Data_encoded.txt")
testy<-X1625Data_encoded$`data11$V2`
testx<-X1625Data_encoded[c(2:length(X1625Data_encoded))]
test_pred<- predict(bestMod,testx)
cm<-table(testy,test_pred)
acc<-sum(diag(cm))/length(testy)
cat("Accuraccy on 764 data:",acc*100,"%")
pred <- prediction(test_pred,testy)
# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf,title="TPR vs. FPR")
auc<- performance(pred, measure = "auc")
cat("AUC ROC:",unlist(auc@y.values)[1])
# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf,title="Precision-Recall")







