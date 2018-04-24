## Application on  ALL gene expression data ##
##------------------------------------------##
source('general.R')
library(kernlab)
library(ROCR)
library(lattice)
# Load the Protein Prediction dataset as text file
fdata <- read.csv('C_elegans_acc_gc.csv', header =FALSE, col.names = c("class","GCIntron(after)", "GC Exon"))

# Inspect them
show(fdata)
summary(fdata)


# Predict type of disease from expression data
y <- fdata$class
x <- fdata[c(2:3)]
x <- as.matrix(sapply(x, as.numeric))
n<-length(y)
# to keep a balanced dataset use first four hundred samples
n<-400
xbal <- x[1:n,]
ybal <- y[1:n]

# Visualize first 400 samples
tindex <- sample(n,round(n*0.8)) # indices of training samples to plot
xtrain <- xbal[tindex,]
xtest <- xbal[-tindex,]
ytrain <- ybal[tindex]
ytest <- ybal[-tindex]
istrain=rep(0,n)
istrain[tindex]=1
plot(xtrain,col=ifelse(ytrain>0,1,2),pch=ifelse(istrain==1,1,2))
legend("topleft",c('Positive Train','Positive Test','Negative Train','Negative Test'),col=c(1,1,2,2),pch=c(1,2,1,2),text.col=c(1,1,2,2))

# Train a linear SVM ------------------------------------------------------

# train the SVM
svp <- ksvm(xtrain,ytrain,type="C-svc",kernel='vanilladot',C=1,scaled=c())


# General summary
svp

# Attributes that you can access
attributes(svp)

# For example, the support vectors
alpha(svp)
alphaindex(svp)
b(svp)

# Use the built-in function to pretty-plot the classifier
plot(svp,data=xtrain)

# Visualize decision boundary with home-made function (available in general.R)
#pdf('linearsvm.pdf',width=8)
plotlinearsvm2D(svp,xtrain)
#dev.off()

# Predict labels on test
ypred = predict(svp,xtest)
# Confusion table
table(ytest,ypred)

pred <- prediction(ypred,ytest)

# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf,title="TPR vs. FPR")

# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf,title="Precision-Recall")

# Effect of C -------------------------------------------------------------

# Check the influence of C
clist <- 2^seq(-1,10)
errlin <- numeric(length(clist))
errrbf <- numeric(length(clist))
for (i in seq(length(clist))) {
  cat('.')
  svp <- ksvm(xbal,ybal,type="C-svc",kernel='vanilladot',C=clist[i],cross=3)
  errlin[i] <- cross(svp)
  svp <- ksvm(xbal,ybal,type="C-svc",kernel='rbf',C=clist[i],cross=3)
  errrbf[i] <- cross(svp)
}

# Plot the CV error as a function of C

plot(clist,errlin,type='l',log="x",ylim=c(0,0.5),xlab="C",ylab="Error rate",col=1,lwd=2,main='Splice Sites GC only')
lines(clist,errrbf,col=2,lwd=2)
grid()
legend("topright",c('Linear','RBF Gaussian'),lwd=2,col=c(1,2))


# Effect of unbalanced data -----------------------------------------------

# Predict type of disease from expression data
y <- fdata$class
x <- fdata[c(2:3)]
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
svp <- ksvm(xtrain,ytrain,type="C-svc",kernel='vanilladot',C=10,scaled=c())

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
table(ytest,ypred)

pred <- prediction(ypred,ytest)

# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf,title="TPR vs. FPR")

# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf,title="Precision-Recall")

