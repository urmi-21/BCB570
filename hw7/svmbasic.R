##--------------------------------------##
## Practical session: initiation to SVM ##
##                                      ##
## Jean-Philippe Vert, 10/11/2010       ##
##--------------------------------------##

# This file is available at
# http://cbio.ensmp.fr/~jvert/svn/tutorials/practical/svmbasic/svmbasic.R
# The corresponding note is available at
# http://cbio.ensmp.fr/~jvert/svn/tutorials/practical/svmbasic/svmbasic_notes.pdf


# Various functions that should be programmed by the students
# Available at http://cbio.ensmp.fr/~jvert/svn/tutorials/practical/svmbasic/general.R
source('general.R')

## Generate data ##
##---------------##

# number of data points
n <- 150
# dimension
p <- 2

sigma <- 2  # variance of the distribution
meanpos <- 0 # centre of the distribution of positive examples
meanneg <- 3 # centre of the distribution of negative examples
npos <- round(n/2) # number of positive examples
nneg <- n-npos # number of negative examples

# Generate the positive and negative examples
xpos <- matrix(rnorm(npos*p,mean=meanpos,sd=sigma),npos,p)
xneg <- matrix(rnorm(nneg*p,mean=meanneg,sd=sigma),npos,p)
x <- rbind(xpos,xneg)

# Generate the labels
y <- matrix(c(rep(1,npos),rep(-1,nneg)))

# Visualize the data
plot(x,col=ifelse(y>0,1,2))
legend("topleft",c('Positive','Negative'),col=seq(2),pch=1,text.col=seq(2))

## Prepare a training and a test set ##
ntrain <- round(n*0.8) # number of training examples
tindex <- sample(n,ntrain) # indices of training samples
xtrain <- x[tindex,]
xtest <- x[-tindex,]
ytrain <- y[tindex]
ytest <- y[-tindex]
istrain=rep(0,n)
istrain[tindex]=1

# Visualize
plot(x,col=ifelse(y>0,1,2),pch=ifelse(istrain==1,1,2))
legend("topleft",c('Positive Train','Positive Test','Negative Train','Negative Test'),col=c(1,1,2,2),pch=c(1,2,1,2),text.col=c(1,1,2,2))


## Train a linear SVM with C=1 ##
##-----------------------------##

# load the kernlab package
library(kernlab)

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

# Compute accuracy
sum(ypred==ytest)/length(ytest)

# Compute at the prediction scores
ypredscore = predict(svp,xtest,type="decision")

# Check that the predicted labels are the signs of the scores
table(ypredscore > 0,ypred)

# Package to compute ROC curve, precision-recall etc...
library(ROCR)

pred <- prediction(ypredscore,ytest)

# Plot ROC curve
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf)

# Plot precision/recall curve
perf <- performance(pred, measure = "prec", x.measure = "rec") 
plot(perf)

# Plot accuracy as function of threshold
perf <- performance(pred, measure = "acc") 
plot(perf)



## Use cross-validation ##
##----------------------##

# We compute the prediction vector by cross-validation
k=5	
ypredscorecv <- cvpred.ksvm(x,y,folds=k,type="C-svc",kernel='vanilladot',C=1,scaled=c(),predtype="decision")

# Check the performance
print(table(ypredscorecv > 0,y))
pred <- prediction(ypredscorecv,y)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf)

# Estimate the CV error with ksvm directly, and compare
svp <- ksvm(x,y,type="C-svc",kernel='vanilladot',C=1,scaled=c(),cross=5)
print(cross(svp))
print(1-sum((ypredscorecv>0)==(y==1))/n)


## Effect of C ##
##-------------##

# Plot SVM for different values of C
# Observe how the margin changes
clist <- 2^seq(-10,10)
err <- numeric(length(clist))
par(ask=T)
for (i in seq(length(clist))) {
	svp <- ksvm(x,y,type="C-svc",kernel='vanilladot',C=clist[i],scaled=c(),cross=3)
	plotlinearsvm2D(svp,x)
	err[i] <- cross(svp)
}
par(ask=F)
# Plot the CV error as a function of C
plot(clist,err,type='l',log="x",ylim=c(0,1),xlab="C",ylab="Error rate")
grid()



## Nonlinear SVM ##
##---------------##

# Make a toy example

nclust <- 30 # number of points in each cluster
n <- 4*nclust # number of data points
p <- 2   # dimension
sigma <- 0.8  # variance of the distribution
meanpos <- 0 # centre of the distribution of positive examples
meanneg <- 3 # centre of the distribution of negative examples

# Generate the positive and negative examples
x1 <- rnorm(nclust,mean= meanpos,sd=sigma)
x2 <- rnorm(nclust,mean= meanpos,sd=sigma)
x = cbind(x1,x2)
x1 <- rnorm(nclust,mean= meanneg,sd=sigma)
x2 <- rnorm(nclust,mean= meanneg,sd=sigma)
x = rbind(x,cbind(x1,x2))
x1 <- rnorm(nclust,mean= meanpos,sd=sigma)
x2 <- rnorm(nclust,mean= meanneg,sd=sigma)
x = rbind(x,cbind(x1,x2))
x1 <- rnorm(nclust,mean= meanneg,sd=sigma)
x2 <- rnorm(nclust,mean= meanpos,sd=sigma)
x = rbind(x,cbind(x1,x2))
y = c(rep(1,2*nclust),rep(-1,2*nclust))

# Visualize
#pdf('toynonlinear.pdf',width=8)
plot(x,col=ifelse(y>0,1,2),pch=ifelse(y>0,1,2))
legend("topleft",c('Positive','Negative'),col=c(1,2),pch=c(1,2),text.col=c(1,2))
grid()
#dev.off()


# Train linear SVM
clist <- 2^seq(-6,10)
errlin <- numeric(length(clist))
par(ask=T)
for (i in seq(length(clist))) {
	svp <- ksvm(x,y,type="C-svc",kernel='vanilladot',C=clist[i],scaled=c(),cross=3)
	plotlinearsvm2D(svp,x)
	errlin[i] <- cross(svp)
}
par(ask=F)
# Plot the CV error as a function of C
plot(clist,errlin,type='l',log="x",ylim=c(0,1),xlab="C",ylab="Error rate")
grid()

# Train a nonlinear SVM with Gaussian RBF kernel and default parameters
svp <- ksvm(x,y,type="C-svc",kernel='rbf')
#pdf('nonlinearsvm.pdf',width=8)
plot(svp,data=x)
#dev.off()


# Effect of C only with automatic choice of sigma
errrbf <- numeric(length(clist))
par(ask=T)
for (i in seq(length(clist))) {
	svp <- ksvm(x,y,type="C-svc",kernel='rbf',C=clist[i],cross=3)
	plot(svp,data=x)
	errrbf[i] <- cross(svp)
}
par(ask=F)
# Plot the CV error as a function of C
plot(clist,errlin,type='l',log="x",ylim=c(0,1),xlab="C",ylab="Error rate",col=1,lwd=2)
lines(clist,errrbf,col=2,lwd=2)
grid()
legend("topright",c('Linear','RBF Gaussian'),lwd=2,col=c(1,2))

# Effect of C and sigma
nc <- length(clist)
sigmalist <- 2^seq(-6,6)
nsigma <- length(sigmalist)
err <- matrix(0,nrow=nc,ncol=nsigma)
for (i in seq(nc)) {
	C <- clist[i]
	cat('C=',C)
	for (j in seq(nsigma)) {
		cat('.')
		svp <- ksvm(x,y,type="C-svc",kernel='rbf',kpar=list(sigma=sigmalist[j]),C=C,cross=3)
		err[i,j] <- cross(svp)
	}
	cat('\n')
}
# Visualize
library(lattice)
dimnames(err) <- list(C=clist,sigma=sigmalist)
levelplot(err,scales=list(x=list(rot=90)),xlab="C",ylab=expression(sigma),main="Error rate")

## Application on  ALL gene expression data ##
##------------------------------------------##

# Load the ALL dataset
source("https://bioconductor.org/biocLite.R")
biocLite("ALL")
library(ALL)
data(ALL)

# Inspect them
?ALL
show(ALL)
print(summary(pData(ALL)))

# Look at type and stage of each patient
print(ALL$BT)

# Predict type of disease from expression data
y <- substr(ALL$BT,1,1)
x <- t(exprs(ALL))
dim(x)
length(y)

# Train a linear SVM
ypred <- cvpred.ksvm(x,y,type="C-svc",kernel='vanilladot',C=1)
print(table(ypred,y))

# Train a nonlinear SVM
ypred <- cvpred.ksvm(x,y,type="C-svc",kernel='rbf',C=1)
print(table(ypred > 0,y))

# Check the influence of C
clist <- 2^seq(-6,10)
errlin <- numeric(length(clist))
errrbf <- numeric(length(clist))
for (i in seq(length(clist))) {
	cat('.')
	svp <- ksvm(x,y,type="C-svc",kernel='vanilladot',C=clist[i],cross=3)
	errlin[i] <- cross(svp)
	svp <- ksvm(x,y,type="C-svc",kernel='rbf',C=clist[i],cross=3)
	errrbf[i] <- cross(svp)
}

# Plot the CV error as a function of C
plot(clist,errlin,type='l',log="x",ylim=c(0,0.5),xlab="C",ylab="Error rate",col=1,lwd=2,main='ALL dataset')
lines(clist,errrbf,col=2,lwd=2)
grid()
legend("topright",c('Linear','RBF Gaussian'),lwd=2,col=c(1,2))

