##-----------------------------------------##
## Practical session: user-defined kernels ##
##                                         ##
## Jean-Philippe Vert, 10/11/2010          ##
##-----------------------------------------##

# This file is available at
# http://cbio.ensmp.fr/~jvert/svn/tutorials/practical/makekernel/makekernel.R
# The corresponding note is available at
# http://cbio.ensmp.fr/~jvert/svn/tutorials/practical/makekernel/makekernel_notes.pdf



## Generate toy data ##
##-------------------##

n <- 100 # number of data points
p <- 2   # dimension
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


## Train a SVM with precomputed kernels ##
##--------------------------------------##

library(kernlab)

# Compute by hands a Gaussian RBF kernel
sigma=1
kk <- tcrossprod(x)
dd <- diag(kk)

# Warning: the sigma in kernlab's RBF kernel is not the usual one
# see help(kernels) for definition
myK <- exp(sigma*(-matrix(dd,n,n)-t(matrix(dd,n,n))+2*kk))

# Train a SVM with the precomputed kernel:
# either we explicitly convert myK to a kernelMatrix, and ksvm understands it
svp1 <- ksvm(as.kernelMatrix(myK),y,type="C-svc")
# or we keep it a matrix and we add the kernel='matrix' argument
svp2 <- ksvm(myK,y,type="C-svc",kernel='matrix')

# Compare with the classical way to do
svp3 <- ksvm(x,y,type="C-svc",kernel='rbf',kpar=list(sigma=1),scale=c())


# Compare the 3 formulations, they should be the same
svp1
svp2
svp3


## Predict with a precomputed kernel ##
##-----------------------------------##

# Split the data into training set and test set
ntrain <- round(n*0.8) # number of training examples
tindex <- sample(n,ntrain) # indices of training samples

# Train the svm with the kernel over the training points
svp <- ksvm(myK[tindex,tindex],y[tindex],type="C-svc",kernel='matrix')

# Then it becomes tricky. We must compute the test-vs-SV kernel matrix
# First the test-vs-train matrix
testK <- myK[-tindex,tindex]
# then we extract the SV from the train
testK <- testK[,SVindex(svp),drop=FALSE]

# Predict with the SVM
# Warning: here we MUST convert the matrix testK to a kernelMatrix
ypred <- predict(svp,as.kernelMatrix(testK))

# Do the same with the usual formulation
svp <- ksvm(x[tindex,],y[tindex],type='C-svc',kernel='rbf',kpar=list(sigma=1),scale=c())
ypred2 <- predict(svp,x[-tindex,])

# Check that the predictions are the same
table(ypred,ypred2)

# Check the performance
table(ypred,y[-tindex])
cat('Error rate = ',sum(ypred!=y[-tindex])/length(ypred))


## Understand the class kernel ##
##-----------------------------##


# A object of class 'kernel' is a function with additional slot for kpar
# Look at how they are created
vanilladot
rbfdot

# Let us create a RBF kernel and look at its attributes
rbf <- rbfdot(sigma=1)
rbf
rbf@.Data # the kernel function itself
rbf@kpar  # the kernel paramters
rbf@class # the class


# Once we have a kernel object such as rbf, we can do several things, eg:

# 1) Compute kernel between two vectors
rbf(x[1,],x[2,])

# 2) Compute a kernel matrix
# between two sets of vectors
K <- kernelMatrix(rbf,x[1:5,],x[6:n,])
dim(K)
# or between a set of vectors and itself
K <- kernelMatrix(rbf,x)
dim(K)

# 3) Train a SVM
m <- ksvm(x,y,kernel=rbf,scale=c())

# 4) Look at the points with kernel PCA
kpc <- kpca(x,kernel=rbf,scale=c())
plot(rotated(kpc),col=ifelse(y>0,1,2))


## Make your own kernel ##
##----------------------##

# To make a kernel let us first define the function.
# Eg, to make our own linear kernel:
kval <- function(x, y = NULL) {
	if (is.null(y)) {
		crossprod(x)
	} else {
		crossprod(x,y)
	}
}

# We then create the kernel object as follows
mylinearK <- new("kernel",.Data=kval,kpar=list())

# Now we can call different function of kernlab
mylinearK(x[1,],x[2,])
kernelMatrix(mylinearK,x[1:5,])
m <- ksvm(x,y,kernel=mylinearK)

# Check that we get the same results as the normal vanilla kernel
linearK <- vanilladot()
linearK(x[1,],x[2,])
kernelMatrix(linearK,x[1:5,])
m <- ksvm(x,y,kernel=mylinearK)


## Example: make a kernel that evaluates a precomputed kernel ##
##------------------------------------------------------------##

# We make a kernel object to handle precomputed kernels.
# Its parameter is a precomputed kernel matrix K
# It is then a function of integers such as preK(i,j)=K[i,j]

mypreK <- function(preK=matrix())
{
	rval <- function(x, y = NULL) {
		## we assume that x and y are just indices to be evaluated
		if (is.null(y)) {
			preK[x,x]
		} else {
			preK[x,y]
		}
	}
	return(new("kernel", .Data=rval, kpar=list(preK = preK)))
}

# Check that it works
mypre <- mypreK(myK)
myK[seq(5),seq(5)]
kernelMatrix(mypre,seq(5))

# Also to train SVM
svp3 <- ksvm(seq(n),y,type="C-svc",kernel=mypre,scale=c())
predic(svp3,seq(10))
svp3
svp1