plotlinearsvm2D=function(svp,xtrain)
## Pretty plot a linear SVM with decision boundary ##
## xtrain should be a 2-dimensional data.
{
	# Define the range of the plot
	# First column is plotted vertically
	yr <- c(min(xtrain[,1]), max(xtrain[,1]))
	# Second column is plotted horizontally
	xr <- c(min(xtrain[,2]), max(xtrain[,2]))

	# Plot the points of xtrain with different signs for positive/negative and SV/non SV
	plot(xr,yr,type='n')
	ymat <- ymatrix(svp)
	points(xtrain[-SVindex(svp),2], xtrain[-SVindex(svp),1], pch = ifelse(ymat[-SVindex(svp)] < 0, 2, 1))
	points(xtrain[SVindex(svp),2], xtrain[SVindex(svp),1], pch = ifelse(ymat[SVindex(svp)] < 0, 17, 16))
	
	# Extract w and b from the model	
	w <- colSums(coef(svp)[[1]] * xtrain[SVindex(svp),])
	b <- b(svp)
	
	# Draw the lines 
	abline(b/w[1],-w[2]/w[1])
	abline((b+1)/w[1],-w[2]/w[1],lty=2)
	abline((b-1)/w[1],-w[2]/w[1],lty=2)
}




cv.folds <- function(n,nfolds=3)
## Randomly split the n samples into folds
## Returns a list of nfolds lists of indices, each corresponding to a fold
{
 	return(split(sample(n),rep(1:nfolds,length=n)))
}
 
cvpred.ksvm <- function(x,y,folds=3,predtype="response",...)
## Return a vector of predictions by cross-validation
## 'predtype' should be one of response (by default), decision or probabilities, depending the prediction we want (SVM label, score or probability, see predict.ksvm())
## Additional parameters are passed to ksvm() to train the SVM
{
 	n <- length(y)
 	ypred <- numeric(n)
 	s <- cv.folds(n,folds)
 	for (i in seq(folds)) {
 		m <- ksvm(x[-s[[i]],],y[-s[[i]]],...)
 		ypred[s[[i]]] <- predict(m,x[s[[i]],],type=predtype)
 		}
 	invisible(ypred)
 }
 	

cvpred.precomp.ksvm <- function(K,y,folds=3,predtype="decision",...)
## Cross-validatin predictions for precomputed kernels
{
 	n <- length(y)
 	ypred <- numeric(n)
 	s <- cv.folds(n,folds)
 	for (i in seq(folds)) {
 		m <- ksvm(as.kernelMatrix(K[-s[[i]],-s[[i]]]),y[-s[[i]]],...)
 		ktest <- as.kernelMatrix(K[s[[i]],-s[[i]]][,SVindex(m),drop=F])
 		ypred[s[[i]]] <- predict(m,ktest,type=predtype)
 		}
 	invisible(ypred)
}