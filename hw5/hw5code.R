setwd("D:/courses/sem4/bcb570/hw/hw5")
library(GENIE3)
n_data_100<-rnorm(100,mean = 2,sd = sqrt(8))
cat("Mean is",mean(n_data_100))
cat("Var is",var(n_data_100))

means<-c()
vars<-c()
N<-20
k<-10
for(i in 1:N){
  this<-sample(n_data_100,k,replace = T)
  #print(this)
  means<-c(means,mean(this))
  vars<-c(vars,var(this))
}

hist(means,col = c("Red"))
hist(vars,col = c("Red"))
mean(means)
mean(vars)


#run GENIE3
library(GENIE3)
library(readr)
library("doParallel", lib.loc="~/R/win-library/3.4")
library("doRNG", lib.loc="~/R/win-library/3.4")
library(plyr)
library(PRROC)
#read true labels
trueEcoli1<-read.csv("NIHW in silico data/NIHW in silico data/Size100/DREAM3 gold standards/DREAM3GoldStandard_InSilicoSize100_Ecoli1.txt",sep = "\t",header = F)
trueEcoli2<-read.csv("NIHW in silico data/NIHW in silico data/Size100/DREAM3 gold standards/DREAM3GoldStandard_InSilicoSize100_Ecoli2.txt",sep = "\t",header = F)
#do for ecoli data 1
adata<-read.csv("NIHW in silico data/NIHW in silico data/Size100/DREAM3 data/InSilicoSize100-Ecoli1-trajectories.tsv",sep = "\t")
adata<-t(adata)
expmat<-as.matrix(adata[2:101,1:966])
rownames(expmat) <- paste("G", 1:100, sep="")
colnames(expmat) <- paste("Samp", 1:966, sep="")

res1<-GENIE3(expmat,nCores = 6,nTrees = 100)
linkList <- getLinkList(res1)
names(linkList)<-c("RG","TG","Weight")
names(linkList_q4)<-c("RG","TG","Weight")
names(trueEcoli1)<-c("RG","TG","Label")

q4thresh<-quantile(linkList$Weight,0.75)
#q4thresh<-0.1
linkList_q4<-linkList[linkList$Weight >= q4thresh,]
joined_data<-plyr::join(linkList,trueEcoli1)
joined_data_q4<-plyr::join(linkList_q4,trueEcoli1)
pr1 <- pr.curve(scores.class0 = joined_data$Weight, weights.class0 = joined_data$Label, curve = T)
plot(pr1,main = "PR Curve for Ecoli dataset 1")
pr1_1 <- pr.curve(scores.class0 = joined_data_q4$Weight, weights.class0 = joined_data_q4$Label, curve = T)
plot(pr1_1,main = "PR Curve for Ecoli dataset 1 (only top quartile edges)")

#do for ecoli data 2
adata2<-read.csv("NIHW in silico data/NIHW in silico data/Size100/DREAM3 data/InSilicoSize100-Ecoli2-trajectories.tsv",sep = "\t")
adata2<-t(adata2)
expmat2<-as.matrix(adata2[2:101,1:966])
rownames(expmat2) <- paste("G", 1:100, sep="")
colnames(expmat2) <- paste("Samp", 1:966, sep="")

res2<-GENIE3(expmat2,nCores = 6,K="all",nTrees = 100)
linkList2 <- getLinkList(res2)
names(linkList2)<-c("RG","TG","Weight")
names(linkList2_q4)<-c("RG","TG","Weight")
names(trueEcoli2)<-c("RG","TG","Label")
q4thresh2<-quantile(linkList2$Weight,0.75)
linkList2_q4<-linkList2[linkList2$Weight >= q4thresh2,]
joined_data2<-plyr::join(linkList2,trueEcoli2)
joined_data2_q4<-plyr::join(linkList2_q4,trueEcoli2)
pr2 <- pr.curve(scores.class0 = joined_data2$Weight, weights.class0 = joined_data2$Label, curve = T)
plot(pr2,main = "PR Curve for Ecoli dataset 2")
pr2_1 <- pr.curve(scores.class0 = joined_data2_q4$Weight, weights.class0 = joined_data2_q4$Label, curve = T)
plot(pr2_1,main = "PR Curve for Ecoli dataset 2 (only top quartile edges)")


#save results to file
write_tsv(linkList,"Ecoli1_results_all.tsv")
write_tsv(linkList_q4,"Ecoli1_results_topq.tsv")
write_tsv(linkList2,"Ecoli2_results_all.tsv")
write_tsv(linkList2_q4,"Ecoli2_results_topq.tsv")

#find best estimate intersection from linkList_q4 and linkList2_q4
temp_linkList_q4<-linkList_q4
names(temp_linkList_q4)<-c("RG","TG","ecoli1")
temp_linkList2_q4<-linkList2_q4
names(temp_linkList2_q4)<-c("RG","TG","ecoli2")
#common edges
common<-plyr::join(temp_linkList_q4,temp_linkList2_q4,type="inner")
#save to file
write_tsv(common,"Ecoli_best_estimate.tsv")
