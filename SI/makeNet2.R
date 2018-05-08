library(readr)
library("igraph", lib.loc="~/R/win-library/3.4")
library("Hmisc", lib.loc="~/R/win-library/3.4")
library("gplots", lib.loc="~/R/win-library/3.4")
library(GENIE3)
library("doParallel", lib.loc="~/R/win-library/3.4")
library("doRNG", lib.loc="~/R/win-library/3.4")
library("MCL", lib.loc="~/R/win-library/3.4")
library("mixOmics", lib.loc="~/R/win-library/3.4")
library("yaml", lib.loc="~/R/win-library/3.4")
library("stringi", lib.loc="~/R/win-library/3.4")
setwd("D:/courses/sem4/bcb570/hw/final_proj")

Urminder <- read_csv("Urminder.txt", col_names = FALSE)
Protein <- read_csv("Protein.txt")
#remove extra prots
Protein<-subset(Protein,(rownames(Protein) %in% rownames(Urminder)))

###gene network using genei3
gene_expMat<-as.matrix(Urminder[,2:length(Urminder)])
rownames(gene_expMat)<-Urminder$X1
res<-GENIE3(gene_expMat,nCores = 8,nTrees = 1000,verbose = T)
res[is.nan(res)] <- 0

gene_linkList <- getLinkList(res, threshold=0)
minw<-min(gene_linkList$weight)
maxw<-max(gene_linkList$weight)
step<-(maxw-minw)/10
gene_densities<-c()
num_edges<-c()
for(t in seq(minw,maxw,step)){
  print(t)
  #gene_linkList <- getLinkList(res, threshold=t)
  afterT<-gene_linkList[which(gene_linkList$weight>=t),]
  gene_elist<-cbind(as.character(afterT$regulatoryGene),as.character(afterT$targetGene))
  #create network
  #gene_network_GENIE<-graph_from_adjacency_matrix(gene_adj_new,diag = F,add.rownames = T,mode = "undirected")
  gene_network_GENIE<-graph_from_edgelist(gene_elist,directed = F)
  plot(degree.distribution(gene_network_GENIE),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "distribution")
  vcount(gene_network_GENIE)
  edgs<-ecount(gene_network_GENIE)
  dens<-graph.density(gene_network_GENIE)
  gene_densities<-c(gene_densities,dens)
  num_edges<-c(num_edges,edgs)
}
plot(gene_densities,x = seq(minw,maxw,step),type = "l")
plot(num_edges,x = seq(minw,maxw,step),type = "l")

#from graph t= 0.0086
afterT<-gene_linkList[which(gene_linkList$weight>=0.0086),]
gene_elist<-cbind(as.character(afterT$regulatoryGene),as.character(afterT$targetGene))
gene_network_GENIE<-graph_from_edgelist(gene_elist,directed = F)
plot(degree.distribution(gene_network_GENIE),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "distribution")
vcount(gene_network_GENIE)
edgs<-ecount(gene_network_GENIE)
dens<-graph.density(gene_network_GENIE)


#####################################################



###prot network using genei3
prot_expMat<-as.matrix(Protein[,2:length(Protein)])
rownames(prot_expMat)<-Protein$Gene
res_prot<-GENIE3(prot_expMat,nCores = 6,nTrees = 1000,verbose = T)
res_prot[is.nan(res_prot)] <- 0

prot_linkList <- getLinkList(res_prot, threshold=0)
minw<-min(prot_linkList$weight)
maxw<-max(prot_linkList$weight)
step<-(maxw-minw)/10
prot_densities<-c()
num_edges_prot<-c()
for(t in seq(minw,maxw,step)){
  print(t)
  #gene_linkList <- getLinkList(res, threshold=t)
  afterT<-prot_linkList[which(prot_linkList$weight>=t),]
  prot_elist<-cbind(as.character(afterT$regulatoryGene),as.character(afterT$targetGene))
  #create network
  #gene_network_GENIE<-graph_from_adjacency_matrix(gene_adj_new,diag = F,add.rownames = T,mode = "undirected")
  prot_network_GENIE<-graph_from_edgelist(prot_elist,directed = F)
  plot(degree.distribution(prot_network_GENIE),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "distribution")
  vcount(prot_network_GENIE)
  edgs<-ecount(prot_network_GENIE)
  dens<-graph.density(prot_network_GENIE)
  prot_densities<-c(prot_densities,dens)
  num_edges_prot<-c(num_edges_prot,edgs)
}
plot(prot_densities,x = seq(minw,maxw,step),type = "l",xlab = "threshold")
plot(num_edges_prot,x = seq(minw,maxw,step),type = "l",xlab = "threshold")

#from graph t= 0.0012
afterT<-prot_linkList[which(prot_linkList$weight>=0.006308359),]
prot_elist<-cbind(as.character(afterT$regulatoryGene),as.character(afterT$targetGene))
prot_network_GENIE<-graph_from_edgelist(prot_elist,directed = F)
plot(degree.distribution(prot_network_GENIE),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "distribution")
vcount(prot_network_GENIE)
edgs<-ecount(prot_network_GENIE)
dens<-graph.density(prot_network_GENIE)

#################################################################
##save as edgelist
gene_el_final<-as_edgelist(gene_network_GENIE)
prot_el_final<-as_edgelist(prot_network_GENIE)
write.csv(gene_el_final,"gedge_final.txt",row.names = F)
write.csv(prot_el_final,"pedge_final.txt",row.names = F)

#find clusters using mcl
gene_mcl <- mcl(x = gene_adj_new, addLoops=TRUE, ESM = TRUE, max.iter = 500)
prot_mcl <- mcl(x = prot_adj_new, addLoops=TRUE, ESM = TRUE, max.iter = 500)
save.image("mclres.Rdata")

#find largest clusters in Y2h components
findlargestcluster<-function(graph,mcldata){
  res<-0
  cres<-V(graph)[ which(mcldata$Cluster %in% c(0))]
  t<-mcldata$K
  for(i in c(0:t))
    #print (i)
    thiscluster<-V(graph)[ which(mcldata$Cluster %in% c(i))]$name
  #print (length(thiscluster))
  if(length(thiscluster) > length(cres)){
    cres <- thiscluster
    res <- i
  }
  return(cres)
}
findlargestcluster(gene_network_GENIE,gene_mcl)
gene_mcl$K

gene_cluster_list<-list()
c<-1
for(i in unique(gene_mcl$Cluster)){
  thiscluster<-V(gene_network_GENIE)[ which(gene_mcl$Cluster == i)]$name
  if(length(thiscluster)>=0){
    gene_cluster_list[[c]]<-thiscluster
    c<-c+1
  }
}

gene_c_lens<-c()
for (l in gene_cluster_list){
  gene_c_lens<-c(gene_c_lens,length(l))
}

prot_cluster_list<-list()
c<-1
for(i in unique(prot_mcl$Cluster)){
  thiscluster<-V(prot_network_GENIE)[ which(prot_mcl$Cluster == i)]$name
  if(length(thiscluster)>=0){
    prot_cluster_list[[c]]<-thiscluster
    c<-c+1
  }
}

prot_c_lens<-c()
for (l in prot_cluster_list){
  prot_c_lens<-c(prot_c_lens,length(l))
}

##cluster lens
hist(gene_c_lens,breaks = 20,col = "red")
plot(sort(gene_c_lens),type="l",log="y",col="red")
hist(prot_c_lens,breaks = 20,col="blue")
plot(sort(prot_c_lens),type="l",log="y",col="blue")


##write cluster to files
# get all clusters
gene_clusters<-cbind(rownames(gene_adj_new),gene_mcl$Cluster)
write.csv(gene_clusters,"gene_clusters.csv")
prot_clusters<-cbind(rownames(prot_adj_new),prot_mcl$Cluster)
write.csv(prot_clusters,"prot_clusters.csv")

##get Network properties
gene_network_GENIE_details<-c(ecount(gene_network_GENIE), vcount(gene_network_GENIE), graph.density(gene_network_GENIE),diameter(gene_network_GENIE), radius(gene_network_GENIE), transitivity(gene_network_GENIE, type = "global"),transitivity(gene_network_GENIE, type = "average"),mean_distance(gene_network_GENIE))
prot_network_GENIE_details<-c(ecount(prot_network_GENIE), vcount(prot_network_GENIE), graph.density(prot_network_GENIE),diameter(prot_network_GENIE), radius(prot_network_GENIE), transitivity(prot_network_GENIE, type = "global"),transitivity(prot_network_GENIE, type = "average"),mean_distance(prot_network_GENIE))

header1<-c("#Edges","#Vertices","Density","Diameter","Radius","Clustering_coeff(global)","Clustering_coeff(avg)","Avg Shortest Path")
tab1<-data.frame(Feature=header1,gene_network_GENIE=gene_network_GENIE_details,prot_network_GENIE=prot_network_GENIE_details)

#scale free
fitA <- power.law.fit(igraph::degree(gene_network_GENIE,mode="all"),implementation="plfit")
fitB <- power.law.fit(igraph::degree(prot_network_GENIE,mode="all"),implementation="plfit")
fitA
fitB
#use powerlaw fit
head3<-c("Gene network","Protein network")
tab3<-data.frame(Graph=head3,Alpha=c(fitA$alpha,fitB$alpha),X_min=c(fitA$xmin,fitB$xmin),KS_stat=c(fitA$KS.stat,fitB$KS.stat),P_Value=c(fitA$KS.p,fitB$KS.p),logLik=c(fitA$logLik,fitB$logLik))
write.csv(tab3,"plaw.csv",row.names = F)

gene_edges<-as_edgelist(gene_network_GENIE)
prot_edges<-as_edgelist(prot_network_GENIE)
write.csv(gene_edges,"gedge.txt",row.names = F,col.names = F)
write.csv(prot_edges,"pedge.txt",row.names = F,col.names = F)

write.csv(Urminder$X1[which(rowSums(Urminder[2:19])<=0)],"lowexp.txt",row.names = F)
write.csv(Urminder$X1[which(rowSums(Urminder[2:19])>=10000)],"highexp.txt",row.names = F)



###DO PCA



