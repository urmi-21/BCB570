library(readr)
library("igraph", lib.loc="~/R/win-library/3.4")
library("Hmisc", lib.loc="~/R/win-library/3.4")
library("gplots", lib.loc="~/R/win-library/3.4")
library("MCL", lib.loc="~/R/win-library/3.4")
setwd("D:/courses/sem4/bcb570/hw/final_proj")
mapRes <- read_delim("mapRes.txt", "\t",escape_double = FALSE, trim_ws = TRUE)
mapRes_hc<-mapRes[which(mapRes$combined_score>=900),]

##create mapres_hc 
protList <- read_csv("protList.txt", col_names = FALSE)

adjMat_PPI<-matrix(nrow = length(protList$X1),ncol=length(protList$X1))
adjMat_PPI[T]<-0
colnames(adjMat_PPI)<-protList$X1
rownames(adjMat_PPI)<-protList$X1

for(k in 1:length(mapRes_hc$combined_score)){
  i<-match(as.character(mapRes_hc$Prot1[k]),rownames(adjMat_PPI))
  j<-match(as.character(mapRes_hc$Prot2[k]),colnames(adjMat_PPI))
  adjMat_PPI[i,j]<-1
  print(k)
}

#create network
network_PPI<-graph_from_adjacency_matrix(adjMat_PPI,diag = F,add.rownames = T,mode = "undirected")
plot(degree.distribution(network_PPI),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "distribution")
vcount(network_PPI)
ecount(network_PPI)
graph.density(network_PPI)

ppi_network_details<-c(ecount(network_PPI), vcount(network_PPI), graph.density(network_PPI),diameter(network_PPI), radius(network_PPI), transitivity(network_PPI, type = "global"),transitivity(network_PPI, type = "average"),mean_distance(network_PPI))

header1<-c("#Edges","#Vertices","Density","Diameter","Radius","Clustering_coeff(global)","Clustering_coeff(avg)","Avg Shortest Path")
tab2<-data.frame(Feature=header1,PPI_network_=ppi_network_details,gene_network_GENIE=gene_network_GENIE_details,prot_network_GENIE=prot_network_GENIE_details)
write.csv(tab2,"ppiinfo.csv",row.names = F)
#find clusters
PPI_mcl <- mcl(x = adjMat_PPI, addLoops=TRUE, ESM = TRUE, max.iter = 500)
save.image("PPI.RData")

# get all clusters
PPI_clusters<-cbind(rownames(adjMat_PPI),PPI_mcl$Cluster)
write.csv(PPI_clusters,"PPI_clusters.csv")
#make list of clusters
PPI_cluster_list<-list()
c<-1
for(i in unique(PPI_mcl$Cluster)){
  thiscluster<-V(network_PPI)[ which(PPI_mcl$Cluster == i)]$name
  if(length(thiscluster)>=0){
    PPI_cluster_list[[c]]<-thiscluster
    c<-c+1
  }
}
ppi_c_lens<-c()
for (l in PPI_cluster_list){
  ppi_c_lens<-c(ppi_c_lens,length(l))
}
##cluster lens
hist(ppi_c_lens,breaks = 20,col = "Green")
plot(sort(ppi_c_lens),type="l",log="y",col="Green")
  



