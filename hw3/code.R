library("igraph")
library(readxl)
setwd("D:/courses/sem4/bcb570/hw/hw3")
##read graph

Y2H_uni<-read.table('Y2H_uniondata.txt')
graph_Y2H_uni<-graph.data.frame(as.matrix(Y2H_uni),directed = FALSE)
graph_Y2H_uni_Details<-c(ecount(graph_Y2H_uni), vcount(graph_Y2H_uni), graph.density(graph_Y2H_uni),diameter(graph_Y2H_uni), radius(graph_Y2H_uni), transitivity(graph_Y2H_uni, type = "global"),transitivity(graph_Y2H_uni, type = "average"),mean_distance(graph_Y2H_uni))

CCSB_YI1<-read.table('CCSB_YI1.txt')
graph_CCSB_YI1<-graph.data.frame(as.matrix(CCSB_YI1),directed = FALSE)
graph_CCSB_YI1_Details<-c(ecount(graph_CCSB_YI1), vcount(graph_CCSB_YI1), graph.density(graph_CCSB_YI1),diameter(graph_CCSB_YI1), radius(graph_CCSB_YI1), transitivity(graph_CCSB_YI1, type = "global"),transitivity(graph_CCSB_YI1, type = "average"),mean_distance(graph_CCSB_YI1))

#wtite to table
header1<-c("#Edges","#Vertices","Density","Diameter","Radius","Clustering_coeff(global)","Clustering_coeff(avg)","Avg Shortest Path")
tab1<-data.frame(Feature=header1,Y2H_uni=graph_Y2H_uni_Details,CCSB_YI1=graph_CCSB_YI1_Details)

#generate 3 random netwroks for each graph
rn1a <- erdos.renyi.game(vcount(graph_Y2H_uni), ecount(graph_Y2H_uni),type='gnm')
rn1b <- erdos.renyi.game(vcount(graph_Y2H_uni), ecount(graph_Y2H_uni),type='gnm')
rn1c <- erdos.renyi.game(vcount(graph_Y2H_uni), ecount(graph_Y2H_uni),type='gnm')
rn2a <- erdos.renyi.game(vcount(graph_CCSB_YI1), ecount(graph_CCSB_YI1),type='gnm')
rn2b <- erdos.renyi.game(vcount(graph_CCSB_YI1), ecount(graph_CCSB_YI1),type='gnm')
rn2c <- erdos.renyi.game(vcount(graph_CCSB_YI1), ecount(graph_CCSB_YI1),type='gnm')

#get cc and mean dist
cc1a <- transitivity(rn1a, type = "global")
cc1b <- transitivity(rn1b, type = "global")
cc1c <- transitivity(rn1c, type = "global")
cc2a <- transitivity(rn2a, type = "global")
cc2b <- transitivity(rn2b, type = "global")
cc2c <- transitivity(rn2c, type = "global")

md1a <- mean_distance(rn1a)
md1b <- mean_distance(rn1b)
md1c <- mean_distance(rn1c)
md2a <- mean_distance(rn2a)
md2b <- mean_distance(rn2b)
md2c <- mean_distance(rn2c)

header2<-c("Clustering_coeff(global)","Avg Shortest Path")
tab2<-data.frame(Network=header2,Y2H_uni_randomized_a=c(cc1a,md1a),Y2H_uni_randomized_b=c(cc1b,md1b),Y2H_uni_randomized_c=c(cc1c,md1c),CCSB_YI1_randomized_a=c(cc2a,md2a),CCSB_YI1_randomized_b=c(cc2b,md2b),CCSB_YI1_randomized_c=c(cc2c,md2c))
tab2<-t(tab2)

#plot degree dist
par(mfrow=c(1,3))
plot(degree.distribution(graph_Y2H_uni),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "a: Degree Distribution Y2H_uni")
plot(degree.distribution(graph_CCSB_YI1),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "b: Degree Distribution CCSB_YI1")
g_powlaw <- barabasi.game(100000)
plot(degree.distribution(g_powlaw),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "c: Powerlaw distribution")
#fit powerlaw
fitA <- power.law.fit(degree(graph_Y2H_uni,mode="all"),implementation="plfit") 
fitB <- power.law.fit(degree(graph_CCSB_YI1,mode="all"),implementation="plfit") 
fitC <- power.law.fit(degree(g_powlaw,mode="all"),implementation="plfit") 
head3<-c("Y2H_uni","CCSB_YI1","Powerlaw graph")
tab3<-data.frame(Graph=head3,Alpha=c(fitA$alpha,fitB$alpha,fitC$alpha),X_min=c(fitA$xmin,fitB$xmin,fitC$xmin),KS_stat=c(fitA$KS.stat,fitB$KS.stat,fitC$KS.stat),P_Value=c(fitA$KS.p,fitB$KS.p,fitC$KS.p),logLik=c(fitA$logLik,fitB$logLik,fitC$logLik))
#find hubs in the networks
#first plot degree dist
plot(degree.distribution(graph_Y2H_uni),type = "l", xlab = "Degree",main = "Degree distribution Y2H_uni", xaxt="n")
axis(side=1, at=c(0:100))
plot(degree.distribution(graph_CCSB_YI1),type = "l", xlab = "Degree",main = "Degree distribution CCSB_YI1", xaxt="n")
axis(side=1, at=c(0:100))
hubs_Y2H_uni<-names(which(degree(graph_Y2H_uni)>11))
hubs_CCSB_YI1<-names(which(degree(graph_CCSB_YI1)>10))
#hubs_common <- intersect(hubs_Y2H_uni,hubs_CCSB_YI1)
#read yeast deletion data
Yeast_deletionProject <- read_excel("Yeast_deletionProject.xlsx")
#remove extra whitespace from first col
ORFs_Yeast_proj <- as.data.frame(apply(Yeast_deletionProject,2,function(x)gsub('\\s+', '',x)))$ORF
#find common ORFs
hubs_common_Y2H_uni <- intersect(hubs_Y2H_uni,ORFs_Yeast_proj)
hubs_common_CCSB_YI1 <- intersect(hubs_CCSB_YI1,ORFs_Yeast_proj)
hubs_common <- intersect(hubs_common_CCSB_YI1,hubs_common_Y2H_uni)
head4<-c("Y2H_uni","CCSB_YI1","Intersection")
tab4<-data.frame(Network=head4,Hubs_common=c(paste(hubs_common_Y2H_uni,collapse = ","), paste(hubs_common_CCSB_YI1,collapse = ","),paste(hubs_common,collapse = ",")))


#get p val
N=10
n1=length(hubs_Y2H_uni)
C1<-length(hubs_common_Y2H_uni)
T1=0
#repeat N times
nodes_Y2H<-names(degree(graph_Y2H_uni))
for(i in 1:N ){
  #print(i)
  #sample n1 genes from Y2H_uni
  thissample <- sample(nodes_Y2H,n1)
  intersect(thissample,ORFs_Yeast_proj)
  if(length(intersect(thissample,ORFs_Yeast_proj)) >= C1){
    T1=T1+1
  }
  
}
print (T1)
print (T1/N)


#read BioGrid2018_uni-2
BioGrid18_data<-read.table('BioGrid2018_uni-2.txt')
graph_BioGrid18<-graph.data.frame(as.matrix(BioGrid18_data),directed = FALSE)
graph_BioGrid18_Details<-c(ecount(graph_BioGrid18), vcount(graph_BioGrid18), graph.density(graph_BioGrid18),diameter(graph_BioGrid18), radius(graph_BioGrid18), transitivity(graph_BioGrid18, type = "global"),transitivity(graph_BioGrid18, type = "average"),mean_distance(graph_BioGrid18))
head5<-c("#Edges","#Vertices","Density","Diameter","Radius","Clustering_coeff(global)","Clustering_coeff(avg)","Avg Shortest Path")
tab5<-data.frame(Feature=head5,BioGrid18=graph_BioGrid18_Details)
plot(degree.distribution(graph_BioGrid18),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "a: Degree Distribution BioGrid18")
fitD <- power.law.fit(degree(graph_BioGrid18,mode="all"),implementation="plfit") 
head6<-c("Y2H_uni","CCSB_YI1","BioGrid18")
tab6<-data.frame(Graph=head3,Alpha=c(fitA$alpha,fitB$alpha,fitD$alpha),X_min=c(fitA$xmin,fitB$xmin,fitD$xmin),KS_stat=c(fitA$KS.stat,fitB$KS.stat,fitD$KS.stat),P_Value=c(fitA$KS.p,fitB$KS.p,fitD$KS.p),logLik=c(fitA$logLik,fitB$logLik,fitD$logLik))
#check for small world
#generate 3 random netwroks for each graph
rn3a <- erdos.renyi.game(vcount(graph_BioGrid18), ecount(graph_BioGrid18),type='gnm')
cc3a <- transitivity(rn3a, type = "global")
md3a <- mean_distance(rn3a)
head7<-c("Clustering_coeff(global)","Avg Shortest Path")
tab7<-data.frame(Feature=head7,BioGrid18=c(transitivity(graph_BioGrid18, type = "global"),mean_distance(graph_BioGrid18)),Random_Net=c(cc3a,md3a))
tab7 %>% knitr::kable(caption = "BioGrid Network vs Random Network")


#find largest connected comp
comps <- components(graph = graph_Y2H_uni)

getlargestCC<-function(g)
  {
  copms<-components(graph = g)
  #return component with max csize
  for(c in comps){
    print(c$csize)
  }
  
}
getlargestCC(graph_Y2H_uni)


getLargestCC<-function(g)
{
  #decompose graph into connected comps
  decomp<-decompose.graph(g)
  #say first component is largest
  largest<-decomp[[1]]
  for(i in 1:length(decomp))
  {
    thisg <- decomp[[i]]
    #if there is largere component
    if( vcount(thisg) > vcount(largest) )
    {
      largest <- thisg
    }
  }
  #return largest
  return(largest)
}

graph_Y2H_uni_lcc <- getLargestCC(graph_Y2H_uni)
graph_CCSB_YI1_lcc <- getLargestCC(graph_CCSB_YI1)
graph_BioGrid18_lcc <- getLargestCC(graph_BioGrid18)
graph_Y2H_uni_lcc_Details<-c(ecount(graph_Y2H_uni_lcc), vcount(graph_Y2H_uni_lcc), graph.density(graph_Y2H_uni_lcc),diameter(graph_Y2H_uni_lcc), radius(graph_Y2H_uni_lcc), transitivity(graph_Y2H_uni_lcc, type = "global"),transitivity(graph_Y2H_uni_lcc, type = "average"),mean_distance(graph_Y2H_uni_lcc))
graph_CCSB_YI1_lcc_Details<-c(ecount(graph_CCSB_YI1_lcc), vcount(graph_CCSB_YI1_lcc), graph.density(graph_CCSB_YI1_lcc),diameter(graph_CCSB_YI1_lcc), radius(graph_CCSB_YI1_lcc), transitivity(graph_CCSB_YI1_lcc, type = "global"),transitivity(graph_CCSB_YI1_lcc, type = "average"),mean_distance(graph_CCSB_YI1_lcc))
graph_BioGrid18_lcc_Details<-c(ecount(graph_BioGrid18_lcc), vcount(graph_BioGrid18_lcc), graph.density(graph_BioGrid18_lcc),diameter(graph_BioGrid18_lcc), radius(graph_BioGrid18_lcc), transitivity(graph_BioGrid18_lcc, type = "global"),transitivity(graph_BioGrid18_lcc, type = "average"),mean_distance(graph_BioGrid18_lcc))
#wtite to table
head8<-c("#Edges","#Vertices","Density","Diameter","Radius","Clustering_coeff(global)","Clustering_coeff(avg)","Avg Shortest Path")
tab8<-data.frame(Feature=head8,Y2H_uni_lcc=graph_Y2H_uni_lcc_Details,CCSB_YI1_lcc=graph_CCSB_YI1_lcc_Details,BioGrid_lcc=graph_BioGrid18_lcc_Details)


#find most essential genes by betweenness
#Assuming Top 5% are essential
bw_Y2H_uni_lcc<-betweenness(graph_Y2H_uni_lcc,directed = F)
bw_CCSB_YI1_lcc<-betweenness(graph_CCSB_YI1_lcc,directed = F)
bw_Y2H_uni_lcc<-names(head(sort(bw_Y2H_uni_lcc,decreasing = T),vcount(graph_Y2H_uni_lcc)*0.05))
bw_CCSB_YI1_lcc<-names(head(sort(bw_CCSB_YI1_lcc,decreasing = T),vcount(graph_CCSB_YI1_lcc)*0.05))
#find common with Yeast_deletion data
bw_common_Y2H_uni <- intersect(bw_Y2H_uni_lcc,ORFs_Yeast_proj)
bw_common_CCSB_YI1 <- intersect(bw_CCSB_YI1_lcc,ORFs_Yeast_proj)
bw_common <- intersect(bw_common_Y2H_uni,bw_common_CCSB_YI1)
head9<-c("Y2H_uni_lcc","CCSB_YI1_lcc","Intersection")
tab9<-data.frame(Network=head9,Hubs_common=c(paste(bw_common_Y2H_uni,collapse = ","), paste(bw_common_CCSB_YI1,collapse = ","),paste(bw_common,collapse = ",")))

#run MCL on networks
adj_Y2H_uni_lcc<- as_adj(graph_Y2H_uni_lcc, type = "both", attr = NULL,edges = FALSE, names = TRUE)
Y2H_uni_lcc_mcl <- mcl(x = adj_Y2H_uni_lcc, addLoops=TRUE, ESM = TRUE)
#randomly remove 10% and 25% edges
graph_Y2H_uni_lcc_10<-delete.edges(graph_Y2H_uni_lcc, head(sample(E(graph_Y2H_uni_lcc)),ecount(graph_Y2H_uni_lcc)*0.10))
graph_Y2H_uni_lcc_25<-delete.edges(graph_Y2H_uni_lcc, head(sample(E(graph_Y2H_uni_lcc)),ecount(graph_Y2H_uni_lcc)*0.25))
adj_Y2H_uni_lcc_10<- as_adj(graph_Y2H_uni_lcc_10, type = "both", attr = NULL,edges = FALSE, names = TRUE)
adj_Y2H_uni_lcc_25<- as_adj(graph_Y2H_uni_lcc_25, type = "both", attr = NULL,edges = FALSE, names = TRUE)
Y2H_uni_lcc_mcl_10 <- mcl(x = adj_Y2H_uni_lcc_10, addLoops=TRUE, ESM = TRUE)
Y2H_uni_lcc_mcl_25 <- mcl(x = adj_Y2H_uni_lcc_25, addLoops=TRUE, ESM = TRUE)
head10<-c("Y2H_uni_lcc","Y2H_uni_lcc_10","Y2H_uni_lcc_25")
tab10<-data.frame(Network=head10,Num_clusters=c(Y2H_uni_lcc_mcl$K,Y2H_uni_lcc_mcl_10$K,Y2H_uni_lcc_mcl_25$K),Num_iter=c(Y2H_uni_lcc_mcl$n.iterations,Y2H_uni_lcc_mcl_10$n.iterations,Y2H_uni_lcc_mcl_25$n.iterations))

#for CCSB network
adj_CCSB_YI1_lcc<- as_adj(graph_CCSB_YI1_lcc, type = "both", attr = NULL,edges = FALSE, names = TRUE)
CCSB_YI1_lcc_mcl <- mcl(x = adj_CCSB_YI1_lcc, addLoops=TRUE, ESM = TRUE)
#randomly remove 10% and 25% edges
graph_CCSB_YI1_lcc_10<-delete.edges(graph_CCSB_YI1_lcc, head(sample(E(graph_CCSB_YI1_lcc)),ecount(graph_CCSB_YI1_lcc)*0.10))
graph_CCSB_YI1_lcc_25<-delete.edges(graph_CCSB_YI1_lcc, head(sample(E(graph_CCSB_YI1_lcc)),ecount(graph_CCSB_YI1_lcc)*0.25))
adj_CCSB_YI1_lcc_10<- as_adj(graph_CCSB_YI1_lcc_10, type = "both", attr = NULL,edges = FALSE, names = TRUE)
adj_CCSB_YI1_lcc_25<- as_adj(graph_CCSB_YI1_lcc_25, type = "both", attr = NULL,edges = FALSE, names = TRUE)
CCSB_YI1_lcc_mcl_10 <- mcl(x = adj_CCSB_YI1_lcc_10, addLoops=TRUE, ESM = TRUE)
CCSB_YI1_lcc_mcl_25 <- mcl(x = adj_CCSB_YI1_lcc_25, addLoops=TRUE, ESM = TRUE)
head11<-c("CCSB_YI1_lcc","CCSB_YI1_lcc_10","CCSB_YI1_lcc_25")
tab11<-data.frame(Network=head11,Num_clusters=c(CCSB_YI1_lcc_mcl$K,CCSB_YI1_lcc_mcl_10$K,CCSB_YI1_lcc_mcl_25$K),Num_iter=c(CCSB_YI1_lcc_mcl$n.iterations,CCSB_YI1_lcc_mcl_10$n.iterations,CCSB_YI1_lcc_mcl_25$n.iterations))

#find largest clusters in Y2h components
findlargestcluster<-function(graph,mcldata){
  res<-0
  cres<-V(graph)[ which(mcldata$Cluster %in% c(0))]
  t<-mcldata$K
  for(i in c(0:t))
    #print (i)
    thiscluster<-V(graph)[ which(mcldata$Cluster %in% c(i))]$name
    #print (length(thiscluster))
    if(length(thiscluster) >= length(cres)){
      cres <- thiscluster
      res <- i
    }
  print("Largest cluster num:")
  print (res)
  return(cres)
}

lc_Y2H_uni_lcc_mcl<-findlargestcluster(graph_Y2H_uni_lcc,Y2H_uni_lcc_mcl)
lc_Y2H_uni_lcc_mcl_10<-findlargestcluster(graph_Y2H_uni_lcc_10,Y2H_uni_lcc_mcl_10)
lc_Y2H_uni_lcc_mcl_25<-findlargestcluster(graph_Y2H_uni_lcc_25,Y2H_uni_lcc_mcl_25)

lc_Y2H<-intersection(lc_Y2H_uni_lcc_mcl,lc_Y2H_uni_lcc_mcl_10,lc_Y2H_uni_lcc_mcl_25)

lc_CCSB_YI1_lcc_mcl<-findlargestcluster(graph_CCSB_YI1_lcc,CCSB_YI1_lcc_mcl)
lc_CCSB_YI1_lcc_mcl_10<-findlargestcluster(graph_CCSB_YI1_lcc_10,CCSB_YI1_lcc_mcl_10)
lc_CCSB_YI1_lcc_mcl_25<-findlargestcluster(graph_CCSB_YI1_lcc_25,CCSB_YI1_lcc_mcl_25)
lc_CCSB<-intersection(lc_CCSB_YI1_lcc_mcl,lc_CCSB_YI1_lcc_mcl_10,lc_CCSB_YI1_lcc_mcl_25)

lc_Y2_CCSB<-intersect(as.array(lc_Y2H_uni_lcc_mcl$name),as.array(lc_CCSB_YI1_lcc_mcl$name))
head12<-c("lc_Y2H_uni_lcc_mcl","lc_Y2H_uni_lcc_mcl_10","lc_Y2H_uni_lcc_mcl_25","lc_CCSB_YI1_lcc_mcl","lc_CCSB_YI1_lcc_mcl_10","lc_CCSB_YI1_lcc_mcl_25")
tab12<-data.frame(ClusterIn=head12,Size_clusters=c(length(lc_Y2H_uni_lcc_mcl),length(lc_Y2H_uni_lcc_mcl_10),length(lc_Y2H_uni_lcc_mcl_25),length(lc_CCSB_YI1_lcc_mcl),length(lc_CCSB_YI1_lcc_mcl_10),length(lc_CCSB_YI1_lcc_mcl_25)))


#get all genes in cluster 0
c1_graph_Y2H_uni_lcc <- as.data.frame(V(graph_Y2H_uni_lcc)[ which(Y2H_uni_lcc_mcl$Cluster %in% c(21))]$name)
#c1_graph_Y2H_uni_lcc_10 <- V(graph_Y2H_uni_lcc_10)[ which(Y2H_uni_lcc_mcl_10$Cluster %in% c(0))]
#c1_graph_Y2H_uni_lcc_25 <- V(graph_Y2H_uni_lcc_25)[ which(Y2H_uni_lcc_mcl_25$Cluster %in% c(0))]


library(readr)
GOres <- read_delim("GOres.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


