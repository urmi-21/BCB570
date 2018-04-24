#HW2
library("igraph")
fivestargraph <- graph(combn(1:5, 2), directed = FALSE)
plot.igraph(fivestargraph)

#q5
library(igraph)
#read GraphA
graphA<-read_graph('graphA.txt',format = "edgelist",directed = "false")
vcA<-vcount(graphA)
ecA<-ecount(graphA)
densityA<-graph.density(graphA)
diaA<-diameter(graphA,directed = "false")
radA<-radius(graphA)
#read GraphB
graphB<-read_graph('graphB.txt',format = "edgelist",directed = "false")
vcB<-vcount(graphB)
ecB<-ecount(graphB)
densityB<-graph.density(graphB)
diaB<-diameter(graphB,directed = "false")
radB<-radius(graphB)
#read GraphC
graphC<-read_graph('graphC.txt',format = "edgelist",directed = "false")
vcC<-vcount(graphC)
ecC<-ecount(graphC)
densityC<-graph.density(graphC)
diaC<-diameter(graphC,directed = "false")
radC<-radius(graphC)
header<-c("#Vertices","#Edges","Density")
tab1<-data.frame(Feature=header,GraphA=c(vcA,ecA,densityA),GraphB=c(vcB,ecB,densityB),GraphC=c(vcC,ecC,densityC))
header2<-c("Radius","Diameter")
tab2<-data.frame(Feature=header2,GraphA=c(radA,diaA),GraphB=c(radB,diaB),GraphC=c(radC,diaC))

#compute clustering coeff
ccA<-transitivity(graphA, type = "average")
ccA_g<-transitivity(graphA, type = "global")
ccB<-transitivity(graphB, type = "average")
ccB_g<-transitivity(graphB, type = "global")
ccC<-transitivity(graphC, type = "average")
ccC_g<-transitivity(graphC, type = "global")
header3<-c(" clustering coefficient(avg)"," clustering coefficient(global)")
tab3<-data.frame(Feature=header3,GraphA=c(ccA,ccA_g),GraphB=c(ccB,ccB_g),GraphC=c(ccC,ccC_g))

#compute avg shortest dist
mdA<-mean_distance(graphA,directed = "false")
mdB<-mean_distance(graphB,directed = "false")
mdC<-mean_distance(graphC,directed = "false")
header4<-c("Average shortest path")
tab4<-data.frame(Feature=header4,GraphA=c(mdA),GraphB=c(mdB),GraphC=c(mdC))

#get central nodes
topDA<-head(order(igraph::degree(graphA),decreasing = TRUE),5)
topDB<-head(order(igraph::degree(graphB),decreasing = TRUE),5)
topDC<-head(order(igraph::degree(graphC),decreasing = TRUE),5)
header5<-c("Vertices with highest degree","Vertices with highest betweenness", "Common vertices")
topBA<-(head(order(igraph::betweenness(graphA),decreasing = TRUE),5))
topBB<-(head(order(igraph::betweenness(graphB),decreasing = TRUE),5))
topBC<-(head(order(igraph::betweenness(graphC),decreasing = TRUE),5))
commA<-intersect(topDA,topBA)
commB<-intersect(topDB,topBC)
commC<-intersect(topDC,topBC)
tab5<-data.frame(Feature=header5,GraphA=c(toString(topDA),toString(topBA),toString(commA)),GraphB=c(toString(topDB),toString(topBB),toString(commB)),GraphC=c(toString(topDC),toString(topBC),toString(commC)))

#fit power law
degA<-degree(graphA,mode = "in")
fit1 <- fit_power_law(degA+1, 10)
fit2 <- fit_power_law(degA+1, 10, implementation="R.mle")
degB<-degree(graphB,mode = "in")
degC<-degree(graphC,mode = "in")

fit1 <- fit_power_law(degC+1)



#7
#plot deg
plot(degree.distribution(graphA),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "Degree Distribution GraphA")
plot(degree.distribution(graphB),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "Degree Distribution GraphB")
plot(degree.distribution(graphC),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "Degree Distribution GraphC")

#fit powerlaw
#find xmin by plotting cdf on log-log scale
plot(cumsum(degree.distribution(graphA)),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "Cumulative Degree Distribution GraphA")
plot(cumsum(degree.distribution(graphB)),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "Cumulative Degree Distribution GraphB")
plot(cumsum(degree.distribution(graphC)),xlab = "log(k)",ylab = "log(P(k))",log="xy",type='o',main = "Cumulative Degree Distribution GraphC")

fitA <- power.law.fit(degree(graphA,mode="all"),6,implementation="plfit") 
fitB <- power.law.fit(degree(graphB,mode="all"),25,implementation="plfit") 
fitC <- power.law.fit(degree(graphC,mode="all"),5,implementation="plfit") 
fitA <- power.law.fit(degree(graphA,mode="all"),implementation="plfit") 
fitB <- power.law.fit(degree(graphB,mode="all"),implementation="plfit") 
fitC <- power.law.fit(degree(graphC,mode="all"),implementation="plfit") 
head6<-c("GraphA","GraphB","GraphC")
tab6<-data.frame(Graph=head6,Alpha=c(fitA$alpha,fitB$alpha,fitC$alpha),X_min=c(fitA$xmin,fitB$xmin,fitC$xmin),KS_stat=c(fitA$KS.stat,fitB$KS.stat,fitC$KS.stat),P_Value=c(fitA$KS.p,fitB$KS.p,fitC$KS.p),logLik=c(fitA$logLik,fitB$logLik,fitC$logLik))
