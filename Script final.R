setwd("C:/Users/asus/Desktop/Resultados finais")

library(vegan)
library(ade4)
library(ggplot2)
library(ellipse)
library(cluster)
library(FD)
library(mclust)

#biological data analysis
#upload of data of biological data with transformation log(data=abundancia)
taxa<-read.table("Bio 1 com log (transformacao).txt",sep="\t",dec=".",header=T)
View(taxa) #checking data frame
taxa$Sites->rownames(taxa)
id<-taxa[,1] #add site to an object
taxa<-taxa[,-1]#remove sites from the original dataset to do an exploratory analysis, calculating
#distances between sites
View(taxa)#checking dataframe
View(id)#checking object id (sites)

#calculating bray-curtis dissimilarity to taxa
#we use bray-curtis dissimillarity index because data has a lot of zeros (abundance data)
bio.d<-vegdist(taxa,method="bray")
bio.d #dissimilarity between sites 
#0 means that sites are the same 
#1 means that sites are diferent

#building a cluster
#k-number of groups selected
#using ward method, the most strict method of classification
hclust(bio.d,method="ward.D")->bio.ward

#statistical inference to choose the number of clusters
#Identifying clusters through mean values of the original variables
library(plyr)
bio.gr <- as.factor(cutree(bio.ward, 6))
ddply(taxa,.(bio.gr),function(x){round(colMeans(x),1)})
###
#Determining optimal number of groups for kmeans clustering

wss.bio1<- (nrow(taxa)-1)*sum(apply(taxa,2,var))
for (i in 2:15) wss.bio1[i] <- sum(kmeans(taxa, centers=i)$withinss)
#k-means does not provide graphic output
#sum of squared error (SSE) scree plot

par(mfrow=c(1,1))
plot(1:15, wss.bio1, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
abline(v=3, lty=2, col="blue")
abline(v=6, lty=2, col="violet")

# Non-hierarchical clustering (centers=number of groups)
#we shouldn't choose a lot of groups, because it's to mutch information to handle
kmeans (bio.d, centers = 6)->env.kmeans.bio1
env.kmeans.bio1$cluster
#Selecting number of trees
bio.gr <- as.factor(cutree(bio.ward, 6))



#plot dendogram
#discriminatory cluster of 1 dimension
#we choose k=6 in the number of groups, because of the communitie knowledge of the sites in the identification
plot(bio.ward)
rect.hclust(bio.ward, k = 6, border = 2:7) # add rectangle


#
#Detrended Correspondence Analysis

bio.dca<-decorana(taxa)
bio.dca
bio.dca$evals->eigenvalues
#downheighting of rare species iweigh=1
bio.dca2<-decorana(taxa,iweigh = 1)
bio.dca2#worts explained variance


#with decorana values from dca
par(mfrow=c(1,2))
barplot(eigenvalues, xlab="Principal Components",ylab="Explained variance",las=1)
barplot(cumsum(eigenvalues), xlab="Principal Components",ylab="Explained variance",las=1)
cumsum(eigenvalues)#variance explained by axis
sum(eigenvalues)#variance explained by the 4 dca axis
eigenvalues*100# % of variance explained by axis


###plot of dca
par(mfrow=c(1,1))
plot(bio.dca, type = "n",xlim =c(-0.6,2.1),ylim = c(-1,1.5),main="DCA")
points(bio.dca, display = "sites",cex = 0.6, pch=21, col="blue", bg="blue")
text(bio.dca, display =  "sites", labels = id, cex=0.8, col="black")
#text(bio.dca, display = "spec", cex= 0.5, col="red")
mtext("a) DCA(38,74%)", line = 2.0, adj = 0.01, cex = 2, font = 2)

##getting symbols and colors for graphical output
sym<-c(rep(16,4),rep(17,3)) # symbols
colour<-c("green3","gold","steelblue3","firebrick1","cyan","magenta3") # colours
s.class(bio.dca$rproj,bio.gr,pch=sym[bio.gr],cpoint=1,col=colour,cellipse=2,cstar=1,axesell = FALSE,
        sub = "DCA", possub = "bottomright",cgrid = 2, csub = 2)

###
####
####
###
####hydrological data analysis
#upload of hydrological data
#hydrology#goal:obtain 2 hydrological variables hydro1  and hydro2 from pca hydrology
hydrology<-read.table("Hydrology Summer 2020.txt",sep="\t",dec=".",header=T)
View(hydrology)
hydrology$site->rownames(hydrology)
hydrology<-hydrology[,-c(2,3)] #remmoving year and campaign
View(hydrology)
hydrology->labels.hydro
id.amb<-hydrology[,1]
View(id.amb)
hydrology<-hydrology[,-1]
View(hydrology)

# Variable transformation
q<-sapply(hydrology,class)=="numeric" | sapply(hydrology,class)=="integer"# selecting quantitative variables

par(mfrow=c(2,3))
for (i in which(q==T)) hist(hydrology[,i], main=names(hydrology)[i])

###data transformation
###performed for all variables
####min.D
par(mfrow=c(1,3))
hist(hydrology$min.D)#raw data with bad normal distribution
hist(sqrt(hydrology$min.D))#use this transformation#better performed
hist(log(hydrology$min.D))#bad transformation

hydrology$min.D<-sqrt(hydrology$min.D) ##best result selected as variable

###max.D
par(mfrow=c(1,3))
hist(hydrology$max.D)#raw data with bad normal distribution
hist(sqrt(hydrology$max.D))#use this transformation#better performed
hist(log(hydrology$max.D))#bad transformation

hydrology$max.D<-sqrt(hydrology$max.D) ##best result selected as variable

####STD.D
par(mfrow=c(1,3))
hist(hydrology$STD.D)#raw data with bad normal distribution
hist(sqrt(hydrology$STD.D))#use this transformation#better performed
hist(log(hydrology$STD.D))#bad transformation

hydrology$STD.D<-sqrt(hydrology$STD.D) ##best result selected as variable

####CV.D
par(mfrow=c(1,3))
hist(hydrology$CV.D)#raw data with bad normal distribution
hist(sqrt(hydrology$CV.D))#use this transformation#better performed
hist(log(hydrology$CV.D))#bad transformation

hydrology$CV.D<-sqrt(hydrology$CV.D) ##best result selected as variable


####min.V
par(mfrow=c(1,3))
hist(hydrology$min.V)#raw data with bad normal distribution
hist(sqrt(hydrology$min.V))#use this transformation#better performed
hist(log(hydrology$min.V+0.0001))#bad transformation

hydrology$min.V<-log(hydrology$min.V+0.0001) ##best result selected as variable
#we use 0.0001 as default value so we can a play a distance matrix
#that we won't be able if we only have (log min.V) because it gives us a (-inf) interval
#that makes impossible to calculate dissimilarity matrix

####max.V
par(mfrow=c(1,3))
hist(hydrology$max.V)#raw data with bad normal distribution
hist(sqrt(hydrology$max.V))#use this transformation#better performed
hist(log(hydrology$max.V))#bad transformation

hydrology$max.V<-sqrt(hydrology$max.V) ##best result selected as variable


####STD.V
par(mfrow=c(1,3))
hist(hydrology$STD.V)#raw data with bad normal distribution
hist(sqrt(hydrology$STD.V))#use this transformation#better performed
hist(log(hydrology$STD.V))#bad transformation

hydrology$STD.V<-sqrt(hydrology$STD.V) ##best result selected as variable


####CV.V
par(mfrow=c(1,3))
hist(hydrology$CV.V)#raw data with bad normal distribution
hist(sqrt(hydrology$CV.V))#use this transformation#better performed
hist(log(hydrology$CV.V))#bad transformation

hydrology$CV.V<-sqrt(hydrology$CV.V) ##best result selected as variable

####mean.Depth
par(mfrow=c(1,3))
hist(hydrology$mean.Depth)#raw data with bad normal distribution
hist(sqrt(hydrology$mean.Depth))#use this transformation#better performed
hist(log(hydrology$mean.Depth))#bad transformation

hydrology$mean.Depth<-sqrt(hydrology$mean.Depth) ##best result selected as variable

####mean.Velocity
par(mfrow=c(1,3))
hist(hydrology$mean.Velocity)#raw data with bad normal distribution
hist(sqrt(hydrology$mean.Velocity))#use this transformation#better performed
hist(log(hydrology$mean.Velocity))#bad transformation

hydrology$mean.Velocity<-sqrt(hydrology$mean.Velocity) ##best result selected as variable

####Width
par(mfrow=c(1,3))
hist(hydrology$Width)#raw data with bad normal distribution
hist(sqrt(hydrology$Width))#use this transformation#better performed
hist(log(hydrology$Width))#bad transformation

hydrology$Width<-log(hydrology$Width) ##best result selected as variable

####Discharge
par(mfrow=c(1,3))
hist(hydrology$Discharge)#raw data with bad normal distribution
hist(sqrt(hydrology$Discharge))#use this transformation#better performed
hist(log(hydrology$Discharge))#bad transformation

hydrology$Discharge<-log(hydrology$Discharge) ##best result selected as variable

###
# Variable transformation
b<-sapply(hydrology,class)=="numeric" | sapply(hydrology,class)=="integer"# selecting quantitative variables

par(mfrow=c(2,3))
for (i in which(b==T)) hist(hydrology[,i], main=names(hydrology)[i])
par(mfrow=c(1,1))

#Calculating the Euclidean distance matrix
#scaling the variables

scale(hydrology)->hydrology.st
View(hydrology)

##calculating euclidean distance ##only min.V has a lot of 0's
vegdist(hydrology.st,method="euclidean")->hydrology.d

###Flexible clustering needs a par.method. 0.7 makes the cluster similar to Ward's method.
hclust(hydrology.d,method="ward.D")->hydrology.ward


#Identifying clusters through mean values of the original variables
library(plyr)
ddply(hydrology,.(hydrology.ward.gr),function(x){round(colMeans(x),1)})
###
#Determining optimal number of groups for kmeans clustering

wss.hydro <- (nrow(hydrology)-1)*sum(apply(hydrology,2,var))
for (i in 2:15) wss.hydro[i] <- sum(kmeans(hydrology,
                                           centers=i)$withinss)
#k-means does not provide graphic output
# Non-hierarchical clustering (centers=number of groups)
kmeans (hydrology.d, centers = 5)->env.kmeans.hydro
env.kmeans.hydro$cluster

#sum of squared error (SSE) scree plot

par(mfrow=c(1,1))
plot(1:15, wss.hydro, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
abline(v=4, lty=2, col="blue")
abline(v=6, lty=2, col="violet")

#Selecting number of trees
#we choose k=5 because of the statistical suport
hydrology.ward.gr <- as.factor(cutree(hydrology.ward, 5))


##ploting clustering
plot(hydrology.ward,hang=-1,main="Ward's method",labels=id.amb)
rect.hclust(hydrology.ward, k = 5, border = 2:7) # add rectangle


###PCA
dudi.pca(hydrology, center=T, scale=T, scann=F, nf=5)->hydro.pca

#How many axes do we have to select?
par(mfrow = c(1,1))
screeplot(princomp(hydrology,scale=T),bstick=T,type="lines")
#we choose 2 dimension, vecause of the broken stick test

#env.pca$eig shows the eigenvalues, it is the importance of each axis
hydro.pca$eig

names(hydro.pca$eig)<-1:length(hydro.pca$eig) # To assign the names of the PC to each eigenvalue
hydro.pca$eig #Note the numbers


#explained variance by axis
barplot(hydro.pca$eig, xlab="Principal Components",ylab="Eigenvalues",las=1)

#How much variance explains each axis?

(hydro.pca$eig*100)/sum(hydro.pca$eig)->var.hydro.pca
barplot(var.hydro.pca, xlab="Principal Components",ylab="Explained variance",las=1)

#How much cummulative variance explains each axis?

cumsum(var.hydro.pca)
barplot(cumsum(var.hydro.pca), xlab="Principal Components",ylab="Explained variance",las=1)


par(mfrow=c(1,2))
barplot(var.hydro.pca, xlab="Principal Components",ylab="Explained variance",las=1)
barplot(cumsum(var.hydro.pca), xlab="Principal Components",ylab="Explained variance",las=1)

#$c1 shows the PCA loadings (the weight of each variable in constructing each retained principal component)

hydro.pca$c1

#$co shows the variable coordinates (arrow coordinates). In the specific case of normalised PCA, these coordinates represent correlation with axes and thus should vary between -1 and 1.

hydro.pca$co

#Graphical representation of PCA Show the pca.biplot (two first axes)

par(mfrow=c(1,1))
#General PCA plot

biplot(hydro.pca)

#Plotting the dendrogram
s.arrow(hydro.pca$c1, lab = names(hydro.pca$tab))

#Plotting cases

s.label(hydro.pca$li, sub = "Environmental PCA", csub = 1.5, possub = "topleft")

##ploting pca with groups
##getting symbols and colors for graphical output
sym<-c(rep(16,4),rep(17,3)) # symbols
colour2<-c("green3","steelblue3","magenta3","firebrick1","cyan") # colours

s.class(hydro.pca$li,hydrology.ward.gr,pch=sym[hydrology.ward.gr],cpoint=1,col=colour2,cellipse=2,cstar=1,axesell = FALSE,
        sub = "PCA", possub = "bottomright",cgrid = 2, csub = 2)

#pca loadings->points ->coordinates in the multivariate space
hydro.pca$li$Axis1->hydro1 #related with min.D,max.D,STD.D, max.V,STD.V,mean.Depth,Width,Discharge
hydro.pca$li$Axis2->hydro2 #related with max.V, mean.Velocity,STD.V
hydrology$mean.Velocity->Mean_velocity
hydrology$Discharge->Discharge


##colinearity (statistical term) and coorelation (assossiation between predictor and variable)
#checking for coorellation
library(GGally)
ggpairs(hydrology)
round(cor(hydrology),3) #discharge correlated with depth, width,velocity because flow (discharge)= A(depth x width) * velocity
#this means that could have autocorrelation
#this means that hydro1 is a good hydrological variable
##


###
####
####
###
####environmental data analysis
#upload of environmental data
env.data<-read.table("var ambientais.txt",sep="\t",dec=".",header=T)
View(env.data)


#add hydro1 and hydro2 to env data
hydro.stress<-cbind(Mean_velocity,Discharge)
View(hydro.stress)

env.data<-cbind.data.frame(env.data,hydro.stress)
View(env.data)
env.data$code->rownames(env.data)
id.env<-env.data[,1]
View(id.env)
env.data<-env.data[,-1]
View(env.data)


#removing lat and lon
env.data<-env.data[,-c(1,2)]
View(env.data)


# Variable transformation
a<-sapply(env.data,class)=="numeric" | sapply(env.data,class)=="integer"# selecting quantitative variables

par(mfrow=c(2,3))
for (i in which(a==T)) hist(env.data[,i], main=names(env.data)[i])

###data transformation
###performed for all variables,and choose the one tha performed better transformations
###altitude
par(mfrow=c(1,3))
hist(env.data$Slope)#raw data with bad normal distribution
hist(sqrt(env.data$Slope))#use this transformation#better performed
hist(log(env.data$Slope+0.0001))#bad transformation

env.data$Slope<-sqrt(env.data$Slope) ##best result selected as variable



###altitude
par(mfrow=c(1,3))
hist(env.data$Altitude)#raw data with bad normal distribution
hist(sqrt(env.data$Altitude))#use this transformation#better performed
hist(log(env.data$Altitude+0.0001))#bad transformation

env.data$Altitude<-sqrt(env.data$Altitude) ##best result selected as variable


###QBR
par(mfrow=c(1,3))
hist(env.data$qbr)#raw data with bad normal distribution
hist(sqrt(env.data$qbr))#use this transformation#better performed
hist(log(env.data$qbr+0.0001))#bad transformation

#better without transformation

###ihf
par(mfrow=c(1,3))
hist(env.data$ihf)#raw data with bad normal distribution
hist(sqrt(env.data$ihf))#use this transformation#better performed
hist(log(env.data$ihf+0.0001))#bad transformation

###better without transformation

###shadow
par(mfrow=c(1,3))
hist(env.data$shadow)#raw data with bad normal distribution
hist(sqrt(env.data$shadow))#use this transformation#better performed
hist(log(env.data$shadow+0.0001))#bad transformation

###better without transformation


###COD
par(mfrow=c(1,3))
hist(env.data$COD)#raw data with bad normal distribution
hist(sqrt(env.data$COD))#use this transformation#better performed
hist(log(env.data$COD+0.0001))#bad transformation

env.data$COD<-log(env.data$COD+0.0001) ##best result selected as variable

###

###DIN
par(mfrow=c(1,3))
hist(env.data$DIN)#raw data with bad normal distribution
hist(sqrt(env.data$DIN))#use this transformation#better performed
hist(log(env.data$DIN+0.0001))#bad transformation

env.data$DIN<-log(env.data$DIN+0.0001) ##best result selected as variable

###P.PO4
par(mfrow=c(1,3))
hist(env.data$P.PO4)#raw data with bad normal distribution
hist(sqrt(env.data$P.PO4))#use this transformation#better performed
hist(log(env.data$P.PO4+0.0001))#bad transformation

env.data$P.PO4<-sqrt(env.data$P.PO4) ##best result selected as variable


#####ph
par(mfrow=c(1,3))
hist(env.data$ph)#raw data with bad normal distribution
hist(sqrt(env.data$ph))#use this transformation#better performed
hist(log(env.data$ph+0.0001))#bad transformation

#better without transformation

##
#####cond
par(mfrow=c(1,3))
hist(env.data$cond)#raw data with bad normal distribution
hist(sqrt(env.data$cond))#use this transformation#better performed
hist(log(env.data$cond+0.0001))#bad transformation

env.data$cond<-log(env.data$cond+0.0001) ##best result selected as variable



############DOmin
par(mfrow=c(1,3))
hist(env.data$DOmin)#raw data with bad normal distribution
hist(sqrt(env.data$DOmin))#use this transformation#better performed
hist(log(env.data$DOmin+0.0001))#bad transformation

#better without transformation

############Tmax
par(mfrow=c(1,3))
hist(env.data$Tmax)#raw data with bad normal distribution
hist(sqrt(env.data$Tmax))#use this transformation#better performed
hist(log(env.data$Tmax+0.0001))#bad transformation

env.data$Tmax<-sqrt(env.data$Tmax) ##best result selected as variable


########X.Artificial.500m.
par(mfrow=c(1,3))
hist(env.data$X.Artificial.500m.)#raw data with bad normal distribution
hist(sqrt(env.data$X.Artificial.500m.))#use this transformation#better performed
hist(log(env.data$X.Artificial.500m.+0.0001))#bad transformation

#better without transformation


#########mean_light
par(mfrow=c(1,3))
hist(env.data$mean_light)#raw data with bad normal distribution
hist(sqrt(env.data$mean_light))#use this transformation#better performed
hist(log(env.data$mean_light+0.0001))#bad transformation

env.data$mean_light<-log(env.data$mean_light+0.0001)##best result selected as variable

####
########X.Agriculture.500m.
par(mfrow=c(1,3))
hist(env.data$X.Agriculture.500m.)#raw data with bad normal distribution
hist(sqrt(env.data$X.Agriculture.500m.))#use this transformation#better performed
hist(log(env.data$X.Agriculture.500m.+0.0001))#bad transformation

env.data$X.Agriculture.500m.<-sqrt(env.data$X.Agriculture.500m.)##best result selected as variable

#####
########X.Pasture.500m.
par(mfrow=c(1,3))
hist(env.data$X.Pasture.500m.)#raw data with bad normal distribution
hist(sqrt(env.data$X.Pasture.500m.))#use this transformation#better performed
hist(log(env.data$X.Pasture.500m.+0.0001))#bad transformation

#better without transformation

########X.Natural.500m.
par(mfrow=c(1,3))
hist(env.data$X.Natural.500m.)#raw data with bad normal distribution
hist(sqrt(env.data$X.Natural.500m.))#use this transformation#better performed
hist(log(env.data$X.Natural.500m.+0.0001))#bad transformation

env.data$X.Natural.500m.<-sqrt(env.data$X.Natural.500m.)##best result selected as variable


######
# Variable transformation
c<-sapply(env.data,class)=="numeric" | sapply(env.data,class)=="integer"# selecting quantitative variables

par(mfrow=c(2,3))
for (i in which(c==T)) hist(env.data[,i], main=names(env.data)[i])
par(mfrow=c(1,1))

# Saving all environmental stressors
write.table(env.data,"environmental_stressors.txt",sep="\t",dec=".")
View(env.data)

#Calculating the Euclidean distance matrix
#scaling the variables
scale(env.data)->env.data.st

#checkin dataframes
View(env.data.st)
View(env.data)

##calculating euclidean distance ##only min.V has a lot of 0's
vegdist(env.data.st,method="euclidean")->env.data.d

###Flexible clustering needs a par.method. 0.7 makes the cluster similar to Ward's method.
hclust(env.data.d,method="ward.D")->env.data.ward

#Selecting number of trees
env.data.ward.gr <- cutree(env.data.ward, 5)

#Identifying clusters through mean values of the original variables
library(plyr)
ddply(env.data,.(env.data.ward.gr),function(x){round(colMeans(x),1)})
###
#Determining optimal number of groups for kmeans clustering

wss.env <- (nrow(env.data)-1)*sum(apply(env.data,2,var))
for (i in 2:15) wss.env[i] <- sum(kmeans(env.data,
                                         centers=i)$withinss)
#k-means does not provide graphic output
# Non-hierarchical clustering (centers=number of groups)
kmeans (env.data.d, centers = 5)->env.kmeans.env
env.kmeans.env$cluster

#sum of squared error (SSE) scree plot

par(mfrow=c(1,1))
plot(1:15, wss.env, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
abline(v=4, lty=2, col="blue")
abline(v=6, lty=2, col="violet")

#Selecting number of trees
##based on statistical and knowledge of the site, we decide to cut de tree with 6 groups
env.data.ward.gr <- as.factor(cutree(env.data.ward, 6))

##ploting clustering
plot(env.data.ward,hang=-1,main="Ward's method",labels=id.env)
rect.hclust(env.data.ward, k = 6, border = 2:7) # add rectangle


###PCA
dudi.pca(env.data.st, center=T, scale=T, scann=F, nf=5)->env.pca

#How many axes do we have to select?
par(mfrow = c(1,1))
screeplot(princomp(env.data.st,scale=T),bstick=T,type="lines")

#env.pca$eig shows the eigenvalues, it is the importance of each axis
env.pca$eig

names(env.pca$eig)<-1:length(env.pca$eig) # To assign the names of the PC to each eigenvalue
env.pca$eig #Note the numbers


#explained variance by axis
barplot(env.pca$eig, xlab="Principal Components",ylab="Eigenvalues",las=1)

#How much variance explains each axis?

(env.pca$eig*100)/sum(env.pca$eig)->var.env.pca

barplot(var.env.pca, xlab="Principal Components",ylab="Explained variance",las=1)

#How much cummulative variance explains each axis?

cumsum(var.env.pca)
barplot(cumsum(var.env.pca), xlab="Principal Components",ylab="Explained variance",las=1)


par(mfrow=c(1,2))
barplot(var.env.pca, xlab="Principal Components",ylab="Explained variance",las=1)
barplot(cumsum(var.env.pca), xlab="Principal Components",ylab="Explained variance",las=1)

#$c1 shows the PCA loadings (the weight of each variable in constructing each retained principal component)

env.pca$c1

#$co shows the variable coordinates (arrow coordinates). In the specific case of normalised PCA, these coordinates represent correlation with axes and thus should vary between -1 and 1.

env.pca$co

#Graphical representation of PCA Show the pca.biplot (two first axes)

par(mfrow=c(1,1))
#General PCA plot

biplot(env.pca)

#Plotting the dendrogram
s.arrow(env.pca$c1, lab = names(env.pca$tab))

#Plotting cases
s.label(env.pca$li, sub = "Environmental PCA", csub = 1.5, possub = "topleft")

##ploting pca with groups
##getting symbols and colors for graphical output
sym<-c(rep(16,4),rep(17,3)) # symbols
colour3<-c("gold","green3","magenta3","firebrick1","cyan","steelblue3") # colours

s.class(env.pca$li,env.data.ward.gr,pch=sym[env.data.ward.gr],cpoint=1,col=colour3,cellipse=2,cstar=1,axesell = FALSE,
        sub = "PCA", possub = "bottomright",cgrid = 2, csub = 2)

env.pca$li$Axis1->stress.gen1
env.pca$li$Axis2->stress.gen2

###
##########################
#   Biodiversity metrics
##########################
####Taxonomic diversity metrics

taxa.diversity<-read.table("Bio 1 sem log (transformacao).txt",sep="\t",dec=".",header=T)
View(taxa.diversity) #checking data frame
taxa.diversity$Sites->rownames(taxa.diversity)
id.diversity<-taxa.diversity[,1] #add site to an object
taxa.diversity<-taxa.diversity[,-1]#remove sites from the original dataset to do an exploratory analysis, calculating
#distances between sites
View(taxa.diversity)#checking dataframe
View(id.diversity) #checking object id (sites)

richness<-specnumber(taxa.diversity) #species richness
shannon<-diversity(taxa.diversity, index = "shannon", MARGIN = 1, base = exp(1))
#Hill's diversity
#Hill number, q = 0 (default) to get species richness, q = 1 to get shannon entropy,
#q = 2 will give inverse Simpson
library(hillR)
Hill.shannon<-hill_taxa(taxa.diversity, q = 1, MARGIN = 1, base = exp(1))
##colinearity (statistical term) and coorelation (assossiation between predictor and variable)
#checking for coorellation

####
####Functional diversity
###loading matrices
bio1.traits<-read.table("Matriz_tracos_murcia_final.txt",sep="\t",dec=".",header=T)
View(bio1.traits)
trait.matrix<-read.table("Trait_matrix_Murriaetall_edited.txt",sep="\t",dec=".",header=T)
View(trait.matrix)
trait.code<-read.table("Traits code.txt",sep="\t",dec=".",header=T)
View(trait.code)

#preparing matrices

# Assigning rownames
trait.matrix$Genus->rownames(trait.matrix)

#arranging biological matrix
bio1.traits$R�tulos.de.Linha->rownames(bio1.traits)
bio1.traits<-bio1.traits[,-1]
View(bio1.traits)

# arranging traits
trait.matrix<-trait.matrix[,-c(1:4)]
View(trait.matrix)

## Trait subset for taxa
trait_sub<-trait.matrix[intersect(colnames(bio1.traits),rownames(trait.matrix)),]
View(trait_sub)

# Checking if all taxa considered in "taxa" matrix are present in the trait matrix
# FALSE means we have no problems. Otherwise we need to check taxon names
any(rownames(trait_sub)==colnames(bio1.traits))==F 

# Checking matrix dimensions
dim(trait_sub)
dim(bio1.traits)
dim(trait.matrix)


# Duplicating trait matrix to preserve original traits
tr<-trait_sub
tr2<-tr

# Loading additional functions
source("0_FD_functions.R")
source("0_quality_funct_space_fromdist.R")
source("1_FD_functions_Mouillot et all.R")
source("FD_functions.R")

# arranging traits
traits.blo<-c(7,2,3,4,8,4,5,4,8,10)
tr <- prep.fuzzy(tr, traits.blo) 

# Combining the traits
tr.ktab<-ktab.list.df(list(tr))
tr.dist <- dist.ktab(tr.ktab, "F") # fuzzy-coding adapted Gower distance
###################################################################################################
# Functional space and groups
################################################################################################### 

# calculate species x species Gower dissimilarity matrix based on functional traits
gowdis(tr2)->tr.dist2
cor(tr.dist,tr.dist2)
cor(tr.dist,tr.dist2)^2

# Estimating the optimum number of dimensions
qual_fs<-quality_funct_space_fromdist(tr.dist, nbdim=55)
qual_fs$meanSD<0.01 # 2D seems to be an appropiate number of dimensions
qual_fs$meanSD
#10 Dimension are aproprietaded

# Functional space
dudi.pco(tr.dist,scannf = F,nf=10)->tr.pco

# checking for negative eigenvalues
length(which(tr.pco$eig<0)) # No negative eigenvalues licenced us to use Ward clustering method

# Axis importance (explained variance)
round(tr.pco$eig[1:50]/sum(tr.pco$eig),2)
sum(tr.pco$eig[1:10]/sum(tr.pco$eig)) #variance explained with 10 D

# Correlation between axes and original quantiative variables
round(cor(tr[,which(sapply(tr,is.numeric))],tr.pco$li, use="pairwise.complete.obs"),2)

##########################################
# Functional groups
##########################################

# Classifying species into functional groups
tr.clust <- hclust(tr.dist, method = "ward.D")

# Plotting cluster
par(mfrow=c(1,1))
plot(as.phylo(tr.clust), cex = 0.9, label.offset= 0.05)

# Five groups were chosen
cut.g <- 5
f.gr <- cutree(tr.clust, k = cut.g)

# Median features for each group
ddply(data.frame(f.gr, tr), .(f.gr), function(x) apply(x[,-1], 2, median))

# Exploring the functional groups
for (i in 1:cut.g) print(summary(tr[which(f.gr==i),]))

# Exploring which functional grouping reduce between-groups differences

# Creating a variable to store r2 results
r2.res <- rep(NA, 14)

for(i in 2:15) {
  
  cut.g <- i
  f.gr <- factor(cutree(tr.clust, k = cut.g))
  
  r2.res[i-1] <-  adonis2(tr.dist~f.gr)$R2[1]
  
}

par(mfrow=c(1,2))

# Grouping results
plot(2:15, r2.res, xlab="Number of groups", ylab="Goodness-of-fit")
lines(loess(r2.res~c(2:15)), col="blue")
abline(v=6, lty=2, col="grey")
abline(h=r2.res[5], lty=2, col="grey")

r2.dif=rep(NA, 13)
names(r2.dif)  <- 3:15

for (i in 1: 13) r2.dif[i]  <- r2.res[i+1]-r2.res[i]

# Barplot
barplot(r2.dif, ylab="Marginal increment in goodness-of-fit", 
        xlab="Number of groups",col=c(rep("blue", 4), rep("black",9)))
abline(v=5, lty=2, col="grey")
##########################################
# Plotting functional cluster and space
##########################################

# Six groups were chosen
cut.g <- 6
f.gr <- cutree(tr.clust, k = cut.g)

# Median features for each group
ddply(data.frame(f.gr, tr), .(f.gr), function(x) apply(x[,-1], 2, median))

# Setting group colours
fgr.col<-c("#FF0000FF", "gold3", "#00FF66FF", "#0066FFFF", "#CC00FFFF", "orange", "dark green")

par(mfrow=c(1,2))
plot(as.phylo(tr.clust), cex = 0.9, label.offset= 0.05, tip.color = fgr.col[f.gr])
s.class(tr.pco$li[,1:2],fac = as.factor(f.gr), col=fgr.col,clabel=1.5, cpoint=1.5,cellipse=0)

#####################################################
### Biodiversity metrics                    
#####################################################

# Community-wide functinal metrics

# Number of dimensions (remember your qual_fs values!)
k=10

# Functional richness
fric_3d(bio1.traits, tr.pco$li, k, prec="Qt")->FRic 

# Functional dispersion
fdisp_k(tr.dist,bio1.traits, k)$FDis->FDis 

# Saving generated metrics
diversity<-data.frame(richness,
                      shannon,
                      Hill.shannon,
                      FRic,
                      FDis)
View(diversity)
# Saving datasets
write.table(diversity,"taxonomic&functional_diversity.txt",sep="\t",dec=".")

# Correlation among functional metrics
as.dist(round(cor(diversity, use="pairwise.complete.obs"),2))

# Plotting scatter plots
par(mfrow=c(2,2), cex=1.25, mar=c(4,4,4,1))

var.sel<-c("richness", "shannon", "Hill.shannon", "FRic", "FDis")

for (i in var.sel) {
  plot(diversity[,i]~richness, xlab="Richness", ylab=i)
  mod <- lm(diversity[,i]~richness)
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
  
}

# PCA on functional metrics: exploring their redundancy

# Replacing NAs by 0 to run PCA (if needed)
diversity$FRic[which(is.na(diversity$FRic)==T)]<-0
dudi.pca(diversity, center = TRUE, scale = TRUE, scannf = F,nf=4)->metric.pca

# PCA Axis importance (explained variance)
round(metric.pca$eig[1:4]/sum(metric.pca$eig),2)
sum(metric.pca$eig[1:4]/sum(metric.pca$eig))

# Interpreting axes
round(cor(diversity,metric.pca$li),2)

# Plotting results
par(mfrow=c(1,1))
s.arrow(metric.pca$c1)

corre<-cbind(env.data,diversity)
View(corre)
cor(corre, method = "spearman")->y
# Saving datasets
write.table(y,"correlation_between_biodiversity_metrics_and_environmental_variables_spearman.txt",sep="\t",dec=".")


####original matrix
####data with transformation

####
##escolha das variaveis finais ponderando de 3 formas
###1Âº correlacoes entre variaveis
####2Âºdirecoes e aparente relacao no espaco funcional (estarem na mesma direcao nao e significado de correlacao)
####3Âºimportancia e impacto na comunidade de macroinvertebrados
####4Âº idealmente selecionar 1 variavel para cada tipo de stressor esperado (natural,vegetacao riparia,hidrologico,eutrofizacao, fisico-quimico,uso do solo,natural)


##colinearity (statistical term) and coorelation (assossiation between predictor and variable)
#checking for coorellation

####oriinal matrix
####data with transformation
library(GGally)
ggpairs(env.data)
round(cor(env.data), 3)->x

View(as.data.frame(round(cor(env.data.st), 3)))


# Saving datasets of correlation
write.table(x,"correlation.txt",sep="\t",dec=".")
View(x)

###env.data; environmental variables selected=stressors
env.data2<-env.data[,-c(1,5,6,9,13,15,16,17,19)]
View(env.data2)
###scaling
scale(env.data2)->env.data2.st

##calculating euclidean distance ##only min.V has a lot of 0's
vegdist(env.data2.st,method="euclidean")->env.data.d2

###Flexible clustering needs a par.method. 0.7 makes the cluster similar to Ward's method.
hclust(env.data.d2,method="ward.D")->env.data.ward2

#Selecting number of trees
env.data.ward.gr2 <- cutree(env.data.ward2, 5)

#Identifying clusters through mean values of the original variables
library(plyr)
ddply(env.data2,.(env.data.ward.gr2),function(x){round(colMeans(x),1)})
###
#Determining optimal number of groups for kmeans clustering

wss.env2 <- (nrow(env.data2)-1)*sum(apply(env.data2,2,var))
for (i in 2:15) wss.env[i] <- sum(kmeans(env.data2,
                                         centers=i)$withinss)
#k-means does not provide graphic output
# Non-hierarchical clustering (centers=number of groups)
kmeans (env.data.d2, centers = 5)->env.kmeans.env2
env.kmeans.env2$cluster

#sum of squared error (SSE) scree plot

par(mfrow=c(1,1))
plot(1:10, wss.env2, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
abline(v=4, lty=2, col="blue")
abline(v=6, lty=2, col="violet")

#Selecting number of trees
##based on statistical and knowledge of the site, we decide to cut de tree with 6 groups
env.data.ward.gr2 <- as.factor(cutree(env.data.ward2, 3))

##ploting clustering
plot(env.data.ward2,hang=-1,main="Ward's method",labels=id.env)
rect.hclust(env.data.ward2, k = 3, border = 2:7) # add rectangle



###PCA
dudi.pca(env.data2.st, center=T, scale=T, scann=F, nf=5)->env.pca2

#How many axes do we have to select?
par(mfrow = c(1,1))
screeplot(princomp(env.data2.st,scale=T),bstick=T,type="lines")

#env.pca$eig shows the eigenvalues, it is the importance of each axis
env.pca2$eig

names(env.pca2$eig)<-1:length(env.pca2$eig) # To assign the names of the PC to each eigenvalue
env.pca2$eig #Note the numbers


#explained variance by axis
barplot(env.pca2$eig, xlab="Principal Components",ylab="Eigenvalues",las=1)

#How much variance explains each axis?

(env.pca2$eig*100)/sum(env.pca2$eig)->var.env.pca2
var.env.pca2 #variance explained by axis
barplot(var.env.pca2, xlab="Principal Components",ylab="Explained variance",las=1)

#How much cummulative variance explains each axis?

cumsum(var.env.pca2)#cummulative variance explained 
barplot(cumsum(var.env.pca2), xlab="Principal Components",ylab="Explained variance",las=1)


par(mfrow=c(1,2))
barplot(var.env.pca2, xlab="Principal Components",ylab="Explained variance",las=1)
barplot(cumsum(var.env.pca2), xlab="Principal Components",ylab="Explained variance",las=1)

#$c1 shows the PCA loadings (the weight of each variable in constructing each retained principal component)

env.pca2$c1 #weigth of each variable in the construction of PCA

#$co shows the variable coordinates (arrow coordinates). In the specific case of normalised PCA, these coordinates represent correlation with axes and thus should vary between -1 and 1.

env.pca2$co#correlation with axis

#Graphical representation of PCA Show the pca.biplot (two first axes)

par(mfrow=c(1,1))
#General PCA plot

biplot(env.pca2)

#Plotting the dendrogram
s.arrow(env.pca2$c1, lab = names(env.pca2$tab))

#Plotting cases
s.label(env.pca2$li, sub = "Environmental PCA", csub = 1.5, possub = "topleft")

##ploting pca with groups
##getting symbols and colors for graphical output
sym<-c(rep(16,4),rep(17,3)) # symbols
colour4<-c("green3","steelblue3","firebrick1","magenta3","gold","cyan") # colours

s.class(env.pca2$li,env.data.ward.gr2,pch=sym[env.data.ward.gr2],cpoint=1,col=colour4,cellipse=2,cstar=1,axesell = FALSE,
        sub = "PCA", possub = "bottomright",cgrid = 2, csub = 2)

####
##escolha das variaveis finais ponderando de 3 formas
###1º correlacoes entre variaveis
####2ºdirecoes e aparente relacao no espaco funcional (estarem na mesma direcao nao e significado de correlacao)
####3ºimportancia e impacto na comunidade de macroinvertebrados
####4º idealmente selecionar 1 variavel para cada tipo de stressor esperado (natural,vegetacao riparia,hidrologico,eutrofizacao, fisico-quimico,uso do solo,natural)


####
####
####
#multiple stressors analysis
####
###
###pearson ->linearity correlation
###spearman ->non-linearity correlation
biodiversity_metrics<-diversity
View(biodiversity_metrics)
stressors<-env.data2
View(stressors)
####correlation between metrics
cor(biodiversity_metrics, method = "pearson")->bio1
cor(biodiversity_metrics,method = "spearman")->bio2
View(bio1)
View(bio2)



library(Hmisc)
###correlation between metrics
pearson<-cor(biodiversity_metrics,stressors, method = "pearson")#linearity correlation
View(pearson)
spearman<-cor(biodiversity_metrics,stressors, method = "spearman")#non-linearity correlation
View(spearman)     
View(bio)
bio$ligth1->pc1_ligth
bio$ligth2->pc2_ligth

# Plotting scatter plots
par(mfrow=c(1,1))
var.sel<-c("richness", "shannon", "Hill.shannon", "FRic","FDis")
var.sel2<-c("Altitude","qbr","ihf","DIN","P.PO4","cond","Tmax","DOmin","X.Artificial.500m.","Mean_velocity")

###data transformation
###performed for all variables

# Variable transformation biodiversity metrics
a<-sapply(biodiversity_metrics,class)=="numeric" | sapply(biodiversity_metrics,class)=="integer"# selecting quantitative variables

par(mfrow=c(2,3))
for (i in which(a==T)) hist(biodiversity_metrics[,i], main=names(biodiversity_metrics)[i])


###data transformation biodiversity metrics
#FRic
par(mfrow=c(1,3))
hist(biodiversity_metrics$FRic)#raw data with bad normal distribution
hist(sqrt(biodiversity_metrics$FRic))#use this transformation#better performed
hist(log(biodiversity_metrics$FRic+0.0001))#bad transformation

biodiversity_metrics$FRic<-sqrt(biodiversity_metrics$FRic) ##best result selected as variable

#FDis
par(mfrow=c(1,3))
hist(biodiversity_metrics$FDis)#raw data with bad normal distribution
hist(sqrt(biodiversity_metrics$FDis))#use this transformation#better performed
hist(log(biodiversity_metrics$FDis+0.0001))#bad transformation

##best without transformation

#Hill Shannon
par(mfrow=c(1,3))
hist(biodiversity_metrics$Hill.shannon)#raw data with bad normal distribution
hist(sqrt(biodiversity_metrics$Hill.shannon))#use this transformation#better performed
hist(log(biodiversity_metrics$Hill.shannon+0.0001))#bad transformation

biodiversity_metrics$Hill.shannon<-log(biodiversity_metrics$Hill.shannon+0.0001) ##best result selected as variable

#richness
par(mfrow=c(1,3))
hist(biodiversity_metrics$richness)#raw data with bad normal distribution
hist(sqrt(biodiversity_metrics$richness))#use this transformation#better performed
hist(log(biodiversity_metrics$richness+0.0001))#bad transformation

biodiversity_metrics$richness<-log(biodiversity_metrics$richness+0.0001) ##best result selected as variable

####shannon
par(mfrow=c(1,3))
hist(biodiversity_metrics$shannon)#raw data with bad normal distribution
hist(sqrt(biodiversity_metrics$shannon))#use this transformation#better performed
hist(log(biodiversity_metrics$shannon+0.0001))#bad transformation
##best without transformation

###stressors descriptive/exploratory analysis

####pairs, correlation
pairs(stressors)

####
View(stressors)
stressors.final<-cbind(stressors,buff2$Buff,env.data$mean_light,hidrology2)


###
library(FactorAssumptions, quietly = T, verbose = F)

kmo_bfi <- kmo_optimal_solution(stressors.final, squared = FALSE)
#e=dataframe das variaveis ambientais
kmo_bfi$removed

kmo_bfi

View(stressors.final)
stressors.final<-stressors.final[,-13]
View(stressors.final)




###best correlation by biodiversity metric
###
###Richness
#Richness
View(round(cor(stressors.final),3))
stressors.final<-cbind(stressors.final,ligth1,ligth2)
stressors.final<-stressors.final[,-c(13,14)]
View(stressors.final)
###data transformation biodiversity metrics
#buff
par(mfrow=c(1,3))
hist(stressors.final$`buff2$Buff`)#raw data with bad normal distribution
hist(sqrt(stressors.final$`buff2$Buff`))#use this transformation#better performed
hist(log(stressors.final$`buff2$Buff`+0.0001))#bad transformation

#mean_ligth
par(mfrow=c(1,3))
hist(stressors.final$`env.data$mean_light`)#raw data with bad normal distribution
hist(sqrt(stressors.final$`env.data$mean_light`))#use this transformation#better performed
hist(log(stressors.final$`env.data$mean_light`+0.0001))#bad transformation

View(stressors.final)
View(round(cor(stressors.final),3))
var.sel2<-c("Altitude","qbr","ihf","DIN","P.PO4","cond","Tmax","DOmin","X.Artificial.500m.","Mean_velocity","buff2$Buff","env.data$mean_light")

for (i in var.sel2) {plot(biodiversity_metrics$richness~stressors.final[,i], xlab=i, ylab = "Richness")
  mod <- lm(biodiversity_metrics$richness~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}
###stress.richness
stress.richness<-c("DOmin","cond","DIN","X.Artificial.500m.","qbr","P.PO4")

par(mfrow=(c(1,1)))
for (i in stress.richness) {plot(biodiversity_metrics$richness~stressors.final[,i], xlab=i, ylab = "Richness")
  mod <- lm(biodiversity_metrics$richness~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}

####Hill Shannon
par(mfrow=c(1,1))
var.sel3<-c("DIN","cond","qbr","DOmin")

for (i in var.sel2) {plot(biodiversity_metrics$Hill.shannon~stressors.final[,i], xlab=i, ylab = "Hill Shannon")
  mod <- lm(biodiversity_metrics$Hill.shannon~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}

###stress.hill.shannon
stress.hill.shannon<-c("qbr","DIN","P.PO4","cond","DOmin","X.Artificial.500m.")

for (i in stress.hill.shannon) {plot(biodiversity_metrics$Hill.shannon~stressors.final[,i], xlab=i, ylab = "Hill Shannon")
  mod <- lm(biodiversity_metrics$Hill.shannon~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}


#####FRic
for (i in var.sel2) {plot(biodiversity_metrics$FRic~stressors.final[,i], xlab=i, ylab = "FRic")
  mod <- lm(biodiversity_metrics$FRic~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}

par(mfrow=c(1,1))
var.sel<-c("DOmin","DIN","cond","qbr")

for (i in var.sel) {plot(biodiversity_metrics$FRic~stressors[,i], xlab=i, ylab = "FRic")
  mod <- lm(biodiversity_metrics$FRic~stressors[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}

ggplot(stressors, aes(DOmin, biodiversity_metrics$FRic) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

ggplot(stressors, aes(qbr, biodiversity_metrics$FRic) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

ggplot(stressors, aes(DIN, biodiversity_metrics$FRic) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

ggplot(stressors, aes(cond, biodiversity_metrics$FRic) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

ggplot(stressors.final, aes(env.data$Tmax, biodiversity_metrics$FRic) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))


###stress.FRIC
stress.FRIC<-c("qbr","ihf","DIN","P.PO4","cond","DOmin","X.Artificial.500m.")
for (i in stress.FRIC) {plot(biodiversity_metrics$FRic~stressors.final[,i], xlab=i, ylab = "FRic")
  mod <- lm(biodiversity_metrics$FRic~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}


#####FDis
for (i in var.sel2) {plot(biodiversity_metrics$FDis~stressors.final[,i], xlab=i, ylab = "FDis")
  mod <- lm(biodiversity_metrics$FDis~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}

par(mfrow=c(1,1))
var.sel4<-c("Mean_velocity","DIN")

ggplot(stressors.final, aes(Mean_velocity, biodiversity_metrics$FDis) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

ggplot(stressors, aes(DIN, biodiversity_metrics$FDis) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

###stress.FDIS
stress.FDIS<-c("Altitude","qbr","DIN","cond","DOmin","X.Artificial.500m.","Mean_velocity")
for (i in stress.FDIS) {plot(biodiversity_metrics$FDis~stressors.final[,i], xlab=i, ylab = "FDis")
  mod <- lm(biodiversity_metrics$FDis~stressors.final[,i])
  abline(mod, col="blue", lwd=3)
  r2 <- bquote(italic(r)^2 == .(round(summary(mod)$r.squared,2)))
  mtext(r2, line = -4.5, adj = 0.9, cex = 1.2, font = 2)
}

ggplot(stressors.final, aes(ihf, biodiversity_metrics$FDis) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

ggplot(stressors.final, aes(Tmax, biodiversity_metrics$FDis) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))

ggplot(stressors.final, aes(DOmin, biodiversity_metrics$FDis) ) +
  geom_point() +
  stat_smooth(method = gam, formula = y ~ s(x))


####
#fazer o boxplot do bmwp ou iptn
#fazer os graficos da hier.part por metrica de diversidade
#fazer os modelos multiplos para as metricas

library(ggplot2)
library(mgcv)

####boxplot IPTN/ IBMWP
###boxplot bmwp, aspt, iptin
library(ggplot2)
library("biotic")
library("vegan")


family.taxa<-read.table("familia taxa.txt",sep="\t",dec=".",header=T)
family.taxa[family.taxa==0]<-NA
View(family.taxa) #checking data frame
family.taxa2<-read.table("familia taxa2.txt",sep="\t",dec=".",header=T)
View(family.taxa2)
id.fam<-family.taxa2[,1] #add site to an object
family.taxa2$Rótulos.de.Linha->row.names(family.taxa2)
family.taxa2<-family.taxa2[,-1]#remove sites from the original dataset to do an exploratory analysis, calculating
#distances between sites
View(family.taxa2)
View(id.fam)
EPT<-read.table("EPT-ETD.txt", sep="\t", dec = ".", header=T)
View(EPT)

Group1<- c("AGR1","AND1","CASM1","EST2","FERRO1","LAM1","MOU1","TOJ1")
Group2<- c("ARG1","COV1","MACI1","RAB1","RPON1","SOU1","VIZ1")
Group3<- c("ARG2","AVE1","BIL1","CAN1","CAV1","CAV2","FER1","GER1","OLI1","RAB2","SAL1","TAVE1","VIL1","VILC1")
Group4<- c("ARG3","AZE1","AZE2","CAB1","CBR1","ERM1","FOR1","VEZ1","VEZ2","VEZ3")
Group5<- c("EST1","LAB1","NES1","PON1","ROD1","SAN1","SEL1","TAVE2","VEI1")
Group6<- c("HOM1","MACE1")

##calcular o BMWP,ASPT e número de taxas, considerando a resolucao taxonomica de familia
calc.ind<-calcindex(family.taxa)
View(calc.ind)
Ntaxa<-calc.ind$Ntaxa
BMWP<-calc.ind$BMWP
ASPT<-calc.ind$ASPT
richness2<-specnumber(family.taxa2) #species richness
# Pielou's evenness (J)
J2 <- diversity(family.taxa2)/log(richness2)
IPtIn<-(richness2*0.25)+(EPT$EPT*0.15)+(J2*0.1)+((ASPT-2)*0.3)+(log(EPT$ETD+1)*0.2)

Grupos<- c("Group1","Group2","Group3","Group4","Group5","Group6")

colour<-c("green3","gold","steelblue3","firebrick1","cyan","magenta3") # colours

data<-data.frame(calc.ind$Sample,EPT$EPT,ASPT,Ntaxa,BMWP,richness2,J2,IPtIn,bio.gr)
View(data)

library(tidyverse)
library(lubridate)

data1<-as.data.frame(data)
View(data1)

# Saving datasets
write.table(data,"boxplot_index.txt",sep="\t",dec=".")

data1<-read.table("boxplot_index.txt",sep="\t",dec=".",header=T)
View(data1)
data1<-data1[,-1]

####tuckey test
# library
library(multcompView)

# What is the effect of the treatment on the value ?
data1$bio.gr<-as.factor(data1$bio.gr)
bmwp.lm<-lm(BMWP~bio.gr, data=data1)
bmwp.aov<-aov(bmwp.lm)
summary(bmwp.aov)

library(agricolae)
tukey.test <- HSD.test(bmwp.aov, trt="bio.gr")
tukey.test

tukey<-tukey.test$groups
bio.gr.tukey<-c("3","4","1","6","2","5")
tukey<-cbind(tukey,bio.gr.tukey)
View(tukey)
boxplot<-boxplot(BMWP ~ bio.gr,data=data1  , ylim=c(min(BMWP) , 1.1*max(BMWP)) ,aes(bio.gr, col=colour[as.factor(tukey$bio.gr.tukey)]), ylab="BMWP" , main="Biotic Index")
boxplot
 
# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=bmwp.aov, "bio.gr", conf.level=0.95)

# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$bio.gr=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$bio.gr) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY , "bio.gr")

# A panel of colors to draw each group with the same color :
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255), 
  rgb(111,145,202,maxColorValue = 255)
)


# Draw the basic boxplot
a <- boxplot(BMWP ~ bio.gr ,data=data1, ylim=c(min(BMWP) , 1.1*max(BMWP)) , col=colour[as.factor(LABELS[,1])] , ylab="BMWP" , main="Biotic Index")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:nlevels(data1$bio.gr)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=colour[as.factor(LABELS[,1])] )

######Richness
# What is the effect of the treatment on the value ?
data1$bio.gr<-as.factor(data1$bio.gr)
richness.lm<-lm(richness~bio.gr, data=data1)
richness.aov<-aov(richness.lm)
summary(richness.aov)

library(agricolae)
tukey.test2 <- HSD.test(richness.aov, trt="bio.gr")
tukey.test2

tukey2<-tukey.test2$groups
bio.gr.tukey2<-c("3","1","4","6","2","5")
tukey2<-cbind(tukey2,bio.gr.tukey)

# Tukey test to study each pair of treatment :
TUKEY1 <- TukeyHSD(x=richness.aov, "bio.gr", conf.level=0.95)

# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$bio.gr=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$bio.gr) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
LABELS <- generate_label_df(TUKEY1 , "bio.gr")

# A panel of colors to draw each group with the same color :
my_colors <- c( 
  rgb(143,199,74,maxColorValue = 255),
  rgb(242,104,34,maxColorValue = 255), 
  rgb(111,145,202,maxColorValue = 255)
)


# Draw the basic boxplot 
a <- boxplot(richness ~ bio.gr ,data=data1, ylim=c(min(richness) , 1.1*max(richness)) , col=colour[as.factor(LABELS[,1])] , ylab="Richness" , main="Biotic Index")

# I want to write the letter over each box. Over is how high I want to write it.
over <- 0.1*max( a$stats[nrow(a$stats),] )

#Add the labels
text( c(1:nlevels(data1$bio.gr)) , a$stats[nrow(a$stats),]+over , LABELS[,1]  , col=colour[as.factor(LABELS[,1])] )

