# directory = "C:/Users/asus/Desktop/Pedro Gonçalves Shotgun Method/Streameco"
setwd("C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco")
# setwd(directory)

# Loading packages
library(vegan)
library(ade4)
library(usdm)
library(mgcv)


# Loading data
filenames <- list.files(directory, pattern="*riqueza.txt")
all_txt <- lapply(filenames, read.csv, sep="\t", dec=".", header=T)
names(all_txt) <- c("bact_genus", "bact_order", "bact_species", "fun_genus", "fun_order", "fun_species")
invisible(lapply(names(all_txt), function(x) assign(x, all_txt[[x]], envir=.GlobalEnv)))

# taxa <- list(species)
# taxa <- lapply(taxa,function(df){
#   df <- df[nrow(df), ]        
#   df <- df[, -c(1, 52)]
#   df <- strtoi(df)
#   df
#   
# })

div <- bact_species
div <- div[6180,-c(1, 52)]
View(div)
# div = strtoi(div)

div = rbind(div, fun_species[999, -c(1, 52)])
div = t(div)
colnames(div) = c("bacteria", "fungi")
div = as.data.frame(div)
div$bacteria = as.numeric(as.character(div$bacteria)) 
div$fungi = as.numeric(as.character(div$fungi))


# par(mar = c(1, 1, 1, 1))
# par(mfrow=c(3,3))
# for (i in 1:length(taxa)) {
#   hist(taxa[[i]])
#   hist(sqrt(taxa[[i]]))
#   hist(log(taxa[[i]]+0.1))
# }


bact_reads = read.csv("Bacteria_species_reads.txt", sep="\t", dec=",", header=T)
rownames(bact_reads) = bact_reads$species
bact_reads = bact_reads[, -1]
View(bact_reads)

fung_reads = read.csv("Fungi_species_reads.txt", sep="\t", dec=",", header=T)
rownames(fung_reads) = fung_reads$species
fung_reads = fung_reads[, -1]
View(fung_reads)



# Índice de Shannnon e Simpson bactérias
shannon_bact = diversity(t(bact_reads), index = "shannon", MARGIN = 1, base = exp(1))
simpson_bact = diversity(t(bact_reads), index = "simpson")
shannon_fung = diversity(t(fung_reads), index = "shannon", MARGIN = 1, base = exp(1))
simpson_fung = diversity(t(fung_reads), index = "simpson")

# Pielou
J_bact = diversity(t(bact_reads)/log(as.numeric(div$bacteria)))
J_fung = diversity(t(fung_reads)/log(as.numeric(div$fungi)))

# bray
S_bact = apply(t(bact_reads)>0, 1, sum)
N_bact = apply(t(bact_reads), 1, sum)
margalef_bact = (S_bact-1)/log(N_bact)

S_fung = apply(t(fung_reads)>0, 1, sum)
N_fung = apply(t(fung_reads), 1, sum)
margalef_fung = (S_fung-1)/log(N_fung)


###bray curtis 1st axis of PCoA
# Functional space
bray_curtis_bact = vegdist(t(bact_reads), method = "bray")
dudi.pco(bray_curtis_bact,scannf = F,nf=10)->tr.pco

Bray_Curtis_bact<-tr.pco$li$A1

###
bray_curtis_fung = vegdist(t(fung_reads), method = "bray")
dudi.pco(bray_curtis_fung,scannf = F,nf=10)->tr.pco2
Bray_Curtis_fung<-tr.pco2$li$A1

div = cbind(div, shannon_bact, simpson_bact, J_bact, margalef_bact, Bray_Curtis_bact, shannon_fung, simpson_fung, J_fung, margalef_fung, Bray_Curtis_fung)
div2=div

par(mfrow=c(2,3))
for (i in 1:ncol(div)) hist(div[,i], main=names(div)[i])

par(mfrow=c(1,3))
hist(div$bacteria)
hist(sqrt(div$bacteria))
hist(log(div$bacteria+0.01))

div$bacteria<-log(div$bacteria+0.01)

hist(div$simpson_bact)
hist(sqrt(div$simpson_bact))
hist(log(div$simpson_bact+0.01))

hist(div$Bray_Curtis_bact)
hist(sqrt(div$Bray_Curtis_bact))
hist(log(div$Bray_Curtis_bact+0.01))

hist(div$Bray_Curtis_fung)
hist(sqrt(div$Bray_Curtis_fung))
hist(log(div$Bray_Curtis_fung+0.01))
shapiro.test(div$bray_fung)
shapiro.test(sqrt(div$bray_fung))
shapiro.test(log(div$bray_fung+0.001))

div<-as.data.frame(scale(div))

env.data <- read.table("var ambientais.txt", sep="\t", dec=".", header=T)
row.names(env.data) <- env.data$code
env.data <- env.data[, -1]
View(env.data)
e.data=env.data

# Selecting key environmental variables
env.data2 <- env.data[, c("DIN", "Tmax", "DOmin", "mean.Velocity", "shadow", "Altitude", "Buff")]

# Transforming environmental variables
a <- sapply(env.data2,class)=="numeric" | sapply(env.data2,class)=="integer"# selecting quantitative variables

par(mfrow=c(3,3))
for (i in which(a==T)) hist(env.data2[,i], main=names(env.data2)[i])

env.data2$DIN <- sqrt(env.data2$DIN)
env.data2$mean.Velocity <- sqrt(env.data2$mean.Velocity)
env.data2$shadow <- car::logit(env.data2$shadow/100)
env.data2$Altitude <- sqrt(env.data2$Altitude)
env.data2$Buff <- log(env.data2$Buff+1)

par(mfrow=c(3,3))
for (i in which(a==T)) hist(env.data2[,i], main=names(env.data2)[i])

# Standardizing environmental predictors and diversity
env.data2 <- data.frame(scale(env.data2))
env.data3<-data.frame(env.data)

env.data$DIN <- sqrt(env.data$DIN)
env.data$mean.Velocity <- sqrt(env.data$mean.Velocity)
env.data$shadow <- car::logit(env.data$shadow/100)
env.data$Altitude <- sqrt(env.data$Altitude)
env.data$Buff <- log(env.data$Buff+1)
env.data$Slope <- log(env.data$Slope+0.1) ##best result selected as variable
env.data$Altitude <- sqrt(env.data$Altitude) ##best result selected as variable
env.data$P.PO4 <- log(env.data$P.PO4+0.01) ##best result selected as variable
env.data$cond <- log(env.data$cond+0.01) ##best result selected as variable
env.data$mean_light <- log(env.data$mean_light+0.01)##best result selected as variable
env.data$X.Agriculture.500m. <- car::logit(env.data$X.Agriculture.500m./100)##best result selected as variable
env.data$LUI_500m <- car::logit(env.data$LUI_500m/100)##best result selected as variable
env.data$mean.Depth <- log(env.data$mean.Depth+0.01)##best result selected as variable
env.data$Width <- log(env.data$Width+1)##best result selected as variable
env.data$Discharge <- log(env.data$Discharge+0.01)##best result selected as variable

env <- env.data
env.data <- data.frame(scale(env.data))

dados = as.data.frame(cbind(div, env.data))

#env-data2/env.data/div
View(cor(div, method="pearson"))
###riqueza e bray_curtis não apresentam colinearidade
pearson<-cor(div,env.data, method = "pearson") #shows linear correlation
spearman<-cor(div,env.data, method = "spearman") #shows linear and non-linear correlation
kendall<-cor(div,env.data, method = "kendall")#also shows non_linear monotonic correlation
write.table(pearson,"pearson.txt",sep="\t",dec=".") #shows linear correlation
write.table(spearman,"spearman.txt",sep="\t",dec=".") #shows linear and non-linear correlation
write.table(kendall,"kendall.txt",sep="\t",dec=".") #also shows non_linear monotonic correlation

#
par(mfrow=c(1,1))
####all linear models
pdf("linear_models.pdf")
for (i in 1:ncol(env.data)) {
  for(j in 1:ncol(div)) {
    
    ##plot(div[,y]~env.data[,i], xlab=colnames(env.data[,i]), ylab=colnames(div[,y]))
    #   }}
    mod <- lm(div[,j]~env.data[,i])
    x<-summary(mod)
    x<-as.numeric(unlist(x$coefficients[2,4]))
    if (x<=0.05) {
      y=if (x<=0.001) {as.character("<0.001***")} else if (x<=0.01) {as.character("<0.01**")
      }else if (x<=0.05) {as.character("<0.05*")} else {as.character("=n.s.")}
      plot(div[,j]~env.data[,i], xlab=colnames(env.data[i]), ylab=colnames(div[j]))
      abline(mod, col="blue", lwd=3)
      y<-as.character(y)
      p<-("p")
      p<-paste(p,y)
      mtext(p, line=-1.5, adj = 1, cex = 1.2, font = 2)
      r2 <-bquote(paste(bold(r^2 == .(round(summary(mod)$r.squared,2)))))
      mtext(r2, line = -3, adj = 1, cex = 1.2, font = 2)
      par(mfrow=c(2,2))
      plot(mod)
      par(mfrow=c(1,1))}
    
    #title(main="Diversity", xlab= colnames(env.data[i]), ylab=colnames(div[i][j]))
    #axis(side=2, labels = colnames(div[y]), las=0)
  }
}


dev.off()

par(mfrow=c(1,1))
library(mgcv)
###non linear models with data transformation 
####all linear models
pdf("non_linear_models_with_data_transformation.pdf")
for (i in 1:ncol(env.data)) {
  for(j in 1:ncol(div)) {
    ##plot(div[,y]~env.data[,i], xlab=colnames(env.data[,i]), ylab=colnames(div[,y]))
    #   }}
    gam<-gam(div[,j]~s(env.data[,i])) 
    x<-summary(gam)
    p=as.character(round(as.numeric(unlist((x$s.table[,4]))),3))
    if (p<=0.05) {
      plot(gam, main='Diversity', xlab=colnames(env.data[i]), ylab=colnames(div[j]), se=TRUE)
      p=if (p<=0.001) {as.character("<0.001***")} else if (p<=0.01) {as.character("<0.01**")
      }else if (p<=0.05) {as.character("<0.05*")} else {as.character("=n.s.")}
      p=paste("p",p)
      r2 <-bquote(paste(bold(r^2 == .(as.character(round(as.numeric(x$r.sq),3))))))
      dev=as.character(round(as.numeric(x$dev.expl),3))
      dev <-(paste("dev", "=",dev))
      mtext(p, line=-1.5, adj = 1, cex = 1.2, font = 2)
      mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
      mtext(dev, line=-3.5, adj = 1, cex = 1.2, font = 2)
      par(mfrow=c(2,2))
      gam.check (gam)
      par(mfrow=c(1,1))
      }

 }
}

dev.off()

pdf("non_linear_models_without_data_transformation.pdf")
for (i in 1:ncol(e.data)) {
  for(j in 1:ncol(div2)) {
    ##plot(div2[,y]~e.data[,i], xlab=colnames(e.data[,i]), ylab=colnames(div2[,y]))
    #   }}
    gam<-gam(div2[,j]~s(e.data[,i])) 
    x<-summary(gam)
    p=as.character(round(as.numeric(unlist((x$s.table[,4]))),3))
    if (p<=0.05) {
      plot(gam, main='Diversity', xlab=colnames(e.data[i]), ylab=colnames(div2[j]), se=TRUE)
      p=if (p<=0.001) {as.character("<0.001***")} else if (p<=0.01) {as.character("<0.01**")
      }else if (p<=0.05) {as.character("<0.05*")} else {as.character("=n.s.")}
      p=paste("p",p)
      r2 <-bquote(paste(bold(r^2 == .(as.character(round(as.numeric(x$r.sq),3))))))
      dev=as.character(round(as.numeric(x$dev.expl),3))
      dev <-(paste("dev", "=",dev))
      mtext(p, line=-1.5, adj = 1, cex = 1.2, font = 2)
      mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
      mtext(dev, line=-3.5, adj = 1, cex = 1.2, font = 2)
      par(mfrow=c(2,2))
      gam.check (gam)
      par(mfrow=c(1,1))
    }
  }
}

dev.off()