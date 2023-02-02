# setwd("C:/Users/asus/Desktop/Pedro Gonçalves Shotgun Method/Streameco")
# directory="C:/Users/asus/Desktop/Pedro Gonçalves Shotgun Method/Streameco"

directory = "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco"
setwd(directory)

# Loading packages
library(vegan)
library(ade4)
library(usdm)
library(mgcv)


# Loading data
filenames <- list.files(directory, pattern="*riqueza.txt")
all_txt <- lapply(filenames, read.csv, sep="\t", dec=",", header=T)
names(all_txt) <- c("bact_genus", "bact_order", "bact_species", "fun_genus", "fun_order", "fun_species")
invisible(lapply(names(all_txt), function(x) assign(x, all_txt[[x]], envir=.GlobalEnv)))

filenames_reads <- list.files(directory, pattern="*reads.txt")
all_txt <- lapply(filenames, read.csv, sep="\t", dec=",", header=T)
names(all_txt) <- c("bact_genus_reads", "bact_order_reads", "bact_species_reads", "fun_genus_reads", "fun_order_reads", "fun_species_reads")
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
shannon_fung = diversity(t(fung_reads), index = "shannon", MARGIN = 1, base = exp(1))



###bray curtis 1st axis of PCoA
# Functional space
bray_curtis_bact = vegdist(t(bact_reads), method = "bray")
dudi.pco(bray_curtis_bact,scannf = F,nf=10)->tr.pco

Bray_Curtis_bact<-tr.pco$li$A1

###
bray_curtis_fung = vegdist(t(fung_reads), method = "bray")
dudi.pco(bray_curtis_fung,scannf = F,nf=10)->tr.pco2
Bray_Curtis_fung<-tr.pco2$li$A1

div = cbind(div, shannon_bact, Bray_Curtis_bact, shannon_fung, Bray_Curtis_fung)
div=as.data.frame(div)


env.data <- read.table("var ambientais.txt", sep="\t", dec=".", header=T)
row.names(env.data) <- env.data$code
env.data <- env.data[, -1]
View(env.data)


library(mgcv)
env.data=as.data.frame(env.data)
env.data2=env.data[,c(12,25,32,37)]
env.data=env.data[,-c(12,25,32,37)]
row.names(env.data2)<-row.names(env.data)



div=as.data.frame(scale(div))
env.data=env.data[,c(1,2,4,5,6,15,17,22,26,28,29,30,31)]
env.data2 <- as.data.frame(scale(env.data))
dudi.pca(env.data2, center = TRUE, scale = FALSE, scannf = F,nf=4)->metric.pca
biplot(metric.pca)
metric.pca$li$Axis1 -> PCA1
metric.pca$li$Axis2 -> PCA2
# PCA Axis importance (explained variance)
round(metric.pca$eig[1:4]/sum(metric.pca$eig),2)
env.data = cbind(env.data, PCA1, PCA2)
env.data = as.data.frame(scale(env.data))

pdf("Gam_models.pdf")
for (i in 1:length(env.data)) {
  for(j in 1:length(div)) {
    ##plot(div2[,y]~e.data[,i], xlab=colnames(e.data[,i]), ylab=colnames(div2[,y]))
    #   }}
    print(i)
    gam<-gam(div[,j]~s(env.data[,i], fx = FALSE, k=-1,  bs = "cr"))
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

biplot(metric.pca)

dev.off()

row.names(fun_order_reads) <- fun_order_reads[, 1]
fun_order_reads <- fun_order_reads[, -1]
fun_order_reads <- as.data.frame(fun_order_reads)

nmds <- metaMDS(fun_order_reads, distance = "bray")
nmds

plot(nmds)
stressplot(nmds)
plot(nmds, type = "n")+ #displays empty ordination space
  ggrepel::geom_text_repel(data = nmds, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)
points(nmds, display = "sites", col = "blue", labels(row.names(nmds$points)))
# displays site points where symbols (pch) are different management options and colour (col) are different land uses
legend("topright", legend = c(levels(dune.env$Management), levels(dune.env$Use)), pch = c(16, 8, 17, 11, 16, 16, 16), col = c("black","black","black","black","blue", "orange", "black"), bty = "n", cex = 1) # displays symbol and colour legend
legend("topleft", "stress = 0.118", bty = "n", cex = 1) # displays legend text of stress value 
ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)

nmds2 <- metaMDS(fun_reads, distance = "bray")
nmds2

plot(div$shannon_bact~env.data$mean.Velocity)
