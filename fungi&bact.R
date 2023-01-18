directory = "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco"
setwd(directory)

# Loading packages
library(vegan)
library(ade4)
library(usdm)
library(xlsx)
library(ggplot2)


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
div <- div[6187,-c(1, 52)]
# div = strtoi(div)

div = rbind(div, fun_species[1000, -c(1, 52)])
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

bact_species_total = bact_species[-6187, ]
fun_species_total = fun_species[-1000, ]


bact_reads = read.csv("Bacteria_species_reads.txt", sep="\t", dec=",", header=T)
bact_reads[is.na(bact_reads)] = 0
rownames(bact_reads) = bact_reads$espécies
bact_reads = bact_reads[, -1]

# Índice de Shannnon e Simpson bactérias
shannon_bact = diversity(t(bact_reads), index = "shannon", MARGIN = 1, base = exp(1))
simpson_bact = diversity(t(bact_reads), index = "simpson")

# Pielou
J_bact = diversity(t(bact_reads)/log(as.numeric(div$bacteria)))

# Margalef
S_bact = apply(t(bact_reads)>0, 1, sum)
N_bact = apply(t(bact_reads), 1, sum)
margalef_bact = (S_bact-1)/log(N_bact)


bray_curtis_bact = vegdist(t(bact_reads), method = "bray")

#building a cluster
#k-number of groups selected
#using ward method, the most strict method of classification
hclust(bray_curtis_bact, method="ward.D") -> bio.ward
bray_bact = bio.ward$height
bray_bact = append(bray_bact, 0)


div = cbind(div, shannon_bact, simpson_bact, J_bact, margalef_bact, bray_bact)


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

hist(div$bray_bact)
hist(sqrt(div$bray_bact))
hist(log(div$bray_bactt+0.01))

# Índice de Shannnon e Simpson fungos
# shannon_fun = diversity(fun_species_total, index = "shannon", MARGIN = 1, base = exp(1))



env.data <- read.table("var ambientais.txt", sep="\t", dec=".", header=T)
row.names(env.data) <- env.data$code
env.data <- env.data[, -1]
View(env.data)

# Selecting key environmental variables
env.data2 <- env.data[, c("DIN", "Tmax", "DOmin", "mean.Velocity", "shadow", "Altitude", "Buff")]

# Transforming environmental variables
a <- sapply(env.data2,class)=="numeric" | sapply(env.data2,class)=="integer"# selecting quantitative variables

par(mar = c(1, 1, 1, 1))
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

dados = cbind(div, env.data)

# No linear correlation
for (i in ncol(dados)) {
  nlcor(dados$i, dados$i, plt = T)
}

# Checking redundancy among diversity indices and creating a file.txt
for (i in 1:length(div)) {
  View(cor(taxa[[i]],env.data, method = "pearson"))
  View(cor(taxa[[i]],env.data, method = "spearman"))
  View(cor(taxa[[i]],env.data, method = "kendall"))
}


#View(cor(taxa[[i]],env.data, method = "pearson"))
#View(cor(taxa[[i]],env.data, method = "spearman"))

write.xlsx(cor(taxa[[1]],env.data, method = "pearson"), "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco/species_pearson.xlsx", sheetName = "Sheet1",
           col.names = TRUE, row.names = TRUE, append = FALSE)

write.xlsx(cor(taxa[[2]],env.data, method = "pearson"), "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco/genus_pearson.xlsx", sheetName = "Sheet1",
           col.names = TRUE, row.names = TRUE, append = FALSE)

write.xlsx(cor(taxa[[3]],env.data, method = "pearson"), "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco/order_pearson.xlsx", sheetName = "Sheet1",
           col.names = TRUE, row.names = TRUE, append = FALSE)

ggplot(cbind(taxa[[1]], env.data), aes(taxa[[1]], P.PO4) ) +
  geom_point() +
  stat_smooth()

ggplot(cbind(taxa[[2]], env.data), aes(taxa[[2]], P.PO4) ) +
  geom_point() +
  stat_smooth()

ggplot(cbind(taxa[[3]], env.data), aes(taxa[[3]], P.PO4) ) +
  geom_point() +
  stat_smooth()

par(mar = c(1, 1, 1, 1))
plot(cbind(taxa[[1]], env.data2))
par(mfrow=c(1,1))
plot(taxa[[1]]~env.data$P.PO4)


modelo_riqueza <- lm(taxa[[1]]~env.data$P.PO4)
summary(modelo_riqueza)
abline(modelo_riqueza, col="blue")

par(mar = c(1, 1, 1, 1))
plot(cbind(taxa[[2]], env.data2))
plot(taxa[[2]]~env.data$P.PO4)
par(mfrow=c(1,1))

modelo_riqueza <- lm(taxa[[2]]~env.data$P.PO4)
summary(modelo_riqueza)
abline(modelo_riqueza, col="blue")

par(mar = c(1, 1, 1, 1))
plot(cbind(taxa[[2]], env.data2))
plot(taxa[[2]]~env.data$P.PO4)
par(mfrow=c(1,1))

modelo_riqueza <- lm(taxa[[2]]~env.data$P.PO4)
summary(modelo_riqueza)
abline(modelo_riqueza, col="blue")

