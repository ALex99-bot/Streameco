# Set the directory and file name
directory <- "/home/pedro/PycharmProjects/Streameco"
# directory <- "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco"
# directory <-"C:/Users/asus/Desktop/Streameco"
setwd(directory)

file <- "total_reads.xlsx"
file2 <- "top10.xlsx"
file3 <- "menor_cutoff_reads.xlsx"
file4 <- "hyphomycetes_aquaticos.xlsx"
path <- file.path(directory, file)
path2 <- file.path(directory, file2)
path3 <- file.path(directory, file3)
path4 <- file.path(directory, file4)


# Load the required packages
library(readxl)
library(mgcv)
library(ade4)
library(usdm)
library(vegan)
library(purrr)
library(ggplot2)
library(ggrepel)
library(factoextra)
library(readxl)
library(abdiv)
library(plyr)
library(dplyr)
library(ggpmisc)
library(aomisc)
library(nlstools)
library(nlshelper)
library(xcms)
library(hillR)
library(nlraa)


read_excel_sheets <- function(path, file){
  my_sheet_names <- excel_sheets(path)
  my_sheets <- lapply(my_sheet_names, function(x) read_excel(file, sheet = x))
  names(my_sheets) <- my_sheet_names
  return(my_sheets)
}


excel2dataframe <- function(x) {
  x <- as.data.frame(x)
  row.names(x) <- x[, 1]
  x <- x[, -1]
}

bact_fung_taxa <- lapply(read_excel_sheets(path, file), excel2dataframe)

# Variáveis ambientais
bacias <- read_excel("STREAMECO database - environment (copy).xlsx", "LandUse")
bacias <- as.data.frame(bacias)
river_basin <- bacias[2:nrow(bacias), 3]
river_basin <- as.data.frame(river_basin)
row.names(river_basin) <- bacias[2:51, 1]

env.data <- read.table("var ambientais.txt", sep="\t", dec=".", header=T)
row.names(env.data) <- env.data$code
env.data <- env.data[, -1]
env_data <- as.data.frame(env.data)
to_remove <- c("mean_light", "Artificial.500m.", "Agriculture.500m.", "Pasture.500m.", "Natural.500m.", "N.NH4", "N.NO2",
               "N.NO3", "Tmean", "Tmin", "TCV", "DOmean", "Domax", "Artificial.100m.", "Agriculture.100m.", "Pasture.100m.",
               "Natural.100m.", "Artificial_subasin", "Agriculture_subasin", "Pasture_subasin", "Natural_subasin", "LUI_100m", "LUI_500m")
env_data <- env_data[ , !(names(env_data) %in% to_remove)]
env_data_nt <- env_data[ , !(names(env_data) %in% to_remove)]

# Transformações
env_data$Altitude <- sqrt(env_data$Altitude)
env_data$shadow <- car::logit(env_data$shadow/100)
env_data$DIN <- sqrt(env_data$DIN)
env_data$P.PO4 <- sqrt(env_data$P.PO4)
env_data$cond <- log(env_data$cond+0.0001)
env_data$mean.Depth <- log(env_data$mean.Depth+0.0001)
env_data$mean.Velocity <- sqrt(env_data$mean.Velocity)
env_data$Discharge <- log(env_data$Discharge+0.0001)

# par(mfrow=c(1, 1))
# hist(env_data$LUI_100m)
# hist(sqrt(env_data$LUI_100m))
# wilcox.test(sqrt(env_data$LUI_100m))

# Apenas para hyphomycetes
# env_data <- env_data[-40,]

# matriz_final <- Reduce(function(df1, df2) cbind(df1, df2), shannon, init = env.data)
# indices <- Reduce(function(df1, df2) cbind(df1, df2), shannon, init = matriz_final)
# indices <- indices[,-1:-31]

diversidade <- lapply(bact_fung_taxa, function(x) specnumber(t(x)))

shannon_index <- lapply(bact_fung_taxa, function(x) diversity(t(x), index = "shannon", base = exp(1)))
hill_shannon <- lapply(bact_fung_taxa, function(x) hill_taxa(t(x), q = 1, MARGIN = 1, base = exp(1)))

pielou_index <- lapply(bact_fung_taxa, function(x) diversity(t(x))/log(specnumber(t(x))))

margalef_index <- lapply(bact_fung_taxa, function(x) apply(x, 2, margalef))

simpson_index <- lapply(bact_fung_taxa, function(x) diversity(t(x), index = "simpson", MARGIN = 1, base = exp(1)))

# Está mal, é preciso corrigir
# bray_index <- lapply(bact_fung_taxa, function(x) {pco <- dudi.pco(vegdist(t(x), method = "bray"), scannf = F, nf = 10)
# x <- return(pco$li$A1)
# row.names(x) <- col.names(bact_fung_taxa)})


# Para hifomicetes, mudar times=1
#env_data <- scale(env_data)
#env_data_nt <- scale(env_data_nt)
#PCA_env <- prcomp(env_data, scale = TRUE)
PCA_env <- dudi.pca(env_data, center = TRUE, scale = TRUE, nf=5, scannf = FALSE)


# PCA_env <- rep(list(dudi.pca(env_data, center = TRUE, scale = TRUE, nf=5, scannf = FALSE)), times = 12)
biplot(PCA_env)
s.arrow(PCA_env$c1, lab=names(PCA_env$tab))
var_exp <- (PCA_env$eig*100)/sum(PCA_env$eig)

env_data <- cbind(env_data, PC1 = PCA_env$li[,1])
env_data_nt <- cbind(env_data_nt, PC1 = PCA_env$li[,1])

indices <- read_excel("bioindices.xlsx")
indices <- as.data.frame(indices)
row.names(indices) <- indices$local
indices <- indices[, -1]

indices <- as.data.frame(indices)
env_data_nt <- as.data.frame(env_data_nt)
env_data <- as.data.frame(env_data)

indices2 <- indices[,c("Bacteria_species_div","Bacteria_species_shannon",
           "Fungi_species_div", "Fungi_species_shannon")]

# cor(indices2, env_data_nt, method="pearson")
# cor(indices2, env_data_nt, method="spearman")
#
# wb <- createWorkbook()
# sheet <- createSheet(wb, "pearson")
# addDataFrame(cor(indices2, env_data_nt, method="pearson"), sheet=sheet,
#              startColumn=1, row.names=T)
#
# sheet <- createSheet(wb, "spearman")
#
# addDataFrame(cor(indices2, env_data_nt, method="spearman"), sheet=sheet,
#              startColumn=1, row.names=T)
#
# saveWorkbook(wb, "correlacoes.xlsx")

modelos <- as.data.frame(cbind(indices2, env_data))
modelos <- cbind(modelos, vel = modelos$mean.Velocity)

modelos_nt <- as.data.frame(cbind(indices2, env_data_nt))
modelos_nt <- cbind(modelos_nt, vel = modelos_nt$mean.Velocity)

# jpeg("bact_shannon_velocidade.jpg")
par(mfrow=c(1,1))
plot(Fungi_species_div ~ cond, data = modelos_nt)
r2 <- bquote(paste(bold(R^2 == .(11.6))))
mtext("p<0.01", line=-1.5, adj = 1, cex = 1.2, font = 2)
mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
mod <- lm(Fungi_species_div ~ Altitude, data = modelos_nt)
abline(mod)
# dev.off()
summary(mod)

# modelos_nt2 <- modelos_nt[!(row.names(modelos_nt) %in% "ROD1"),]

########## Fungos #############
# Fungi_species_div ~ PC1
Fungi_div_PC1 <- nls(Fungi_species_div ~ NLS.expoDecay(PC1, a, k), data = modelos_nt)

# Fungi_species_div ~ alt
fungi_div_alt <- nls(Fungi_species_div ~ NLS.expoGrowth(Altitude, a, k),
             data = modelos_nt)

# Fungi_species_div ~ LUI_subasin
fungi_div_subasin <- nls(Fungi_species_div ~ NLS.expoGrowth(LUI_subasin, a, k), data = modelos_nt)

# Fungi_species_div ~ DOmin
modelos_nt_retirado <- modelos_nt[!(row.names(modelos_nt) %in% "SEL1"),]
fungi_div_DOmin <- nls(Fungi_species_div ~ NLS.expoGrowth(DOmin, a, k), data = modelos_nt_retirado)

# Fungi_species_div ~ cond
fung_div_cond <- nls(Fungi_species_div ~ NLS.expoDecay(cond, a, k), data = modelos_nt)

# Fungi_species_div ~ DIN
fung_div_din <- nls(Fungi_species_div ~ NLS.lorentz.3(DIN, b, d, e), data = modelos_nt)
fung_div_din <- nls(Fungi_species_div ~ SSgauss(DIN, mu, sigma, h), data = modelos_nt)
plot(Fungi_species_div ~ DIN, data = modelos_nt)
par(mfrow = c(1, 1))
plot_nls(fung_div_din)
r2 <- bquote(paste("R"^2 == .(format(R2nls(fung_div_din)$PseudoR2, digits = 4))))
pval <- bquote(paste(bold("p-value: " == .(format(summary(fung_div_din)$coefficients[2,4], digits = 5)))))
mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

par(mfrow = c(2, 2))
plot(nlsResiduals(fung_div_din), which = 0)

# Fungi_species_div ~ P.PO4
fung_div_po4 <- nls(Fungi_species_div ~ NLS.bragg.3(P.PO4, b, d, e), data = modelos_nt)
fung_div_po4 <- nls(Fungi_species_div ~ SSgauss(P.PO4, mu, sigma, h), data = modelos_nt)
plot(Fungi_species_div ~ P.PO4, data = modelos_nt)
par(mfrow = c(1, 1))
plot_nls(fung_div_po4)
r2 <- bquote(paste("R"^2 == .(format(R2nls(fung_div_po4)$PseudoR2, digits = 4))))
pval <- bquote(paste(bold("p-value: " == .(format(summary(fung_div_po4)$coefficients[2,4], digits = 5)))))
mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

par(mfrow = c(2, 2))
plot(nlsResiduals(fung_div_po4), which = 0)

# Fungi_species_div ~ Tmax
fung_div_tmax <- nls(Fungi_species_div ~ NLS.bragg.3(Tmax, b, d, e), data = modelos_nt)
fung_div_tmax <- nls(Fungi_species_div ~ SSgauss(Tmax, mu, sigma, h), data = modelos_nt)


par(mfrow = c(1, 1))
plot_nls(fung_div_tmax)
r2 <- bquote(paste("R"^2 == .(format(R2nls(fung_div_tmax)$PseudoR2, digits = 4))))
pval <- bquote(paste(bold("p-value: " == .(format(summary(fung_div_tmax)$coefficients[2,4], digits = 5)))))
mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

par(mfrow = c(2, 2))
plot(nlsResiduals(fung_div_tmax), which = 0)

# Fungi_species_shannon ~ ihf
fung_shannon_ihf <- nls(Fungi_species_shannon ~ NLS.expoDecay(ihf, a, k), data = modelos)


pdf("modelos_nao_lineares_fungos.pdf")
  # Fungi_species_div ~ PC1
  par(mfrow = c(4, 4))
  plot_nls(Fungi_div_PC1)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(Fungi_div_PC1)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(Fungi_div_PC1)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(Fungi_div_PC1), which = 0)

  # Fungi_species_div ~ alt
  par(mfrow = c(1, 1))
  plot_nls(fungi_div_alt)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_alt)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_alt)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fungi_div_alt), which = 0)

  # Fungi_species_div ~ LUI_subbasin
  par(mfrow = c(1, 1))
  plot_nls(fungi_div_subasin)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_subasin)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_subasin)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fungi_div_subasin), which = 0)

  # Fungi_species_div ~ LUI_100m
  par(mfrow = c(1, 1))
  plot_nls(fungi_div_100m)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_100m)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_100m)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fungi_div_100m), which = 0)

  # Fungi_species_div ~ LUI_500m
  par(mfrow = c(1, 1))
  plot_nls(fungi_div_500m)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_500m)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_500m)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fungi_div_500m), which = 0)

  # Fungi_species_div ~ DOmin
  par(mfrow = c(1, 1))
  plot_nls(fungi_div_DOmin)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_DOmin)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_DOmin)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fungi_div_DOmin), which = 0)

  # Fungi_species_div ~ cond
  par(mfrow = c(1, 1))
  plot_nls(fung_div_cond)
  # Obtain R-squared value
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fung_div_cond)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fung_div_cond)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fung_div_cond), which = 0)

  # Fungi_species_shannon ~ ihf
  par(mfrow = c(1, 1))
  plot_nls(fung_shannon_ihf)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fung_shannon_ihf)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fung_shannon_ihf)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fung_shannon_ihf), which = 0)

dev.off()

########## Bactérias #############
# Bacteria_species_div ~ qbr
bact_div_qbr <- nls(Bacteria_species_div ~ NLS.expoDecay(qbr, a, k), data = modelos)

# Bacteria_species_div ~ alt
bact_div_alt <- nls(Bacteria_species_div ~ SSgauss(Altitude, mu, sigma, h), data = modelos)

# Bacteria_species_div ~ LUI_500m
bact_div_500m <- nls(Bacteria_species_div ~ NLS.expoGrowth(LUI_500m, a, k), data = modelos)

# Bacteria_species_div ~ DOmin
bact_div_DOmin <- nls(Bacteria_species_div ~ NLS.expoGrowth(DOmin, a, k), data = modelos)

# Bacteria_species_div ~ cond
bact_div_cond <- nls(Bacteria_species_div ~ SSgauss(cond, mu, sigma, h), data = modelos)

# Bacteria_species_div ~ DIN
bact_div_din <- nls(Bacteria_species_div ~ SSgauss(DIN, mu, sigma, h), data = modelos)

# Bacteria_species_div ~ P.PO4
modelos_nt_retirado <- modelos_nt[!(row.names(modelos_nt) %in% "SEL1"),]
# bact_div_po4 <- nls(Bacteria_species_div ~ NLS.bragg.3(P.PO4, b, d, e), data = modelos)
bact_div_po4 <- nls(Bacteria_species_div ~ SSgauss(P.PO4, mu, sigma, h), data = modelos_nt_retirado)

# Bacteria_species_div ~ Tmax
# bact_div_tmax <- nls(Bacteria_species_div ~ NLS.bragg.3(Tmax, b, d, e), data = modelos)
bact_div_tmax <- nls(Bacteria_species_div ~ SSgauss(Tmax, mu, sigma, h), data = modelos)
# bact_div_tmax <- nls(Bacteria_species_div ~ SSquadp(Tmax, a, b, c, xs),
#                      data = modelos)

# Bacteria_species_shannon ~ Discharge
bact_shannon_discharge <- nls(Bacteria_species_shannon ~ NLS.expoDecay(Discharge, a, k), data = modelos)

# Bacteria_species_shannon ~ Tmax
bact_shannon_tmax <- nls(Bacteria_species_shannon ~ NLS.expoDecay(Tmax, a, k), data = modelos)

# Bacteria_species_shannon ~ subasin
bact_shannon_subasin <- nls(Bacteria_species_shannon ~ SSgauss(subasin, mu, sigma, h), data = modelos)

# Bacteria_species_shannon ~ mean.Velocity
bact_shannon_tmax <- nls(Bacteria_species_shannon ~ SSgauss(Tmax, mu, sigma, h), data = modelos)


pdf("modelos_nao_lineares_bacterias.pdf")
  # Bacteria_species_div ~ qbr
  par(mfrow = c(4, 4))
  plot_nls(bact_div_qbr)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_qbr)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_qbr)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_qbr), which = 0)

  # Bacteria_species_div ~ alt
  par(mfrow = c(1, 1))
  plot_nls(bact_div_alt)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_alt)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_alt)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_alt), which = 0)

  # Bacteria_species_div ~ LUI_500m
  par(mfrow = c(1, 1))
  plot_nls(bact_div_500m)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_500m)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_500m)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_500m), which = 0)

  # Bacteria_species_div ~ DOmin
  par(mfrow = c(1, 1))
  plot_nls(bact_div_DOmin)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_DOmin)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_DOmin)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_DOmin), which = 0)

  # Bacteria_species_div ~ cond
  par(mfrow = c(1, 1))
  plot_nls(bact_div_cond)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_cond)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_cond)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_cond), which = 0)

  # Bacteria_species_div ~ DIN
  par(mfrow = c(1, 1))
  plot_nls(bact_div_din)
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_din)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_din)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_din), which = 0)

  # Bacteria_species_div ~ P.PO4
  par(mfrow = c(1, 1))
  plot_nls(bact_div_po4)
  # Obtain R-squared value
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_po4)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_po4)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_po4), which = 0)

  # Bacteria_species_div ~ Tmax
  par(mfrow = c(1, 1))
  plot_nls(bact_div_tmax)
  # Obtain R-squared value
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_div_tmax)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_div_tmax)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_div_tmax), which = 0)

  # Bacteria_species_shannon ~ Discharge
  par(mfrow = c(1, 1))
  plot_nls(bact_shannon_discharge)
  # Obtain R-squared value
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_shannon_discharge)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_shannon_discharge)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_shannon_discharge), which = 0)

  # Bacteria_species_shannon ~ Tmax
  par(mfrow = c(1, 1))
  plot_nls(bact_shannon_tmax)
  # Obtain R-squared value
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_shannon_tmax)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_shannon_tmax)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_shannon_tmax), which = 0)

  # Bacteria_species_shannon ~ subasin
  par(mfrow = c(1, 1))
  plot_nls(bact_shannon_subasin)
  # Obtain R-squared value
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_shannon_subasin)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_shannon_subasin)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_shannon_subasin), which = 0)


  # Bacteria_species_shannon ~ mean.Velocity
  par(mfrow = c(1, 1))
  plot_nls(bact_shannon_tmax)
  # Obtain R-squared value
  r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_shannon_tmax)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_shannon_tmax)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(bact_shannon_tmax), which = 0)



dev.off()

# nls fit
Fung_div_qbr <- nls(Fungi_species_div ~ NLS.expoGrowth(qbr, a, k),
             data = modelos_nt)
summary(Fung_div_qbr)
plot_nls(Fung_div_qbr)
R2nls(Fung_div_qbr)$PseudoR2
residuos_expodecay <- nlsResiduals(Fung_div_qbr)
par(mfrow = c(2, 2))
plot(residuos_expodecay, which = 0)


pdf("modelos_lineares_hill.pdf")
for (i in 1:length(hill_shannon)) {
  for(j in colnames(env_data)) {
    dados <- as.data.frame(cbind(hill_shannon[[i]], env_data[[j]]))
    #print(dados)
    par(mfrow=c(1,1))
    plot(hill_shannon[[i]] ~ env_data[[j]], data = dados,
         xlab = colnames(env_data[j]), ylab = names(hill_shannon)[i])
    mod <- lm(hill_shannon[[i]] ~ env_data[[j]], data = dados)
    r2 <- bquote(paste(bold(R^2 == .(summary(mod)$r.squared))))
    x <- as.numeric(summary(mod)$coefficients[2,4])
    pval <- bquote(paste(bold("p-value: " == .(round(x, 4)))))
    mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
    mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)
    abline(mod)
    par(mfrow=c(2,2))
    plot(mod)
  }
}
dev.off()


