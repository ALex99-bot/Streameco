# directory <- "/home/pedro/PycharmProjects/Streameco"
directory <- "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco"
# directory <-"C:/Users/asus/Desktop/Streameco"
setwd(directory)

file <- "total_reads.xlsx"
top10 <- "top10.xlsx"
less10 <- "top10inv.xlsx"
file3 <- "menor_cutoff_reads.xlsx"
hypho_file <- "hyphomycetes_aquaticos.xlsx"
path <- file.path(directory, file)
path2 <- file.path(directory, top10)
path3 <- file.path(directory, file3)
hypho_path <- file.path(directory, hypho_file)
less10_path <- file.path(directory, less10)

# Load the required packages
library(readxl)
library(mgcv)
library(ade4)
library(usdm)
library(vegan)
library(ggrepel)
library(factoextra)
library(readxl)
library(abdiv)
library(plyr)
library(ggpmisc)
library(aomisc)
library(nlstools)
library(nlshelper)
library(xcms)
library(hillR)
library(nlraa)
library(tidyverse)
library(gridExtra)

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
# bact_fung_taxa_less10 <- lapply(read_excel_sheets(less10_path, less10), excel2dataframe)

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
               "Natural.100m.", "Artificial_subasin", "Agriculture_subasin", "Pasture_subasin", "Natural_subasin", "LUI_100m", "LUI_subasin")
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

# Índices
diversidade <- lapply(bact_fung_taxa[-c(2, 5, 6, 8, 11, 12)], function(x) specnumber(t(x)))

shannon_index <- lapply(bact_fung_taxa[-c(2, 5, 6, 8, 11, 12)], function(x) diversity(t(x), index = "shannon", base = exp(1)))
hill_shannon <- lapply(bact_fung_taxa[-c(2, 5, 6, 8, 11, 12)], function(x) hill_taxa(t(x), q = 1, MARGIN = 1, base = exp(1)))

pielou_index <- lapply(bact_fung_taxa[-c(2, 5, 6, 8, 11, 12)], function(x) diversity(t(x))/log(specnumber(t(x))))

margalef_index <- lapply(bact_fung_taxa[-c(2, 5, 6, 8, 11, 12)], function(x) apply(x, 2, margalef))

simpson_index <- lapply(bact_fung_taxa[-c(2, 5, 6, 8, 11, 12)], function(x) diversity(t(x), index = "simpson", MARGIN = 1, base = exp(1)))


PCA_env <- dudi.pca(env_data, center = TRUE, scale = TRUE, nf=5, scannf = FALSE)

# PCA_env <- rep(list(dudi.pca(env_data, center = TRUE, scale = TRUE, nf=5, scannf = FALSE)), times = 12)
biplot(PCA_env)
s.arrow(PCA_env$c1, lab=names(PCA_env$tab))
var_exp <- (PCA_env$eig*100)/sum(PCA_env$eig)

env_data <- cbind(env_data, PC1 = PCA_env$li[,1])
env_data_nt <- cbind(env_data_nt, PC1 = PCA_env$li[,1])

all_lists <- list(
  diversidade = diversidade,
  hill_shannon = hill_shannon,
  shannon_index = shannon_index,
  pielou_index = pielou_index,
  margalef_index = margalef_index,
  simpson_index = simpson_index
)


indices <- lapply(names(all_lists), function(list_name) {
  lista <- all_lists[[list_name]]
  repla <- gsub("_?index$", "", list_name)
  names(lista) <- lapply(names(lista), function(df) gsub("reads", repla, df))
  as.data.frame(lista)
})

indices <- as.data.frame(do.call(cbind, indices))

env_data_nt <- as.data.frame(env_data_nt)
env_data <- as.data.frame(env_data)

indices2 <- indices[,c("Bacteria_species_diversidade","Bacteria_species_shannon",
          "Fungi_species_diversidade", "Fungi_species_shannon",
          "Bacteria_genus_diversidade","Bacteria_genus_shannon",
          "Fungi_genus_diversidade", "Fungi_genus_shannon",
          "Bacteria_family_diversidade","Bacteria_family_shannon",
          "Fungi_family_diversidade", "Fungi_family_shannon")]

modelos <- as.data.frame(cbind(indices2, env_data))
modelos <- cbind(modelos, vel = modelos$mean.Velocity)

modelos_nt <- as.data.frame(cbind(indices2, env_data_nt))
modelos_nt <- cbind(modelos_nt, vel = modelos_nt$mean.Velocity)

########## Fungos #############
# Fungi_species_div ~ PC1
Fungi_div_PC1 <- nls(Fungi_species_diversidade ~ NLS.expoDecay(PC1, a, k), data = modelos)

# Fungi_species_div ~ alt
fungi_div_alt <- nls(Fungi_species_diversidade ~ NLS.expoGrowth(Altitude, a, k),
              data = modelos)

# Fungi_species_div ~ LUI_500m
fungi_div_500m <- nls(Fungi_species_diversidade ~ NLS.expoDecay(LUI_500m, a, k), data = modelos)

# Fungi_species_div ~ cond
fung_div_cond <- nls(Fungi_species_diversidade ~ NLS.expoDecay(cond, a, k), data = modelos)

# Generate predicted values from the models
modelos$predicted_div_PC1 <- predict(Fungi_div_PC1)
modelos$predicted_div_alt <- predict(fungi_div_alt)
modelos$predicted_div_500m <- predict(fungi_div_500m)
modelos$predicted_div_cond <- predict(fung_div_cond)

# Create individual plots for each model
plot_PC1 <- ggplot(modelos, aes(x = PC1, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_PC1)) +
  labs(x = "PC1", y = "S") # +
  # ggtitle("Fungi_species_div ~ PC1")

plot_alt <- ggplot(modelos, aes(x = Altitude, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_alt)) +
  labs(x = "Altitude (m)", y = "S") # +
  # ggtitle("Fungi_species_div ~ alt")

plot_500m <- ggplot(modelos, aes(x = LUI_500m, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_500m)) +
  labs(x = "LUI", y = "S") # +
  # ggtitle("Fungi_species_div ~ LUI_500m")

plot_cond <- ggplot(modelos, aes(x = cond, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_cond)) +
  labs(x = "Conductivity (\u03BCS/cm)", y = "S") # +
  # ggtitle("Fungi_species_div ~ cond")

# Combine all the plots into a single image
combined_plot <- grid.arrange(plot_PC1, plot_alt, plot_500m, plot_cond, ncol = 2)

# Save the combined plot as an image
ggsave("modelos_fungos.png", combined_plot, width = 10, height = 8, dpi = 300)


# png("modelos_fungos.png")
#   # Fungi_species_div ~ PC1
#   par(mfrow = c(2, 2))
#   plot_nls(Fungi_div_PC1, ylab = "S")
#   r2 <- bquote(paste("R"^2 == .(format(R2nls(Fungi_div_PC1)$PseudoR2, digits = 4))))
#   pval <- bquote(paste(bold("p-value: " == .(format(summary(Fungi_div_PC1)$coefficients[2,4], digits = 5)))))
#   mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
#   mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)
#
#   # par(mfrow = c(2, 2))
#   # plot(nlsResiduals(Fungi_div_PC1), which = 0)
#
#
#   # Fungi_species_div ~ alt
#   # par(mfrow = c(1, 1))
#   plot_nls(fungi_div_alt, ylab = "S")
#   r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_alt)$PseudoR2, digits = 4))))
#   pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_alt)$coefficients[2,4], digits = 5)))))
#   mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
#   mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)
#
#   # par(mfrow = c(2, 2))
#   # plot(nlsResiduals(fungi_div_alt), which = 0)
#
  # Fungi_species_div ~ LUI_500m
  # par(mfrow = c(1, 1))
  plot_nls(fungi_div_500m, ylab = "S")
  r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_500m)$PseudoR2, digits = 4))))
  pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_500m)$coefficients[2,4], digits = 5)))))
  mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)

  par(mfrow = c(2, 2))
  plot(nlsResiduals(fungi_div_500m), which = 0)
#
#   # Fungi_species_div ~ cond
#   # par(mfrow = c(1, 2))
#   plot_nls(fung_div_cond, ylab = "S")
#   r2 <- bquote(paste("R"^2 == .(format(R2nls(fung_div_cond)$PseudoR2, digits = 4))))
#   pval <- bquote(paste(bold("p-value: " == .(format(summary(fung_div_cond)$coefficients[2,4], digits = 5)))))
#   mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
#   mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)
#
#   # par(mfrow = c(2, 2))
#   # plot(nlsResiduals(fung_div_cond), which = 0)
#
# dev.off()


########## Bactérias #############
# Bacteria_species_shannon ~ mean.Velocity
bact_shannon_mv <- nls(Bacteria_species_shannon ~ SSgauss(mean.Velocity, mu, sigma, h), data = modelos)

# modelos_retirado |>
#   ggplot(aes(mean.Velocity, Bacteria_species_shannon)) +
#   geom_point() +
#   geom_text_repel(aes(label = rownames(modelos_retirado)))

# Create the ggplot
bacterias_modelo <- ggplot(modelos, aes(mean.Velocity, Bacteria_species_shannon)) +
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ SSgauss(x, mu, sigma, h), method.args = list(start = coef(bact_shannon_mv)), se = FALSE) +
  labs(x = "Mean velocity (m s^-1)", y = "H'") #+
  #theme_minimal()

ggsave("modelo_bacterias.png", bacterias_modelo, width = 8, height = 6, dpi = 300)


# modelos_retirado |>
#   ggplot(aes(mean.Velocity, Bacteria_species_shannon)) +
#   geom_point() +
#   geom_text_repel(aes(label = rownames(modelos_retirado)))
#
# par(mfrow = c(1, 1))
# plot_nls(bact_shannon_mv)
# r2 <- bquote(paste("R"^2 == .(format(R2nls(bact_shannon_mv)$PseudoR2, digits = 4))))
# pval <- bquote(paste(bold("p-value: " == .(format(summary(bact_shannon_mv)$coefficients[2,4], digits = 5)))))
# mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
# mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)
#
# par(mfrow = c(2, 2))
# plot(nlsResiduals(bact_shannon_mv), which = 0)


# NMDS
nmds <- lapply(bact_fung_taxa, metaMDS, distance = "bray")

scores(nmds)



