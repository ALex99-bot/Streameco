# Set the directory and file name
# directory <- "/home/pedro/PycharmProjects/Streameco"
directory <- "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco"
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
library(xlsx)
library(abdiv)
library(dplyr)
library(ggpmisc)

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

# list2env(bact_fung_taxa, envir=.GlobalEnv)

# nmds <- lapply(bact_fung_taxa, metaMDS, distance = "bray")

# # organise4plot <- function (nmds) {
# #   Site <- as.data.frame(nmds$points)
# #   Site <- as.data.frame(cbind(Site, name="Species"))
# #   Species <- as.data.frame(nmds$species)
# #   Species <- as.data.frame(cbind(Species, name="Site"))
# #   df <- rbind(Site, Species)
# # }
#
# prepare_nmds_data <- function(nmds) {
#   if (!is.numeric(nmds$points) || !is.numeric(nmds$species)) {
#     stop("Input must be a list with 'points' and 'species' as numeric data frames")
#   }
#
#   Site <- data.frame(nmds$points, name="Species")
#   Species <- data.frame(nmds$species, name="Site")
#   df <- rbind(Site, Species)
#   return(df)
# }
#
# # all_nmds <- lapply(nmds, prepare_nmds_data)
# all_nmds <- sapply(nmds, FUN = prepare_nmds_data, simplify = FALSE, USE.NAMES = TRUE)
#
# plot_text = function (data, df_name) {
#     ggplot(data, aes(x=MDS1, y=MDS2)) +
#     geom_text_repel(aes(color=name, label=row.names(data)))+labs(col="Legend")+
#     guides(fill="name") +
#     theme_gray() +
#     labs(title=df_name)+
#     scale_color_manual(labels = c("Site", strsplit(df_name, split = "_")[[1]][2]), values = c("red", "blue"))
# }
#
# scatter_plot = function (data, df_name) {
#   ggplot() +
#   geom_point(data = data, aes(x = MDS1, y = MDS2, color = name)) +
#   labs(col = "Legend") +
#   guides(fill = "name") +
#   theme_gray() +
#   labs(title = df_name) +
#   scale_color_manual(labels = c("Site", strsplit(df_name, split = "_")[[1]][2]), values = c("red", "blue"))
# }
#
# myplots_text <- lapply(names(all_nmds), function(df_name) {
#   plot_text(all_nmds[[df_name]], df_name)
# })
#
# myplots_scatter <- lapply(names(all_nmds), function(df_name) {
#   scatter_plot(all_nmds[[df_name]], df_name)
# })
#
# pdf("top10_nomes.pdf")
# for (i in 1:length(myplots_text)) {
#   print(myplots_text[[i]])
# }
# dev.off()
#
# pdf("top10_pontos.pdf")
# for (i in 1:length(myplots_scatter)) {
#   print(myplots_scatter[[i]])
# }
# dev.off()

# Índice de Shannon
# shannon_by_column <- function(data, name) {
#   # Create an empty data frame to store the results
#   shannon_df <- data.frame(shannon = numeric(), row.names = character(), stringsAsFactors = FALSE)
#   for(i in 1:ncol(data)) {       # for-loop over columns
#   shannon <- data.frame(shannon = diversity(t(data[,i]), index = "shannon", MARGIN = 1, base = exp(1)),
#                          row.names = colnames(data)[i])
#   # Add the current result to the shannon_df data frame
#   shannon_df <- rbind(shannon_df, shannon)
#   }
#   # Name of dataframe
#   name_col <- paste("shannon_index", sub("_([^_]*)$", "", name), sep="_")
#   colnames(shannon_df)[1] <- name_col
#   return(shannon_df)
# }
#
#
# shannon <- lapply(names(bact_fung_taxa), function(df_name) {
#   shannon_by_column(bact_fung_taxa[[df_name]], df_name)
# })
#
# # Índice de Pielou
# pielou_by_column <- function(data, name) {
#   # Create an empty data frame to store the results
#   pielou_df <- data.frame(pielou = numeric(), row.names = character(), stringsAsFactors = FALSE)
#   for(i in 1:ncol(data)) {       # for-loop over columns
#   pielou <- data.frame(pielou = diversity(t(data[,i]))/log(specnumber(t(data[,i]))),
#                          row.names = colnames(data)[i])
#   # Add the current result to the shannon_df data frame
#   pielou_df <- rbind(pielou_df, pielou)
#   }
#   # Name of dataframe
#   name_col <- paste("pielou_index", sub("_([^_]*)$", "", name), sep="_")
#   colnames(pielou_df)[1] <- name_col
#   return(pielou_df)
# }
#
# pielou <- lapply(names(bact_fung_taxa), function(df_name) {
#   pielou_by_column(bact_fung_taxa[[df_name]], df_name)
# })
#
# # Bray-Curtis
# bray <- lapply(bact_fung_taxa, function(x) vegdist(t(x)))
#
#
# # Margalef
# margalef_by_column <- function(data, name) {
#   # Create an empty data frame to store the results
#   margalef_df <- data.frame(margalef = numeric(), row.names = character(), stringsAsFactors = FALSE)
#   for(i in 1:ncol(data)) {       # for-loop over columns
#   margalef <- data.frame(margalef = margalef(data[, i]),
#                          row.names = colnames(data)[i])
#   # Add the current result to the shannon_df data frame
#   margalef_df <- rbind(margalef_df, margalef)
#   }
#   # Name of dataframe
#   name_col <- paste("margalef_index", sub("_([^_]*)$", "", name), sep="_")
#   colnames(margalef_df)[1] <- name_col
#   return(margalef_df)
# }
#
# margalef_ind <- lapply(names(bact_fung_taxa), function(df_name) {
#   margalef_by_column(bact_fung_taxa[[df_name]], df_name)
# })
#
# # Simpson
# simpson_by_column <- function(data, name) {
#   # Create an empty data frame to store the results
#   simpson_df <- data.frame(simpson = numeric(), row.names = character(), stringsAsFactors = FALSE)
#   for(i in 1:ncol(data)) {       # for-loop over columns
#   simpson <- data.frame(simpson = diversity(t(data[,i]), index = "simpson", MARGIN = 1, base = exp(1)),
#                          row.names = colnames(data)[i])
#   # Add the current result to the shannon_df data frame
#   simpson_df <- rbind(simpson_df, simpson)
#   }
#   # Name of dataframe
#   name_col <- paste("simpson_index", sub("_([^_]*)$", "", name), sep="_")
#   colnames(simpson_df)[1] <- name_col
#   return(simpson_df)
# }
#
#
# simpson_ind <- lapply(names(bact_fung_taxa), function(df_name) {
#   simpson_by_column(bact_fung_taxa[[df_name]], df_name)
# })
#
# # # Guardar excel dos índices
# wb <- createWorkbook()
# #
# sheet <- createSheet(wb, "shannon")
#
# addDataFrame(shannon_index, sheet=sheet,
#              startColumn=1, row.names=FALSE)
#
# sheet <- createSheet(wb, "pielou")
#
# addDataFrame(pielou_index, sheet=sheet,
#              startColumn=1, row.names=FALSE)
#
# sheet <- createSheet(wb, "margalef")
#
# addDataFrame(margalef_index, sheet=sheet,
#              startColumn=1, row.names=FALSE)
#
# sheet <- createSheet(wb, "simpson")
#
# addDataFrame(simpson_index, sheet=sheet,
#              startColumn=1, row.names=FALSE)
#
# sheet <- createSheet(wb,"div")
#
# addDataFrame(diversidade, sheet=sheet,
#              startColumn=1, row.names=TRUE)
#
# saveWorkbook(wb, "bioindices2.xlsx")

# Variáveis ambientais
bacias <- read_excel("STREAMECO database - environment (copy).xlsx", "LandUse")
bacias <- as.data.frame(bacias)
river_basin <- bacias[2:nrow(bacias), 3]
river_basin <- as.data.frame(river_basin)
row.names(river_basin) <- bacias[2:51, 1]

# ggplot(bacias, aes(lon, lat)) +
#   geom_point() +
#   geom_text_repel(aes(label = code))

env.data <- read.table("var ambientais.txt", sep="\t", dec=".", header=T)
row.names(env.data) <- env.data$code
env.data <- env.data[, -1]
env_data <- as.data.frame(env.data)
to_remove <- c("mean_light", "Artificial.500m.", "Agriculture.500m.", "Pasture.500m.", "Natural.500m.", "N.NH4", "N.NO2",
               "N.NO3", "Tmean", "Tmin", "TCV", "DOmean", "Domax", "Artificial.100m.", "Agriculture.100m.", "Pasture.100m.",
               "Natural.100m.", "Artificial_subasin", "Agriculture_subasin", "Pasture_subasin", "Natural_subasin")
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

river_basin$alt <- cut(env_data_nt$Altitude,
              breaks=c(0, 200, 400, 600, 800, 1000),
              labels=c('<200', '200-400', '400-600', '600-800', ">800"))

river_basin$river_basin <- as.factor(river_basin$river_basin)
river_basin$alt <- as.factor(river_basin$alt)

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

pielou_index <- lapply(bact_fung_taxa, function(x) diversity(t(x))/log(specnumber(t(x))))

margalef_index <- lapply(bact_fung_taxa, function(x) apply(x, 2, margalef))

simpson_index <- lapply(bact_fung_taxa, function(x) diversity(t(x), index = "simpson", MARGIN = 1, base = exp(1)))

bray_index <- lapply(bact_fung_taxa, function(x) {pco <- dudi.pco(vegdist(t(x), method = "bray"), scannf = F, nf = 10)
x <- return(pco$li$A1)
row.names(x) <- col.names(bact_fung_taxa)})


# Para hifomicetes, mudar times=1
#env_data <- scale(env_data)
#env_data_nt <- scale(env_data_nt)
#PCA_env <- prcomp(env_data, scale = TRUE)
PCA_env <- dudi.pca(env_data, center = TRUE, scale = TRUE, nf=5, scannf = FALSE)


# PCA_env <- rep(list(dudi.pca(env_data, center = TRUE, scale = TRUE, nf=5, scannf = FALSE)), times = 12)
biplot(PCA_env)
s.arrow(PCA_env$c1, lab=names(PCA_env$tab))
var_exp <- (PCA_env$eig*100)/sum(PCA_env$eig)

env_data <- cbind(env_data, PCA_env$li[,1])
env_data_nt <- cbind(env_data_nt, PCA_env$li[,1])

# PCA_dataframe <- map2(diversidade, PCA_env, function(ind, env) {
#   data.frame(bioindex = ind, environmental_data = env$li[,1])
# })

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

modelos_nt <- as.data.frame(cbind(indices2, env_data_nt))

modelos_nt<-cbind(modelos_nt, vel = modelos_nt$mean.Velocity)

jpeg("bact_shannon_velocidade.jpg")
par(mfrow=c(1,1))
plot(Bacteria_species_shannon~vel, data = modelos_nt, xlab = "Mean velocity (m2/s)")
r2 <- bquote(paste(bold(R^2 == .(11.6))))
mtext("p<0.01", line=-1.5, adj = 1, cex = 1.2, font = 2)
mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
mod <- lm(Bacteria_species_shannon~vel, data = modelos_nt)
abline(mod)
dev.off()
summary(mod)


#jpeg("bact_riqueza_QBR.jpg")
par(mfrow=c(1,1))
plot(Bacteria_species_div~qbr, data = modelos_nt, xlab = "QBR", ylab = "Bacteria species richness")
# r2 <- bquote(paste(bold(R^2 == .(11.6))))
# mtext("p<0.01", line=-1.5, adj = 1, cex = 1.2, font = 2)
# mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
# mod2 <- glm(Bacteria_species_div~I(1/qbr), data = modelos_nt, family = gaussian("inverse"))
# mod2 <- nls(Bacteria_species_div~Vm*qbr/(K+qbr), data=modelos_nt,start=list(K=1, Vm=95), trace = TRUE)
# mod2 <- nls(Bacteria_species_div~SSlogis(log(qbr+0.001)), data=modelos_nt)
mod2 <- loess(Bacteria_species_div~qbr, data=modelos_nt)
smooth1 <- predict(mod2, data.frame(qbr = seq(0.1,95)))
line(smooth1)
summary(mod2)

par(mfrow=c(1,1))

jpeg("fung_riqueza_lui100.jpg")
plot(Fungi_species_div~LUI_100m, data = modelos_nt,  ylab = "Fungi species richness")
# r2 <- bquote(paste(bold(R^2 == .(9.4))))
# mtext("p<0.05", line=-1.5, adj = 1, cex = 1.2, font = 2)
# mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
mod <- lm(Fungi_species_div~LUI_500m, data = modelos_nt)
abline(mod)
dev.off()
summary(mod)

jpeg("bact_riqueza_LUI500.jpg")
par(mfrow=c(2,2))
plot(mod)
dev.off()
# jpeg("pca_env.jpg")
# fviz_eig(PCA_env)
# fviz_pca_ind(PCA_env,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
#              )
#
# fviz_pca_var(PCA_env,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
#              )
#
# fviz_pca_biplot(PCA_env, repel = TRUE,
#                 col.var = "#2E9FDF", # Variables color
#                 col.ind = "#696969"  # Individuals color
#                 )

# dev.off()


# all_graphs <- lapply(seq_along(PCA_dataframe), function(i) {
#   all_dataframe <- PCA_dataframe[[i]]
#   ggplot(all_dataframe, aes(x=environmental_data, y=bioindex)) +
#     geom_point() +
#     # geom_text(aes(label = rownames(all_dataframe))) +
#     #geom_text(aes(label = rownames(all_dataframe)), size = 2, nudge_x = 1, nudge_y = 0.03) +
#     ggtitle(names(PCA_dataframe)[i]) +  # add the dataframe name as the plot title
#     labs(x="Environmental variables", y="Diversity index") +
#     theme(legend.position = "none")
# })

# pdf("PCA_div.pdf")
# for (i in 1:length(all_graphs)) {
#   print(all_graphs[[i]])
# }
# dev.off()

pdf("modelos_lineares.pdf")
for (i in 1:ncol(modelos_nt)) {
  for(j in 1:ncol(indices2)) {
    dados = as.data.frame(cbind(env_data_nt[i], indices2[j]))
    plot(Fungi_species_div~LUI_100m, data = modelos_nt,  ylab = "Fungi species richness")
    # r2 <- bquote(paste(bold(R^2 == .(9.4))))
    # mtext("p<0.05", line=-1.5, adj = 1, cex = 1.2, font = 2)
    # mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
    mod <- lm(Fungi_species_div~LUI_500m, data = modelos_nt)
    abline(mod)

  }
}

dev.off()


# div <- as.data.frame(scale(div))
# env.data <- env.data[, c(1, 2, 4, 5, 6, 15, 17, 22, 26, 28, 29, 30, 31)]
# env.data2 <- as.data.frame(scale(env.data))
# dudi.pca(env.data2, center = TRUE, scale = FALSE, scannf = F,nf=4)->metric.pca
# biplot(metric.pca)
# metric.pca$li$Axis1 -> PCA1
# metric.pca$li$Axis2 -> PCA2
# # PCA Axis importance (explained variance)
# round(metric.pca$eig[1:4]/sum(metric.pca$eig),2)
# env.data <- cbind(env.data, PCA1, PCA2)
# env.data <- as.data.frame(scale(env.data))

# pdf("Gam_models.pdf")
# for (i in 1:length(env.data)) {
#   for(j in 1:length(div)) {
#     ##plot(div2[,y]~e.data[,i], xlab=colnames(e.data[,i]), ylab=colnames(div2[,y]))
#     #   }}
#     print(i)
#     gam<-gam(div[,j]~s(env.data[,i], fx = FALSE, k=-1,  bs = "cr"))
#     x<-summary(gam)
#     p <- as.character(round(as.numeric(unlist((x$s.table[, 4]))), 3))
#     if (p<=0.05) {
#       plot(gam, main='Diversity', xlab=colnames(env.data[i]), ylab=colnames(div[j]), se=TRUE)
#       p <- if (p<=0.001) {as.character("<0.001***")} else if (p<=0.01) {as.character("<0.01**")
#       }else if (p<=0.05) {as.character("<0.05*")} else {as.character("=n.s.")}
#       p <- paste("p", p)
#       r2 <-bquote(paste(bold(r^2 == .(as.character(round(as.numeric(x$r.sq),3))))))
#       dev <- as.character(round(as.numeric(x$dev.expl), 3))
#       dev <-(paste("dev", "=",dev))
#       mtext(p, line=-1.5, adj = 1, cex = 1.2, font = 2)
#       mtext(r2, line=-2.5, adj = 1, cex = 1.2, font = 2)
#       mtext(dev, line=-3.5, adj = 1, cex = 1.2, font = 2)
#       par(mfrow=c(2,2))
#       gam.check (gam)
#       par(mfrow=c(1,1))
#     }
#   }
# }
#
# biplot(metric.pca)

# dev.off()

# row.names(fun_order_reads) <- fun_order_reads[, 1]
# fun_order_reads <- fun_order_reads[, -1]
# fun_order_reads <- as.data.frame(fun_order_reads)
#
# nmds <- metaMDS(fun_order_reads, distance = "bray")
# nmds

# plot(nmds)
# stressplot(nmds)
# plot(nmds, type = "n")+ #displays empty ordination space
#   ggrepel::geom_text_repel(data = nmds, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)
# points(nmds, display = "sites", col = "blue", labels(row.names(nmds$points)))
# # displays site points where symbols (pch) are different management options and colour (col) are different land uses
# legend("topright", legend = c(levels(dune.env$Management), levels(dune.env$Use)), pch = c(16, 8, 17, 11, 16, 16, 16), col = c("black","black","black","black","blue", "orange", "black"), bty = "n", cex = 1) # displays symbol and colour legend
# legend("topleft", "stress = 0.118", bty = "n", cex = 1) # displays legend text of stress value
# ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)
#
# plot(div$shannon_bact~env.data$mean.Velocity)

plot(indices$Bacteria_phylum_pielou~env_data_nt$V20)
mod<-lm(indices$Bacteria_phylum_pielou~env_data_nt$V20)
abline(mod)
summary(mod)

plot(data$V1~data$V2,data = data )
mod <- lm(data$V1~data$V2,data = data )
abline(mod)
summary(mod)
par(mfrow=c(2,2))
plot(mod)

plot(data2$V1~data2$V2, data = data2)
mod <- lm(data2$V1~data2$V2, data = data2)
abline(mod)
summary(mod)

data <- as.data.frame(cbind(indices$Bacteria_phylum_pielou, env_data_nt$`PCA_env$li[, 1]`))
data2 <- data[!row.names(data) %in% "VEZ1",]
row.names(data) <- row.names(indices)

jpeg("Pielou_bacteria_phylum.jpg")

ggplot(data, aes(x=V2, y=V1)) +
   geom_point() +
      geom_text(aes(label = rownames(data)), nudge_y = 0.005) +
   geom_smooth(method='lm') +
     #geom_text(aes(label = rownames(all_dataframe)), size = 2, nudge_x = 1, nudge_y = 0.03) +
     ggtitle(names("Bacteria Phylum")) +  # add the dataframe name as the plot title
     labs(x="Environmental variables", y="Pielou index") +
     theme(legend.position = "none")

dev.off()

