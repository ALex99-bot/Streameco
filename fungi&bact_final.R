# Set the directory and file name
#directory <- "/home/pedro/PycharmProjects/Streameco"
directory <- "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco"
setwd(directory)

file <- "total_reads.xlsx"
file2 <- "top10.xlsx"
file3 <- "menor_cutoff_reads.xlsx"
path <- file.path(directory, file)
path2 <- file.path(directory, file2)
path3 <- file.path(directory, file3)


# Load the required packages
library(readxl)
library(mgcv)
library(ade4)
library(usdm)
library(vegan)
library(ggplot2)
library(ggrepel)
library(factoextra)

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

list2env(bact_fung_taxa, envir=.GlobalEnv)

# nmds <- lapply(new_list, metaMDS, distance = "bray")
#
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


shannon_by_column <- function(data, name) {
  # Create an empty data frame to store the results
  shannon_df <- data.frame(shannon = numeric(), row.names = character(), stringsAsFactors = FALSE)
  for(i in 1:ncol(data)) {       # for-loop over columns
  shannon <- data.frame(shannon = diversity(t(data[,i]), index = "shannon", MARGIN = 1, base = exp(1)),
                         row.names = colnames(data)[i])
  # Add the current result to the shannon_df data frame
  shannon_df <- rbind(shannon_df, shannon)
  }
  # Name of dataframe
  name_col <- paste("shannon_index", sub("_([^_]*)$", "", name), sep="_")
  colnames(shannon_df)[1] <- name_col
  return(shannon_df)
}

# Índice de Shannnon
shannon <- lapply(names(bact_fung_taxa), function(df_name) {
  shannon_by_column(bact_fung_taxa[[df_name]], df_name)
})
# shannon_bact <- diversity(t(bact_reads), index = "shannon", MARGIN = 1, base = exp(1))
# shannon_fung <- diversity(t(fung_reads), index = "shannon", MARGIN = 1, base = exp(1))


# bray curtis 1st axis of PCoA
# Functional space
bray_curtis <- lapply(t(bact_fung_taxa), vegdist, method = "bray")
# _bact <- vegdist(t(bact_reads), method = "bray")
# dudi.pco(bray_curtis_bact,scannf = F,nf=10)-> tr.pco
#
# Bray_Curtis_bact <- tr.pco$li$A1
#
# bray_curtis_fung <- vegdist(t(fung_reads), method = "bray")
# dudi.pco(bray_curtis_fung,scannf = F,nf=10) -> tr.pco2
# Bray_Curtis_fung<-tr.pco2$li$A1
#
# div <- cbind(div, shannon_bact, Bray_Curtis_bact, shannon_fung, Bray_Curtis_fung)
# div <- as.data.frame(div)
#
#
env.data <- read.table("var ambientais.txt", sep="\t", dec=".", header=T)
row.names(env.data) <- env.data$code
env.data <- env.data[, -1]
View(env.data)

env.data <- as.data.frame(env.data)
env.data2 <- env.data[, c(12, 25, 32, 37)]
env.data <- env.data[, -c(12, 25, 32, 37)]
row.names(env.data2) <- row.names(env.data)

PCA_teste1 <- prcomp(env.data, scale = FALSE)
PCA_teste3 <- prcomp(t(Bacteria_species_reads))
a <- data.frame(especies_bact = PCA_teste3$x[,1], abiotico = PCA_teste1$x[,1])
pdf("PCA_espécies_bactérias.pdf")

fviz_eig(PCA_teste3)

fviz_pca_ind(PCA_teste3,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(PCA_teste3,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

ggplot(a, aes(x=abiotico, y=especies_bact)) + geom_point()
dev.off()


fviz_eig(PCA_teste1)
fviz_pca_ind(PCA_teste1,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(PCA_teste1,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

PCA_teste2 <- princomp(env.data, cor = FALSE, scores = TRUE)

summary(PCA_teste3)




matriz_final <- Reduce(function(df1, df2) cbind(df1, df2), shannon, init = env.data)


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

