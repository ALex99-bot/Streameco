# directory <- "/home/pedro/PycharmProjects/Streameco"
directory <- "C:/Users/pedro/OneDrive/Ambiente de Trabalho/Streameco"
# directory <-"C:/Users/asus/Desktop/Streameco"
setwd(directory)

file <- "total_reads.xlsx"
top10 <- "top10.xlsx"
less10 <- "top10inv.xlsx"
file3 <- "menor_cutoff_reads.xlsx"
hypho_file <- "hyphomycetes_aquaticos.xlsx"
sem_cutoff <- "todos.xlsx"
top10_percent <- "top10_percent.xlsx"
path <- file.path(directory, file)
path2 <- file.path(directory, top10)
path3 <- file.path(directory, file3)
hypho_path <- file.path(directory, hypho_file)
less10_path <- file.path(directory, less10)
sem_cutoff_path <- file.path(directory, sem_cutoff)
top10_percent_path <- file.path(directory, top10_percent)

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
library(tools)
library(scales)
library(cowplot)


set.seed(2023)

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

excel2dataframe_total <- function(x) {
  x <- as.data.frame(x)
}

bact_fung_taxa_total <- lapply(read_excel_sheets(path, file), excel2dataframe_total)
bact_fung_taxa <- lapply(read_excel_sheets(path, file), excel2dataframe)
# bact_fung_taxa_less10 <- lapply(read_excel_sheets(less10_path, less10), excel2dataframe)
sem_cutoff_taxa <- lapply(read_excel_sheets(sem_cutoff_path, sem_cutoff), excel2dataframe_total)
top10_taxa <- lapply(read_excel_sheets(path2, top10), excel2dataframe_total)
top10_percent_taxa <- lapply(read_excel_sheets(top10_percent_path, top10_percent), excel2dataframe_total)
top10_taxa_d <- lapply(read_excel_sheets(path2, top10), excel2dataframe)

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
png("biplot.png", width = 800, height = 600)
fviz_pca_biplot(PCA_env, title = NULL)
dev.off()


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

# Fungi_species_div ~ DOmin
fung_div_domin <- nls(Fungi_species_diversidade ~ NLS.expoGrowth(DOmin, a, k), data = modelos)

# Fungi_species_div ~ mean.Velocity
fung_div_vel <- nls(Fungi_species_diversidade ~ NLS.expoDecay(mean.Velocity, a, k), data = modelos)

# Fungi_species_div ~ P.PO4
fung_shannon_ppo4 <- nls(Fungi_species_shannon ~ NLS.expoDecay(P.PO4, a, k), data = modelos)

# Fungi_species_div ~ Tmax
fung_div_tmax <- nls(Fungi_species_diversidade ~ NLS.expoDecay(Tmax, a, k), data = modelos)

# Fungi_species_div ~ DIN
fung_div_din <- nls(Fungi_species_diversidade ~ NLS.expoDecay(DIN, a, k), data = modelos)

# Generate predicted values from the models
modelos$predicted_div_PC1 <- predict(Fungi_div_PC1)
modelos$predicted_div_alt <- predict(fungi_div_alt)
modelos$predicted_div_500m <- predict(fungi_div_500m)
modelos$predicted_div_cond <- predict(fung_div_cond)
modelos$predicted_div_domin <- predict(fung_div_domin)
modelos$predicted_div_vel <- predict(fung_div_vel)
modelos$predicted_shannon_ppo4  <- predict(fung_shannon_ppo4)
modelos$predicted_div_tmax <- predict(fung_div_tmax)
modelos$predicted_div_din <- predict(fung_div_din)


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

plot_domin <- ggplot(modelos, aes(x = DOmin, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_domin)) +
  labs(x = expression(paste("Minimum dissolved oxygen (mg ", "L"^-1, ")")), y = "S") # +
  # ggtitle("Fungi_species_div ~ cond")

plot_vel <- ggplot(modelos, aes(x = mean.Velocity, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_vel)) +
  labs(x = expression(paste("Mean velocity (m ", "s"^-1, ")")), y = "S") # +
  # ggtitle("Fungi_species_div ~ cond")

plot_ppo4 <- ggplot(modelos, aes(x = P.PO4, y = Fungi_species_shannon)) +
  geom_point() +
  geom_line(aes(y = predicted_shannon_ppo4)) +
  labs(x = expression(paste("Phosphorus (mg ", "L"^-1, ")")), y = "H'") # +
  # ggtitle("Fungi_species_div ~ cond")

plot_tmax <- ggplot(modelos, aes(x = Tmax, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_tmax)) +
  labs(x = "Thermal stress (\u00B0C)", y = "S") # +
  # ggtitle("Fungi_species_div ~ cond")

plot_din <- ggplot(modelos, aes(x = DIN, y = Fungi_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_din)) +
  labs(x = expression(paste("Nutrient enrichment (DIN) (mg ", "L"^-1, ")")), y = "S") # +
  # ggtitle("Fungi_species_div ~ cond")

# Combine all the plots into a single image
combined_plot <- grid.arrange(plot_PC1, plot_alt, plot_500m, plot_cond, plot_vel, plot_domin, plot_tmax, plot_din, ncol = 2)
# combined_plot2 <- grid.arrange(plot_vel, plot_domin, plot_ppo4, plot_tmax, ncol = 2)

# Save the combined plot as an image
ggsave("modelos_fungos_div.png", combined_plot, width = 10, height = 8, dpi = 300)
# ggsave("modelos_fungos2.png", combined_plot2, width = 10, height = 8, dpi = 300)
ggsave("fung_shannon_ppo4.png", plot_ppo4, width = 10, height = 8, dpi = 300)

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
  # plot_nls(fungi_div_500m, ylab = "S")
  # r2 <- bquote(paste("R"^2 == .(format(R2nls(fungi_div_500m)$PseudoR2, digits = 4))))
  # pval <- bquote(paste(bold("p-value: " == .(format(summary(fungi_div_500m)$coefficients[2,4], digits = 5)))))
  # mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
  # mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)
  #
  # par(mfrow = c(2, 2))
  # plot(nlsResiduals(fungi_div_500m), which = 0)
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

# Bacteria_species_shannon ~ Tmax
bact_div_tmax <- nls(Bacteria_species_diversidade ~ NLS.expoDecay(Tmax, a, k), data = modelos)

modelos$predicted_div_tmax_bact <- predict(bact_div_tmax)
# modelos |>
#   ggplot(aes(Tmax, Bacteria_species_diversidade)) +
#   geom_point()
#
# par(mfrow = c(1, 1))
# plot_nls(fung_shannon_tmax, ylab = "S")
# r2 <- bquote(paste("R"^2 == .(format(R2nls(fung_shannon_tmax)$PseudoR2, digits = 4))))
# pval <- bquote(paste(bold("p-value: " == .(format(summary(fung_shannon_tmax)$coefficients[2,4], digits = 5)))))
# mtext(r2, line=-2.5, adj = 0.9, cex = 1.2, font = 2)
# mtext(pval, line=-3.5, adj = 0.9, cex = 1.2, font = 2)
#
# par(mfrow = c(2, 2))
# plot(nlsResiduals(fung_shannon_tmax), which = 0)

# Create the ggplot
bacterias_vel_shannon <- ggplot(modelos, aes(mean.Velocity, Bacteria_species_shannon)) +
  geom_point() +
  stat_smooth(method = "nls", formula = y ~ SSgauss(x, mu, sigma, h), method.args = list(start = coef(bact_shannon_mv)), se = FALSE) +
  labs(x = expression(paste("Mean velocity (m ", "s"^-1, ")")), y = "H'")
  #theme_minimal()

plot_tmax_bact <- ggplot(modelos, aes(x = Tmax, y = Bacteria_species_diversidade)) +
  geom_point() +
  geom_line(aes(y = predicted_div_tmax_bact)) +
  labs(x = "Thermal stress (\u00B0C)", y = "S")

ggsave("bacterias_vel_shannon.png", bacterias_vel_shannon, width = 8, height = 6, dpi = 300)
ggsave("bact_div_tmax.png", plot_tmax_bact, width = 8, height = 6, dpi = 300)

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
# nmds <- lapply(top10_taxa_d, metaMDS, distance = "bray")
#
# scores(nmds)
# ordiplot(example_NMDS,type="n")
# orditorp(example_NMDS,display="species",col="red",air=0.01)
# orditorp(example_NMDS,display="sites",cex=1.25,air=0.01)

# Create a list of NMDS results for each data frame in bact_fung_taxa
# Set the output file format and name (PNG)
# Change "nmds_plots" to your desired file name (without the extension)

nmds_list <- lapply(top10_taxa_d, function(df) {
  metaMDS(df, distance = "bray")
})

lapply(nmds_list, stressplot)
# result_list <- lapply(nmds_list, function(nmds_data) {
#   envfit_result <- envfit(nmds_data, env_data)
#   return(envfit_result)
# })

# Access and plot the scores for each NMDS result
scores_list <- lapply(nmds_list, scores)

for (x in seq_along(nmds_list)) {
  # Open a new PNG graphics device for each plot
  png(paste0(directory, "nmds_plot_", names(nmds_list)[x], ".png"), width = 800, height = 600)  # Adjust width and height as needed

  ordiplot(nmds_list[[x]])
  orditorp(nmds_list[[x]], display = "species", col = "red", air = 0.01)
  orditorp(nmds_list[[x]], display = "sites", cex = 1.25, air = 0.01)

  # Add a title to each plot (optional)
  # title(main = paste("NMDS Plot", names(nmds_list)[x]))

  # Close the PNG graphics device
  dev.off()
}


# Create empty plots for each NMDS result
plot_list <- lapply(nmds_list, function(nmds) {
  ordiplot(nmds, type = "n")
})

# Add species labels to each plot
plot_list <- lapply(seq_along(nmds_list), function(i) {
  orditorp(nmds_list[[i]], display = "species", col = "red", air = 0.01)
})

# Add site labels to each plot
plot_list <- lapply(seq_along(nmds_list), function(i) {
  orditorp(nmds_list[[i]], display = "sites", cex = 1.25, air = 0.01)
})

# Combine the individual plots into a single composite plot
composite_plot <- plot_grid(plotlist = plot_list, ncol = 2)  # Adjust ncol as needed



# Save the composite plot to an image file (e.g., a PNG file)
ggsave("output_plot.png", composite_plot, width = 10, height = 6)

# You can specify the width and height as per your preference

# Create a mapping of singular to plural forms
plural_mapping <- c("genus" = "genera", "order" = "orders", "species" = "species",
                    "family" = "families", "phylum" = "phyla", "class" = "classes")

# Function to update the names of list elements
update_names_in_list <- function(lst, mapping) {
  new_names <- lapply(names(lst), function(name) {
    parts <- unlist(strsplit(name, "_"))
    singular <- parts[length(parts)]
    plural <- mapping[singular]
    if (!is.na(plural)) {
      parts[length(parts)] <- plural
      return(paste(parts, collapse = "_"))
    } else {
      return(name)
    }
  })

  # Loop through list elements and update the first column name
  for (i in seq_along(lst)) {
    singular <- unlist(strsplit(names(lst)[i], "_"))[length(unlist(strsplit(names(lst)[i], "_")))]
    plural <- mapping[singular]
    if (!is.na(plural)) {
      colnames(lst[[i]])[1] <- plural
    }
  }

  names(lst) <- new_names
  return(lst)
}

# Apply the function to your list
top10_taxa <- update_names_in_list(top10_taxa, plural_mapping)


# Diversidade
tidy_data <- setNames(
  lapply(names(top10_taxa), function(list_name) {
    lista <- top10_taxa[[list_name]]
    gathered_data <- pivot_longer(lista, cols = -1, names_to = "Sample", values_to = "Percent")
    attr(gathered_data, "original_name") <- list_name  # Store the original name as an attribute
    return(gathered_data)
  }),
  names(top10_taxa)
)

barplot_list <- setNames(lapply(names(tidy_data), function(list_name) {
  df <- tidy_data[[list_name]]
  repla <- gsub(".*_", "", list_name)
  df_filtered <- df %>% filter(Percent != 0)
  ggplot(df_filtered, aes(x = Sample, y = Percent, fill = .data[[repla]])) +
    geom_bar(position="fill", stat="identity") +
    labs(x = "Samples", y = "Average percentage of reads", fill = paste(toTitleCase(repla), sep="")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5),
      #legend.position = "bottom",
      legend.text = element_text(size = 8,       # Adjust the font size of the legend text (set to 8 or adjust as needed)
                               hjust = 0.5,     # Center the text horizontally in the legend key
                               vjust = 0.5)) +
    guides(fill = guide_legend(label.hjust = 0.5, label.vjust = 0.5)) +
    scale_y_continuous(labels = percent_format(scale = 100)) # Set y-axis labels as percentages
}), names(tidy_data))

# Print the bar plots
# print(barplot_list[[1]]) # Print the first bar plot
# print(barplot_list[[2]]) # Print the second bar plot


# Save each plot with a unique filename based on the original name of the data frame
for(i in seq_along(barplot_list)) {
  filename <- paste0(directory, attr(tidy_data[[i]], "original_name"), ".png")
  ggsave(filename, plot = barplot_list[[i]], width = 8, height = 5.2, dpi = 300)  # Adjust width and height as needed
}


# Diversidade - total
tidy_data <- setNames(
  lapply(names(bact_fung_taxa_total), function(list_name) {
    lista <- bact_fung_taxa_total[[list_name]]
    gathered_data <- pivot_longer(lista, cols = -1, names_to = "Sample", values_to = "Percent")
    attr(gathered_data, "original_name") <- list_name  # Store the original name as an attribute
    return(gathered_data)
  }),
  names(bact_fung_taxa_total)
)

tidy_data2 <- tidy_data[c(3, 9)]
names(tidy_data2) <- c("Bacteria_species", "Fungi_species")

barplot_list <- setNames(lapply(names(tidy_data2), function(list_name) {
  df <- tidy_data2[[list_name]]
  repla <- gsub(".*_", "", list_name)
  df_filtered <- df %>% filter(Percent != 0)
  ggplot(df_filtered, aes(x = Sample, y = Percent, fill = .data[[repla]])) +
    geom_bar(position="fill", stat="identity") +
    labs(x = "Samples", y = "Average percentage of reads", fill = paste(toTitleCase(repla), "es", sep="")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5),
      #legend.position = "bottom",
      legend.text = element_text(size = 8,       # Adjust the font size of the legend text (set to 8 or adjust as needed)
                               hjust = 0.5,     # Center the text horizontally in the legend key
                               vjust = 0.5)) +
    guides(fill = guide_legend(label.hjust = 0.5, label.vjust = 0.5)) +
    scale_y_continuous(labels = percent_format(scale = 100)) # Set y-axis labels as percentages
}), names(tidy_data2))


barplot_list <- setNames(lapply(seq_along(tidy_data2), function(i) {
  df <- tidy_data2[[i]]
  repla <- gsub(".*_", "", names(tidy_data2)[i])
  df_filtered <- df %>% filter(Percent != 0)
  ggplot(df_filtered, aes(x = Sample, y = Percent, fill = .data[[repla]])) +
    geom_bar(position="fill", stat="identity") +
    labs(x = "Samples", y = "Average percentage of reads", fill = paste(toTitleCase(repla), "es", sep="")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 4.5),
      #legend.position = "bottom",
      legend.text = element_text(size = 8,       # Adjust the font size of the legend text (set to 8 or adjust as needed)
                               hjust = 0.5,     # Center the text horizontally in the legend key
                               vjust = 0.5)) +
    guides(fill = guide_legend(label.hjust = 0.5, label.vjust = 0.5)) +
    scale_y_continuous(labels = percent_format(scale = 100)) # Set y-axis labels as percentages
}), seq_along(tidy_data2))



# Print the bar plots
print(barplot_list[[1]]) # Print the first bar plot
print(barplot_list[[2]]) # Print the second bar plot


# Save each plot with a unique filename based on the original name of the data frame
for(i in seq_along(barplot_list)) {
  filename <- paste0(directory, attr(tidy_data[[i]], "original_name"), ".png")
  ggsave(filename, plot = barplot_list[[i]], width = 8, height = 5.2, dpi = 300)  # Adjust width and height as needed
}