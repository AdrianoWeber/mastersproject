## This is a script to plot the output of a Plink pca run.
## An info file must also be given to ensure population coloration. The file must be a 3 columns tab (or space) separated values table (Sample ID, Continent, Population) 
## Usage: Rscript plot_plinkpca.R [pca_res_dir] [info_pop_dir] [TRUE | FALSE] 
## The third parameter is set to False and should be "T" if samples should be named on the plot.

## This script is part of the master thesis of Adriano Weber at the University of Geneva (09/2024-09/2025).
## This script is licensed under the GNU General Public License v3.0 or later.
## See https://www.gnu.org/licenses/gpl-3.0.en.html for more details.

#Library calling
library(ggplot2)
library(ggrepel)

#Argument passing
arg<-commandArgs(trailingOnly=T)
input_dir<-arg[1]
eigenvec_files <- list.files(path = input_dir, pattern = "\\.eigenvec$", full.names = TRUE)
eigenval_files <- list.files(path = input_dir, pattern = "\\.eigenval$", full.names = TRUE)
info_pop<- arg[2]

plot_name <- if (length(arg) >= 3 && !is.na(arg[3])) arg[3] else "F"

#Reading inputs
pca_dat<-read.table(eigenvec_files)
pc_val <- read.table(eigenval_files)
colnames(pca_dat) <- c("FID", "IID", paste0("PC", 1:(ncol(pca_dat) - 2)))
samples<-read.table(info_pop,header=T,fill=T)
colnames(samples) <- c("IID","Continent", "Population")
samples$Continent <- as.factor(samples$Continent)
samples$Population <- as.factor(samples$Population)
merged_data <- merge(pca_dat, samples, by = "IID", all.x = TRUE)

#Making nice plots
lim <- max(abs(c(merged_data$PC1, merged_data$PC2)))+0.02
p1 <- ggplot(merged_data, aes(x = PC1, y = PC2, colour = Population)) +
  geom_point() +
  geom_text_repel(aes(label = IID), max.overlaps = 100, box.padding = 0.5, point.padding = 0.3) +
  labs(title = "PCA Plot", x = paste0("PC1 (", pc_val$V1[1], "%)"), y = paste0("PC2 (", pc_val$V1[2], "%)")) +
  theme_minimal() +
  coord_fixed() +
  scale_x_continuous(limits = c(-lim, lim), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-lim, lim), expand = expansion(mult = 0)) +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))

p2 <- ggplot(merged_data, aes(x = PC1, y = PC2, colour = Continent)) +
  geom_point() +
  labs(title = "PCA Plot", x = paste0("PC1 (", pc_val$V1[1], "%)"), y = paste0("PC2 (", pc_val$V1[2], "%)")) +
  theme_minimal() +
  coord_fixed() +
  scale_x_continuous(limits = c(-lim, lim), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-lim, lim), expand = expansion(mult = 0)) +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))
p3 <- ggplot(merged_data, aes(x = PC1, y = PC2, colour = Population)) +
  geom_point() +
  labs(title = "PCA Plot", x = paste0("PC1 (", pc_val$V1[1], "%)"), y = paste0("PC2 (", pc_val$V1[2], "%)")) +
  theme_minimal() +
  coord_fixed() +
  scale_x_continuous(limits = c(-lim, lim), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-lim, lim), expand = expansion(mult = 0)) +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"))

#Saving outputs
time <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
ggsave(filename=paste("pca_plot_continent", time, ".jpeg", sep="_"), plot=p2, device=jpeg)
ggsave(filename=paste("pca_plot_population", time, ".jpeg", sep="_"), plot=p3, device=jpeg)
if (plot_name=="T"){
  ggsave(filename=paste("pca_plot_pop_named", time, ".jpeg", sep="_"), plot=p1, device=jpeg)
} else message("Argument 'T' non provided or incorrect. Skipping named samples plot.")