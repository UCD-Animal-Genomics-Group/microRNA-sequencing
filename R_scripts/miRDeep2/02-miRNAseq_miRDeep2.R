##########################################################################
# miRNA-seq analysis of a time course experimental infection in cattle   #
#                                                                        #
#            --- R workflow for the Novo-feature approach ---            #
#                                 Part 2                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 11/01/2018

############################################
# 14 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(biobroom)
library(Cairo)
library(extrafont)
library(forcats)
library(ggjoy)
library(PerformanceAnalytics)
library(ggrepel)
library(cowplot)

# Uncomment functions below to install packages in case you don't have them

#install.packages("ggjoy")
#install.packages("PerformanceAnalytics")
#install.packages("ggrepel")
#install.packages("cowplot")

##################################
# 15 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_miRDeep2.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

##########################################################
# 16 Tidy DGElist for exploratory data analysis plotting #
##########################################################

tidy_dgelist <- tidy(dgelist_filt, addSamples = TRUE)

tidy_dgelist

# Clean animal IDs
tidy_dgelist$animal %<>%
  stringr::str_replace("A", "") %>%
  fct_inorder()

# Change time point info for labels
tidy_dgelist$time.point %<>%
  str_replace("pre1", "-1 wk") %>%
  str_replace("pre2", "-2 wk") %>%
  str_replace("W1$", "+1 wk") %>%
  str_replace("W2", "+2 wk") %>%
  str_replace("W6", "+6 wk") %>%
  str_replace("W10", "+10 wk") %>%
  str_replace("W12", "+12 wk") %>%
  factor(levels = c("-2 wk", "-1 wk", "+1 wk", "+2 wk",
                    "+6 wk", "+10 wk", "+12 wk"))

# Combine animal and time point info for
# plotting labels
tidy_dgelist %<>%
  dplyr::mutate(labels = paste0(time.point, "_", animal))

tidy_dgelist$labels %<>%
  factor(levels = c("-2 wk_6511", "-2 wk_6514", "-2 wk_6520", "-2 wk_6522", "-2 wk_6526",
                    "-2 wk_6635", "-2 wk_6636", "-2 wk_6637", "-2 wk_6644", "-2 wk_6698",
                    "-1 wk_6511", "-1 wk_6514", "-1 wk_6520", "-1 wk_6522", "-1 wk_6526",
                    "-1 wk_6635", "-1 wk_6636", "-1 wk_6637", "-1 wk_6644", "-1 wk_6698",
                    "+1 wk_6511", "+1 wk_6514", "+1 wk_6520", "+1 wk_6522", "+1 wk_6526",
                    "+1 wk_6635", "+1 wk_6636", "+1 wk_6637", "+1 wk_6644", "+1 wk_6698",
                    "+2 wk_6511", "+2 wk_6514", "+2 wk_6520", "+2 wk_6522", "+2 wk_6526",
                    "+2 wk_6635", "+2 wk_6636", "+2 wk_6637", "+2 wk_6644", "+2 wk_6698",
                    "+6 wk_6511", "+6 wk_6514", "+6 wk_6520", "+6 wk_6522", "+6 wk_6526",
                    "+6 wk_6635", "+6 wk_6636", "+6 wk_6637", "+6 wk_6644", "+6 wk_6698",
                    "+10 wk_6511", "+10 wk_6514", "+10 wk_6520", "+10 wk_6522", "+10 wk_6526",
                    "+10 wk_6635", "+10 wk_6636", "+10 wk_6637", "+10 wk_6644", "+10 wk_6698",
                    "+12 wk_6511", "+12 wk_6514", "+12 wk_6520", "+12 wk_6522", "+12 wk_6526",
                    "+12 wk_6635", "+12 wk_6636", "+12 wk_6637", "+12 wk_6644", "+12 wk_6698"))

# Check factors
levels(tidy_dgelist$labels)

# Check data frame
tidy_dgelist

########################################################
# 17 Plot: density of filtered gene counts per library #
########################################################

ggplot(tidy_dgelist, aes(x = log10(count + 1),
                         y = labels)) +
  scale_y_discrete(limits = rev(levels(tidy_dgelist$labels))) +
  geom_joy(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Time point",
                    values = c("#b2b2b2", rep("#e06377", 6))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("Density of filtered gene counts per sample (Novo-feature)") +
  ylab("Time point_Animal number") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_filt


density_filt

# Export high quality image
ggsave(paste(method, "_density-filt.pdf", sep = ""),
       plot      = density_filt,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

##########################################################
# 18 Plot: Gene expression correlation between libraries #
##########################################################

# Check gene expression correlation between libraries at each time point to
# identify potential outliers (CPM values not required with Spearman)

for (val in c("pre2", "pre1", "1", "2", "6", "10", "12")) {
  png(file = file.path(paste0(imgDir, method, "_", val, "_cor.png", sep = "")),
      width = 1500, height = 1500, units = "px")
  chart.Correlation(R = log(x = (dgelist_norm$counts[, grep(
    pattern = paste("_", val, "$", sep = ""), x = colnames(dgelist_norm$counts),
    perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
    main = paste(val, " time point expression correlation", sep = ""))
  dev.off()
}

##################################################
# 19 Define function for getting MDS coordinates #
##################################################

getBCVcoord <- function(dgelst, time_pattrn) {

  mds <- plotMDS.DGEList(x = dgelst[ , grepl(paste(time_pattrn,
                                                   collapse = "|"),
                                             x = colnames(dgelst))],
                         plot = FALSE,
                         method = "bcv")

  mds_coord <- mds$cmdscale.out # Get coords to plot with ggplot2

  mds_coord %<>% # Tidy coords
    tidy() %>%
    dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
    dplyr::mutate(animal = sample, group = sample)

  mds_coord$animal %<>% # Clean animal IDs for plotting
    str_replace(paste(time_pattrn,
                      collapse = "|"), "") %>%
    factor()

  mds_coord$group %<>% # Clean group info for plotting
    str_replace(".*_", "") %>%
    fct_rev()

  return(mds_coord)

}

##########################
# 20 Get MDS coordinates #
##########################

# Get tidy MDS coordinates for all time points
all_coord <- getBCVcoord(dgelist_norm, c("_pre2", "_pre1", "_1$",
                                         "_2", "_6", "_10", "_12"))

# Get tidy MDS coordinates for each time point post-infection
# with both controls pre-infection
## +1 wk and controls
W1_coord <- getBCVcoord(dgelist_norm, c("_pre2", "_pre1", "_1$"))

## +2 wk and controls
W2_coord <- getBCVcoord(dgelist_norm, c("_pre2", "_pre1", "_2"))

## +6 wk and controls
W6_coord <- getBCVcoord(dgelist_norm, c("_pre2", "_pre1", "_6"))

## +10 wk and controls
W10_coord <- getBCVcoord(dgelist_norm, c("_pre2", "_pre1", "_10"))

## +12 wk and controls
W12_coord <- getBCVcoord(dgelist_norm, c("_pre2", "_pre1", "_12"))

#####################################
# 21 Plot: MDS with all time points #
#####################################

MDS_all <- ggplot(all_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c(rep("#b2b2b2", 2), rep("#e06377", 5))) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 0, 1, 2, 5, 6)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("All time points") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

MDS_all

# Export high quality image
ggsave(filename  = paste(method, "_MDS_all.pdf", sep = ""),
       plot      = MDS_all,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 10,
       units     = "in")

###################################################
# 22 Plots: MDS at each time point post-infection #
###################################################

# Plot +1 wk MDS
MDS_W1 <- ggplot(W1_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+1 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

MDS_W1

# Plot +2 wk MDS
MDS_W2 <- ggplot(W2_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+2 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

MDS_W2

# Plot +6 wk MDS
MDS_W6 <- ggplot(W6_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+6 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

MDS_W6

# Plot +10 wk MDS
MDS_W10 <- ggplot(W10_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+10 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

MDS_W10

# Plot +12 wk MDS
MDS_W12 <- ggplot(W12_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+12 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

MDS_W12

### Combine all MDS plots into single figure

# Set grid
MDS_grid <- plot_grid(MDS_W1, MDS_W2, MDS_W6,
                      MDS_W10, MDS_W12,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 2)

# Check plot
MDS_grid

# Export high quality image
ggsave(filename  = paste(method, "_MDS_grid.pdf", sep = ""),
       plot      = MDS_grid,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 15,
       width     = 20,
       units     = "in")

#################################################
# 23 Create a design matrix for paired analysis #
#################################################

# Create a design matrix with animal as a blocking factor
block_animal <- model.matrix(~ animal + group,
                             data = dgelist_norm$samples)

dim(block_animal)
dim(dgelist_norm$samples)
block_animal

# Rename design matrix columns for simplicity
colnames(block_animal) %<>%
  str_replace("animal", "") %>%
  str_replace("group", "")

head(block_animal)

# Output the design matrix info
write_csv(as.data.frame(block_animal),
          path = file.path(paste0(tablesDir, method, "_design-matrix.csv")),
          col_names = TRUE)

#########################################
# 24 Estimate the dispersion parameters #
#########################################

# Common and trended dispersions are estimated with the
# Cox-Reid method and tagwise dispersions with the
# empirical Bayes method
dgelist_disp <- estimateDisp.DGEList(y       = dgelist_norm,
                                     design  = block_animal,
                                     robust  = TRUE,
                                     verbose = TRUE)

names(dgelist_disp)

# Check the calculated dispersion
dgelist_disp$common.dispersion

# Check the calculated dispersion's square root,
# which corresponds to the biological coefficient of variation (BCV)
sqrt(dgelist_disp$common.dispersion)
sqrt(dgelist_disp$tagwise.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information
Tagwisedisp <- as.data.frame(cbind(dgelist_disp$genes,
                                   dgelist_disp$tagwise.dispersion))
head(Tagwisedisp)
dim(Tagwisedisp)

# Output tagwise dispersion values with gene info
write_csv(Tagwisedisp,
          path = file.path(paste0(tablesDir, method, "_Tagwise_dispersion.csv")),
          col_names = TRUE)

################################
# 25 Plot: BCV and dispersions #
################################

# Create a dataframe with the dispersion values
names(dgelist_disp)

Disp <- as.data.frame(cbind(dgelist_disp$genes,
                            dgelist_disp$tagwise.dispersion,
                            dgelist_disp$common.dispersion,
                            dgelist_disp$trended.dispersion,
                            dgelist_disp$AveLogCPM))

colnames(Disp) %<>%
  str_replace("dgelist_disp\\$", "")

Disp %<>%
  dplyr::mutate(type_point = "Tagwise dispersion") %>%
  dplyr::mutate(type_hline = "Common dispersion") %>%
  dplyr::mutate(type_smooth = "Trended dispersion")

# Plot all dispersions
ggplot(Disp) +
  geom_point(aes(x = AveLogCPM,
                 y = sqrt(tagwise.dispersion),
                 fill = type_point),
             alpha = 0.5) +
  geom_hline(aes(yintercept = sqrt(common.dispersion),
                 colour = type_hline)) +
  geom_smooth(aes(x = AveLogCPM,
                  y = sqrt(trended.dispersion),
                  colour = type_smooth),
              linetype = 2) +
  scale_fill_manual("", values = c("black")) +
  scale_colour_manual("", values = c("red", "blue")) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("Estimated dispersions") +
  xlab(expression(paste("Average ", log[2],"CPM"))) +
  ylab("Biological Coefficient of Variation") -> dgelist_BCV

dgelist_BCV

# Output high resolution plot
ggsave(paste(method, "_BCV.pdf", sep = ""),
       plot      = dgelist_BCV,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 10,
       units     = "in")

#######################
# 26 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 27 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 3 of this analysis #
######################################

# File: 03-miRNAseq_Novo-feature.R