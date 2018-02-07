############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#     --- R workflow for the Novo-feature approach ---     #
#                       Part 2                             #
############################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 07/02/2018

############################################
# 13 Load and/or install required packages #
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
# 14 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_Bio-Val_Novo-feature.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

##########################################################
# 15 Tidy DGElist for exploratory data analysis plotting #
##########################################################

tidy_dgelist <- tidy(dgelist_filt, addSamples = TRUE)

# Sample as factor
tidy_dgelist$sample %<>%
  factor(levels = c("E1", "E2", "E3", "E4", "E5", "E6",
                    "E7", "E8", "E9", "E10", "E11", "E12",
                    "E13", "E14", "E15", "E16", "E17", "E18",
                    "E19", "E20", "E21", "E22", "E23", "E24"))

# Check factors
levels(tidy_dgelist$sample)

# Animal as factor
tidy_dgelist$animal %<>%
  factor()

# Check factors
levels(tidy_dgelist$animal)

# Check data frame
tidy_dgelist

########################################################
# 16 Plot: density of filtered gene counts per library #
########################################################

ggplot(tidy_dgelist, aes(x = log10(count + 1),
                         y = sample)) +
  scale_y_discrete(limits = rev(levels(tidy_dgelist$sample))) +
  geom_joy(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Time point",
                    values = c("#b2b2b2", rep("#e06377", 6))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("Density of filtered gene counts per sample (", method, ")")) +
  ylab("Sample") +
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
# 17 Plot: Gene expression correlation between libraries #
##########################################################

# Check gene expression correlation between libraries at each time point to
# identify potential outliers (CPM values not required with Spearman)

png(file = file.path(paste0(imgDir, method, "_Control", "_cor.png", sep = "")),
    width = 1500, height = 1500, units = "px")
chart.Correlation(R = log(x = (dgelist_norm$counts[, grep(
  pattern = "E(1|2|3|4|5|6|7|8|9|10|11|12)$", x = colnames(dgelist_norm$counts),
  perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
  main = " Control group expression correlation")
dev.off()

png(file = file.path(paste0(imgDir, method, "_Infected", "_cor.png", sep = "")),
    width = 1500, height = 1500, units = "px")
chart.Correlation(R = log(x = (dgelist_norm$counts[, grep(
  pattern = "E(13|14||15|16|17|18|19|20|21|22|23|24)$", x = colnames(dgelist_norm$counts),
  perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
  main = " Infected group expression correlation")
dev.off()

##################################################
# 18 Define function for getting MDS coordinates #
##################################################

getBCVcoord <- function(dgelst, pattrn) {

  mds <- plotMDS.DGEList(x = dgelst[ , grepl(paste(pattrn,
                                                   collapse = "|"),
                                             x = colnames(dgelst))],
                         plot = FALSE,
                         method = "bcv")

  mds_coord <- mds$cmdscale.out # Get coords to plot with ggplot2

  mds_coord %<>% # Tidy coords
    tidy() %>%
    dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
    dplyr::mutate(group = sample)

  mds_coord$group %<>% # Clean group info for plotting
    str_replace("E(1|2|3|4|5|6|7|8|9|10|11|12)$", "Control") %>%
    str_replace("E(13|14||15|16|17|18|19|20|21|22|23|24)$", "Infected")

  return(mds_coord)

}

##########################
# 19 Get MDS coordinates #
##########################

# Get tidy MDS coordinates for all time points
all_coord <- getBCVcoord(dgelist_norm, colnames(dgelist_norm$counts))

# Add animal info
Samples %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(sample, animal) -> animal_id

all_coord %<>%
  dplyr::inner_join(animal_id,
                    by = "sample")

#####################################
# 20 Plot: MDS with all time points #
#####################################

MDS_all <- ggplot(all_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(method) +
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
       height    = 5,
       width     = 6,
       units     = "in")

######################################################
# 21 Create a design matrix for independent analysis #
######################################################

# Create a design matrix with animal as a blocking factor
design <- model.matrix(~ group,
                       data = dgelist_norm$samples)

dim(design)
dim(dgelist_norm$samples)
design

# Rename design matrix columns for simplicity
colnames(design) %<>%
  str_replace("group", "")

head(design)

# Output the design matrix info
design %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  write_csv(path = file.path(paste0(tablesDir, method, "_design-matrix.csv")),
            col_names = TRUE)

#########################################
# 22 Estimate the dispersion parameters #
#########################################

# Common and trended dispersions are estimated with the
# Cox-Reid method and tagwise dispersions with the
# empirical Bayes method
dgelist_disp <- estimateDisp.DGEList(y       = dgelist_norm,
                                     design  = design,
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
Tagwisedisp %>%
  rownames_to_column(var = "miRBase_ID") %>%
  write_csv(path = file.path(paste0(tablesDir, method, "_Tagwise_dispersion.csv")),
            col_names = TRUE)

################################
# 23 Plot: BCV and dispersions #
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
# 24 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 25 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 3 of this analysis #
######################################

# File: 03-Bio_validation_miRNAseq_Novo-feature.R
