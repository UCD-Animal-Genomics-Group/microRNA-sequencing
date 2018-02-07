##########################################################################
# miRNA-seq analysis of a time course experimental infection in cattle   #
#                                                                        #
#            --- R workflow for the Novo-feature approach ---            #
#                                 Part 4                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 02/02/2018

############################################
# 01 Load and/or install required packages #
############################################

library(here)
library(tidyverse)
library(magrittr)
library(stringr)
library(forcats)
library(Cairo)
library(extrafont)
library(waffle)

##################################
# 02 Working directory and RData #
##################################

# Check working directory
here()

# Define variables for subdirectories
imgDir <- here("Figures/")

# Define the method used for miRNA identification
method <- "Novo-feature"

# Set time zone
Sys.setenv(TZ = "Europe/London")

##########################################
# 03 Import cutadapt and Novoalign stats #
##########################################

# Import stats and log ids
stats <- read_csv("miRNAseq_stats.csv", col_names = TRUE)

# Check data frame
stats

##############################################
# 04 Plots: Waffle charts of miRNA-seq stats #
##############################################

# Cutadapt
stats %<>%
  dplyr::mutate(removed_reads = `Processed reads` - `Output reads`) %>%
  dplyr::mutate(perc_removed_reads = (removed_reads / `Processed reads`) * 100)

  filt_reads <- c(`Reads post-filtering (84%)` =
                    round(mean(stats$`Percentage output reads`)),
                  `Removed reads (16%)` =
                    round(mean(stats$perc_removed_reads)))

filt_waffle <- waffle(filt_reads,
                      rows = 5,
                      size = 0.5,
                      title = "Filtering of reads",
                      colors = c("#5ab4ac", "#d8b365"))

filt_waffle <- filt_waffle + theme(text=element_text(size = 16, family = "Calibri"))

# Check plot
filt_waffle

# Save plot
ggsave("cutadapt_waffle.pdf",
       plot      = filt_waffle,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 6,
       width     = 9,
       units     = "in")

# Novoalign
align_reads <- c(`Percentage uniquely aligned reads (17%)` =
                   round(mean(stats$`Percentage uniquely aligned reads`)),
                 `Percentage multi mapped reads (82%)` =
                   round(mean(stats$`Percentage multi mapped reads`)),
                 `Unmapped reads (0.35%)` = mean(stats$`Percentage unmapped reads`))

align_waffle <- waffle(align_reads,
                       rows = 5,
                       size = 0.5,
                       title = paste0("Alignment of reads (", method, ")"),
                       colors = c("#5e3c99", "#b2abd2", "#e66101"))

align_waffle <- align_waffle +
                theme(text = element_text(size = 16, family = "Calibri"))

# Check plot
align_waffle

# Save plot
ggsave(paste(method, "_waffle.pdf", sep = ""),
       plot      = align_waffle,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 6,
       width     = 9,
       units     = "in")

#######################
# 06 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_stats-cutadapt-and-", method, ".RData", sep = ""))

##########################
# 07 Save R session info #
##########################

devtools::session_info()

#######
# End #
#######











