############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#      --- R workflow for the miRDeep2 approach ---        #
#                       Part 4                             #
############################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 21/02/2018

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
method <- "Bio-Val-miRDeep2"

# Set time zone
Sys.setenv(TZ = "Europe/London")

#############################
# 03 Import mapper.pl stats #
#############################

# Import stats and log ids
stats <- read_table2("mapper_stats.txt", col_names = TRUE)
log_ids <- read_tsv("id_sample.txt", col_names = FALSE)

# Check data frames
stats
log_ids

# Merge stats and log ids
log_ids %>%
  dplyr::rename(Log_Id = X1, Sample = X2) %>%
  dplyr::inner_join(stats, by = "Log_Id") %>%
  dplyr::select(-Log_Id) -> stats2

# Check data frame
stats2

#######################
# 04 Clean data frame #
#######################

# Sample name
stats2$Sample %<>%
  str_replace("_collapsed.fa", "")

# Merge with sample info to get animal numbers
sample_info <- read_csv("Bio-Val-sample-info.csv", col_names = TRUE)

stats2 %<>%
  dplyr::inner_join(sample_info, by = c("Sample" = "Sample code name")) %>%
  dplyr::rename(Animal = `Animal id`) %>%
  dplyr::select(Sample, Animal, Treatment, everything())

# Sample as factor
stats2$Sample %<>%
  factor(levels = c("E1", "E2", "E3", "E4", "E5", "E6",
                    "E7", "E8", "E9", "E10", "E11", "E12",
                    "E13", "E14", "E15", "E16", "E17", "E18",
                    "E19", "E20", "E21", "E22", "E23", "E24"))

levels(stats2$Sample)

# Export data
stats2 %>%
  dplyr::arrange(Sample) %>%
  write_csv("Bio-Val-mirdeep2_stats.csv", col_names = TRUE)

#############################################
# 05 Plot: Waffle charts of miRNA-seq stats #
#############################################

stats2 %<>%
  dplyr::mutate(perc_mapped = (Total_Mapped_Reads / Input_reads) * 100) %>%
  dplyr::mutate(perc_unmapp = (Total_Unmapped_Reads / Input_reads) * 100)

align_reads <- c(`Mapped reads (77%)` = round(mean(stats2$perc_mapped)),
                 `Unmapped reads (23%)` = round(mean(stats2$perc_unmapp)))
align_waffle <- waffle(align_reads,
                       rows = 5,
                       size = 0.5,
                       title = paste0("Alignment of reads (", method, ")"),
                       colors = c("#5e3c99", "#b2abd2"))

align_waffle <- align_waffle +
                theme(text = element_text(size = 16, family = "Calibri"))

# Check plot
align_waffle

# Save plot
ggsave(paste(method, "_waffle.pdf", sep = ""),
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

save.image(file = paste0("miRNAseq_stats", method, ".RData", sep = ""))

##########################
# 07 Save R session info #
##########################

devtools::session_info()

#######
# End #
#######











