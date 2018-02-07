##########################################################################
# miRNA-seq analysis of a time course experimental infection in cattle   #
#                                                                        #
#              --- R workflow for the miRDeep2 approach ---              #
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
method <- "miRDeep2"

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
  dplyr::select(-Log_Id) %>%
  dplyr::mutate(Animal = Sample, Time_point = Sample) %>%
  dplyr::select(Sample, Animal, Time_point, everything()) -> stats2

# Check data frame
stats2

#######################
# 04 Clean data frame #
#######################

# Animal factor
stats2$Animal %<>%
  str_replace("_\\w*\\d+", "") %>%
  fct_inorder()

levels(stats2$Animal)

# Time point factor
stats2$Time_point %<>%
  str_replace("\\d\\d\\d\\d_", "") %>%
  str_replace("pre2", "-2 wk") %>%
  str_replace("pre1", "-1 wk") %>%
  str_replace("1$", "+1 wk") %>%
  str_replace("12", "+12 wk") %>%
  str_replace("2$", "+2 wk") %>%
  str_replace("6", "+6 wk") %>%
  str_replace("10", "+10 wk") %>%
  factor(levels = c("-2 wk", "-1 wk", "+1 wk",
                    "+2 wk", "+6 wk", "+10 wk", "+12 wk"))

levels(stats2$Time_point)

# Export data
stats2 %>%
  dplyr::arrange(Animal, Time_point) %>%
  write_csv("mirdeep2_stats.csv", col_names = TRUE)

#############################################
# 05 Plot: Waffle charts of miRNA-seq stats #
#############################################

stats2 %<>%
  dplyr::mutate(perc_mapped = (Total_Mapped_Reads / Input_reads) * 100) %>%
  dplyr::mutate(perc_unmapp = (Total_Unmapped_Reads / Input_reads) * 100)

align_reads <- c(`Mapped reads (10%)` = round(mean(stats2$perc_mapped)),
                 `Unmapped reads (90%)` = round(mean(stats2$perc_unmapp)))
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











