############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#     --- R workflow for prediction of novel miRNAs ---    #
############################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 22/02/2018

############################################
# 01 Load and/or install required packages #
############################################

library(here)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(edgeR)
library(Cairo)
library(extrafont)

# Uncomment functions below to install packages in case you don't have them
# Bioconductor packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")


# CRAN packages
#install.packages("here")
#install.packages("tidyverse")
#install.packages("biobroom")
#install.packages("Cairo")
#install.packages("extrafont")

##############################################
# 02 Working directory, fonts, and time zone #
##############################################

# Check working directory
here()

# Define variables for subdirectories
countsDir <- here("clean_novel/")

# Define the method used for miRNA identification
method <- "Bio-Val-Novel"

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

############################################
# 03 Import novel miRNA info from miRDeep2 #
############################################

# Create a vector of file names
files <- list.files(path        = countsDir,
                    pattern     = "^E",
                    all.files   = TRUE,
                    full.names  = FALSE,
                    recursive   = FALSE,
                    ignore.case = FALSE)

files
names(files) <- files
length(files)

# Create a dataframe with novel miRNAs for all samples
files %>%
  map_df(~ read_tsv(paste0(countsDir, .x),
                    col_names = TRUE),
         .id = "filename") %>%
  dplyr::filter(`miRDeep2 score` >= 4,
                `mature read count` > 200,
                `rfam alert` != "rRNA/tRNA",
                `rfam alert` != "rRNA",
                `significant randfold p-value`
                == "yes") -> novel_info

View(novel_info)

# Checks
unique(novel_info$filename)
unique(novel_info$`provisional id`)
unique(novel_info$`consensus mature sequence`)
unique(novel_info$`consensus precursor sequence`)
unique(novel_info$`precursor coordinate`)

# summary
novel_info %>%
  dplyr::group_by(`provisional id`) %>%
  dplyr::summarise(mean_mat_count =
                     mean(`mature read count`),
                   tot_samples = n_distinct(filename)) %>%
  View()

write_csv(novel_info, "novel_miRNAs.csv", col_names = TRUE)


















