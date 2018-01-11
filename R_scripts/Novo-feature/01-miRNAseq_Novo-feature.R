##########################################################################
# miRNA-seq analysis of a time course experimental infection in cattle   #
#                                                                        #
#            --- R workflow for the Novo-feature approach ---            #
#                                 Part 1                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 09/01/2018

############################################
# 01 Load and/or install required packages #
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
countsDir <- here("featC_mature_counts")
imgDir <- here("Figures/")
tablesDir <- here("Tables/")

# Define the method used for miRNA identification
method <- "Novo-feature"

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

#########################################
# 03 Import featureCounts sense counts  #
#########################################

# Create a vector of file names
files <- list.files(path        = countsDir,
                    pattern     = "^6",
                    all.files   = TRUE,
                    full.names  = FALSE,
                    recursive   = FALSE,
                    ignore.case = FALSE)

files
length(files)

# Create a dataframe with raw counts for all samples
rawCounts <- readDGE(path         = countsDir,
                     files        = files,
                     header       = TRUE,
                     comment.char = "#",
                     columns      = c(1, 7))
names(rawCounts)
head(rawCounts$samples)
head(rawCounts$counts)

# Add a prefix to avoid sample names starting with numbers
colnames(rawCounts$counts) %<>%
  str_replace("65", "A65") %>%
  str_replace("66", "A66")

rownames(rawCounts$samples) %<>%
  str_replace("65", "A65") %>%
  str_replace("66", "A66")

rawCounts$samples$files %<>%
  str_replace("65", "A65") %>%
  str_replace("66", "A66")

# Check data frames
head(rawCounts$samples)
head(rawCounts$counts)

###############################################################
# 04 Gene annotation using information obtained from GTF file #
###############################################################

# Read in the annotation information
miRNA_info <- read.table(file   = "miRNA_Btaurus.txt",
                         header = TRUE,
                         sep    = "\t",
                         quote  = "")

head(miRNA_info)
dim(miRNA_info)

# Determine which miRNAs have identical mature sequence but
# originate from different precursors and add it to annotation
miRNA_info %<>%
  dplyr::group_by(sequence) %>%
  dplyr::summarise(identical_sequence = paste(gene_id, collapse = ",")) %>%
  dplyr::inner_join(miRNA_info, by = "sequence") %>%
  dplyr::select(starts_with("gene"), chromosome, contains("position"),
                strand, sequence, starts_with("precursor"),
                identical_sequence)

head(miRNA_info)
dim(miRNA_info)

# Merge gene annotation with the counts and output data
rawCounts$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id") %>%
  dplyr::inner_join(x = miRNA_info, y = ., by = "gene_id") %>%
  write_csv(file.path(paste0(tablesDir, method, "_Raw-counts.csv",
                             sep = "")),
            col_names = TRUE)

# Add genes as rownames
miRNA_info %<>%
  data.frame() %>%
  column_to_rownames(var = "gene_id")

head(miRNA_info)
dim(miRNA_info)

#########################################
# 05 Add sample information for DGElist #
#########################################

# Treatment group (avoid using underscores)
rawCounts$samples$group <- rownames(rawCounts$samples)
rawCounts$samples$group %<>%
  str_replace("A\\d+_", "") %>%
  str_replace("pre(1|2)", "Control") %>%
  str_replace("^1", "W1") %>%
  str_replace("^2", "W2") %>%
  str_replace("^6", "W6") %>%
  factor(levels = c("Control", "W1", "W2", "W6", "W10", "W12"))

# Animal ID (cannot start with numbers)
rawCounts$samples$animal <- rownames(rawCounts$samples)
rawCounts$samples$animal %<>%
  str_replace("_.+", "") %>%
  factor()

# Time points (avoid using underscores)
rawCounts$samples$time.point <- rownames(rawCounts$samples)
rawCounts$samples$time.point %<>%
  str_replace("A\\d+_", "") %>%
  str_replace("^1", "W1") %>%
  str_replace("^2", "W2") %>%
  str_replace("^6", "W6") %>%
  factor(levels = c("pre2", "pre1", "W1", "W2", "W6", "W10", "W12"))

# Check data frame
head(rawCounts$samples)

#####################
# 06 Create DGElist #
#####################

# Create DGElist
dgelist <- DGEList(counts       = as.data.frame(rawCounts$counts),
                   group        = rawCounts$samples$group,
                   genes        = miRNA_info,
                   lib.size     = NULL,
                   norm.factors = NULL,
                   remove.zeros = FALSE)

names(dgelist)
dim(dgelist)
head(dgelist$counts)
head(dgelist$samples)
head(dgelist$genes)

# Include addtional sample information to DGElist
identical(rownames(rawCounts$samples), rownames(dgelist$samples))

dgelist$samples$animal <- rawCounts$samples$animal
dgelist$samples$time.point <- rawCounts$samples$time.point

head(dgelist$samples)

# Output sample information
dgelist$samples %>%
write_csv(file.path(paste0(tablesDir, method, "_samples-info.csv", sep = "")),
          col_names = TRUE)

################################################
# 07 Density plot: raw gene counts per library #
################################################

# Tidy DGElist and plot data
dgelist %>%
  tidy() %>%
  ggplot() +
    geom_density(aes(x     = log10(count + 1),
                     group = sample), size = 0.1) +
    theme_bw(base_size = 14, base_family = "Calibri") +
    ggtitle(method) +
    ylab("Density of raw gene counts per sample") +
    xlab(expression(paste(log[10], "(counts + 1)"))) -> density_raw


density_raw

# Export image
ggsave(paste0(method, "_Raw-density.pdf", sep = ""),
       plot      = density_raw,
       device    = cairo_pdf,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 10,
       path      = imgDir)

###########################################
# 08 Remove zero and lowly expressed tags #
###########################################

# Filter non-expressed tags (all genes that have zero counts in all samples)
dgelist_no_zeros <- dgelist[rowSums(dgelist$counts) > 0, ]
dim(dgelist_no_zeros$counts)
head(dgelist_no_zeros$counts)
colnames(dgelist_no_zeros$counts)

# Filter lowly expressed tags, retaining only tags with
# more than 50 counts per million in 10 or more libraries
# (10 libraries correspond to 10 biological replicates and represent
# one time point)
dgelist_filt <- dgelist_no_zeros[rowSums(cpm(dgelist_no_zeros) > 50) >= 10, ]
dim(dgelist_filt$counts)
head(dgelist_filt$counts)

# Ouptut filtered counts
dgelist_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "miRBaseID") %>%
  write_csv(file.path(paste0(tablesDir, method, "_Filt_counts.csv",
                             sep = "")),
            col_names = TRUE)

##############################
# 09 Recompute library sizes #
##############################

dgelist_filt$samples$lib.size <- colSums(dgelist_filt$counts)
head(dgelist_filt$samples)
head(dgelist$samples)

###########################################################################
# 10 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after
# calculating normalisation factors
dgelist_norm <- calcNormFactors(dgelist_filt, method = "TMM")
head(dgelist_norm$samples)

#######################
# 11 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 12 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 02-miRNAseq_Novo-feature.R













