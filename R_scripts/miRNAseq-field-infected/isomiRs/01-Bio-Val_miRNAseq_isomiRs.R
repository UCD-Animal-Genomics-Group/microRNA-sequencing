############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#             --- R workflow for isomiRs ---               #
#                          Part 1                          #
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
countsDir <- here("iso_counts/")
imgDir <- here("Figures/")
tablesDir <- here("Tables/")

# Define the method used for miRNA identification
method <- "Bio-Val-isomiRs"

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

############################################
# 03 Import isomiR counts from miRDeep2.pl #
############################################

# Create a vector of file names
files <- list.files(path        = countsDir,
                    pattern     = "^E",
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
                     columns      = c(2, 4))

names(rawCounts)
head(rawCounts$samples)
head(rawCounts$counts)
dim(rawCounts$counts)

# Add a prefix to avoid sample names starting with numbers
colnames(rawCounts$counts) %<>%
  str_replace("_isomiR", "")

rownames(rawCounts$samples) %<>%
  str_replace("_isomiR", "")

# Convert sequences to upper case
rownames(rawCounts$counts) %<>%
  toupper()

# Check data frames
head(rawCounts$samples)
head(rawCounts$counts)

#################################################
# 04 Gene information from miRDeep2 count files #
#################################################

# Name the vector of files for purrr::map_df
names(files) <- files

# Import counts from miRDeep2.pl again to get the
# gene information, convert miRNA sequences to uppercase so
# they match the counts imported with edgeR,
# and split gene information into multiple columns
files %>%
  map_df(~ read_tsv(paste0(countsDir, .x), col_names = TRUE),
         .id = "filename") %>%
  dplyr::select(-c(filename, IsomiR_ID, Count)) %>%
  dplyr::rename(sequence = Sequence) %>%
  dplyr::mutate(sequence = toupper(sequence)) %>%
  dplyr::distinct(sequence, .keep_all = TRUE) %>%
  separate(`Matching_primiRNA_(isomiR_position,pri_count,mismatch[Matching_miRNA,miRNA_position,count])`,
           c("primiRNA_ID", "isomiR_position", "pcount",
             "mismatch", "canonical_miRNA_ID", "canon_position",
             "ccount"), sep = " ", extra = "merge") %>%
  separate(isomiR_position,
           c("isomiR_position_start", "isomiR_position_end"),
           sep = "-", fill = "right") %>%
  separate(canon_position,
           c("canon_position_start", "canon_position_end"),
           sep = "-", fill = "right") %>%
  dplyr::select(-c(pcount, ccount)) -> miRNA_info

# Check data frame
miRNA_info

#################################
# 05 Clean gene annotation info #
#################################

# Clean metadata columns
unique(miRNA_info$isomiR_position_start)
miRNA_info$isomiR_position_start %<>%
  str_replace("\\(", "") %>%
  str_replace(",", "") %>%
  as.numeric() # Multihits will be converted to NA

unique(miRNA_info$isomiR_position_end)
miRNA_info$isomiR_position_end %<>%
  str_replace(",", "") %>%
  as.numeric()

unique(miRNA_info$mismatch)
miRNA_info$mismatch %<>%
  str_replace("mm", "") %>%
  as.numeric()

unique(miRNA_info$canonical_miRNA_ID)
miRNA_info$canonical_miRNA_ID %<>%
  str_replace("\\[", "") %>%
  str_replace(",", "")

unique(miRNA_info$canon_position_start)
miRNA_info$canon_position_start %<>%
  str_replace(",", "") %>%
  as.numeric()

unique(miRNA_info$canon_position_end)
miRNA_info$canon_position_end %<>%
  str_replace(",", "") %>%
  as.numeric()

# Check data frame
miRNA_info

######################################
# 06 Merge annotation and raw counts #
######################################

rawCounts$counts %>%
  data.frame() %>%
  rownames_to_column(var = "sequence") %>%
  as.tibble() %>%
  dplyr::inner_join(x = miRNA_info, by = "sequence") -> annotCounts

# Check data frame
annotCounts

# Check which columns have duplicated values
anyDuplicated(annotCounts$sequence) # no duplicated values
anyDuplicated(annotCounts$primiRNA_ID) # has duplicated

# How many NAs in the data frame
annotCounts %>%
  dplyr::filter_all(any_vars(is.na(.)))

# How many isomiRs don't have a canonical match
annotCounts %>%
  dplyr::filter(canonical_miRNA_ID == "NA")

# Output raw counts with gene information
annotCounts %>%
  write_csv(file.path(paste0(tablesDir, method, "_Raw-counts.csv",
                             sep = "")),
            col_names = TRUE)

##############################################
# 07 Create groups and variables for DGElist #
##############################################

# Select counts only and add sequence as rownames
annotCounts %>%
  dplyr::select(sequence, starts_with("E")) %>%
  data.frame() %>%
  column_to_rownames(var = "sequence") -> counts

head(counts)
dim(counts)

# Read sample info file
sample_info <- read_csv("Bio-Val-sample-info.csv", col_names = TRUE)
sample_info

# Define control and treatment groups
colnames(counts) %>%
  as.tibble() %>%
  dplyr::select(sample = value) %>%
  dplyr::inner_join(sample_info, by = c("sample" = "Sample code name")) -> group_info

group_info
identical(group_info$sample, colnames(counts))

# Select annotation only
annotCounts %>%
  dplyr::select(-starts_with("E")) %>%
  data.frame() %>%
  column_to_rownames(var = "sequence") -> gene_annot

head(gene_annot)
dim(gene_annot)

#####################
# 08 Create DGElist #
#####################

# Create DGElist
dgelist <- DGEList(counts       = counts,
                   group        = group_info$Treatment,
                   genes        = gene_annot,
                   lib.size     = NULL,
                   norm.factors = NULL,
                   remove.zeros = FALSE)

names(dgelist)
dim(dgelist)
head(dgelist$counts)
head(dgelist$samples)
head(dgelist$genes)

################################################
# 09 Additional sample information for DGElist #
################################################

# Include Animal ID (cannot start with numbers)
dgelist$samples$animal <- group_info$`Animal id`
dgelist$samples$animal %<>%
  str_c("A", .) %>%
  factor()

# Check order of factor levels
levels(dgelist$samples$animal)

# Convert treatment group to factor
dgelist$samples$group %<>%
  factor(levels = c("Control", "Infected"))

# Check order of factor levels
levels(dgelist$samples$group)

# Check data frame
head(dgelist$samples)

# Output sample information
dgelist$samples %>%
  rownames_to_column(var = "Sample") %>%
  write_csv(file.path(paste0(tablesDir, method, "_samples-info.csv", sep = "")),
            col_names = TRUE)

################################################
# 10 Density plot: raw gene counts per library #
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
# 11 Remove zero and lowly expressed tags #
###########################################

# Filter non-expressed tags (all genes that have zero counts in all samples)
dgelist_no_zeros <- dgelist[rowSums(dgelist$counts) > 0, ]
dim(dgelist_no_zeros$counts)
head(dgelist_no_zeros$counts)
colnames(dgelist_no_zeros$counts)

# Filter lowly expressed tags, retaining only tags with
# more than 50 counts per million in 12 or more libraries
# (12 libraries correspond to any one condition group)
dgelist_filt <- dgelist_no_zeros[rowSums(cpm(dgelist_no_zeros) > 50) >= 12, ]
dim(dgelist_filt$counts)
head(dgelist_filt$counts)

# Ouptut filtered counts
dgelist_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "IsomiR Sequence") %>%
  write_csv(file.path(paste0(tablesDir, method, "_Filt_counts.csv",
                             sep = "")),
            col_names = TRUE)

##############################
# 12 Recompute library sizes #
##############################

dgelist_filt$samples$lib.size <- colSums(dgelist_filt$counts)
head(dgelist_filt$samples)
head(dgelist$samples)

###########################################################################
# 13 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after
# calculating normalisation factors
dgelist_norm <- calcNormFactors(dgelist_filt, method = "TMM")
head(dgelist_norm$samples)

#######################
# 14 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 15 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 02-Bio-Val_miRNAseq_isomiRs.R













