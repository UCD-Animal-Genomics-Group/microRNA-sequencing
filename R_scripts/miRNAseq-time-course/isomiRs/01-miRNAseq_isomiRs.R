##########################################################################
# miRNA-seq analysis of a time course experimental infection in cattle   #
#                                                                        #
#                     --- R workflow for isomiRs ---                     #
#                                 Part 1                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 26/02/2018

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
countsDir <- here("isomiR_counts/")
imgDir <- here("Figures/")
tablesDir <- here("Tables/")

# Define the method used for miRNA identification
method <- "isomiRs"

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

############################################
# 03 Import isomiR counts from miRDeep2.pl #
############################################

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
                     columns      = c(2, 4))

names(rawCounts)
head(rawCounts$samples)
head(rawCounts$counts)
dim(rawCounts$counts)

# Add a prefix to avoid sample names starting with numbers
colnames(rawCounts$counts) %<>%
  str_replace("_full_isomiR", "") %>%
  str_replace("65", "A65") %>%
  str_replace("66", "A66")

rownames(rawCounts$samples) %<>%
  str_replace("_full_isomiR", "") %>%
  str_replace("65", "A65") %>%
  str_replace("66", "A66")

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

# Add canonical gene name
read_tsv("miRNA_Btaurus.txt", col_names = TRUE) %>%
  dplyr::select(gene_id, gene_name, chromosome, strand) %>%
  dplyr::rename(canonical_miRNA_ID = gene_id,
                canonical_miRNA_name = gene_name,
                canon_chromosome = chromosome,
                canon_strand = strand) %>%
  dplyr::right_join(miRNA_info, by = c("canonical_miRNA_ID")) %>%
  dplyr::select(sequence, isomiR_position_start,
                isomiR_position_end, mismatch,
                everything()) -> miRNA_info

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
  dplyr::select(sequence, starts_with("A6")) %>%
  data.frame() %>%
  column_to_rownames(var = "sequence") -> counts

head(counts)
dim(counts)

# Define control and treatment groups
colnames(counts) %>%
  as.tibble() %>%
  dplyr::select(sample = value) %>%
  dplyr::mutate(group = sample) -> group

group

group$group %<>%
  str_replace("A\\d\\d\\d\\d", "") %>%
  str_replace("_pre(1|2)", "Control") %>%
  str_replace("_", "W")

group

group %<>%
  data.frame() %>%
  column_to_rownames(var = "sample")

head(group)
dim(group)
identical(rownames(group), colnames(counts))

# Select annotation only
annotCounts %>%
  dplyr::select(-starts_with("A6")) %>%
  data.frame() %>%
  column_to_rownames(var = "sequence") -> gene_annot

head(gene_annot)
dim(gene_annot)

#####################
# 08 Create DGElist #
#####################

# Create DGElist
dgelist <- DGEList(counts       = counts,
                   group        = group$group,
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
dgelist$samples$animal <- rownames(dgelist$samples)
dgelist$samples$animal %<>%
  str_replace("_.+", "") %>%
  factor()

# Check order of factor levels
levels(dgelist$samples$animal)

# Include time points (avoid using underscores)
dgelist$samples$time.point <- rownames(dgelist$samples)
dgelist$samples$time.point %<>%
  str_replace("A\\d+_", "") %>%
  str_replace("^1", "W1") %>%
  str_replace("^2", "W2") %>%
  str_replace("^6", "W6") %>%
  factor(levels = c("pre2", "pre1", "W1", "W2", "W6", "W10", "W12"))

# Check order of factor levels
levels(dgelist$samples$time.point)

# Convert treatment group to factor
dgelist$samples$group %<>%
  factor(levels = c("Control", "W1", "W2", "W6", "W10", "W12"))

# Check order of factor levels
levels(dgelist$samples$group)

# Check data frame
head(dgelist$samples)

# Output sample information
dgelist$samples %>%
  rownames_to_column(var = "sample") %>%
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

# File: 02-miRNAseq_isomiRs.R













