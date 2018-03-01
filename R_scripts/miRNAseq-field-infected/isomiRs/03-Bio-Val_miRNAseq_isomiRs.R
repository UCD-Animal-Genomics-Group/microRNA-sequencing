############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#             --- R workflow for isomiRs ---               #
#                          Part 3                          #
############################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 01/03/2018

############################################
# 29 Load and/or install required packages #
############################################

library(here)
library(edgeR)
library(tidyverse)
library(devtools)
library(stringr)
library(magrittr)
library(Cairo)
library(extrafont)
library(VennDiagram)
library(forcats)
library(treemap)

# Uncomment functions below to install packages in case you don't have them

#install.packages("VennDiagram")
#install.packages("treemap")

##################################
# 30 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_Bio-Val-isomiRs.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

####################
# 31 Fit GLM model #
####################

# Fit a negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersions
dgelist_fit <- glmFit(y = dgelist_disp,
                      design = design)

names(dgelist_fit)
colnames(dgelist_fit$design)

#####################################################
# 32 Determine differential expression by fitting a #
# negative binomial GLM with Likelihood Ratio Tests #
#####################################################

# Test for differential expression between each of the post-infection
# time points and both controls pre-infection,
# using the coefficients from degelist_fit$design

Infect <- glmLRT(dgelist_fit, coef = "Infected")
testDE <- topTags(object        = Infect,
                  n             = "inf",
                  adjust.method = "BH",
                  sort.by       = "PValue")

dim(testDE$table)
head(testDE$table)

#####################################
# 33 Output all genes tested for DE #
#####################################

testDE$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_AllGenes.csv")),
            col_names = TRUE)

##############################################
# 34 Filter genes considered DE (FDR < 0.05) #
##############################################

testDE$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> DE_FDR_05

# Check data frame
DE_FDR_05

# Output list of DE genes
DE_FDR_05 %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_DE_Genes.csv")),
            col_names = TRUE)

##############################################
# 35 Plot: Barplots of DE genes (FDR < 0.05) #
##############################################

# Get numbers of up and down regulated genes
# at each time point
DE_FDR_05 %>%
  dplyr::count(up = sum(logFC > 0),
               down = abs(sum(logFC < 0)),
               zero = sum(logFC == 0)) -> Up_Down

# Check data frame
Up_Down

# Gather up and down, and get rid of unnecessary columns
Up_Down %<>%
  gather(key = "DE_genes", value = "number", c(up, down)) %>%
  dplyr::select(-c(zero, n))

# Check data frame
Up_Down

# Up and down as factor
Up_Down$DE_genes %<>%
  factor() %>%
  fct_inorder()

levels(Up_Down$DE_genes)

# Plot chart
ggplot(Up_Down, aes(x = DE_genes,
                    y = number,
                    fill = DE_genes)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(label = Up_Down$number) +
  scale_fill_manual("Expression",
                    values = alpha(c("#762A83", "#1B7837"), 0.6)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("Canonical DE miRNAs (", method, ")")) +
  ylab(NULL) +
  xlab(NULL) -> DE_barplot

DE_barplot

ggsave(paste(method, "_DE_barplot.pdf", sep = ""),
       plot      = DE_barplot,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 7,
       units     = "in")

#########################################################
# 36 Get number of DE isomiR sequences per canonical ID #
#########################################################

head(testDE$table)
dim(testDE$table)

testDE$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::group_by(canonical_miRNA_ID) %>%
  dplyr::mutate(total_isomirs_per_canonical =
                  length(unique(`IsomiR sequence`))) %>%
  dplyr::select(canonical_miRNA_ID, canonical_miRNA_name,
                total_isomirs_per_canonical) %>%
  dplyr::arrange(desc(total_isomirs_per_canonical)) %>%
  dplyr::distinct() -> iso_canon

iso_canon

# Output data
iso_canon %>%
  write_csv(file.path(paste0(tablesDir, method, "-per-canonical.csv")),
            col_names = TRUE)

###########################################################
# 37 Plot: number of DE isomiR sequences per canonical ID #
###########################################################

# Plotting labels
iso_canon %<>%
  unite("label",
        c(canonical_miRNA_name, canonical_miRNA_ID),
        sep = " (")

iso_canon$label %<>%
  str_c(")") %>%
  factor() %>%
  fct_inorder()

# Check data frame
iso_canon

# Plot
ggplot(iso_canon) +
  geom_point(aes(x = total_isomirs_per_canonical,
                 y = label)) +
  scale_y_discrete(limits = rev(levels(iso_canon$label))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("IsomiRs per canonical DE miRNAs (", method, ")")) +
  ylab(NULL) +
  xlab("Total isomiR number") -> iso_canon_plot

iso_canon_plot

ggsave(paste(method, "_DE_iso_canon_plot.pdf", sep = ""),
       plot      = iso_canon_plot,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 15,
       width     = 8,
       units     = "in")

#######################
# 38 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 39 Save R session info #
##########################

devtools::session_info()

#######
# END #
#######
