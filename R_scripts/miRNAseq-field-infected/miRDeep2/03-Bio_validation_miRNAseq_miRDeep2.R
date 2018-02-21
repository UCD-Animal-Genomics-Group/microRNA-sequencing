############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#      --- R workflow for the miRDeep2 approach ---        #
#                       Part 3                             #
############################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 21/02/2018

############################################
# 27 Load and/or install required packages #
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

# Uncomment functions below to install packages in case you don't have them

#install.packages("VennDiagram")
#install.packages("treemap")

##################################
# 28 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_Bio-Val-miRDeep2.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

####################
# 29 Fit GLM model #
####################

# Fit a negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersions
dgelist_fit <- glmFit(y = dgelist_disp,
                      design = design)

names(dgelist_fit)
colnames(dgelist_fit$design)

#####################################################
# 30 Determine differential expression by fitting a #
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
# 31 Output all genes tested for DE #
#####################################

testDE$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_AllGenes.csv")),
            col_names = TRUE)

##############################################
# 32 Filter genes considered DE (FDR < 0.05) #
##############################################

testDE$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
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
# 33 Plot: Barplots of DE genes (FDR < 0.05) #
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

#######################
# 34 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 35 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 4 of this analysis #
######################################

# File: 04-Bio_validation_miRNAseq-miRDeep2-stats.R

