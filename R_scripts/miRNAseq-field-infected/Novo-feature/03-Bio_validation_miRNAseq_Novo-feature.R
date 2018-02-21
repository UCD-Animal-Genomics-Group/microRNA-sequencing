############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#     --- R workflow for the Novo-feature approach ---     #
#                       Part 3                             #
############################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 07/02/2018

############################################
# 26 Load and/or install required packages #
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
# 27 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_Bio-Val_Novo-feature.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

####################
# 28 Fit GLM model #
####################

# Fit a negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersions
dgelist_fit <- glmFit(y = dgelist_disp,
                      design = design)

names(dgelist_fit)
colnames(dgelist_fit$design)

#####################################################
# 29 Determine differential expression by fitting a #
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
# 30 Output all genes tested for DE #
#####################################

testDE$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_AllGenes.csv")),
            col_names = TRUE)

##############################################
# 31 Filter genes considered DE (FDR < 0.05) #
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
# 32 Plot: Barplots of DE genes (FDR < 0.05) #
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

#################################################
# 33 Tidy data for plotting most abundant genes #
#################################################

# Get CPM values normalised by library size
norm_cpm <- cpm(dgelist_fit$counts,
                normalized.lib.sizes = TRUE)

head(norm_cpm)

# Tidy data frame
norm_cpm %<>%
  data.frame() %>%
  rownames_to_column(var = "miRBaseID") %>%
  gather(key = "sample", value = "cpm", starts_with("E"))

head(norm_cpm)

# Merge with gene information
dgelist_fit$genes %>%
  rownames_to_column(var = "miRBaseID") %>%
  dplyr::select(miRBaseID, gene_name) %>%
  dplyr::inner_join(x = norm_cpm,
                    y = ., by = "miRBaseID") %>%
  unite(col = "label",
        c("gene_name", "miRBaseID"),
        sep = " (") -> norm_cpm

norm_cpm$label %<>%
  str_c(., ")")

head(norm_cpm)

#################################################
# 34 Plot: Boxplot of most abundant miRNA genes #
#################################################

# Get top 30 genes
norm_cpm %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(sum_cpm = sum(cpm)) %>%
  dplyr::arrange(desc(sum_cpm)) %>%
  dplyr::filter(!str_detect(label, "MIMAT0003844_1")) %>%
  dplyr::filter(!str_detect(label, "MIMAT0024579_5")) %>%
  dplyr::filter(!str_detect(label, "MIMAT0009225_1")) %>%
  dplyr::filter(!str_detect(label, "MIMAT0025559_1")) %>%
  dplyr::filter(!str_detect(label, "MIMAT0025559_6")) %>%
  head(30) -> top30


# Filter and plot
norm_cpm %>%
  dplyr::filter(label %in% top30$label) %>%
  ggplot() +
  geom_boxplot(aes(label, log10(cpm)),
               fill = "#efedf5") +
  scale_x_discrete(limits = top30$label) +
  theme_bw() +
  ggtitle(method) +
  xlab(NULL) +
  ylab(expression(paste(log[10], "CPM", sep = ""))) +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        text = element_text(size = 14,
                            family = "Calibri")) -> boxplot_top

ggsave(paste0(method, "-boxplot", ".pdf"),
       plot      = boxplot_top,
       device    = cairo_pdf,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 15,
       units     = "in",
       path      = imgDir)

#######################
# 35 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 36 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 4 of this analysis #
######################################

# File: 04-Bio_validation_miRNAseq-Novo-feature-stats.R
