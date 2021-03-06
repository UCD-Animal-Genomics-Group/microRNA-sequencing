##########################################################################
# miRNA-seq analysis of a time course experimental infection in cattle   #
#                                                                        #
#              --- R workflow for the miRDeep2 approach ---              #
#                                 Part 3                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 05/02/2018

############################################
# 28 Load and/or install required packages #
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
# 29 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_miRDeep2.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

####################
# 30 Fit GLM model #
####################

# Fit a negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersions
dgelist_fit <- glmFit(y = dgelist_disp,
                      design = block_animal)

names(dgelist_fit)
colnames(dgelist_fit$design)

#####################################################
# 31 Determine differential expression by fitting a #
# negative binomial GLM with Likelihood Ratio Tests #
#####################################################

# Test for differential expression between each of the post-infection
# time points and both controls pre-infection,
# using the coefficients from degelist_fit$design

# +1 wk
W1 <- glmLRT(dgelist_fit, coef = "W1")
testDE.W1 <- topTags(object        = W1,
                     n             = "inf",
                     adjust.method = "BH",
                     sort.by       = "PValue")

dim(testDE.W1$table)
head(testDE.W1$table)


# +2 wk
W2 <- glmLRT(dgelist_fit, coef = "W2")
testDE.W2 <- topTags(object        = W2,
                     n             = "inf",
                     adjust.method = "BH",
                     sort.by       = "PValue")

dim(testDE.W2$table)
head(testDE.W2$table)

# +6 wk
W6 <- glmLRT(dgelist_fit, coef = "W6")
testDE.W6 <- topTags(object        = W6,
                     n             = "inf",
                     adjust.method = "BH",
                     sort.by       = "PValue")

dim(testDE.W6$table)
head(testDE.W6$table)

# +10 wk
W10 <- glmLRT(dgelist_fit, coef = "W10")
testDE.W10 <- topTags(object        = W10,
                      n             = "inf",
                      adjust.method = "BH",
                      sort.by       = "PValue")

dim(testDE.W10$table)
head(testDE.W10$table)

# +12 wk
W12 <- glmLRT(dgelist_fit, coef = "W12")
testDE.W12 <- topTags(object        = W12,
                      n             = "inf",
                      adjust.method = "BH",
                      sort.by       = "PValue")

dim(testDE.W12$table)
head(testDE.W12$table)

#####################################
# 32 Output all genes tested for DE #
#####################################

# +1 wk
testDE.W1$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W1_AllGenes.csv")),
            col_names = TRUE)

# +2 wk
testDE.W2$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W2_AllGenes.csv")),
            col_names = TRUE)

# +6 wk
testDE.W6$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W6_AllGenes.csv")),
            col_names = TRUE)

# +10 wk
testDE.W10$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W10_AllGenes.csv")),
            col_names = TRUE)

# +12 wk
testDE.W12$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W12_AllGenes.csv")),
            col_names = TRUE)

##############################################
# 33 Filter genes considered DE (FDR < 0.05) #
##############################################

# +1 wk
testDE.W1$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W1_FDR_05

W1_FDR_05

# +2 wk
testDE.W2$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W2_FDR_05

W2_FDR_05

# +6 wk
testDE.W6$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W6_FDR_05

W6_FDR_05

# +10 wk
testDE.W10$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W10_FDR_05

W10_FDR_05

# +12 wk
testDE.W12$table %>%
  rownames_to_column(var = "miRBase_ID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W12_FDR_05

W12_FDR_05

### Output lists of genes considered DE (FDR < 0.05)
DElists <- list(W1_FDR_05, W2_FDR_05, W6_FDR_05,
                W10_FDR_05, W12_FDR_05)
DEfiles <- c(paste0(c("W1_FDR_05", "W2_FDR_05", "W6_FDR_05",
                      "W10_FDR_05", "W12_FDR_05"), "_genes.csv"))
DEpaths <- file.path(paste0(tablesDir, method, "_", DEfiles))

pwalk(list(DElists, DEpaths),
      write_csv,
      col_names = TRUE)

###################################################
# 34 Plot: Venn diagram of DE genes (FDR < 0.05) #
###################################################

# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
venn.diagram(list(A = as.vector(W1_FDR_05$miRBase_ID),
                  B = as.vector(W12_FDR_05$miRBase_ID),
                  C = as.vector(W10_FDR_05$miRBase_ID),
                  D = as.vector(W6_FDR_05$miRBase_ID),
                  E = as.vector(W2_FDR_05$miRBase_ID)),
             filename        = file.path(paste0(imgDir, method,
                                                "_Venn_DE_FDR_05.png")),
             imagetype       = "png",
             col             = "Black",
             fill            = c("#ffffcc", "#225ea8", "#ab9abc",
                                 "#a1dab4", "#41b6c4"),
             alpha           = 0.50,
             label.col       = "#003333",
             cex             = 3,
             fontfamily      = "Calibri",
             category.names  = c("+1 wk", "+12 wk", "+10 wk",
                                 "+6 wk", "+2 wk"),
             cat.col         = "black",
             cat.cex         = 3,
             cat.pos         = c(45, -20, 70, 125, 30),
             cat.dist        = c(0.20, 0.19, -0.23, 0.20, 0.19),
             cat.fontfamily  = "Calibri",
             margin          = 0,
             height          = 12,
             width           = 12,
             units           = 'in',
             compression     = 'lzw',
             resolution      = 300)

#############################################
# 35 Plot: Barplot of DE genes (FDR < 0.05) #
#############################################

# Get numbers of up and down regulated genes
# at each time point
list_DE <- list(W1_FDR_05, W2_FDR_05, W6_FDR_05,
                W10_FDR_05, W12_FDR_05)
names(list_DE) <- c("+1 wk", "+2 wk", "+6 wk",
                    "+10 wk", "+12 wk")

Up_Down <- map_df(list_DE,
                  ~ dplyr::count(.x,
                                 up = sum(logFC > 0),
                                 down = -abs(sum(logFC < 0)),
                                 zero = sum(logFC == 0)),
                  .id = "time_point")

# Check data frame
Up_Down

# Add +12 wk
Up_Down %<>%
  tibble::add_row(time_point = "+12 wk",
                  up = 0,
                  down = 0,
                  zero = 0,
                  n = 0)

# Time point as factor
Up_Down$time_point %<>%
  factor() %>%
  fct_inorder()

levels(Up_Down$time_point)

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
ggplot(Up_Down, aes(x = time_point,
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
# 36 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 37 Save R session info #
##########################

devtools::session_info()

######################################
# Proceed to Part 4 of this analysis #
######################################

# File: 04-miRNAseq-miRDeep2-stats.R
