##########################################################################
# miRNA-seq analysis of a time course experimental infection in cattle   #
#                                                                        #
#                     --- R workflow for isomiRs ---                     #
#                                 Part 3                                 #
##########################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 04/03/2018

############################################
# 30 Load and/or install required packages #
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
# 31 Working directory and RData #
##################################

# Check working directory
here()

# Load previously saved data
load("miRNAseq_isomiRs.RData")

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Load fonts
loadfonts()

####################
# 32 Fit GLM model #
####################

# Fit a negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersions
dgelist_fit <- glmFit(y = dgelist_disp,
                      design = block_animal)

names(dgelist_fit)
colnames(dgelist_fit$design)

#####################################################
# 33 Determine differential expression by fitting a #
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
# 34 Output all genes tested for DE #
#####################################

# +1 wk
testDE.W1$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W1_AllGenes.csv")),
            col_names = TRUE)

# +2 wk
testDE.W2$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W2_AllGenes.csv")),
            col_names = TRUE)

# +6 wk
testDE.W6$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W6_AllGenes.csv")),
            col_names = TRUE)

# +10 wk
testDE.W10$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W10_AllGenes.csv")),
            col_names = TRUE)

# +12 wk
testDE.W12$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "TRUE", "FALSE")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(tablesDir, method, "_W12_AllGenes.csv")),
            col_names = TRUE)

##############################################
# 35 Filter genes considered DE (FDR < 0.05) #
##############################################

# +1 wk
testDE.W1$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W1_FDR_05

W1_FDR_05

# +2 wk
testDE.W2$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W2_FDR_05

W2_FDR_05

# +6 wk
testDE.W6$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W6_FDR_05

W6_FDR_05

# +10 wk
testDE.W10$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(FDR) %>%
  as.tibble() -> W10_FDR_05

W10_FDR_05

# +12 wk
testDE.W12$table %>%
  rownames_to_column(var = "IsomiR sequence") %>%
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
# 36 Plot: Venn diagram of DE genes (FDR < 0.05) #
###################################################

# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
venn.diagram(list(A = as.vector(W1_FDR_05$`IsomiR sequence`),
                  B = as.vector(W12_FDR_05$`IsomiR sequence`),
                  C = as.vector(W10_FDR_05$`IsomiR sequence`),
                  D = as.vector(W6_FDR_05$`IsomiR sequence`),
                  E = as.vector(W2_FDR_05$`IsomiR sequence`)),
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

##############################################
# 37 Plot: Treemaps of DE genes (FDR < 0.05) #
##############################################

# Get numbers of up and down regulated genes
# at each time point
list_DE <- list(W1_FDR_05, W2_FDR_05, W6_FDR_05,
                W10_FDR_05, W12_FDR_05)
names(list_DE) <- c("+1 wk", "+2 wk", "+6 wk",
                    "+10 wk", "+12 wk")

Up_Down <- map_df(list_DE,
                  ~ dplyr::count(.x,
                                 up = sum(logFC > 0),
                                 down = sum(logFC < 0),
                                 zero = sum(logFC == 0)),
                  .id = "time_point")

# Check data frame
Up_Down

# Time point as factor
Up_Down$time_point %<>%
  factor() %>%
  fct_inorder()

# Plotting labels
Up_Down %<>% dplyr::mutate(labelsUp = paste(time_point, up, sep = ' '))
Up_Down %<>% dplyr::mutate(labelsDown = paste(time_point, down, sep = ' '))

# Check data frame
Up_Down

# Plot chart increased expression
# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, method, "_tree_up.pdf")),
          width    = 8,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
treemap(Up_Down,
        index             = "labelsUp",
        vSize             = "up",
        type              = "index",
        palette           = "PRGn",
        title             = "Increased expression",
        fontsize.title    = 14,
        fontfamily.title  = "Calibri",
        fontfamily.labels = "Calibri",
        fontsize.labels   = 16)
dev.off()

# Plot chart decreased expression
# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, method, "_tree_down.pdf")),
          width    = 8,
          height   = 4,
          family   = "Calibri",
          fallback_resolution = 300)
treemap(Up_Down,
        index             = "labelsDown",
        vSize             = "down",
        type              = "index",
        palette           = "-PRGn",
        title             = "Decreased expression",
        fontsize.title    = 14,
        fontfamily.title  = "Calibri",
        fontfamily.labels = "Calibri",
        fontsize.labels   = 16)
dev.off()

###############################################
# 38 Number of DE isomiRs per canonical miRNA #
###############################################

# +1 wk
W1_FDR_05 %>%
  dplyr::group_by(canonical_miRNA_ID) %>%
  dplyr::mutate(total_isomirs_per_canonical =
                  length(unique(`IsomiR sequence`))) %>%
  dplyr::select(canonical_miRNA_ID, canonical_miRNA_name,
                total_isomirs_per_canonical) %>%
  dplyr::arrange(desc(total_isomirs_per_canonical)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_iso_per_canon =
                  (total_isomirs_per_canonical /
                     sum(.$total_isomirs_per_canonical)
                   * 100)) %>%
  na.omit() -> iso_W1

iso_W1

# +2 wk
W2_FDR_05 %>%
  dplyr::group_by(canonical_miRNA_ID) %>%
  dplyr::mutate(total_isomirs_per_canonical =
                  length(unique(`IsomiR sequence`))) %>%
  dplyr::select(canonical_miRNA_ID, canonical_miRNA_name,
                total_isomirs_per_canonical) %>%
  dplyr::arrange(desc(total_isomirs_per_canonical)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_iso_per_canon =
                  (total_isomirs_per_canonical /
                     sum(.$total_isomirs_per_canonical)
                   * 100)) %>%
  na.omit() -> iso_W2

iso_W2

# +6 wk
W6_FDR_05 %>%
  dplyr::group_by(canonical_miRNA_ID) %>%
  dplyr::mutate(total_isomirs_per_canonical =
                  length(unique(`IsomiR sequence`))) %>%
  dplyr::select(canonical_miRNA_ID, canonical_miRNA_name,
                total_isomirs_per_canonical) %>%
  dplyr::arrange(desc(total_isomirs_per_canonical)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_iso_per_canon =
                  (total_isomirs_per_canonical /
                     sum(.$total_isomirs_per_canonical)
                   * 100)) %>%
  na.omit() -> iso_W6

iso_W6

# +10 wk
W10_FDR_05 %>%
  dplyr::group_by(canonical_miRNA_ID) %>%
  dplyr::mutate(total_isomirs_per_canonical =
                  length(unique(`IsomiR sequence`))) %>%
  dplyr::select(canonical_miRNA_ID, canonical_miRNA_name,
                total_isomirs_per_canonical) %>%
  dplyr::arrange(desc(total_isomirs_per_canonical)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_iso_per_canon =
                  (total_isomirs_per_canonical /
                     sum(.$total_isomirs_per_canonical)
                   * 100)) %>%
  na.omit() -> iso_W10

iso_W10

# +12 wk
W12_FDR_05 %>%
  dplyr::group_by(canonical_miRNA_ID) %>%
  dplyr::mutate(total_isomirs_per_canonical =
                  length(unique(`IsomiR sequence`))) %>%
  dplyr::select(canonical_miRNA_ID, canonical_miRNA_name,
                total_isomirs_per_canonical) %>%
  dplyr::arrange(desc(total_isomirs_per_canonical)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(percent_iso_per_canon =
                  (total_isomirs_per_canonical /
                     sum(.$total_isomirs_per_canonical)
                   * 100)) %>%
  na.omit() -> iso_W12

iso_W12

# Output data
isoDE <- list(iso_W1, iso_W2, iso_W6, iso_W10, iso_W12)
isofiles <- c(paste0(c("DE-iso_W1", "DE-iso_W2", "DE-iso_W6",
                       "DE-iso_W10", "DE-iso_W12"),
                     "_per-canonical.csv"))
isopaths <- file.path(paste0(tablesDir, method, "_", isofiles))

pwalk(list(isoDE, isopaths),
      write_csv,
      col_names = TRUE)

###########################################################
# 39 Plot: number of DE isomiR sequences per canonical ID #
###########################################################

# +1 wk
iso_W1 %<>% # Plotting labels
  unite("label",
        c(canonical_miRNA_name, canonical_miRNA_ID),
        sep = " (")

iso_W1$label %<>%
  str_c(")") %>%
  factor() %>%
  fct_inorder()

# Plot
ggplot(iso_W1) +
  geom_point(aes(x = total_isomirs_per_canonical,
                 y = label)) +
  scale_y_discrete(limits = rev(levels(iso_W1$label))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("DE isomiRs per canonical miRNA at +1 wk")) +
  ylab(NULL) +
  xlab("Total isomiR number") -> iso_W1_plot

iso_W1_plot


# +2 wk
iso_W2 %<>% # Plotting labels
  unite("label",
        c(canonical_miRNA_name, canonical_miRNA_ID),
        sep = " (")

iso_W2$label %<>%
  str_c(")") %>%
  factor() %>%
  fct_inorder()

# Plot
ggplot(iso_W2) +
  geom_point(aes(x = total_isomirs_per_canonical,
                 y = label)) +
  scale_y_discrete(limits = rev(levels(iso_W2$label))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("DE isomiRs per canonical miRNA at +2 wk")) +
  ylab(NULL) +
  xlab("Total isomiR number") -> iso_W2_plot

iso_W2_plot


# +6 wk
iso_W6 %<>% # Plotting labels
  unite("label",
        c(canonical_miRNA_name, canonical_miRNA_ID),
        sep = " (")

iso_W6$label %<>%
  str_c(")") %>%
  factor() %>%
  fct_inorder()

# Plot
ggplot(iso_W6) +
  geom_point(aes(x = total_isomirs_per_canonical,
                 y = label)) +
  scale_y_discrete(limits = rev(levels(iso_W6$label))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("DE isomiRs per canonical miRNA at +6 wk")) +
  ylab(NULL) +
  xlab("Total isomiR number") -> iso_W6_plot

iso_W6_plot


# +10 wk
iso_W10 %<>% # Plotting labels
  unite("label",
        c(canonical_miRNA_name, canonical_miRNA_ID),
        sep = " (")

iso_W10$label %<>%
  str_c(")") %>%
  factor() %>%
  fct_inorder()

# Plot
ggplot(iso_W10) +
  geom_point(aes(x = total_isomirs_per_canonical,
                 y = label)) +
  scale_y_discrete(limits = rev(levels(iso_W10$label))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("DE isomiRs per canonical miRNA at +10 wk")) +
  ylab(NULL) +
  xlab("Total isomiR number") -> iso_W10_plot

iso_W10_plot


# +12 wk
iso_W12 %<>% # Plotting labels
  unite("label",
        c(canonical_miRNA_name, canonical_miRNA_ID),
        sep = " (")

iso_W12$label %<>%
  str_c(")") %>%
  factor() %>%
  fct_inorder()

# Plot
ggplot(iso_W12) +
  geom_point(aes(x = total_isomirs_per_canonical,
                 y = label)) +
  scale_y_discrete(limits = rev(levels(iso_W12$label))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("DE isomiRs per canonical miRNA at +12 wk")) +
  ylab(NULL) +
  xlab("Total isomiR number") -> iso_W12_plot

iso_W12_plot



# Export all plots
ggsave(paste(method, "_DE_iso_canon_plot.pdf", sep = ""),
       plot      = iso_canon_plot)

Plotfiles <- c(paste0(c("DE-iso_W1", "DE-iso_W2", "DE-iso_W6",
                       "DE-iso_W10", "DE-iso_W12"),
                     "_per-canonical.pdf"))
Plotpaths <- file.path(paste0(imgDir, Plotfiles))
isoPlots <- list(iso_W1_plot, iso_W2_plot, iso_W6_plot,
                 iso_W10_plot, iso_W12_plot)

pwalk(list(Plotpaths, isoPlots),
      ggsave,
      device    = cairo_pdf,
      limitsize = FALSE,
      dpi       = 300,
      height    = 15,
      width     = 8,
      units     = "in")

############################################
# 40 Plot: Common DE ismoiRs per canonical #
############################################

iso_W1 %>%
  dplyr::inner_join(iso_W2,
                    by = "label",
                    suffix =c("_W1", "_W2")) %>%
  dplyr::inner_join(iso_W6,
                    by = "label") %>%
  dplyr::inner_join(iso_W10,
                    by = "label",
                    suffix =c("_W6", "_W10")) %>%
  dplyr::inner_join(iso_W12,
                    by = "label") %>%
  dplyr::select(-starts_with("percent")) %>%
  gather(key = "time_point", value = "total_isomirs",
         starts_with("total_isomirs")) %>%
  ggplot() +
  geom_point(aes(x = total_isomirs,
                 y = label,
                 colour = time_point),
             size = 3) +
  scale_colour_brewer(palette = "BrBG") +
  scale_y_discrete(limits = rev(levels("label"))) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("Common DE isomiRs per canonical miRNA")) +
  ylab(NULL) +
  xlab("Total isomiR number") -> common_iso_plot

ggsave(filename  = paste0(imgDir, "common_iso_plot.pdf"),
       plot      = common_iso_plot,
       device    = cairo_pdf,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 9,
       units     = "in")

#######################
# 41 Save .RData file #
#######################

save.image(file = paste0("miRNAseq_", method, ".RData", sep = ""))

##########################
# 42 Save R session info #
##########################

devtools::session_info()

#######
# END #
#######
