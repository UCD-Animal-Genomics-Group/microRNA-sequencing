############################################################
#   miRNA-seq analysis of positive retested reactor cattle #
#                   from Irish herds                       #
#                                                          #
#       --- Comparisons between Novo-feature and           #
#                   miRDeep2 approaches ---                #
############################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 21/02/2018

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
library(VennDiagram)
library(ggcorrplot)

# Uncomment line below to install package:
#install.packages("ggcorrplot")

##################################
# 02 Working directory and RData #
##################################

# Check working directory
here()

# Define variables for subdirectories
imgDir <- here("Figures_comparison/")
filesDir <- here("Docs_comparison/")


# Set time zone
Sys.setenv(TZ = "Europe/London")

###################################################################
# 03 Plot: Correlation of common miRNAs Novo-feature and miRDeep2 #
###################################################################

Novo_filt <- read_csv(paste0(filesDir, "Bio-Val_Novo-feature_Filt_counts.csv"),
                      col_names = TRUE)

mirdeep_filt <- read_csv(paste0(filesDir, "Bio-Val-miRDeep2_Filt_counts.csv"),
                         col_names = TRUE)

Novo_filt %>%
  dplyr::inner_join(mirdeep_filt,
                    by = "miRBaseID",
                    suffix = c("_novo", "_mirdeep")) -> common_filt

common_filt

correlation <- cor(x = log2((common_filt[, grep(pattern = "_novo$",
                                                x = colnames(common_filt),
                                                perl = TRUE)] + 1)),
                   y = log2((common_filt[, grep(pattern = "_mirdeep$",
                                                x = colnames(common_filt),
                                                perl = TRUE)] + 1)),
                   method = "spearman")

cor_plot <- ggcorrplot(correlation, hc.order = TRUE, type = "lower",
                       outline.col = "white",
                       colors = c("#6D9EC1", "white", "#E46726")) +
            theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1),
            text = element_text(size = 13, family = "Calibri"))

cor_plot

ggsave("correlation_common.pdf",
       plot      = cor_plot,
       device    = cairo_pdf,
       limitsize = FALSE,
       dpi       = 300,
       height    = 11,
       width     = 10,
       units     = "in",
       path      = imgDir)

####################################################################
# 04 Plot: Venn diagram of Common miRNAs Novo-feature and miRDeep2 #
####################################################################

Novo_filt %<>%
  dplyr::select(miRBaseID)

mirdeep_filt %<>%
  dplyr::select(miRBaseID)

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.diagram(list(A = as.vector(Novo_filt$miRBaseID),
                  B = as.vector(mirdeep_filt$miRBaseID)),
             filename        = "merged_Venn_filt.png",
             imagetype       = "png",
             col             = "Black",
             fill            = c("#fc8d59", "#91bfdb"),
             alpha           = 0.50,
             label.col       = "#003333",
             cex             = 3,
             fontfamily      = "Calibri",
             category.names  = c("Novo-feature", "miRDeep2"),
             main            = "Canonical miRNAs post-filtering",
             main.fontfamily = "Calibri",
             main.cex        = 3,
             cat.col         = "black",
             cat.cex         = 2.5,
             cat.fontfamily  = "Calibri",
             margin          = 0,
             height          = 12,
             width           = 12,
             units           = 'in',
             compression     = 'lzw',
             resolution      = 300)

#########################################################
# 05 Plot: Barplots DE miRNAs Novo-feature and miRDeep2 #
#########################################################

novo_DE <- read_csv(paste0(filesDir, "Bio-Val_Novo-feature_DE_Genes.csv"),
                    col_names = TRUE)
mirdeep_DE <- read_csv(paste0(filesDir, "Bio-Val-miRDeep2_DE_Genes.csv"),
                       col_names = TRUE)

# Get numbers of up and down regulated genes
# at each time point
novo_DE %>%
  dplyr::count(up = sum(logFC > 0),
               down = -abs(sum(logFC < 0)),
               zero = sum(logFC == 0)) -> Up_Down_Novo

mirdeep_DE %>%
  dplyr::count(up = sum(logFC > 0),
               down = -abs(sum(logFC < 0)),
               zero = sum(logFC == 0)) -> Up_Down_mirdeep

# Check data frames
Up_Down_Novo
Up_Down_mirdeep

# Gather up and down, and get rid of unnecessary columns
Up_Down_Novo %<>%
  gather(key = "DE_genes", value = "number", c(up, down)) %>%
  dplyr::select(-c(zero, n)) %>%
  dplyr::mutate(method = "Novo-feature")

Up_Down_mirdeep %<>%
  gather(key = "DE_genes", value = "number", c(up, down)) %>%
  dplyr::select(-c(zero, n)) %>%
  dplyr::mutate(method = "miRDeep2")

# Check data frames
Up_Down_Novo
Up_Down_mirdeep

# Bind data frames by row
Up_Down_Novo %>%
  dplyr::bind_rows(Up_Down_mirdeep) -> Up_Down

# Regulation as factor
Up_Down$DE_genes %<>%
  factor() %>%
  fct_inorder()

levels(Up_Down$DE_genes)

# Method as factor
Up_Down$method %<>%
  factor() %>%
  fct_inorder()

levels(Up_Down$method)

# Plot
ggplot(Up_Down, aes(x = method,
                    y = number,
                    fill = DE_genes)) +
  geom_bar(stat = "identity", position = "identity") +
  geom_text(label = Up_Down$number) +
  scale_fill_manual("Expression",
                    values = alpha(c("#762A83", "#1B7837"), 0.6)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle(paste0("Canonical DE miRNAs")) +
  ylab(NULL) +
  xlab(NULL) -> DE_barplot

DE_barplot

ggsave("common_DE_barplot.pdf",
       plot      = DE_barplot,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 5,
       width     = 7,
       units     = "in")

#################################################
# 06 Merge DE miRNAs Novo-feature and miRDeep2 #
#################################################

novo_DE %>%
  dplyr::full_join(mirdeep_DE,
                   by = c("miRBase_ID",
                          "gene_name",
                          "chromosome",
                          "start_position",
                          "end_position",
                          "strand",
                          "sequence",
                          "precursor_id",
                          "precursor_name",
                          "precursor_start",
                          "precursor_end",
                          "identical_sequence"),
                   suffix = c(" Novo-feature", " miRDeep2")) %>%
  dplyr::arrange(`FDR Novo-feature`) %>%
  write_csv(paste0(filesDir, "Bio-Val-novo-mirdeep2_DE-merged.csv"),
            col_names = TRUE)

#############################################
# 07 Plot: Venn diagram of common DE miRNAs #
#############################################

novo_DE %<>%
  dplyr::select(miRBase_ID)

mirdeep_DE %<>%
  dplyr::select(miRBase_ID)


futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn.diagram(list(A = as.vector(novo_DE$miRBase_ID),
                  B = as.vector(mirdeep_DE$miRBase_ID)),
             filename        = "Bio-Val-merged_Venn_DE.png",
             imagetype       = "png",
             col             = "Black",
             fill            = c("#fc8d59", "#91bfdb"),
             alpha           = 0.50,
             label.col       = "#003333",
             cex             = 3,
             fontfamily      = "Calibri",
             category.names  = c("Novo-feature", "miRDeep2"),
             main            = "DE Canonical miRNAs",
             main.fontfamily = "Calibri",
             main.cex        = 3,
             cat.col         = "black",
             cat.cex         = 2.5,
             cat.fontfamily  = "Calibri",
             margin          = 0,
             height          = 12,
             width           = 12,
             units           = 'in',
             compression     = 'lzw',
             resolution      = 300)


##############################
# 08 Filter common DE miRNAs #
##############################

novo_DE2 <- read_csv(paste0(filesDir, "Bio-Val_Novo-feature_DE_Genes.csv"),
                    col_names = TRUE)
mirdeep_DE2 <- read_csv(paste0(filesDir, "Bio-Val-miRDeep2_DE_Genes.csv"),
                       col_names = TRUE)

# Keep only common DE miRNAs to both methods
novo_DE2 %>%
  dplyr::inner_join(mirdeep_DE2,
                   by = c("miRBase_ID",
                          "gene_name",
                          "chromosome",
                          "start_position",
                          "end_position",
                          "strand",
                          "sequence",
                          "precursor_id",
                          "precursor_name",
                          "precursor_start",
                          "precursor_end",
                          "identical_sequence"),
                   suffix = c(" Novo-feature", " miRDeep2")) -> common_DE

# Check data frame
common_DE

# Output data
common_DE %>%
  write_csv(paste0(filesDir, "Bio-Val-commom-DE_novo-mirdeep.csv"), col_names = TRUE)

########################################
# 09 Plot: Heatmap of common DE miRNAs #
########################################

# Subset logFC
common_DE %<>%
  dplyr::select(miRBase_ID, gene_name,
                starts_with("logFC"))

common_DE

# Gather logFC data
common_DE %<>%
  gather(key = "method",
         value = "logFC",
         na.rm = FALSE,
         starts_with("logFC"))

common_DE

# Merge gene name and ID
common_DE %<>%
  unite(col = "label",
        c("gene_name", "miRBase_ID"),
        sep = " (")

# Clean labels
common_DE$label %<>%
  str_c(., ")")

# Clean method and convert to factor
common_DE$method %>%
  str_replace("logFC ", "") %>%
  factor(levels = c("Novo-feature", "miRDeep2"))

levels(common_DE$method)

# Check data frame
common_DE

# Plot heatmap
purpleGreen <- c("#A6DBA0", "#008837", "#78ba93", "gray",
                 "white", "#e7b9f7", "#C2A5CF", "#7B3294")

ggplot(common_DE) +
  geom_tile(aes(x = method, y = label,
                fill = logFC),
            colour = "black",
            size = 0.2) +
  coord_fixed() +
  xlab(NULL) +
  ylab(NULL) +
  ggtitle("Canonical DE miRNAs") +
  scale_fill_gradientn(name = expression(paste(log[2], "FC")),
                       colours = purpleGreen) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA,
                                    colour = "black", size = 1),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        axis.text.y = element_text(hjust = 1),
        text = element_text(size = 14, family = "Calibri"),
        plot.title = element_text(hjust = 0.5)) -> heatmap_common

heatmap_common

ggsave("Bio-Val-heatmap_common_DE.pdf",
       plot      = heatmap_common,
       device    = cairo_pdf,
       limitsize = FALSE,
       dpi       = 300,
       height    = 8,
       width     = 8,
       units     = "in",
       path      = imgDir)

#######################
# 10 Save .RData file #
#######################

save.image(file = "Bio-VAl-miRNAseq_comparison.RData")

##########################
# 11 Save R session info #
##########################

devtools::session_info()

#######
# End #
#######