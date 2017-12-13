#############################################################
#############################################################
##                                                         ##
##  miRNA-seq data analysis of sense counts (paired data)  ##
##                                                         ##
#############################################################
#############################################################

# Analysis of 10 animals infected over a 15 weeks time course for which
# miRNA-seq libraries where prepared from serum samples

#############################
# List of required packages #
#############################

# Source the common functions used across this script
source(file = paste("General_function.R", sep = ""))

# Define variables for working and file directories
workDir <- getwd()
workDir
imgDir <- paste0(workDir, "/Figures")
tablesDir <- paste0(workDir, "/Tables")

# Load the required packages
library(edgeR)
library(tidyverse)
library(devtools)
library(magrittr)
library(stringr)
library(forcats)
library(MASS)
library(biobroom)
library(ggjoy)
library(PerformanceAnalytics)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(reshape2)
library(gplots)
library(ggrepel)
library(cowplot)
library(Cairo)
library(VennDiagram)


install.packages("PerformanceAnalytics")
install.packages("gplots")

#######################################################################
# Analysis of counts data obtained via different pipelines (pick one) #
#######################################################################

# Define the method used for count generation ("Novoalign", "miRdeep2"
# or "miRdeepstar")
method <- "Novoalign"
#method <- "miRdeep2"


############################################################
# Gene annotation using information obtained from GTF file #
############################################################

# Read in the annotation information
miRNA.info <- read.table(file = "miRNA_Btaurus.txt", header = TRUE,
                         sep = "\t", quote = "")

# Determine which miRNA have identical mature sequence
miRNA.duplicate <- aggregate(gene_id ~ sequence, FUN = "as.vector",
                             data = miRNA.info, na.action = "as.vector")

# Combine new info with miRNA annotation
miRNA.info <- merge(x = miRNA.info, y = miRNA.duplicate, by.x = "sequence",
                    by.y = "sequence", all = TRUE)
miRNA.info <- miRNA.info[, c(2:7, 1, 8:ncol(miRNA.info))]
colnames(miRNA.info)[c(1,7,12)] <- c("gene_id", "sequence",
                                     "identical_sequence")
head(miRNA.info)
dim(miRNA.info)

#####################################################################
# Read in and concatenate input files within R for method Novoalign #
#####################################################################

# Reads and merges a set of files containing counts according to the method
{
  if (method == "Novoalign") {
    Count <- readDGE(files = list.files(
      path = "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/TIDA/Bioinformatics/Carol_analyses/miRNAseq/featC_mature_counts",
      pattern = "6"),
      path = "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/TIDA/Bioinformatics/Carol_analyses/miRNAseq/featC_mature_counts",
      columns = c(1,7),
      skip = 1, header = TRUE)
    Count <- data.frame(gene_id = rownames(Count), Count$counts)
    colnames(Count) <- gsub(pattern = "^X", replacement = "",
                            x = colnames(Count), perl = TRUE)
  }
  else if (method == "miRdeep2") {
    files <- list.files(path = "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/TIDA/Bioinformatics/Carol_analyses/miRNAseq/quant_mature_counts",
                        pattern = "6")
    for (i in 1:length(files)) {
      Dat <- read.table(file = paste("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/TIDA/Bioinformatics/Carol_analyses/miRNAseq/quant_mature_counts",
                                     files[i], sep = "/"), quote = "")
      Dat <- Dat[, c(1:3)]
      sample <- gsub(pattern = "_expressed.csv", replacement = "",
                     x = files[i])
      colnames(Dat) <- c("gene_name", sample, "precursor_name")
      Dat <- merge(x = Dat, y = miRNA.info, by = c("gene_name",
                                                   "precursor_name"))
      Dat <- Dat[, c("gene_id", sample)]
      if (i == 1) {
        Count <- Dat
      }
      else {
        Count <- merge(x = Count, y = Dat, by = "gene_id")
      }
    }
  }
  else {
    print("Error: The provided method for count generation is unknown!")
  }
}

# Merge the annotation information with the count table
Count <- merge(x = miRNA.info, y = Count, by = "gene_id", all = TRUE)
head(Count)

# Check the merged matrix in terms of size and content
table(duplicated(Count[, "gene_id"]))
length(unique(Count[, "gene_id"], incomparables = FALSE))

# Check the raw data
head(Count)
tail(Count)
dim(Count)

# Ouptut samples data
write.matrix(x = Count, file = paste(method, "_counts.txt", sep = ""),
             sep = "\t")

###############################
# Create groups and a DGElist #
###############################

# Create a target matrix and a experimental group vector and animal block
target <- colnames(Count[, (ncol(Count)-69):ncol(Count)])
target <- cbind(target, matrix(data = unlist(strsplit(
  x = target, split = "_", fixed = TRUE)), nrow = 70, ncol = 2, byrow = TRUE))
colnames(target) = c("sample", "animal", "time_point")
group <- factor(x = target[, "time_point"])
animal <- factor(x = target[, "animal"])

# Create a DGElist containing the group information
dgelist <- DGEList(counts = Count[, (ncol(Count)-69):ncol(Count)],
                   lib.size = NULL, norm.factors = NULL, group = group,
                   genes = Count[, 1:(ncol(Count)-70)], remove.zeros = FALSE)
rownames(dgelist$counts) <- rownames(
  dgelist$genes) <- dgelist$genes[, "gene_id"]
dgelist$genes["gene_id"] <- NULL
names(dgelist)
head(dgelist$samples)
head(dgelist$counts)
head(dgelist$genes)

# Create a variable containing the current DGElist
assign(x = paste(method, ".dgelist", sep = ""), value = dgelist)

# Ouptut samples data
write.table(x = dgelist$samples, file = paste(method, "_sample.txt", sep = ""),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Tidy DGElist and plot data
dgelist %>%
  tidy() %>%
  ggplot() +
  geom_density(aes(x     = log10(count + 1),
                   group = sample)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ylab("Density of raw gene counts per sample") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_raw


density_raw

# Export image
ggsave(paste(method, "_density_plot_raw_counts.png", sep = ""),
       plot      = density_raw,
       device    = "png",
       limitsize = FALSE,
       dpi       = 300,
       path      = imgDir)

#####################################
# Filtering of lowly expressed tags #
#####################################

# Identify genes with zero counts across all samples
dim(dgelist$counts[rowSums(dgelist$counts) == 0, ])
head(dgelist$counts[rowSums(dgelist$counts) == 0, ])

# Filter lowly expressed tags, retaining only tags with at least 50 counts per
# million in 10 or more libraries (10 libraries correspond to one time point)
dgelist.filt <- dgelist[rowSums(
  cpm(dgelist$counts) > 50) >= median(summary(group)), ]
dim(dgelist.filt$counts)

# Compute the library size
dgelist.filt$samples$lib.size <- colSums(dgelist.filt$counts)
head(dgelist.filt$samples)
head(dgelist$samples)

#######################################################
# Tidy DGElist for exploratory data analysis plotting #
#######################################################

tidy_dgelist <- tidy(dgelist.filt, addSamples = TRUE)

# Correct PPDb stimulation info
tidy_dgelist$group %<>%
  factor(levels = c("pre2", "pre1", "1", "2", "6", "10", "12"))

# Check factors
levels(tidy_dgelist$group)

# Add animal IDs
tidy_dgelist$animal <- tidy_dgelist$sample
tidy_dgelist$animal %<>%
  stringr::str_replace("_.*", "") %>%
  fct_inorder()

# Check factors
levels(tidy_dgelist$animal)

# Combine animal and time point (group) info for
# plotting labels
tidy_dgelist %<>%
  dplyr::mutate(labels = paste0(group, "_", animal))

tidy_dgelist$labels %<>%
  factor(levels = c("pre2_6511", "pre2_6514", "pre2_6520", "pre2_6522", "pre2_6526",
                    "pre2_6635", "pre2_6636", "pre2_6637", "pre2_6644", "pre2_6698",
                    "pre1_6511", "pre1_6514", "pre1_6520", "pre1_6522", "pre1_6526",
                    "pre1_6635", "pre1_6636", "pre1_6637", "pre1_6644", "pre1_6698",
                    "1_6511", "1_6514", "1_6520", "1_6522", "1_6526",
                    "1_6635", "1_6636", "1_6637", "1_6644", "1_6698",
                    "2_6511", "2_6514", "2_6520", "2_6522", "2_6526",
                    "2_6635", "2_6636", "2_6637", "2_6644", "2_6698",
                    "6_6511", "6_6514", "6_6520", "6_6522", "6_6526",
                    "6_6635", "6_6636", "6_6637", "6_6644", "6_6698",
                    "10_6511", "10_6514", "10_6520", "10_6522", "10_6526",
                    "10_6635", "10_6636", "10_6637", "10_6644", "10_6698",
                    "12_6511", "12_6514", "12_6520", "12_6522", "12_6526",
                    "12_6635", "12_6636", "12_6637", "12_6644", "12_6698"))

# Check factors
levels(tidy_dgelist$labels)

###########################################################################
# Quality check of libraries by plotting density of count after filtering #
###########################################################################

ggplot(tidy_dgelist, aes(x = log10(count + 1),
                          y = labels)) +
  scale_y_discrete(limits = rev(levels(tidy_dgelist$labels))) +
  geom_joy(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Time point",
                    values = c(rep("#b2b2b2", 2), rep("#e06377", 5))) +
  theme_bw(base_size = 12, base_family = "Calibri") +
  ggtitle("Density of filtered gene counts per sample") +
  ylab("Time point_Animal number") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> density_filt


density_filt

# Export high quality image
ggsave(paste(method, "_density-filt.pdf", sep = ""),
       plot      = density_filt,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

#################################################
# Gene expression correlation between libraries #
#################################################

# Perform gene expression correlation between libraries at each time points to
# identify potential outlier library (CPM values not required with Spearman)
for (val in levels(group)) {
  png(filename = paste(method, val, "cor.png", sep = "_"),
      width = 1500, height = 1500, units = "px")
  chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
    pattern = paste("_", val, "$", sep = ""), x = colnames(dgelist.filt$counts),
    perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
    main = paste(val, " time point expression correlation", sep = ""))
  dev.off()
}

########################################################
# Normalization of data using trimmed mean of M-values #
########################################################

# Calculate normalisation factor for our DGElist, note that with edgeR
# the counts are not transformed in any way after normalization
dgelist.norm <- calcNormFactors(dgelist.filt)
dgelist.norm$samples

#################################
# Plots: MDS at each time point #
#################################

# Define function for getting MDS coordinates
getBCVcoord <- function(dgelst, time_pattrn) {

  mds <- plotMDS.DGEList(x = dgelst[ , grepl(paste(time_pattrn,
                                                   collapse = "|"),
                                             x = colnames(dgelst))],
                         plot = FALSE,
                         method = "bcv")

  mds_coord <- mds$cmdscale.out # Get coords to plot with ggplot2

  mds_coord %<>% # Tidy coords
    tidy() %>%
    dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
    dplyr::mutate(animal = sample, group = sample)

  mds_coord$animal %<>% # Clean animal IDs for plotting
    str_replace(paste(time_pattrn,
                      collapse = "|"), "") %>%
    factor()

  mds_coord$group %<>% # Clean group info for plotting
    str_replace(".*_", "") %>%
    fct_rev()

  return(mds_coord)

}

# Get tidy MDS coordinates for each time point with both controls

## +1 wk and controls
W1_coord <- getBCVcoord(dgelist.norm, c("_pre2", "_pre1", "_1$"))

## +2 wk and controls
W2_coord <- getBCVcoord(dgelist.norm, c("_pre2", "_pre1", "_2"))

## +6 wk and controls
W6_coord <- getBCVcoord(dgelist.norm, c("_pre2", "_pre1", "_6"))

## +10 wk and controls
W10_coord <- getBCVcoord(dgelist.norm, c("_pre2", "_pre1", "_10"))

## +12 wk and controls
W12_coord <- getBCVcoord(dgelist.norm, c("_pre2", "_pre1", "_12"))


# Plot +1 wk MDS
MDS_W1 <- ggplot(W1_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+1 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Plot +2 wk MDS
MDS_W2 <- ggplot(W2_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+2 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Plot +6 wk MDS
MDS_W6 <- ggplot(W6_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+6 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Plot +10 wk MDS
MDS_W10 <- ggplot(W10_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+10 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Plot +12 wk MDS
MDS_W12 <- ggplot(W12_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  scale_colour_manual("Treatment",
                      values = c("#b2b2b2", "#b2b2b2", "#e06377")) +
  scale_shape_manual("Treatment",
                     values = c(19, 17, 15)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 14, base_family = "Calibri") +
  ggtitle("+12 wk") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

### Combine all MDS plots into single figure

# Set grid
MDS_grid <- plot_grid(MDS_W1, MDS_W2, MDS_W6,
                      MDS_W10, MDS_W12,
                      labels = c("A", "B", "C", "D", "E"),
                      ncol = 2,
                      scale = .96)

# Check plot
MDS_grid

# Export high quality image
ggsave(filename  = "Novoalign_MDS_grid.pdf",
       plot      = MDS_grid,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 20,
       width     = 18,
       units     = "in")

##############################################
# Create a design matrix for paired analysis #
##############################################

# Create a design matrix
design <- model.matrix(~ animal + group)
rownames(design) <- rownames(dgelist.norm$samples)
colnames(design) <- gsub(pattern = "(animal)|(group)", replacement = "",
                         x = colnames(design), perl = TRUE)
design

########################################################################
# Estimate the dispersion parameter for each tag using Cox-Reid method #
########################################################################

# Calculate the dispersion (common, trended and tagwise)
dgelist.disp <- estimateGLMCommonDisp(y = dgelist.norm, design = design,
                                      verbose = TRUE)
dgelist.disp <- estimateGLMTrendedDisp(y = dgelist.disp, design = design)
dgelist.disp <- estimateGLMTagwiseDisp(y = dgelist.disp, design = design)
names(dgelist.disp)

# Plot the dispersion
png(filename = paste(method, "_BCV.png", sep = ""),
    width = 1366, height = 768, units = "px")
plotBCV(dgelist.disp)
dev.off()

# Show the calculated dispersion and the coefficient of biological variation
dgelist.disp$common.dispersion
sqrt(dgelist.disp$common.dispersion)

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
dgeglm.fit <- glmFit(y = dgelist.disp, design = design)
names(dgeglm.fit)

#####################################
# Main differential expression call #
#####################################

# Test for differential expression for -1 week versus -2 week
multi.DE("pre1", 1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         dedata = paste(method, ".de.pre1w", sep = ""),
         smearfile = paste(method, "_smear_pre1w", sep = ""))

# Test for differential expression for 1 week versus -2 week and -1 week
multi.DE(1, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         dedata = paste(method, ".de.1w", sep = ""),
         smearfile = paste(method, "_smear_1w", sep = ""))

# Test for differential expression for 2 week versus -2 week and -1 week
multi.DE(2, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         dedata = paste(method, ".de.2w", sep = ""),
         smearfile = paste(method, "_smear_2w", sep = ""))

# Test for differential expression for 6 week versus -2 week and -1 week
multi.DE(6, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         dedata = paste(method, ".de.6w", sep = ""),
         smearfile = paste(method, "_smear_6w", sep = ""))

# Test for differential expression for 10 week versus -2 week and -1 week
multi.DE(10, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         dedata = paste(method, ".de.10w", sep = ""),
         smearfile = paste(method, "_smear_10w", sep = ""))

# Test for differential expression for 12 week versus -2 week and -1 week
multi.DE(12, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         dedata = paste(method, ".de.12w", sep = ""),
         smearfile = paste(method, "_smear_12w", sep = ""))

#########################################################
# Merge all DE call data from the different time points #
#########################################################

# Create a vector of DE table to merge
DEtable.name <- c("1w", "2w", "6w", "10w", "12w") %>% lapply(
  X = ., FUN = function(x) paste(method, "de", x, sep = '.')) %>% unlist(.)

# Merge all DE table for the different time points into a single dataframe
DEtable.merge(DEtable.name, output = paste(method, ".DE", sep = ""),
              pattern = paste("^", method, '.de.', sep = ""))

# Create a variable containing the current DEtable
DEtable <- eval(parse(text = paste(method, '.DE', sep = "")))
head(DEtable)

# Write into a table the full DE call data
write.matrix(x = DEtable, file = paste(method, "_full_DE.txt", sep = ""),
             sep = "\t")

##############################################
# Comparison of DE genes between time points #
##############################################

# Identify as a vector list the significant DE genes per time point
venn.de(data = DEtable,
        comparison = c("1w", "2w", "6w", "10w", "12w"),
        picname = paste(method, "_Venn_All", sep = ""),
        overlapname = paste(method, ".overlap", sep = ""),
        lwd = 0.7, cex = 1, cat.cex = 1)

# Create a variable containing the current overlap genes
overlap <- eval(parse(text = paste(method, ".overlap", sep = "")))
overlap <- grep(pattern = "_", x = overlap, value = TRUE, invert = TRUE)

##########################################################
# Identification of reference genes for RT-qPCR analysis #
##########################################################

# Merge verage CPM and Standard deviation for each gene together
gene.stable <- merge(
  x = as.matrix(apply(X = CPM, MARGIN = 1, FUN = function(x) mean(x))),
  y = as.matrix(apply(X = CPM, MARGIN = 1, FUN = function(x) sd(x))),
  by = "row.names", all = TRUE)
gene.stable <- merge(x = miRNA.info[, 1:2], y = gene.stable, by.x = "gene_id",
                     by.y = "Row.names")
rownames(gene.stable) <- gene.stable[, "gene_id"]
gene.stable <- gene.stable[, -1]
colnames(gene.stable)[2:3] <- c("meanCPM", "SdevCPM")

# Calculate percentage of Standard deviation in order to compare all genes
gene.stable <- cbind(gene.stable, Perc_Sdev = (
  gene.stable[, "SdevCPM"]*100/gene.stable[, "meanCPM"]))
gene.stable <- gene.stable[order(gene.stable$Perc_Sdev),]
head(gene.stable)

# Plot the mean CPM and Standard deviation CPM of all genes
ggplot(data = gene.stable, aes(x = log(x = gene.stable$meanCPM, base = 2),
                               y = log(x = gene.stable$SdevCPM,
                                       base = 2))) + geom_point(size = 4.5)

# Write into a table the analysis for reference gene discovery
write.matrix(x = gene.stable, file = paste(
  method, "_potential_refgene.txt", sep = ""), sep = "\t")

###########################################################
# Comparison of Exiqon PCR results with miRNA-seq results #
###########################################################

# Read in the csv file containing PCR results from Kirsten
PCR.Kirsten <- read.table(file = paste(
  user, "/Home_work_sync/Work/Colleagues shared work/Kirsten McLoughlin/",
  "miRNA_work/PCR_data/PCR_to_compare.txt", sep = ""), header = TRUE,
  sep = "\t", quote = "", na.strings = '\'N/A\'')

# Reformat Kirsten's PCR data
PCR.Kirsten <- PCR.Kirsten[, c(1,2,4,7,12)]
PCR.Kirsten <- cbind(PCR.Kirsten, as.matrix(as.numeric(lapply(X = PCR.Kirsten[
  , 'Geometric.mean.fold.change..relative.to..time.point.1.'],
  FUN = function(i) {if (i < 0){-log(x = abs(i), base = 2)} else {log(x = abs(
    i), base = 2)}}))))
colnames(PCR.Kirsten) <- c("hsapiens_gene_name", "sequence", "time_point",
                   "fold_change", "FDR", "logFC")
PCR.Kirsten$time_point %<>% gsub(pattern = "\\+(\\d*) wk", replacement = "\\1",
                                 x = .) %>% gsub(pattern = "(-\\d*) wk",
                                                 replacement = "\\1", x = .)
PCR.Kirsten <- cbind(PCR.Kirsten, Methods = rep(
  x = "Kirsten_PCR", times = length(rownames(PCR.Kirsten))))
head(PCR.Kirsten)

# Read the file containing biomarkers RT-qPCR results
RT.qPCR <- read.table(file = paste(user, "/Home_work_sync/Work/TIDA/miRNA_PCR",
                                   "/Diff_expr/RT-qPCR_DE.txt", sep = ""),
                      sep = "\t", header = TRUE)
rownames(RT.qPCR) <- RT.qPCR$gene_name %>% paste(
  ., RT.qPCR$time_point, sep = "_") %>% gsub(
    pattern = "-(\\d*)$", replacement = "pre\\1", x = .)
RT.qPCR <- RT.qPCR[, -c(1)]
head(RT.qPCR)

# Merge the miRNA-seq, RT-qPCR and Kirsten PCR data
Comp <- NULL
for (time in c(1, 2, 6, 10, 12)) {
  data1 <- RT.qPCR[RT.qPCR$time_point == time, c(
    "sequence", "meanlogFC", "se", "final.Pvalue", "Methods")]
  colnames(data1) <- c("sequence", "logFC_N.C.PCR", "se_N.C.PCR",
                       "Pvalue_N.C.PCR", "Methods_N.C.PCR")
  data2 <- PCR.Kirsten[PCR.Kirsten$time_point == time & !is.na(
    PCR.Kirsten$sequence), c("sequence", "hsapiens_gene_name", "logFC", "FDR",
                             "Methods")]
  colnames(data2) <- c("sequence", "hsapiens_gene_name", "logFC_K.PCR",
                       "FDR_K.PCR", "Methods_K.PCR")
  method.DE <- eval(parse(text = paste(method, ".DE", sep = "")))
  pattern.time <- colnames(method.DE) %>% grep(pattern = paste(
    time, "w", sep = ""), x = ., value = TRUE) %>% grep(
      pattern = "logFC|FDR", x = ., value = TRUE)
  pattern.id <- grep(pattern = "_", x = rownames(method.DE),
                     invert = TRUE)
  data3 <- method.DE[pattern.id, c("gene_id", "gene_name", "sequence",
                                        pattern.time[1], pattern.time[2])]
  colnames(data3) <- c("gene_id", "gene_name", "sequence", "logFC_RNAseq",
                       "FDR_RNAseq")
  if (is.null(Comp)) {
    Comp <- merge(x = data1, y = data2, by = "sequence", all = TRUE) %>% merge(
      x = data3, y = ., by = "sequence", all.y = TRUE)
    Comp <- cbind(Comp, time_point = rep(x = time, times = length(
      Comp$sequence)))
  }
  else {
    tmp <- merge(x = data1, y = data2, by = "sequence", all = TRUE) %>% merge(
      x = data3, y = ., by = "sequence", all.y = TRUE)
    tmp <- cbind(tmp, time_point = rep(x = time, times = length(tmp$sequence)))
    Comp <- rbind(Comp, tmp)
  }
}
Comp <- Comp[!is.na(Comp$sequence) & !is.na(Comp$gene_id),]
Comp <- Comp[!(is.na(Comp$logFC_N.C.PCR & is.na(Comp$logFC_RNAseq))), ]
Comp

# Perform overall correlation of log fold-change between PCR and RNA-seq data
cor.gene <- list()
cor.gene["logFC_RNAseq_vs_logFC_N.C.PCR"] <- list(cor.test(
  x = Comp$logFC_RNAseq, y = Comp$logFC_N.C.PCR, alternative = "two.sided",
  method = "spearman", na.action = na.omit()))
cor.gene["logFC_RNAseq_vs_logFC_K.PCR"] <- list(cor.test(
  x = Comp$logFC_RNAseq, y = Comp$logFC_K.PCR, alternative = "two.sided",
  method = "spearman", na.action = na.omit()))
cor.gene["logFC_N.C.PCR_vs_logFC_K.PCR"] <- list(cor.test(
  x = Comp$logFC_N.C.PCR, y = Comp$logFC_K.PCR, alternative = "two.sided",
  method = "spearman", na.action = na.omit()))

# Perform correlation of log fold-change between PCR and RNA-seq data per gene
for (g in unique(Comp[!is.na(Comp$logFC_N.C.PCR), "gene_id"])) {
  dat1 <- Comp[Comp$gene_id == g, "logFC_N.C.PCR"]
  dat2 <- Comp[Comp$gene_id == g, "logFC_RNAseq"]
  res <- cor.test(x = dat1, y = dat2, alternative = "two.sided",
                  method = "spearman", na.action = na.omit())
  cor.gene[g] <- list(res)
}

# Create a variable containing the current correlation results
assign(paste(method, ".cor.gene", sep = ""), cor.gene)

# Determine concordance of direction of expression between techniques
# and concordance of P-value significance between techniques
concordance <- list()
concor.fc <- 0
concor.Pval <- 0
counter <- 0
for (x in 1:nrow(Comp)) {
  if (!is.na(Comp$logFC_N.C.PCR[x])) {
    counter <- counter+1
    if (Comp[x, "logFC_N.C.PCR"] > 0 & Comp[x, "logFC_RNAseq"] > 0) {
      concor.fc <- concor.fc+1
    }
    else if (Comp[x, "logFC_N.C.PCR"] < 0 & Comp[x, "logFC_RNAseq"] < 0) {
      concor.fc <- concor.fc+1
    }
    if (Comp[x, "Pvalue_N.C.PCR"] > 0.05 & Comp[x, "FDR_RNAseq"] > 0.05) {
      concor.Pval <- concor.Pval+1
    }
    else if (Comp[x, "Pvalue_N.C.PCR"] < 0 & Comp[x, "FDR_RNAseq"] < 0) {
      concor.Pval <- concor.Pval+1
    }
  }
  if (x == nrow(Comp)) {
    concor.fc <- (concor.fc*100/counter)
    concor.Pval <- (concor.Pval*100/length(Comp$gene_id))
  }
}
concordance["N.C.PCR-RNAseq"] <- list(list(
  "logFC_concordance" = concor.fc, "Pvalue_concordance" = concor.Pval))

concor.fc <- 0
concor.Pval <- 0
counter <- 0
for (x in 1:nrow(Comp)) {
  if (!is.na(Comp$logFC_K.PCR[x])) {
    counter <- counter+1
    if (Comp[x, "logFC_K.PCR"] > 0 & Comp[x, "logFC_RNAseq"] > 0) {
      concor.fc <- concor.fc+1
    }
    else if (Comp[x, "logFC_K.PCR"] < 0 & Comp[x, "logFC_RNAseq"] < 0) {
      concor.fc <- concor.fc+1
    }
    if (Comp[x, "FDR_K.PCR"] > 0.05 & Comp[x, "FDR_RNAseq"] > 0.05) {
      concor.Pval <- concor.Pval+1
    }
    else if (Comp[x, "FDR_K.PCR"] < 0 & Comp[x, "FDR_RNAseq"] < 0) {
      concor.Pval <- concor.Pval+1
    }
  }
  if (x == nrow(Comp)) {
    concor.fc <- (concor.fc*100/counter)
    concor.Pval <- (concor.Pval*100/length(Comp$gene_id))
  }
}
concordance["K.PCR-RNAseq"] <- list(list(
  "logFC_concordance" = concor.fc, "Pvalue_concordance" = concor.Pval))

# Create a variable containing the current concordance results
assign(paste(method, ".concordance", sep = ""), concordance)

############################################
# Additional differential expression calls #
############################################

# Perform all differential expression analysis within a for loop
time <- c("pre2", "pre1", 1, 2, 6, 10, 12)
for (i in time) {
  print("Performing comparison versus:")
  print(i)
  DEtable.name <- c()
  for (j in 1:length(time)) {
    if (time[j] != i) {
      # Test for differential expression for comparison
      name <- paste(time[j], "wvs", i, "w", sep = "")
      DEtable.name <- c(DEtable.name, name)
      multi.DE(time[j], 1, i, -1, data = dgeglm.fit, design = design,
               group = group, adjpvalue = 0.05, method = "BH",
               dedata = paste(method, ".de.", name, sep = ""))
    }
    else {
      next
    }
  }
  # Merge all DE call data from the different time points
  # Create a vector of DE table to merge
  DEtable.method <- DEtable.name %>% lapply(X = ., FUN = function(x) paste(
    method, "de", x, sep = '.')) %>% unlist(.)
  # Merge all DE table for the different time points into a single dataframe
  DEtable.merge(DEtable.method, output = paste(method, ".DE.vs", i, "w",
                                             sep = ""),
                pattern = paste("^", method, '.de.', sep = ""))
  # Create a variable containing the current DEtable
  DEtable <- eval(parse(text = paste(method, ".DE.vs", i, "w", sep = "")),
                  envir = .GlobalEnv)
  colnames(DEtable)
  # Write into a table the full DE call data
  write.matrix(x = DEtable, file = paste(method, "_DE_vs", i, "w.txt",
                                         sep = ""),
               sep = "\t")
  # Plot the number of differentially expressed genes per comparisons
  plot.numb.DE(data = DEtable,
               comparison = DEtable.name,
               pattern = "vs", replace = " vs ",
               filename = paste(method, "_DE_vs", i, "w.tif", sep = ""))
  # Comparison of DE genes between time points
  # Identify as a vector list the significant DE genes per time point
  venn.de(data = DEtable,
          comparison = DEtable.name[-1],
          pattern = "vs", replace = " vs ",
          picname = paste(method, "_Venn_vs", i, "w", sep = ""),
          overlapname = paste(method, ".overlap.vs", i, "w", sep = ""),
          lwd = 0.7, cex = 1, cat.cex = 1)
}


##############################################################################
# Comparison of time overlapping DE genes between various analysis pipelines #
##############################################################################

# Define a set of colors
colours <- brewer.pal(n = 3, name = "Set1")

# Create the Venn diagram of overlapping DE genes between time points
venn.diagram(x = list("Novoalign" = Novoalign.overlap,
                      "miRdeep2" = miRdeep2.overlap,
                      "miRdeepstar" = miRdeepstar.overlap),
             filename = "Venn_3pipelines.tiff", na="remove",
             res = 600, fill = colours, cat.col = colours)

# Get all the gene expression information for the candidate biomarkers miRNA
biomarker <- merge(x = miRNA.info[(miRNA.info$gene_id %in% unique(c(
  Novoalign.overlap, miRdeep2.overlap, miRdeepstar.overlap))),],
  y = Novoalign.DE[, (ncol(Novoalign.DE)-24):ncol(Novoalign.DE)],
  by.x = "gene_id", by.y = "row.names", all.x = TRUE)
biomarker <- merge(x = biomarker,
                   y = miRdeep2.DE[, (ncol(miRdeep2.DE)-24):ncol(miRdeep2.DE)],
                   by.x = "gene_id", by.y = "row.names", all.x = TRUE)
biomarker <- merge(x = biomarker, y = miRdeepstar.DE[, (ncol(
  miRdeepstar.DE)-24):ncol(miRdeepstar.DE)], by.x = "gene_id",
  by.y = "row.names", all.x = TRUE)

# Format the variable columns
colnames(biomarker) %<>% gsub(pattern = "w$", replacement = "_miRdeepstar",
                              x = ., perl = TRUE) %>% gsub(
                                pattern = "\\.y$", replacement = "_miRdeep2",
                                x = ., perl = TRUE) %>% gsub(
                                  pattern = "\\.x$", replacement = "_novoalign",
                                  x = ., perl = TRUE)

# Write into a table all the time overlapping DE genes
write.matrix(x = biomarker, file = "Candidate_biomarkers.txt", sep = "\t")

####################
# Save .RData file #
####################

save.image(file = paste(method, ".RData", sep = ""))

#######################
# Save R session info #
#######################

devtools::session_info()

#######
# END #
#######