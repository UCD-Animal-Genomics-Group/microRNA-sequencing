################################################################
#     miRNA RT-qPCR statistical analysis of a time course      #
#             experimental infection in cattle                 #
#                      Part 2 of 2:                            #
#  -- all post-infection time points vs -2 wk pre-infection -- #
################################################################

# Based on the workflow created by Nalpas, Nicolas and Correia, Carol (2015)
# DOI badge: https://doi.org/10.5281/zenodo.16164

# Authors of current version (2.0.0): Correia, C.N. and Nalpas, N.C.
# DOI badge of current version:
# Last updated on 29/01/2018

############################################
# 01 Load and/or install required packages #
############################################

# Source the common functions used across this script
source(file = "General_function.R")

# Define method
method <- "pre2"

# Load the required packages
library(ggplot2)
library(grid)
library(RColorBrewer)
library(magrittr)
library(psych)
library(readxl)
library(readr)
library(extrafont)

# Set time zone
Sys.setenv(TZ = "Europe/London")

# Register fonts with R for the PDF output device
loadfonts()

###################################
# 02 Read in input files within R #
###################################

# Read in the excel RT-qPCR expression file

PCR <- data.frame(read_excel("cnrq_mean.xlsx",
                             col_names = TRUE,
                             na = "NaN"))

# Clean up the column name in the qPCR dataset
colnames(PCR) %<>% gsub(
  pattern = "(^X|^miR\\.|^miR)(\\d)", replacement = "miR\\2", x = .,
  perl = TRUE) %>% gsub(pattern = "\\.SE.*", replacement = "_SE", x = .,
                        perl = TRUE) %>% gsub(
                          pattern = "\\.CNRQ.*", replacement = "_CNRQ", x = .,
                          perl = TRUE) %>% gsub(
                            pattern = "^(hs|bt)a\\.", replacement = "", x = .,
                            perl = TRUE) %>% gsub(
                              pattern = "\\.", replacement = "_", x = .,
                              perl = TRUE)
head(PCR)

# Read in and format the excel sample code name file
sample <- data.frame(read_excel("Sample_code.xlsx", col_names = TRUE))
colnames(sample) <- c("animal", "time_point", "Samples")
sample$time_point <- gsub(pattern = "-", replacement = "pre",
                          x = sample$time_point)
head(sample)

# Merge the sample information with the PCR data
PCR.data <- as.data.frame(merge(x = sample, y = PCR, by.x = "Samples",
                                by.y = "Samples", all.x = TRUE))
head(PCR.data)
rownames(PCR.data) <- paste(PCR.data$animal, PCR.data$time_point, sep = "_")
head(PCR.data)

##################################################
# 03 Calculate the log fold-change in expression #
##################################################

# Log 2 transform the CNRQ value
PCR.data <- data.frame(PCR.data[, 1:(ncol(sample))], log(
  x = PCR.data[, grep(pattern = "_CNRQ", x = colnames(PCR.data))],
  base = 2), row.names = row.names(PCR.data))
colnames(PCR.data) <- gsub(pattern = "_(CNRQ)", replacement = "_log\\1",
                           x = colnames(PCR.data))
head(PCR.data)
dim(PCR.data)

# Define the variables required to compute fold-change in expression
targets.animal <- unique(PCR.data$animal)
targets.time <- unique(PCR.data$time_point)
animal.rownames <- lapply(X = unique(PCR.data$animal), FUN = function(x) rep(
  x = x, length(unique(PCR.data$time_point)))) %>% unlist()
time.rownames <- rep(targets.time, length(targets.animal))
full.rownames <- paste(animal.rownames, time.rownames, sep = "_")
logFC <- data.frame(row.names = full.rownames)
logFC <- cbind(logFC, matrix(data = unlist(x = strsplit(
  x = row.names(logFC), split = "_", fixed = TRUE)), nrow = length(rownames(
    logFC)), ncol = 2, byrow = TRUE))
colnames(logFC) <- c("animal", "time_point")
head(logFC)
dim(logFC)

# Compute the log fold-change in expression versus pre2 week
for (gene in colnames(PCR.data)[(ncol(sample)+1):ncol(PCR.data)]){
  col.data <- c()
  for (animal in targets.animal){
    for (t.target in targets.time){
      data.ref <- PCR.data[PCR.data$time_point == 'pre2' &
                             PCR.data$animal == animal, gene]
      data.target <- PCR.data[PCR.data$time_point == t.target &
                                PCR.data$animal == animal, gene]
      col.data <- c(col.data, data.target-data.ref)
    }
  }
  logFC[, gsub(pattern = "CNRQ", replacement = "FC", x = gene)] <- col.data
}
head(logFC)

###########################################################
# 04 Assess normal distribution of PCR data for each gene #
###########################################################

# Use the Shapiro-Wilk test on overall data
shapiro <- apply(X = logFC[, ncol(sample):ncol(logFC)],
                 MARGIN = 2, FUN = function(x) shapiro.test(x = x))

# Plot the Q-Q plots for the overall data
qqplot <- apply(X = logFC[, ncol(sample):ncol(logFC)], MARGIN = 2,
      FUN = function(x) qqnorm(y = x, main = colnames(x)))

###############################################################
# 05 Compute significance values of fold-change in expression #
###############################################################

# Prepare dataframe to include all fold-change and significance evaluation
gene.list <- colnames(logFC)[ncol(sample):ncol(logFC)] %>%
lapply(X = ., FUN = function(x) rep(
  x = x, length(targets.time))) %>% unlist()
full.rownames <- paste(gene.list, time.rownames, sep = ".") %>% gsub(
  pattern = "_logFC.", replacement = ".", x = .)
full.colnames <- c("gene", "time_point", "shapiro", "n", "meanlogFC", "median",
                   "sd", "se", "min", "max", "t.Pvalue", "w.Pvalue",
                   "final.Pvalue")
sig <- matrix(nrow = length(full.rownames), ncol = length(full.colnames))
sig <- data.frame(x = sig, row.names = full.rownames)
colnames(sig) <- full.colnames
sig[, c("gene", "time_point")] <- matrix(data = unlist(x = strsplit(
  x = row.names(sig), split = ".", fixed = TRUE)), nrow = length(rownames(
    sig)), ncol = 2, byrow = TRUE)

# Compute the descriptive statistics and Pvalues of the log Fold-changes
shapiro.time <- list()
t.value <- list()
w.value <- list()

targets.time <- c("1", "2", "6", "10", "12", "pre2", "pre1")

for (gene in gene.list){
  for (t.target in targets.time){
    x.target <- logFC[logFC$time_point == t.target  &
                        logFC$animal %in% targets.animal, gene]
    ref <- logFC[logFC$time_point == "pre2"  &
                   logFC$animal %in% targets.animal, gene]
    stat.value <- describe(x.target)
    if (t.target == "pre2") {
      shap.eval$p.value <- list("NaN")
    }
    else {
      shap.eval <- shapiro.test(x = x.target)
    }
    shapiro.time[t.target] <- list(shap.eval)
    t.eval <- t.test(x = x.target, y = ref, alternative = "two.sided",
                     paired = TRUE, conf.level = 0.95, na.action = omit())
    t.value[t.target] <- list(t.eval)
    w.eval <- wilcox.test(x = x.target, y = ref,alternative = "two.sided",
                          paired=TRUE, na.action = omit())
    w.value[t.target] <- list(w.eval)
    if (shap.eval == "NaN") {
      final.pvalue <- "NaN"
    }
    else if (shap.eval < 0.1) {
      final.pvalue <- w.eval$p.value
    }
    else {
      final.pvalue <- t.eval$p.value
    }
    gene.val <- paste(gene, t.target, sep = ".") %>% gsub(
      pattern = "_logFC.", replacement = ".", x = .)
    sig[gene.val, c("shapiro", "n", "meanlogFC", "median", "sd", "se", "min",
                    "max", "t.Pvalue", "w.Pvalue", "final.Pvalue")] <- c(
                      shap.eval$p.value, stat.value$n, stat.value$mean,
                      stat.value$median, stat.value$sd, stat.value$se,
                      stat.value$min, stat.value$max, t.eval$p.value,
                      w.eval$p.value, final.pvalue)
  }
}

head(sig)
dim(sig)

#####################################
# 06 Plot fold-change in expression #
#####################################

# Add significance label
sig <- sig_label(arg1 = sig, arg2 = "final.Pvalue")
sig <- cbind(sig, Methods = rep(x = "RT-qPCR", times = length(rownames(sig))))
sig$time_point <- as.numeric(gsub(pattern = "pre", replacement = "-",
                                  x = sig$time_point))
sig$meanlogFC <- as.numeric(sig$meanlogFC)
head(sig)

# Plot the expression data
for (i in unique(sig$gene)) {
  file <- paste(i, "pdf", sep = ".")
  dat <- sig[sig$gene == i, c("time_point", "meanlogFC", "se",
                              "Significance_label", "Methods")]
  plot1 <- ggplot(data = dat, aes(x = dat$time_point, y = dat$meanlogFC,
                                  colour = dat$Methods)) +
    geom_point(size = 3) + geom_line(size = 0.5) + geom_text(
      aes(x = dat$time_point, y = (dat$meanlogFC + dat$se),
          label = dat$Significance_label), size = 8) +
    scale_x_discrete(limits = dat$time_point) +
    geom_errorbar(aes(
        x = dat$time_point, ymin = (dat$meanlogFC - dat$se),
        ymax = (dat$meanlogFC + dat$se)),
        width = 0.4, size = 0.5) +
    theme_bw(base_size = 14, base_family = "Calibri") +
    ggtitle(label = i) +
    xlab("Time point (weeks) vs. -2 wk") +
    ylab(expression(paste(log[10], "fold-change"))) +
    scale_colour_discrete(name = "Methods")
  cairo_pdf(filename = file, width = 6, height = 6)
  print(plot1)
  dev.off()
}

# Read in the excel file containing miRNA gene information
gene.info <- data.frame(read_excel("miRNA_RTqPCR_gene.xlsx",
                        col_names = TRUE,
                        na = "NA"))
head(gene.info)

# Include the gene information with the differential expression results
sig <- merge(x = gene.info, y = sig,
             by.x = "Primer_name", by.y = "gene",
             all.y = TRUE)
head(sig)

# Output the RT-qPCR differential expression results
write_csv(sig,
          "RT-qPCR_DE_vs_pre2.csv",
          col_names = TRUE,
          na = "NA")

#######################
# 07 Save .RData file #
#######################

save.image(file = paste0("RT-qPCR_", method, ".RData", sep = ""))

##########################
# 08 Save R session info #
##########################

devtools::session_info()

#######
# END #
#######