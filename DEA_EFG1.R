# Differential expression analysis for gene EFG1 in Candida Albicans
# R version 3.6.1


# import the necessary packages
knitr::opts_chunk$set(echo = TRUE)
auto_install.cvec <- c("easypackages") # All of the packages listed below could go here for automatic install instead
to.install <- auto_install.cvec[!(auto_install.cvec %in% installed.packages()[,"Package"])]
if(length(to.install)) install.packages(to.install)
library(easypackages)
request_install.cvec <- c("devtools","eulerr","tidyverse","lemon","colorspace")
packages(request_install.cvec,prompt=FALSE)
library("limma")
library("edgeR")

# Load the input dataset
DEAs_Morgan_Limma <- read.delim("~/Desktop/Aaron_Lab/DEAs_Morgan_Limma.txt", row.names=1)
DEA_Sample <- DEAs_Morgan_Limma[, c("WT_O1", "WT_O2", "WT_W1", "WT_W2", "EFG1_O1", "EFG1_O2", "EFG1_W1", "EFG1_W2")]
d0 <- DGEList(DEA_Sample)
d0 <- calcNormFactors(d0)

# Filtering genes with low expression
cutoff <- 2
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 

# Create the design matrix
sample<-c(rep("WT_O",2), rep("WT_W",2), rep("EFG1_O",2), rep("EFG1_W",2))
sample<-factor(sample)
design.mat <- model.matrix(~0+sample)
colnames(design.mat)<-levels(sample)

# Voom Plot
y <- voom(d, design.mat, plot = T)

# Fit the data to LIMMA
fit <- lmFit(y, design.mat)

# Create the contrast matrix
contrast <- makeContrasts(EFG1W = EFG1_W - WT_W, EFG1O = EFG1_O - WT_O, Diff = (EFG1_W - WT_W) - (EFG1_O - WT_O), WT = WT_W - WT_O, levels = coef(fit))

# Fit the contrast table and eBayes
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

# Call the topTable
deg.WT <- topTable(fit2, coef = "WT", number = nrow(DEA_Sample), adjust.method = "fdr", p.value = 0.05, lfc = log2(1.5))
WT <- dim(deg.WT)[1]
dim(deg.WT)[1]

deg.EFG1W <- topTable(fit2, coef = "EFG1W", number = nrow(DEA_Sample), adjust.method = "fdr", p.value = 0.05, lfc = log2(1.5))
EFG1W <- dim(deg.EFG1W)[1]
dim(deg.EFG1W)[1]

deg.EFG1O <- topTable(fit2, coef = "EFG1O", number = nrow(DEA_Sample), adjust.method = "fdr", p.value = 0.05, lfc = log2(1.5))
EFG1O <- dim(deg.EFG1O)[1]
dim(deg.EFG1O)[1]

deg.Diff <- topTable(fit2, coef = "Diff", number = nrow(DEA_Sample), adjust.method = "fdr", p.value = 0.05, lfc = log2(1.5))
Diff <- dim(deg.Diff)[1]
dim(deg.Diff)[1]

# Calculate intersects among the comparison groups
WT_EFG1W <- length(intersect(rownames(deg.WT), rownames(deg.EFG1W)))
WT_EFG1O <- length(intersect(rownames(deg.WT), rownames(deg.EFG1O)))
WT_Diff <- length(intersect(rownames(deg.WT), rownames(deg.Diff)))
EFG1O_EFG1W <- length(intersect(rownames(deg.EFG1O), rownames(deg.EFG1W)))
Diff_EFG1W <- length(intersect(rownames(deg.Diff), rownames(deg.EFG1W)))
Diff_EFG1O <- length(intersect(rownames(deg.Diff), rownames(deg.EFG1O)))

WT_EFG1O_EFG1W <- length(intersect(intersect(rownames(deg.WT), rownames(deg.EFG1W)), rownames(deg.EFG1O)))
EFG1O_EFG1W_Diff <- length(intersect(intersect(rownames(deg.Diff), rownames(deg.EFG1W)), rownames(deg.EFG1O)))
WT_EFG1O_Diff <- length(intersect(intersect(rownames(deg.WT), rownames(deg.EFG1O)), rownames(deg.Diff)))
WT_EFG1W_Diff <- length(intersect(intersect(rownames(deg.WT), rownames(deg.EFG1W)), rownames(deg.Diff)))

WT_EFG1W_Diff_EFG1O <- length(intersect(intersect(intersect(rownames(deg.WT), rownames(deg.EFG1W)), rownames(deg.Diff)), rownames(deg.EFG1O)))


# Euler plot
combo <- c("WT" = WT, "EFG1W" = EFG1W, "EFG1O" = EFG1O, "Diff" = Diff, "WT&EFG1W" = WT_EFG1W, "WT&EFG1O" = WT_EFG1O, "WT&Diff" = WT_Diff, "EFG1O&EFG1W" = EFG1O_EFG1W, "Diff&EFG1W" = Diff_EFG1W, "Diff&EFG1O" = Diff_EFG1O, "WT&EFG1O&EFG1W" = WT_EFG1O_EFG1W, "EFG1O&EFG1W&Diff" = EFG1O_EFG1W_Diff, "WT&EFG1O&Diff" = WT_EFG1O_Diff, "WT&EFG1W&Diff" = WT_EFG1W_Diff, "WT&EFG1O&EFG1W&Diff" = WT_EFG1W_Diff_EFG1O)
plot(euler(combo), quantities = TRUE, main = "Differentially Expressed Genes")
