#########################################
# Differential Gene Expression Analysis #
#########################################

## Loading required libraries
library(DESeq2)
library(tidyverse)
library(gprofiler2)

## Data import
countdata <- read.csv("data/exon.csv")
exclude_genes <- read.csv("data/exclude_genes.csv")

## Remove rRNA/tRNA/noncoding/pseudogenes
countdata <- countdata %>%
  anti_join(exclude_genes, by = "Geneid")
rownames(countdata) <- countdata[,1]
countdata[,1] <- NULL

## Data preprocessing
countdata <- countdata %>%
  select(patr1_T1:esgTS_T3) %>%
  as.matrix()

## Assign condition
(genotype <- factor(c(rep("patr1_i", 3), rep("esgTS", 3))))

## Coldata dataframe
(coldata <- data.frame(row.names=colnames(countdata), genotype))

## Instantiate the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=countdata, 
  colData=coldata, 
  design=~genotype)

## Pre-filtering based on read counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

## Retrieve results
res <- results(dds, contrast = c("genotype", "patr1_i", "esgTS"))
resdata <- merge(as.data.frame(res), 
                   as.data.frame(counts(dds, normalized=TRUE)), 
                   by="row.names", 
                   sort=FALSE)
names(resdata)[1] <- "Symbol"

## Results arrangement and export
resdata <- resdata %>%
  as_tibble() %>%
  arrange(padj) %>% write_csv(file = "data/DESeq2_Results.csv")

## Export significantly upregulated and downregulated genes
up <- resdata %>%
  filter(log2FoldChange > 1 & padj < 0.05) %>%
  select(Symbol, baseMean, log2FoldChange, padj) %>%
  write_csv(file = "data/Upregulated_genes.csv")

down <- resdata %>%
  filter(log2FoldChange < -1 & padj < 0.05) %>%
  select(Symbol, baseMean, log2FoldChange, padj) %>%
  write_csv(file = "data/Downregulated_genes.csv")

###########################
# PCA Plot viasualization #
###########################

## Transform count data using variance stabilizing transformation (VST)
vsd <- varianceStabilizingTransformation(dds)

## Generate first two PCs, PC1 and PC2
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

## PCA plot visualization with ggplot2
pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, color=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme(axis.title.x = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position="right", legend.text = element_text(size = 14), 
        legend.title = element_text(face = "bold", size = 16))

pdf("data/PCA_plot.pdf", width=10, height = 12)
pcaplot
dev.off()

#########################
# Intersample variances #
#########################

## Transform count data using regularized logarithm (rlog)
expmatrix_DESeq <- rlog(dds, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)
head(expmatrix)

## Control correlogram and export
ctrlCorr <- expmatrix[,grep("esgTS", colnames(expmatrix))]
pdf("data/Ctrl_Correlogram.pdf", width=5, height=5)
corrgram::corrgram(ctrlCorr, order=TRUE, lower.panel=corrgram::panel.pie,
                   upper.panel=corrgram::panel.pts, text.panel=corrgram::panel.txt,
                   main="Correlogram of Controls",
                   col.regions=colorRampPalette("#1E90FF"))
dev.off()

## Experimental correlogram and export
expCorr <- expmatrix[,grep("patr1", colnames(expmatrix))]
pdf("data/Exp_Correlogram.pdf", width=5, height=5)
corrgram::corrgram(expCorr, order=TRUE, lower.panel=corrgram::panel.pie,
                   upper.panel=corrgram::panel.pts, text.panel=corrgram::panel.txt,
                   main="Correlogram of Patr-1 RNAi",
                   col.regions=colorRampPalette("#DC143C"))
dev.off()

##########################################
# Gene Ontology Analysis with gProfiler2 #
##########################################

go_up <- gost(up$Symbol, organism = "dmelanogaster", evcodes = TRUE,
              correction_method = "fdr") %>% .$result %>%
  select(c(source, term_id, term_name, p_value, 
           intersection_size, intersection)) %>%
  filter(p_value < 0.05) %>%
  write_csv("data/go_upregulated.csv")

go_down <- gost(down$Symbol, organism = "dmelanogaster", evcodes = TRUE,
                correction_method = "fdr") %>% .$result %>%
  select(c(source, term_id, term_name, p_value, 
           intersection_size, intersection)) %>%
  filter(p_value < 0.05) %>%
  write_csv("data/go_downregulated.csv")
