######################
# Data visualization #
######################

## This script must be run after the deseq2.R script

## Loading required libraries
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)

## Data import
dedata <- read_csv("data/DESeq2_Results.csv")

##############################
# Volcano plot visualization #
##############################

## Enhanced Volcano Function
EVolcano_Plotter <- function(input_data, p_cutOFF, fc_cutOFF, 
                             cus_labels, plt_title){
  EnhancedVolcano(input_data,
                  lab = input_data$Symbol,
                  x = "log2FoldChange",
                  y = "padj",
                  xlab = bquote(~Log[2]~ 'Fold Change'),
                  ylab = bquote(~-Log[10]~adjusted~italic(P)),
                  pCutoff = p_cutOFF,
                  FCcutoff = fc_cutOFF,
                  selectLab = cus_labels,
                  drawConnectors = TRUE,
                  colConnectors = "grey30",
                  labCol = 'black',
                  labSize = 8,
                  colAlpha = 0.7,
                  legendPosition = "None",
                  gridlines.major = FALSE,
                  gridlines.minor = FALSE,
                  pointSize = 3,
                  title = plt_title,
                  subtitle = "Volcano Plot",
                  caption = paste("adjP cutoff = ", p_cutOFF, 
                                  " and ", "FC CutOff = ", fc_cutOFF),
                  border = "full",
                  col = c("#808080", "#808080", "#808080", "#DC143C")
  )
}

## Enhanced Volcano visualization and export
cus_labels <- c("esg", "pros", "nub", "bbg", "Myo61F", "Patr-1", "ths", "vtd")
plot_title <- "Patr-1 RNAi vs Control"
vol <- EVolcano_Plotter(dedata, 0.05, 1.0, cus_labels, plot_title)

pdf("data/volcanoPlot.pdf", width=6, height = 9)
vol
dev.off()

##########################
# Heat Map visualization #
##########################

## Create tibbles including row names
sample_names <- coldata %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()

## Recover normalized read counts
normalized_counts <- counts(dds,normalized=TRUE) %>% 
  data.frame() %>% 
  tibble::rownames_to_column(var="Gene") %>% 
  as_tibble()

## Subset significantly changed genes
sigGenes <- dedata %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 1) %>%
  filter(threshold == TRUE)

## Extract normalized expression for significant genes
norm_OEsig <- normalized_counts[,c(1:7)] %>% 
  filter(Gene %in% sigGenes$Symbol) %>% 
  data.frame() %>% 
  column_to_rownames(var = "Gene") 
norm_OEsig <- norm_OEsig[, c(4:6, 1:3)]

### Set a color palette
heat_colors <- brewer.pal(9, "RdYlBu")
annoCol <- data.frame(ID = factor(rep(c("esgTS","Patr1 RNAi"), each=3)))
rownames(annoCol) <- colnames(norm_OEsig[1:6])
annoCols <- list(ID = c("esgTS" = "darkgreen", "Patr1 RNAi"="purple3"))

### Annotate the heatma
annotation <- sample_names %>% 
  dplyr::select(samplename, genotype) %>% 
  data.frame(row.names = "samplename")

### Run pheatmap
heatMap <- pheatmap(norm_OEsig, 
                    color = heat_colors, 
                    gaps_col =  3, 
                    cluster_rows = TRUE, 
                    cluster_cols = FALSE,
                    show_rownames = FALSE, 
                    show_colnames = TRUE, 
                    annotation = annotation, 
                    annotation_col = annoCol, 
                    annotation_colors = annoCols[1],
                    fontsize = 12, 
                    scale = "row", 
                    fontsize_row = 14, 
                    border_color = "black")

## Export the heatmap
save_pheatmap_pdf <- function(x, filename, width=6, height=9) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(heatMap, "Heat_Map.pdf")



