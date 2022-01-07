library(pheatmap)
library(dplyr)
library(mosaic)
library(gplots)
library(ggplot2)
library(Kendall)
library(ggpubr)
library(tidyverse)
library(dendsort)
library(DESeq2)
library(ggrepel)

##### import data #####

# import raw nanostring counts, meta data, gene annotations
rawcounts_path <- "melanomaIL2_rawData.xlsx"
metadata_path <- "sample_annotations.txt"
probe_path <- "ProbeAnnotations.txt"
pathways_path <- "genes.txt"
ihc_path <- "ihc_scores.txt"
rawcounts <- read.xlsx(rawcounts_path, colNames = TRUE, rowNames = TRUE)
metadata <- read.table(metadata_path, header=TRUE, row.names=1)
probe <- read.table(probe_path, header=TRUE, row.names=1, sep = "\t")
pathways <- read.table(pathways_path, sep = "\t", header = TRUE, row.names = NULL, stringsAsFactors = F)
ihc <- read.table(ihc_path, sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = F)
rawcounts <- as.matrix(rawcounts)
identical(colnames(rawcounts), rownames(metadata))

##### differential expression #####

# run deseq2 - lesion response comparisons
dds <- DESeqDataSetFromMatrix(countData=rawcounts, colData=metadata, design=~Lesion_response)
dds$Lesion_response <- relevel(dds$Lesion_response, "UT")
dds <- DESeq(dds)
res_nonCR <- results(dds, contrast = c("Lesion_response", "non_CR", "UT"))
res_CR <- results(dds, contrast = c("Lesion_response", "CR", "UT"))
res_nonCR_Ordered <- res_nonCR[order(res_nonCR$pvalue),]
res_CR_Ordered <- res_CR[order(res_CR$pvalue),]
res_nonCR_Ordered <- data.frame(GeneID = rownames(res_nonCR_Ordered), res_nonCR_Ordered)
res_CR_Ordered <- data.frame(GeneID = rownames(res_CR_Ordered), res_CR_Ordered)
#write.csv(res_nonCR_Ordered, file = "res_nonCR_Ordered.csv")
#write.csv(res_CR_Ordered, file = "res_CR_Ordered.csv")

# run deseq2 - patient comparisons UT only
untreated <- c("0_1", "1_1", "2_1", "3_1", "4_1", "4_2", "5_1", "6_1", "6_4")
rawcounts_UT <- rawcounts[, untreated]
metadata_UT <- data.frame(t((t(metadata))[, untreated]))
dds <- DESeqDataSetFromMatrix(countData=rawcounts_UT, colData=metadata_UT, design=~Patient_response)
dds$Patient_response <- relevel(dds$Patient_response, "ER")
dds <- DESeq(dds)
res <- results(dds, contrast = c("Patient_response", "NR", "ER"))
res_Ordered <- res[order(res$padj),]
res_Ordered <- data.frame(GeneID = rownames(res_Ordered), res_Ordered)
#write.csv(res_Ordered, file = "res_Ordered.csv")

##### gene expression heatmaps #####

# make heatmap of DE genes - lesion response comparisons
housekeeping <- filter(probe, Is_Control == 1)
rawcounts_filtered <- rawcounts[!(row.names(rawcounts) %in% housekeeping$Official_Gene_Name), ]
# filter significant DE genes, reorder by log2FoldChange (includes housekeeping), add significant non-CR genes to end
DE_CR_filter <- filter(res_CR_Ordered, as.numeric(res_CR_Ordered$padj) < 0.0001)
DE_nonCR_filter <- filter(res_nonCR_Ordered, as.numeric(res_nonCR_Ordered$padj) < 0.0001)
DE_CR_filter <- DE_CR_filter[order(DE_CR_filter$log2FoldChange, decreasing=TRUE),]
DE_genes <- rbind(DE_CR_filter, DE_nonCR_filter)
# subset DE genes
all_subset <- subset(rawcounts_filtered, rownames(rawcounts_filtered) %in% rownames(DE_genes))
# reorder subset by log2FoldChange
all_subset <- All_subset[order(match(rownames(All_subset), rownames(DE_genes))), ]
# take log2, zscore of counts
all_log2_counts <- log2(all_subset)
all_z_score <- apply(all_log2_counts, 1, zscore)
all_z_score <- data.frame(t(all_z_score))
# prep column annotations
annot_col <- data.frame(Lesion_response = metadata[,4], row.names = row.names(metadata))
annot_color <- list(Lesion_response = c(UT="grey", CR="green", non_CR="orange"))
# color, column labels
my_breaks=seq(-2,2, by=0.1)
my_color <- colorpanel(n=length(my_breaks)-1,low="blue",mid="grey88",high="red")
# create heatmap
#pdf("lesion_response_heatmap.pdf", 10, 10)
pheatmap(all_z_score, legend = T, color = my_color, breaks = my_breaks, show_colnames = T, show_rownames = T, 
         cluster_rows = F, cluster_cols = T, annotation_col = annot_col, 
         annotation_colors = annot_color, annotation_legend = T, angle_col = 90, border_color = "white", 
         main = "0.0001_L2FC")
#dev.off()

# make heatmap of DE genes - patient comparisons UT only
DE_filter <- filter(res_Ordered, padj < 0.05)
DE_filter <- DE_filter[order(DE_filter$log2FoldChange, decreasing=FALSE),]
# remove housekeeping genes
counts_filtered <- rawcounts_UT[!(row.names(rawcounts_UT) %in% housekeeping$Official_Gene_Name), ]
# calculate z-score
log_counts <- log(counts_filtered)
z_score <- apply(log_counts, 1, zscore)
z_score <- data.frame(t(z_score))
#subset DE genes (FDR<0.05), reorder by log2FC
subset <- subset(z_score, rownames(rawcounts_UT) %in% rownames(DE_filter))
subset <- subset[order(match(rownames(subset), rownames(DE_filter))), ]
subset <- subset[1:96,]
# prep row annotations
row_annot <- data.frame(matrix(0, ncol=6, nrow=770))
row.names(row_annot) <- row.names(rawcounts_UT)
row_annot <- row_annot %>% rename(gamma = X1, alpha = X2, il2 = X3, antigen = X4, dysfunction = X5, TLS = X6)
row_annot$gamma <- row.names(row_annot) %in% pathways$gamma_genes
row_annot$alpha <- row.names(row_annot) %in% pathways$alpha_genes
row_annot$il2 <- row.names(row_annot) %in% pathways$il2_genes
row_annot$antigen <- row.names(row_annot) %in% pathways$antigen_presentation
row_annot$dysfunction <- row.names(row_annot) %in% pathways$dysfunction
row_annot$TLS <- row.names(row_annot) %in% pathways$TLS
cols <- sapply(row_annot, is.logical)
row_annot[,cols] <- lapply(row_annot[,cols], as.numeric)
# prep column annotations
annot_col <- data.frame(Patient_response = metadata[,3], row.names = row.names(metadata))
annot_color <- list(Patient_response = c(ER="green", NR="red"))
# color, column labels
breaks=seq(-2,2, by=0.09)
my_color <- colorpanel(n=length(breaks)-1,low="blue",mid="grey88",high="red")
labels_col <- c("0_1", "1_1", "2_1", "3_1", "4_1", "4_2", "5_1", "6_1", "6_4")
# create heatmap
#pdf("patient_response_UT_heatmap.pdf", 10, 10)
pheatmap(subset, legend = T, color = my_color, show_colnames = T, show_rownames = F, labels_col = labels_col,
         cluster_rows = F, cluster_cols = T, annotation_col = annot_col, annotation_row = row_annot, 
         annotation_colors = annot_color, annotation_legend = T, angle_col = 90, border_color = NA)
#dev.off()

##### B2M/MHCI UT analysis #####

#subset specific HLA/B2M genes
names <- c("HLA-A", "HLA-B", "HLA-C", "B2M")
hla <- subset(z_score, rownames(z_score) %in% names)
col_order <- c("X4_1", "X4_2", "X6_1", "X6_4", "X0_1", "X1_1", "X2_1", "X3_1", "X5_1")
hla_reorder <- hla[, col_order]
breaks=seq(-2,2, by=0.09)
my_color <- colorpanel(n=length(breaks)-1,low="blue",mid="grey88",high="red")
labels_col <- c("0_1", "1_1", "2_1", "3_1", "4_1", "4_2", "5_1", "6_1", "6_4")
#pdf("HLAexpression.pdf", 10, 10)
pheatmap(hla_reorder, color = my_color, show_colnames = T, show_rownames = T, 
         cluster_rows = F, cluster_cols = F, angle_col = 90, border_color = NA, scale = "row")
#dev.off()

#B2M gene, MHCI IHC correlation
x <- data.frame(cbind(ihc$MHCI, t(hla[1,])))
x$MHCI_summary <- c("low", "high", "low", "low", "high", "high", "low", "high", "high")
p <- ggscatter(x, x="V1", y="B2M",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE, 
               xlab = "%Tumor cells membrane MHCI pos", 
               ylab = "B2M Expression Level (RNA)") # Add confidence interval 
p <- p + stat_cor(method = "kendall", cor.coef.name	="tau")
ggsave(plot=p, height=6, width=6, dpi=300, filename="210203_B2Mgene_MHCIihc_correl.pdf", useDingbats=FALSE)
# boxplot stats by MHCI_summary
stat <- wilcox.test(x[which(x$MHCI_summary=='high'), 2], x[which(x$MHCI_summary=='low'), 2], alternative = "two.sided")
p <- ggplot(x, aes(x=MHCI_summary, y=B2M)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + ggtitle(paste('p=', stat$p.value))
ggsave(plot=p, height=6, width=6, dpi=300, filename="210126_B2Mgene_MHCIihc_correl_MHCIsummary.pdf", useDingbats=FALSE)
# boxplot stats by patient response (HLAs and B2M)
y <- data.frame(cbind(ihc$MHCI, t(hla)))
y$patient_response <- c("non", "non", "non", "non", "eres", "eres", "non", "eres", "eres")

stat <- wilcox.test(y[which(y$patient_response=='eres'), 2], y[which(y$patient_response=='non'), 2], alternative = "two.sided")
p <- ggplot(y, aes(x=patient_response, y=B2M)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + ggtitle(paste('p=', stat$p.value))
ggsave(plot=p, height=6, width=6, dpi=300, filename="210225_B2Mgene_patientresponse.pdf", useDingbats=FALSE)

stat <- wilcox.test(y[which(y$patient_response=='eres'), 3], y[which(y$patient_response=='non'), 3], alternative = "two.sided")
p <- ggplot(y, aes(x=patient_response, y=HLA.A)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + ggtitle(paste('p=', stat$p.value))
ggsave(plot=p, height=6, width=6, dpi=300, filename="210225_HLAAgene_patientresponse.pdf", useDingbats=FALSE)

stat <- wilcox.test(y[which(y$patient_response=='eres'), 4], y[which(y$patient_response=='non'), 4], alternative = "two.sided")
p <- ggplot(y, aes(x=patient_response, y=HLA.B)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + ggtitle(paste('p=', stat$p.value))
ggsave(plot=p, height=6, width=6, dpi=300, filename="210225_HLABgene_patientresponse.pdf", useDingbats=FALSE)

stat <- wilcox.test(y[which(y$patient_response=='eres'), 5], y[which(y$patient_response=='non'), 5], alternative = "two.sided")
p <- ggplot(y, aes(x=patient_response, y=HLA.C)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1) + ggtitle(paste('p=', stat$p.value))
ggsave(plot=p, height=6, width=6, dpi=300, filename="210225_HLACgene_patientresponse.pdf", useDingbats=FALSE)


