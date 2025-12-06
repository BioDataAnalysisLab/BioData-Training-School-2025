# # Install CRAN packages
# install.packages(c(
#   "tidyverse",
#   "ggpubr",
#   "pheatmap",
#   "data.table"
# ))
# 
# # Install Bioconductor Packages
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c(
#   "Biobase",
#   "GEOquery",
#   "DESeq2",
#   "org.Hs.eg.db",
#   "ashr",
#   "clusterProfiler",
#   "msigdbr",
#   "fgsea"
# ))

# Call them
# Tidyverse is the MOA Tidyverse Ecosystem
library(tidyverse)
# ggpubr is for cool plots
library(ggpubr)
# to make publication-ready heatmaps
library(pheatmap)
# for reading big and compressed files directly
library(data.table)
# Biobabse and GEOquery are necessary for downloading datasets from GEO
library(Biobase)
library(GEOquery)
# DESeq2 if sor differential expression
library(DESeq2)
# You can annotate your data easily with org.Hs.eg.db
library(org.Hs.eg.db)
# ashr is for log-foldshrinkage
library(ashr)
# The rest are for biological pathway enrichment
library(clusterProfiler)
library(msigdbr)
library(fgsea)

# Set some working directory of your choice
getwd()
# Then Create two directories
## One to save csv tables
out_dir   <- file.path(getwd(), "DE_outputs")
## Second, to save the figures
plot_dir  <- file.path(getwd(), "DE_plots")

# R handles everything by itself
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# We will download this dataset from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31084
# Get raw data all at once
# options(timeout=100) if internet is slow :)
filePaths <- getGEOSuppFiles("GSE225537")

################################################################################
## Part 1
# Load Raw read count matrix in R: GSE225537 is created by defult
raw_counts <- fread("GSE225537/GSE225537_Raw_data_matrix.txt.gz")
head(raw_counts)

# Remove the decimal after ensemble ID: simple regex
raw_counts[, ID := sub("\\..*$", "", ID)]

# A bit of a headache here!
# Set Ensemble as rownames to have a clean numerical matrix: You have duplicates
rownames(raw_counts) <- raw_counts$ID

# You can check for fun, its a recurring theme with STAR alignment
dups <- raw_counts[duplicated(raw_counts$ID) | duplicated(raw_counts$ID, fromLast = TRUE), ]
dups

# Lets clean it up
raw_counts <- raw_counts %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  group_by(ID) %>%
  slice_max(total, n = 1, with_ties = FALSE) %>%
  ungroup() %>% 
  select(-total) %>% 
  as.data.frame()

# Now set row names and remove the first column
rownames(raw_counts) <- raw_counts$ID
raw_counts <- raw_counts[,-1]
head(raw_counts)

# Save it for future use: next time you can automatically set rownames = 1
write.csv(raw_counts, file = paste0(out_dir, "/clean_raw_read_counts.csv"))
# Now the matrix is clean, ready for DESeq2

################################################################################
## Part 2: Metadata
# Get MetaData
gse <- getGEO("GSE225537", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])

# Create a meatdata file for differential expression
sample_info <- pheno %>% dplyr::select(title, characteristics_ch1.4)
sample_info <- sample_info %>%
  mutate(Treatment = characteristics_ch1.4 %>%
           str_remove("treatment:\\s*") %>%
           str_extract("^\\w+")) %>% 
  select(title, Treatment)

sample_info$Treatment <- as.factor(sample_info$Treatment)
rownames(sample_info) <- sample_info$title

##################################################################################
## Part 3: Differential Expression
# You should definitely read this documment and run it on your free time
# Michael is also cool guy, will answer any question on bioconductor
# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = sample_info,
                              design = ~ Treatment)

# Remove zero counts across rows, or you can check documentation to make it more stringent
# Its an iterative process 
dds <- dds[rowSums(counts(dds)) > 0, ]
dds <- DESeq(dds)

# Get normalized reads for PCA
rlog_eda <- rlog(dds, blind = TRUE)

# Plot PCA Quickly with deseq's default
plotPCA(rlog_eda, intgroup="Treatment")

# Customize it
pcaData <- plotPCA(rlog_eda, intgroup= "Treatment", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = paste0(plot_dir, "/PCA_Plot.pdf"), width = 4, height = 3)
ggplot(pcaData, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + 
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_classic()
dev.off()

# Expression results
diff_exp <- results(dds)
summary(diff_exp, alpha = 0.05)
resultsNames(dds)

# Save it for later
write.csv(diff_exp, file = paste0(out_dir, "/sorefanib_vs_control_diff_exp.csv"))
# Reverse if you want control - treatment
# diff_exp <- results(dds, contrast=c("Treatment", "Control", "Sorafenib"))

pdf(file = paste0(plot_dir, "/MA_Plot.pdf"), width = 5, height = 4)
plotMA(diff_exp)
dev.off()

# Too many low counts, lets shrink it
res <- lfcShrink(dds, coef="Treatment_Sorafenib_vs_Control", type="apeglm")
plotMA(res, main = "ashr shrunk")

# You can save if you want
write.csv(res, file = paste0(out_dir, "/sorefanib_vs_control_diff_exp.csv"))

# Lets Annotate Ensemble IDs into Symbols
ensemble_ids = mapIds(org.Hs.eg.db,
                      keys= rownames(diff_exp),
                      column="SYMBOL",
                      keytype="ENSEMBL",
                      multiVals="first")

diff_exp$Symbol <- ensemble_ids

# You wanna make it cool
ma_plot <- ggmaplot(diff_exp, main = "Sorefanib vs Control" ,
                    fdr = 0.05, fc = 1, size = 0.4,
                    palette = c("#B31B21", "#1465AC", "darkgray"),
                    genenames = as.vector(diff_exp$Symbol),
                    label.rectangle = TRUE,
                    ylim = c(-5,5),
                    legend = "top", 
                    top = 20,
                    font.label = c("bold", 5),
                    font.legend = "bold",
                    font.main = "bold",
                    ggtheme = ggplot2::theme_minimal())

pdf(file = paste0(plot_dir, "/MAPlot.pdf"), width = 7, height = 5)
print(ma_plot)
dev.off()

# You can do the same for shrunk log-fold
res$Symbol <- ensemble_ids
pdf(file = paste0(plot_dir, "/ASHR_MAPlot.pdf"), width = 7, height = 5)
ggmaplot(res, main = "Sorefanib vs Control" ,
         fdr = 0.05, fc = 1, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res$Symbol),
         label.rectangle = TRUE,
         ylim = c(-5,5),
         legend = "top",
         top = 20,
         font.label = c("bold", 5),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())
dev.off()

# Let's build volcano plots
# Here use your imagination, check EnhancedVolcano from bioconductor, color, customize ...
# Increase dot size, color by biological pathways, annotate genes of interest ...

# Convert DESeq results into a data frame
soraf_cont <- dplyr::mutate(
  data.frame(diff_exp),
  Sig = padj < 0.05 & abs(log2FoldChange) > 1) %>%
  na.omit(Sig)

pdf(file = paste0(plot_dir, "/Volcano_Plot.pdf"), width = 5, height = 6)
ggplot(soraf_cont, aes(x = log2FoldChange, y = -log10(padj), color = Sig)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano: Treatment vs Control",
       x = "log2 fold-change", 
       y = "-log10(FDR)") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  # Comment them out and see what happens
  ylim(c(0, 100)) +
  xlim(c(-5, 5))
dev.off()

# Do the same for log-fold shrunk
soraf_cont <- dplyr::mutate(
  data.frame(res),
  Sig = padj < 0.05 & abs(log2FoldChange) > 1) %>%
  na.omit(Sig)

pdf(file = paste0(plot_dir, "/ASHR_Volcano.pdf"), width = 5, height = 6)
ggplot(soraf_cont, aes(x = log2FoldChange, y = -log10(padj), color = Sig)) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(title = "Volcano: Treatment vs Control",
       x = "log2 fold-change", 
       y = "-log10(FDR)") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(c(0, 100)) +
  xlim(c(-5, 5))
dev.off()

# Create a heatmap
# blind = FALSE is important for downstream analysis
rlog_downstream <- assay(rlog(dds, blind = FALSE))

# Define the cut-ff
FDR_cutoff <- 0.01
lfc_cutoff <- 1
de_select <- diff_exp[diff_exp$padj < FDR_cutoff & !is.na(diff_exp$padj) & abs(diff_exp$log2FoldChange) > lfc_cutoff,]

# Select normalized genes according to the cutoff
subset_rnaseq_data <- subset(rlog_downstream, rownames(rlog_downstream) %in% rownames(de_select))

# Column annotations for better visualization
annotation_col <- data.frame(
  Condition = ifelse(grepl("SF", colnames(subset_rnaseq_data)), "Treat", "Cont"),
  row.names = colnames(subset_rnaseq_data))

# Give them colors you like
cond_colors <- c(
  Treat = "#D95F02",
  Cont  = "#1B9E77")

ann_colors <- list(
  Condition = cond_colors)

breaks <- seq(-1, 1, length.out = 100)
pheatmap(subset_rnaseq_data, 
         scale="row", 
         color = colorRampPalette(c("#1465AC", "white", "#B31B21"))(100),
         breaks = breaks,
         # These parameters can change
         clustering_distance_rows="correlation", 
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames =F, 
         show_colnames =F, 
         cluster_cols=T,
         cluster_rows = T,
         treeheight_row = 20,
         treeheight_col = 20,
         main = "Significant Genes Heatmap",
         filename = paste0(plot_dir, "/Significant_Genes_heatmap.pdf"),
         width = 4,
         height = 6)

# I will select 50 most significant genes: You can select more
# Once you get the GO terms down, you can come back here and select genes of your interest
top_genes <- diff_exp %>% as.data.frame() %>% 
  dplyr::arrange(padj) %>% slice_head(n=50)
mat <- rlog_downstream[rownames(top_genes),]

# Make sure rownames are in order: They generally are
identical(rownames(top_genes), rownames(mat))
rownames(mat) <- top_genes$Symbol

# I can use the same annotation defined above
pheatmap(mat, 
         scale="row", 
         color = colorRampPalette(c("#1465AC", "white", "#B31B21"))(100),
         breaks = breaks,
         clustering_distance_rows="correlation", 
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames =T, 
         show_colnames =F, 
         cluster_cols=T, 
         fontsize_row=6,
         border_color = "grey100",
         # angle_col = 45,
         treeheight_row = 20,
         treeheight_col = 20,
         main = "Top 50 Genes",
         filename = paste0(plot_dir, "/50_Significant_Genes_heatmap.pdf"),
         width = 4,
         height = 6)

##################################################################################
# Part 4: Let's do some Systems Biology Now
# Convert ENSEMBLE into Entrez IDs
entrez_id = mapIds(org.Hs.eg.db,
                   keys= rownames(res),
                   column="ENTREZID",
                   keytype="ENSEMBL",
                   multiVals="first")

# For GSEA we will use shrunk logFold to avoid noise
shrunk_lfc <- as.data.frame(res)
shrunk_lfc$Entrez <- entrez_id

# You can check here
# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/geseca-tutorial.html
library(msigdbr)
library(fgsea)

# If you are into biology, you can explore this further
msigdbr_df <- msigdbr(species = "human", collection = "H")

# Here we need a named vector...
## Remove any entrez NA values
res_entrez <- shrunk_lfc %>% dplyr::filter(Entrez != "NA")
## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$Entrez) == F), ]
## Extract the gene statistic
gsea.gene_stat <- res_entrez$log2FoldChange
## Name gene statistic with the corresponding Entrez ID
names(gsea.gene_stat) <- res_entrez$Entrez
## Sort gene statistic in decreasing order
gsea.gene_stat <- sort(gsea.gene_stat, decreasing = TRUE)

set.seed(2020)
gsea_results <- GSEA(
  geneList = gsea.gene_stat,
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    msigdbr_df,
    gs_name,
    ncbi_gene
  )
)

head(gsea_results@result)
# Remove "HALLMARK" from all names
gsea_results@result$Description <- gsub("^HALLMARK[_]?", "", gsea_results$Description)
pdf(file = paste0(plot_dir, "/GSEA_Hallmark_Ridgeplot.pdf"), width = 10, height = 8)
ridgeplot(gsea_results, showCategory = 15) +
  ggtitle("GSEA Hallmarks") +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    plot.title = element_text(size = 11, hjust = 0.5)
  )
dev.off()

# You can save if you want
write.csv(gsea_results, file = paste0(out_dir, "/gsea_hallmarks_sorefanib_positive.csv"))

###
# GSEA GO: Almost the same story, but now for gene ontology, much wider
# Create named vector
gene_list <- shrunk_lfc$log2FoldChange
names(gene_list) <- rownames(shrunk_lfc)
gene_list <- gene_list[!is.na(names(gene_list)) & names(gene_list) != ""]
gene_list <- sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             verbose = FALSE, 
             eps = 0,
             OrgDb = org.Hs.eg.db,
             pAdjustMethod = "BH")

ridgeplot(gse, showCategory = 20)

pdf(file = paste0(plot_dir, "/GSEA_GO.pdf"), width = 10, height = 10)
dotplot(gse, showCategory=15, split=".sign") + 
  facet_grid(.~.sign, 
             labeller = labeller(.sign = c("activated" = "Sorefanib", 
                                           "suppressed" = "Control")))
dev.off()
write.csv(gse, file = paste0(my_dir, "/GSEA_go.csv"))

pdf(file = paste0(plot_dir, "/GSEA_GO_Ridgeplot.pdf"), width = 10, height = 12)
ridgeplot(gse, showCategory = 20)
dev.off()

write.csv(gse, file = paste0(out_dir, "/gsea_go_sorefanib_positive.csv"))

################################################################################
# Part 5: More classical system biology stuff
# We can also select by cut-off (some people dont like this anymore)
signif_upregulated <- diff_exp %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange > 1)
head(signif_upregulated)

signif_downregulated <- diff_exp %>% 
  as.data.frame() %>% 
  dplyr::filter(padj < 0.01 & log2FoldChange < -1)
head(signif_downregulated)

data_universe <- diff_exp %>% as.data.frame()

sorefanib_go <- clusterProfiler::enrichGO(gene = rownames(signif_upregulated), 
                                         universe = rownames(data_universe),
                                         keyType = "ENSEMBL",
                                         OrgDb = org.Hs.eg.db,
                                         ont = "BP", 
                                         pAdjustMethod = "BH",
                                         readable = TRUE)

dotplot(sorefanib_go, showCategory=20, font.size=9)

sorefanib_go <- clusterProfiler::simplify(sorefanib_go, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(sorefanib_go, showCategory=20, font.size=9)

pdf(file = paste0(plot_dir, "/Sorefanib_GO_Upregulated.pdf"), width = 6, height = 7)
dotplot(sorefanib_go, showCategory=20, font.size=9)
dev.off()
write.csv(sorefanib_go, file = paste0(out_dir, "/sorefanib_upregulated_goterms.csv"))

control_go <- clusterProfiler::enrichGO(gene = rownames(signif_downregulated), 
                                          universe = rownames(data_universe),
                                          keyType = "ENSEMBL",
                                          OrgDb = org.Hs.eg.db,
                                          ont = "BP", 
                                          pAdjustMethod = "BH",
                                          readable = TRUE)

dotplot(control_go, showCategory=20, font.size=9)

control_go <- clusterProfiler::simplify(control_go, cutoff=0.7, by="p.adjust", select_fun=min)
dotplot(control_go, showCategory=20, font.size=9)

pdf(file = paste0(plot_dir, "/Control_GO_Upregulated.pdf"), width = 7, height = 3)
dotplot(control_go, showCategory=20, font.size=9)
dev.off()
write.csv(control_go, file = paste0(out_dir, "/Control_upregulated_goterms.csv"))

# Save whatever you want for future
save(dds, diff_exp, gsea_results, gse, file = "sorefanib_analysis_data.RData")
