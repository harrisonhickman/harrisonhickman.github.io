if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("rhdf5")
BiocManager::install("DESeq2")
install.packages("pheatmap")
BiocManager::install("EnhancedVolcano")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler", force = TRUE)


library(EnhancedVolcano)
library(pheatmap)
library(tximport)
library(rhdf5)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)


files <- c(
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 1/Kallisto_results/SRR4423266/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 2/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 3/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 4/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 5/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 6/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 7/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 8/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 9/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 10/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 11/kallisto_output/abundance.tsv",
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/sample 12/kallisto_output/abundance.tsv"
)


length(files)
all(file.exists(files))

names(files) <- c(
  "Melanocytes_1","Melanocytes_2","Melanocytes_3",
  "SK-Mel-28_1","SK-Mel-28_2","SK-Mel-28_3",
  "SK-Mel-147_1","SK-Mel-147_2","SK-Mel-147_3",
  "UACC-62_1","UACC-62_2","UACC-62_3"
)

txi <- tximport(files,
                type = "kallisto",
                txOut = TRUE)

#Now let's verify what you've got:
dim(txi$counts)

#Then look at the first few rows:
head(txi$counts)

#create metadata
sample_info <- data.frame(
  row.names = c(
    "Melanocytes_1","Melanocytes_2","Melanocytes_3",
    "SK-Mel-28_1","SK-Mel-28_2","SK-Mel-28_3",
    "SK-Mel-147_1","SK-Mel-147_2","SK-Mel-147_3",
    "UACC-62_1","UACC-62_2","UACC-62_3"
  ),
  condition = c(
    "Melanocytes","Melanocytes","Melanocytes",
    "SK-Mel-28","SK-Mel-28","SK-Mel-28",
    "SK-Mel-147","SK-Mel-147","SK-Mel-147",
    "UACC-62","UACC-62","UACC-62"
  )
)

sample_info

#verify
colnames(txi$counts)
rownames(sample_info)

#create DESeq2 object
?DESeqDataSetFromTximport

dds <- DESeqDataSetFromTximport(
  txi,
  colData = sample_info,
  design = ~ condition
)

dds <- DESeq(dds)

resultsNames(dds)

#Extract results for each comparison
#Results are currently at the transcript level (RefSeq transcript IDs such as NM_000090.3 and XM_006726136.1). Additional annotation may be required to convert transcript IDs to gene symbols for biological interpretation.
res_SK28 <- results(
  dds,
  contrast = c("condition", "SK-Mel-28", "Melanocytes")
)

res_SK147 <- results(
  dds,
  contrast = c("condition", "SK-Mel-147", "Melanocytes")
)

res_UACC62 <- results(
  dds,
  contrast = c("condition", "UACC-62", "Melanocytes")
)

#inspection
head(res_SK28)
head(res_SK147)
head(res_UACC62)

#find DEGs
sig_SK28 <- subset(res_SK28, padj < 0.05)
nrow(sig_SK28)

sig_SK147 <- subset(res_SK147, padj < 0.05)
nrow(sig_SK147)

sig_UACC62 <- subset(res_UACC62, padj < 0.05)
nrow(sig_UACC62)

cat("SK-Mel-28:", nrow(sig_SK28), "\n")
cat("SK-Mel-147:", nrow(sig_SK147), "\n")
cat("UACC-62:", nrow(sig_UACC62), "\n")

#PCA
vsd <- vst(dds)

plotPCA(vsd, intgroup = "condition")

#heatmap
vsd <- vst(dds, blind = FALSE)

topVarGenes <- head(
  order(rowVars(assay(vsd)), decreasing = TRUE),
  50
)

pheatmap(
  assay(vsd)[topVarGenes, ],
  scale = "row",
  annotation_col = sample_info,
  show_rownames = FALSE
)
#make it into pdf
pdf("Heatmap.pdf", width = 8, height = 6)

pheatmap(
  assay(vsd)[topVarGenes, ],
  scale = "row",
  annotation_col = sample_info,
  show_rownames = FALSE
)

dev.off()

#locate
list.files(pattern = "\\.pdf$")

normalizePath("Heatmap.pdf")

shell.exec(getwd())


#volcano plot
#check package is loaded and correct pdf changes
?EnhancedVolcano

graphics.off()
dev.list()
#for SK-Mel-28
EnhancedVolcano(
  res_SK28,
  lab = rownames(res_SK28),
  x = "log2FoldChange",
  y = "padj",
  title = "SK-Mel-28 vs Melanocytes",
  pCutoff = 0.05,
  FCcutoff = 1
)

#for SK-mel-147
EnhancedVolcano(
  res_SK147,
  lab = rownames(res_SK147),
  x = "log2FoldChange",
  y = "padj",
  title = "SK-Mel-147 vs Melanocytes",
  pCutoff = 0.05,
  FCcutoff = 1
)

#for UACC-62
EnhancedVolcano(
  res_UACC62,
  lab = rownames(res_UACC62),
  x = "log2FoldChange",
  y = "padj",
  title = "UACC-62 vs Melanocytes",
  pCutoff = 0.05,
  FCcutoff = 1
)

#Significant transcripts
sig_SK28 <- subset(
  res_SK28,
  !is.na(padj) & padj < 0.05
)

#Annotate all significant transcripts
transcripts <- sub("\\..*$", "", rownames(sig_SK28))

gene_map <- select(
  org.Hs.eg.db,
  keys = transcripts,
  keytype = "REFSEQ",
  columns = c("SYMBOL", "GENENAME")
)

#Convert DESeq2 object to dataframe
sig_SK28_df <- as.data.frame(sig_SK28)

#Add transcript IDs as a column
sig_SK28_df$REFSEQ <- sub("\\..*$", "", rownames(sig_SK28_df))

#Merge annotation with DE results
sig_SK28_annotated <- merge(
  sig_SK28_df,
  gene_map,
  by = "REFSEQ"
)

#top upregulated genes in SK28
top_up <- sig_SK28_annotated[
  order(sig_SK28_annotated$log2FoldChange,
        decreasing = TRUE),
]
head(top_up, 20)
#make this a table
top_up <- top_up[, c(
  "SYMBOL",
  "GENENAME",
  "log2FoldChange",
  "padj"
)]

head(top_up, 20)

#for KEGG enrichment

dir.exists("C:/Users/harri/AppData/Local/clusterProfiler")


top20_up <- head(top_up, 20)
View(top20_up)

gene_df <- bitr(
  top20_up$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#GO 
up_genes <- sig_SK28_annotated$SYMBOL[
  sig_SK28_annotated$log2FoldChange > 0
]

up_genes <- na.omit(up_genes)

gene_df <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

#lets rerun GO
go_up <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP"
)

head(as.data.frame(go_up))
dotplot(go_up, showCategory = 10)
#it worked yay


#top downregulated genes in SK28
top_down <- sig_SK28_annotated[
  order(sig_SK28_annotated$log2FoldChange),
]
head(top_down, 20)
#turn into table
top_down <- top_down[, c(
  "SYMBOL",
  "GENENAME",
  "log2FoldChange",
  "padj"
)]

head(top_down, 20)

#GO enrichment
down_genes <- sig_SK28_annotated$SYMBOL[
  sig_SK28_annotated$log2FoldChange < 0
]

down_genes <- na.omit(down_genes)


gene_df_down <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_down <- enrichGO(
  gene = gene_df_down$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP"
)

#view
head(as.data.frame(go_down))

#dotplot
dotplot(go_down, showCategory = 10)


#this sections separate
top_up_SK28_clean <- subset(
  top_up,
  !is.na(SYMBOL)
)

top_down_SK28_clean <- subset(
  top_down,
  !is.na(SYMBOL)
)

#SK-Mel-147 lets repeat the process

sig_SK147 <- subset(
  res_SK147,
  !is.na(padj) & padj < 0.05
)

transcripts <- sub("\\..*$", "", rownames(sig_SK147))

gene_map <- select(
  org.Hs.eg.db,
  keys = transcripts,
  keytype = "REFSEQ",
  columns = c("SYMBOL", "GENENAME")
)

sig_SK147_df <- as.data.frame(sig_SK147)

sig_SK147_df$REFSEQ <- sub("\\..*$", "", rownames(sig_SK147_df))

sig_SK147_annotated <- merge(
  sig_SK147_df,
  gene_map,
  by = "REFSEQ"
)

#GO upreg
up_genes <- sig_SK147_annotated$SYMBOL[
  sig_SK147_annotated$log2FoldChange > 0
]

up_genes <- na.omit(up_genes)

gene_df <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_up <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP"
)

head(as.data.frame(go_up))
dotplot(go_up, showCategory = 10)

#GO downreg
down_genes <- sig_SK147_annotated$SYMBOL[
  sig_SK147_annotated$log2FoldChange < 0
]

down_genes <- na.omit(down_genes)

gene_df <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_down <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP"
)
head(as.data.frame(go_down))
dotplot(go_down, showCategory = 10)




top_up_SK147 <- sig_SK147_annotated[
  order(sig_SK147_annotated$log2FoldChange,
        decreasing = TRUE),
]

top_down_SK147 <- sig_SK147_annotated[
  order(sig_SK147_annotated$log2FoldChange),
]


top_up_SK147_clean <- subset(
  top_up_SK147,
  !is.na(SYMBOL)
)

top_up_SK147_clean <- top_up_SK147_clean[, c(
  "SYMBOL",
  "GENENAME",
  "log2FoldChange",
  "padj"
)]

top_down_SK147_clean <- subset(
  top_down_SK147,
  !is.na(SYMBOL)
)

top_down_SK147_clean <- top_down_SK147_clean[, c(
  "SYMBOL",
  "GENENAME",
  "log2FoldChange",
  "padj"
)]


#UAC-62 lets repeat the process

sig_UACC62 <- subset(
  res_UACC62,
  !is.na(padj) & padj < 0.05
)

transcripts <- sub("\\..*$", "", rownames(sig_UACC62))

gene_map <- select(
  org.Hs.eg.db,
  keys = transcripts,
  keytype = "REFSEQ",
  columns = c("SYMBOL", "GENENAME")
)

sig_UACC62_df <- as.data.frame(sig_UACC62)

sig_UACC62_df$REFSEQ <- sub("\\..*$", "", rownames(sig_UACC62_df))

sig_UACC62_annotated <- merge(
  sig_UACC62_df,
  gene_map,
  by = "REFSEQ"
)

#GO UPREG
up_genes <- sig_UACC62_annotated$SYMBOL[
  sig_UACC62_annotated$log2FoldChange > 0
]

up_genes <- na.omit(up_genes)

gene_df <- bitr(
  up_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_up <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP"
)

head(as.data.frame(go_up))
dotplot(go_up, showCategory = 10)

#GO DOWNREG
down_genes <- sig_UACC62_annotated$SYMBOL[
  sig_UACC62_annotated$log2FoldChange < 0
]

down_genes <- na.omit(down_genes)

gene_df <- bitr(
  down_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

go_down <- enrichGO(
  gene = gene_df$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP"
)

head(as.data.frame(go_down))
dotplot(go_down, showCategory = 10)









top_up_UACC62 <- sig_UACC62_annotated[
  order(sig_UACC62_annotated$log2FoldChange,
        decreasing = TRUE),
]

top_down_UACC62 <- sig_UACC62_annotated[
  order(sig_UACC62_annotated$log2FoldChange),
]

top_up_UACC62_clean <- subset(
  top_up_UACC62,
  !is.na(SYMBOL)
)

top_up_UACC62_clean <- top_up_UACC62_clean[, c(
  "SYMBOL",
  "GENENAME",
  "log2FoldChange",
  "padj"
)]

top_down_UACC62_clean <- subset(
  top_down_UACC62,
  !is.na(SYMBOL)
)

top_down_UACC62_clean <- top_down_UACC62_clean[, c(
  "SYMBOL",
  "GENENAME",
  "log2FoldChange",
  "padj"
)]


#which genes are consistently altered across all three melanoma cell lines
common_up <- Reduce(
  intersect,
  list(
    top_up$SYMBOL,
    top_up_SK147$SYMBOL,
    top_up_UACC62$SYMBOL
  )
)

common_down <- Reduce(
  intersect,
  list(
    top_down$SYMBOL,
    top_down_SK147$SYMBOL,
    top_down_UACC62$SYMBOL
  )
)

#debugging/fixing
common_up_top20 <- Reduce(
  intersect,
  list(
    head(top_up_SK28_clean$SYMBOL, 20),
    head(top_up_SK147_clean$SYMBOL, 20),
    head(top_up_UACC62_clean$SYMBOL, 20)
  )
)

#convert to dataframe view top20 commonly shared gened
common_up_df <- data.frame(
  SYMBOL = common_up_top20
)

View(common_up_df)


#less stringent
common_up <- Reduce(
  intersect,
  list(
    top_up_SK28_clean$SYMBOL,
    top_up_SK147_clean$SYMBOL,
    top_up_UACC62_clean$SYMBOL
  )
)

View(data.frame(SYMBOL = common_up))
length(common_up)
head(common_up, 20)









#lets save all the tables as csvs
write.csv(
  head(top_up, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Upregulated_SKMel28.csv",
  row.names = FALSE
)
write.csv(
  head(top_down, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Downregulated_SKMel28.csv",
  row.names = FALSE
)


write.csv(
  head(top_up_clean, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Upregulated_SK147.csv",
  row.names = FALSE
)
write.csv(
  head(top_down_SK147, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Downregulated_SK147.csv",
  row.names = FALSE
)

write.csv(
  head(top_up_UACC62, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Upregulated_UACC62.csv",
  row.names = FALSE
)
write.csv(
  head(top_down_UACC62, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Downregulated_UACC62.csv",
  row.names = FALSE
)

#SK28 and SK147 cleanup
write.csv(
  head(top_up_SK147_clean, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Upregulated_SK147.csv",
  row.names = FALSE
)

write.csv(
  head(top_down_SK147_clean, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Downregulated_SK147.csv",
  row.names = FALSE
)

write.csv(
  head(top_up_SK28_clean, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Upregulated_SK28.csv",
  row.names = FALSE
)

write.csv(
  head(top_down_SK28_clean, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Downregulated_SK28.csv",
  row.names = FALSE
)


#UACC-62 cleanup
write.csv(
  head(top_up_UACC62_clean, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Upregulated_UACC62.csv",
  row.names = FALSE
)

write.csv(
  head(top_down_UACC62_clean, 20),
  "C:/Users/harri/OneDrive/Documents/A Masters/Side project/Melanoma rna-seq project/Results/Top20_Downregulated_UACC62.csv",
  row.names = FALSE
)

#debugging
head(res_SK28$log2FoldChange)
head(res_SK147$log2FoldChange)

head(top_up_SK28_clean$SYMBOL, 20)

head(top_up_SK147_clean$SYMBOL, 20)

identical(top_up_SK28_clean, top_up_SK147_clean)
 
