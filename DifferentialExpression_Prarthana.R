BiocManager::install("tximportData")
library("tximportData")
BiocManager::install("tximport", force = TRUE)
library("tximport")
install.packages("readr")
library("readr")
BiocManager::install("DESeq2", force = TRUE)
library("DESeq2")
install.packages("dplyr")
library("dplyr")
BiocManager::install("GenomicFeatures")
library("GenomicFeatures")
install.packages("ggplot2")
library("ggplot2")
BiocManager::install("apeglm")
library("apeglm")
library("readr")
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
install.packages("magrittr")
install.packages("readr")
library(magrittr)
library(readr)
install.packages("tximport")
library(tximport)
install.packages("DESeq2")
library(DESeq2)





#creates transcript database, gets transcript Ids with corresponding geneIDs, gets rid of rows with N/A
txdb <- makeTxDbFromGFF("/storage/ice-shared/biol6150/Data/DiffrentialExpression/Reference/GRCh38_latest_genomic.gff.gz")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene = tx2gene %>% filter(!is.na(TXNAME))
head(tx2gene)
write.table(tx2gene, "~/biol6150/ProjectSubmissions/Group24/Project56/tx2gene.RefSeq.All.tsv", 
            sep = ",", 
            row.names = FALSE)

# Read tx2gene if already created before.
tx2gene <- read_csv("~/biol6150/ProjectSubmissions/Group24/Project56/tx2gene.RefSeq.All.tsv") %>% as.data.frame()
head(tx2gene)

# Read the samples map file.
samples <- c ("SRR24360141", "SRR24360142", "SRR24360143", "SRR24360144", "SRR24360145", "SRR24360146", "SRR24360147", "SRR24360148", "SRR24360149", "SRR24360150", "SRR24360151", "SRR24360152", "SRR24360153", "SRR24360154", "SRR24360155", "SRR24360156", "SRR24360157", "SRR24360158", "SRR24360159", "SRR24360160", "SRR24360161", "SRR24360162", "SRR24360163", "SRR24360164" )
print(samples)

file_paths = paste0("~/biol6150/ProjectSubmissions/Group24/Project56/rnaseqanalysis/Quantificationfin/", samples, "_quant/quant.sf")
print(file_paths)

# Create a txi object. Read Salmon SF files.
txi.salmon <- tximport(file_paths, type = "salmon", txOut = TRUE, tx2gene = tx2gene)
print(head(txi.salmon$counts))

# Deseq needs specific format for the map file.
sampleTable <- data.frame(condition = factor(c(rep("Treated",12), rep("Untreated",12))
))
rownames(sampleTable) <- colnames(txi.salmon$counts)
sampleTable

ncol(counts)
nrow(colData)

#Create the dds object. This prepares the data for DESeq2 to run diffrential expression.
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)

## Run the DESeq pipeline
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
head(res)
write.csv(res, "~/biol6150/ProjectSubmissions/Group24/Project56/res.csv", row.names = FALSE)

# Check how many transcripts pass the filter by adjusted p value.
table(res$padj<0.05)

# Order by adjusted p-value
res <- res[order(res$padj), ]
head(res)

#Visualize the results.
res = res %>% as.data.frame()

ggplot(res, aes(x = log2FoldChange,
                y = -log10(padj))) + 
  geom_point() +
  theme_bw(base_size = 24)

#Get significant results.
res_sig = res %>% filter(padj < 0.05) 
dim(res_sig)
head(res_sig)
write.csv(res_sig, "~/biol6150/ProjectSubmissions/Group24/Project56/res_sig.csv", row.names = FALSE)


#Volcano plot for just the significant results.
ggplot(res_sig, aes(x = log2FoldChange,
                    y = -log10(padj))) + 
  geom_point(pch = 21, fill = "magenta") +
  ylim(c(0,100)) +
  xlim(c(-15,15)) +
  theme_bw(base_size = 24)



