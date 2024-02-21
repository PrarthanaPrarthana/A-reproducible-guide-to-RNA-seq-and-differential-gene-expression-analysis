# A-reproducible-guide-to-RNA-seq-and-differential-gene-expression-analysis

## Objective
The objective of this guide is to provide a detailed walkthrough of my experience conducting RNA-seq data analysis and differential gene expression analysis in a Jupyter Notebook environment. My goal was to identify up- and down-regulated transcripts between cases and controls, as defined by a randomized controlled trial.

## Requirements
- **RNA-seq Tools:** RNA-seq analysis can be performed using Salmon, Kallisto, or Star.
- **Differential Gene Expression Tools:** Differential gene expression analysis can be performed using DESeq2, Sleuth, or EdgeR. DESeq2 is available in both Python and R.

## Procedure

### 1. Data Acquisition
I obtained RNA-seq data from a randomized controlled trial comprising cases and controls. The dataset consisted of 30 samples, with 15 cases and 15 controls.

### 2. Data Preprocessing
#### Quality Control
To start, I assessed the quality of the raw sequencing reads using FastQC. Upon analysis, I noticed that some reads had low quality scores at the ends. Therefore, I performed quality trimming using Trimmomatic to improve the overall quality of the data.

#### Read Alignment or Quantification
For read quantification, I decided to use Salmon. After preprocessing the data and removing low-quality reads, I quantified the remaining reads against the reference transcriptome using Salmon.

### 3. Differential Gene Expression Analysis

#### Differential Expression Analysis
To identify differentially expressed genes between cases and controls, I utilized DESeq2. After setting up the analysis, I performed differential expression analysis, comparing the expression levels of genes in cases versus controls. Significant differentially expressed genes were determined based on adjusted p-values (< 0.05) and fold change thresholds (|log2FC| > 1).

### 4. Interpretation and Visualization
To better understand the results, I created various plots, including volcano plots, heatmaps, and gene expression profiles. These visualizations helped me interpret the findings and identify the biological relevance of the differentially expressed genes.


## Conclusion

This guide provides a comprehensive approach to conducting such analyses in a Jupyter Notebook environment. By following the steps outlined in this guide, others should be able to analyze their own RNA-seq datasets and identify differentially expressed genes between cases and controls in a reproducible manner.
