RNA-seq Differential Expression Analysis
================

``` r
library(tximport)
library(readr)
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

``` r
library(GenomicFeatures)
```

    ## Loading required package: AnnotationDbi

``` r
library(AnnotationDbi)
```

Upload RNA sample vs Time Condition

``` r
setwd("/Users/hariniramanan/Desktop/DOP/RNAseq dataset")
samples <- read.csv("samples.csv", header=TRUE) 
```

Upload Quant files from Salmon Pseudo alignment

``` r
files <- file.path("quants", samples$sample, "quant.sf")
names(files) <- samples$sample
```

Build a TxDb object from GTF – getting transcript data from public

``` r
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.110.gtf.gz")
```

    ## Warning in call_fun_in_txdbmaker("makeTxDbFromGFF", ...): makeTxDbFromGFF() has moved from GenomicFeatures to the txdbmaker package,
    ##   and is formally deprecated in GenomicFeatures >= 1.59.1. Please call
    ##   txdbmaker::makeTxDbFromGFF() to get rid of this warning.

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ...

    ## Warning in .get_cds_IDX(mcols0$type, mcols0$phase): The "phase" metadata column contains non-NA values for features of type
    ##   stop_codon. This information was ignored.

    ## Warning in .makeTxDb_normarg_chrominfo(chrominfo): genome version information
    ## is not available for this TxDb object

    ## OK

Extract transcript-to-gene mapping

``` r
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
write.csv(tx2gene, file="tx2gene.csv", row.names=FALSE)
tx2gene <- read.csv("tx2gene.csv")
```

Now txi\$counts will be gene-wise counts.

``` r
txdb <- GenomicFeatures::makeTxDbFromGFF("Homo_sapiens.GRCh38.110.gtf.gz")
```

    ## Warning in call_fun_in_txdbmaker("makeTxDbFromGFF", ...): makeTxDbFromGFF() has moved from GenomicFeatures to the txdbmaker package,
    ##   and is formally deprecated in GenomicFeatures >= 1.59.1. Please call
    ##   txdbmaker::makeTxDbFromGFF() to get rid of this warning.

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ...

    ## Warning in .get_cds_IDX(mcols0$type, mcols0$phase): The "phase" metadata column contains non-NA values for features of type
    ##   stop_codon. This information was ignored.

    ## Warning in .makeTxDb_normarg_chrominfo(chrominfo): genome version information
    ## is not available for this TxDb object

    ## OK

``` r
library(txdbmaker)
```

    ## 
    ## Attaching package: 'txdbmaker'
    ## 
    ## The following objects are masked from 'package:GenomicFeatures':
    ## 
    ##     browseUCSCtrack, getChromInfoFromBiomart, makeFDbPackageFromUCSC,
    ##     makeFeatureDbFromUCSC, makePackageName, makeTxDb,
    ##     makeTxDbFromBiomart, makeTxDbFromEnsembl, makeTxDbFromGFF,
    ##     makeTxDbFromGRanges, makeTxDbFromUCSC, makeTxDbPackage,
    ##     makeTxDbPackageFromBiomart, makeTxDbPackageFromUCSC,
    ##     supportedMiRBaseBuildValues, supportedUCSCFeatureDbTables,
    ##     supportedUCSCFeatureDbTracks, supportedUCSCtables,
    ##     UCSCFeatureDbTableSchema

``` r
txdb <- txdbmaker::makeTxDbFromGFF("Homo_sapiens.GRCh38.110.gtf.gz")
```

    ## Import genomic features from the file as a GRanges object ... OK
    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ...

    ## Warning in .get_cds_IDX(mcols0$type, mcols0$phase): The "phase" metadata column contains non-NA values for features of type
    ##   stop_codon. This information was ignored.
    ## Warning in .get_cds_IDX(mcols0$type, mcols0$phase): genome version information is not available for this TxDb object

    ## OK

Extract mapping with annotation file - create table with transcript and
to its according gene

``` r
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb,
                                 keys = k,
                                 keytype = "TXNAME",
                                 columns = "GENEID")
```

    ## 'select()' returned 1:1 mapping between keys and columns

Clean transcript IDs (remove version numbers like .1, .2)

``` r
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)
head(tx2gene)
```

    ##            TXNAME          GENEID
    ## 1 ENST00000456328 ENSG00000290825
    ## 2 ENST00000450305 ENSG00000223972
    ## 3 ENST00000473358 ENSG00000243485
    ## 4 ENST00000469289 ENSG00000243485
    ## 5 ENST00000607096 ENSG00000284332
    ## 6 ENST00000606857 ENSG00000268020

``` r
write.csv(tx2gene, "tx2gene.csv", row.names = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
```

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 
    ## transcripts missing from tx2gene: 14313
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
samples <- data.frame(
  SampleID = c("RMC192","RMC196","RMC200",
               "RMC193","RMC197","RMC201",
               "RMC194","RMC198","RMC202",
               "RMC195","RMC199","RMC203"),
  Timepoint = c("3hr","3hr","3hr",
                "3day","3day","3day",
                "6day","6day","6day",
                "10day","10day","10day"),
  stringsAsFactors = FALSE
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ Timepoint
)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## using 'avgTxLength' from assays(dds), correcting for library size

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Gene Annotation

``` r
library(ggplot2)
library(org.Hs.eg.db)
```

    ## 

``` r
library(clusterProfiler)
```

    ## 

    ## clusterProfiler v4.16.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan,
    ## X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal
    ## enrichment tool for interpreting omics data. The Innovation. 2021,
    ## 2(3):100141

    ## 
    ## Attaching package: 'clusterProfiler'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename

    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(AnnotationDbi)
library(UniProt.ws)
```

``` r
res_3hr_vs_3day <- results(dds, contrast = c("Timepoint", "3hr", "3day"))
head(res_3hr_vs_3day)
```

    ## log2 fold change (MLE): Timepoint 3hr vs 3day 
    ## Wald test p-value: Timepoint 3hr vs 3day 
    ## DataFrame with 6 rows and 6 columns
    ##                  baseMean log2FoldChange     lfcSE       stat    pvalue
    ##                 <numeric>      <numeric> <numeric>  <numeric> <numeric>
    ## ENSG00000000003   8.41577     -1.1464825  0.769869 -1.4891926 0.1364367
    ## ENSG00000000005   0.00000             NA        NA         NA        NA
    ## ENSG00000000419 288.46365      0.0783251  0.181342  0.4319189 0.6658003
    ## ENSG00000000457  30.90924     -0.0200842  0.439435 -0.0457047 0.9635456
    ## ENSG00000000460  14.10520      1.2925515  0.661199  1.9548585 0.0505998
    ## ENSG00000000938   6.47810      0.0000000  2.082011  0.0000000 1.0000000
    ##                      padj
    ##                 <numeric>
    ## ENSG00000000003  0.597790
    ## ENSG00000000005        NA
    ## ENSG00000000419  0.960490
    ## ENSG00000000457  1.000000
    ## ENSG00000000460  0.407676
    ## ENSG00000000938  1.000000

``` r
res_3day_vs_6day <- results(dds, contrast = c("Timepoint", "3day", "6day"))
res_6day_vs_10day <- results(dds, contrast = c("Timepoint", "6day", "10day"))

# Define the contrasts you want to check
contrasts <- list(
  c("Timepoint", "3hr", "3day"),
  c("Timepoint", "3day", "6day"),
  c("Timepoint", "6day", "10day"),
  c("Timepoint", "3hr", "10day"),
  c("Timepoint", "3hr", "6day"),
  c("Timepoint", "3day", "10day")
)
```

``` r
###Data Visualization 
library(ggplot2)  # nice package for volcano plots

# Define the contrasts
contrasts <- list(
  c("Timepoint", "3hr", "3day"),
  c("Timepoint", "3day", "6day"),
  c("Timepoint", "6day", "10day"),
  c("Timepoint", "3hr", "10day"),
  c("Timepoint", "3hr", "6day"),
  c("Timepoint", "3day", "10day")
)

library(ggplot2)
```

``` r
  # 3hr vs 3day results
  res <- results(dds, contrast = c("Timepoint", "3hr", "3day"))
  res_df <- as.data.frame(res)
  
  # Add a column for significance (padj < 0.05 & |log2FC| > 1)
  res_df$Significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")
  
  # Volcano plot
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.6, size = 1.8) +
    scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot: 3hr vs 3day",
         x = "Log2 Fold Change",
         y = "-log10 adjusted p-value")
```

    ## Warning: Removed 23941 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
  ## --- MA Plot ---
# MA plot
comp_name <- "3hr vs 3day"
plotMA(res, main = paste("MA Plot:", comp_name), ylim = c(-5, 5))
```

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# PCA
vsd <- vst(dds)
plotPCA(vsd, intgroup = "Timepoint")
```

    ## using ntop=500 top features by variance

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation idioms with `aes()`.
    ## ℹ See also `vignette("ggplot2-in-packages")` for more information.
    ## ℹ The deprecated feature was likely used in the DESeq2 package.
    ##   Please report the issue to the authors.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-14-2.png)<!-- -->

``` r
library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main="Sample Distance Heatmap")
```

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
topGenes <- head(order(res_3hr_vs_3day$padj), 50)
mat <- assay(vsd)[topGenes, ]
pheatmap(mat, scale="row", show_rownames=FALSE,
         annotation_col = as.data.frame(colData(dds)["Timepoint"]))
```

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# Gene Annotation
library(org.Hs.eg.db)
# Extract Ensembl IDs from results
ens_ids <- rownames(res_3hr_vs_3day)

# Map to SYMBOLS and ENTREZ IDs
gene_map <- mapIds(org.Hs.eg.db,
                   keys = ens_ids,
                   column = c("SYMBOL"),
                   keytype = "ENSEMBL",
                   multiVals = "first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
##GO ANALYSIS
library(clusterProfiler)
library(org.Hs.eg.db)

# Example: 3hr vs 3day results
res <- results(dds, contrast = c("Timepoint", "3hr", "3day"))
res_df <- as.data.frame(res)

# Filter significant DEGs
sig_res <- subset(res_df, padj < 0.05)

# Map Ensembl → Entrez
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = rownames(sig_res),
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
# Remove NAs
entrez_ids <- na.omit(entrez_ids)


ego <- enrichGO(gene          = entrez_ids,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",      # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

# View top results
head(ego)
```

    ##                    ID                 Description GeneRatio   BgRatio
    ## GO:0016049 GO:0016049                 cell growth    19/226 492/18805
    ## GO:0001558 GO:0001558   regulation of cell growth    16/226 410/18805
    ## GO:0042692 GO:0042692 muscle cell differentiation    16/226 429/18805
    ## GO:0010506 GO:0010506     regulation of autophagy    14/226 359/18805
    ## GO:0001503 GO:0001503                ossification    16/226 458/18805
    ## GO:0007160 GO:0007160        cell-matrix adhesion    11/226 240/18805
    ##            RichFactor FoldEnrichment   zScore       pvalue   p.adjust    qvalue
    ## GO:0016049 0.03861789       3.213316 5.486740 8.333580e-06 0.02776749 0.0240358
    ## GO:0001558 0.03902439       3.247140 5.073882 3.867707e-05 0.04380116 0.0379147
    ## GO:0042692 0.03729604       3.103327 4.860469 6.617168e-05 0.04380116 0.0379147
    ## GO:0010506 0.03899721       3.244879 4.736496 1.192311e-04 0.04380116 0.0379147
    ## GO:0001503 0.03493450       2.906828 4.556486 1.411992e-04 0.04380116 0.0379147
    ## GO:0007160 0.04583333       3.813698 4.838423 1.614460e-04 0.04380116 0.0379147
    ##                                                                                                             geneID
    ## GO:0016049 2273/23327/51447/79960/51199/7040/59284/51147/1981/8835/23767/23025/10439/91/9472/55364/3020/3490/91584
    ## GO:0001558                   2273/23327/51447/79960/7040/59284/51147/1981/8835/23025/10439/91/9472/3020/3490/91584
    ## GO:0042692                        51177/87/8516/4208/7040/3280/3398/7139/1948/84466/27063/9472/70/3020/255743/1466
    ## GO:0010506                               7107/23256/8140/8992/1981/25782/529/9550/10769/9638/1915/3099/94241/79837
    ## GO:0001503                     58476/4208/7040/3199/1490/91/51176/51232/4041/3020/5364/255743/7059/2194/1435/55366
    ## GO:0007160                                                 3675/395/87/8516/10580/1490/54453/7057/255743/7059/1435
    ##            Count
    ## GO:0016049    19
    ## GO:0001558    16
    ## GO:0042692    16
    ## GO:0010506    14
    ## GO:0001503    16
    ## GO:0007160    11

``` r
# Plot
barplot(ego, showCategory=20, title="GO Biological Process")
```

    ## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • by = x
    ## ℹ Did you misspell an argument name?

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
dotplot(ego, showCategory=20, title="GO Enrichment")
```

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-18-2.png)<!-- -->

``` r
#KEGG
ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = "hsa",   # hsa = Homo sapiens
                    pvalueCutoff = 0.05)
```

    ## Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...

    ## Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...

``` r
barplot(ego, showCategory=20, title="GO Biological Process")
```

    ## Warning in fortify(object, showCategory = showCategory, by = x, ...): Arguments in `...` must be used.
    ## ✖ Problematic argument:
    ## • by = x
    ## ℹ Did you misspell an argument name?

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

``` r
##Venn Diagram
# Get DEGs for each timepoint vs 3hr
DEG_3day <- subset(results(dds, contrast=c("Timepoint","3day","3hr")), padj < 0.05)
DEG_6day <- subset(results(dds, contrast=c("Timepoint","6day","3hr")), padj < 0.05)
DEG_10day <- subset(results(dds, contrast=c("Timepoint","10day","3hr")), padj < 0.05)

# Extract gene IDs
genes_3hr  <- rownames(DEG_3day)  # DEGs appear here when 3day != 3hr
genes_3day <- rownames(DEG_3day)
genes_6day <- rownames(DEG_6day)
genes_10day <- rownames(DEG_10day)

library(VennDiagram)
```

    ## Loading required package: grid

    ## Loading required package: futile.logger

``` r
venn.plot <- venn.diagram(
  x = list(
    "3hr"  = genes_3hr,
    "3day" = genes_3day,
    "6day" = genes_6day,
    "10day"= genes_10day
  ),
  filename = NULL,
  fill = c("red", "blue", "green", "purple"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Overlap of DEGs across Timepoints"
)

grid::grid.draw(venn.plot)
```

![](RNAseq-Analysis-Notebook_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->
