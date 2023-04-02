---
title: "DEG_IFN-L_Signature_12h"
author: "Candice SAKREF"
date: "25/10/2021"
output: 
  html_document:
    toc: true
    toc_float: true
    toc_depth: 5
---
This markdown will be fully understood only taken with the previous one "V2_Markdown_Microbulk-RNA seq", and with "Verification_Stim_ISG"  


```{r}
#Definition de l'environnement de travail 
setwd("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG")

#Overview of the first elements of the df
look <- function(x){
  head(x)[1:5, 1:5]
}
```

```{r message = FALSE}
library(dplyr)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(RColorBrewer)
library(pheatmap)
library(forcats)
library(reshape)
library(viridis)
```

#I. Processing: 

Import the R object:    
```{r}
dds12_2<- readRDS("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/1_Processing/dds12_2.rds")
dds12_2

#first run DESeq2
dds12_2 <- DESeq(dds12_2)
```
**22 542 genes left!** x 35 samples => OK  

## a) Dispersion: 

DESeq2 as the *Dispertion* to study the mean and variance of each genes. 
```{r}
#Then plot
plotDispEsts(dds12_2)
```
Good fit, since the dispersions decrease with increasing mean and cluster around the maximum likelihood (ML) line.  

#II. DEG Medium vs IFNs 
## A) Lambda vs Medium
### 1. Process the data
In order to create DEg we have to pass by the function "results". However, this function only enable us to compare conditions 2 by 2...  


Creation of the **results** table.  
There is 6 columns to the table **results** 
$~~~~~~$baseMean  
$~~~~~~$log2FoldChange  
$~~~~~~$lfcSE  
$~~~~~~$stat  
$~~~~~~$pvalue  
$~~~~~~$padj  

Results tables are generated using the function results, which extracts a results table with log2 fold changes, p values and adjusted p values (padj). With no additional arguments to results, the log2 fold change and Wald test p value will be for the last variable in the design formula, and if this is a factor (we are), the comparison will be the **last level** of this variable over the **reference level*.   
This is why we will use the argument "contrast" to chose the 2 activation that we want to compare.   

Extraction of the result table : lambda vs Medium

```{r}
#For now I would like to compare "IFN-L" with "medium" as the base level 
res_unshrunken <- results(dds12_2, 
                         contrast = c("Stim", "lambda", "medium"), 
                         alpha = 0.05)
look(res_unshrunken)

plotMA(res_unshrunken, ylim=c(-8, 8), xlab="mean expression", colSig = "red3", colLine = "#ff000080", main = "DEG - Medium vs Lambda 12h")
```
We can see that we have **much more DEG** at 12h compared to 6h as below:  

```{r}
#We should technically use this formula below, but I don't know why it does not work. So we have to use another formula. 
res1_shrunken <- lfcShrink(dds12_2,
                  contrast = c("Stim", "lambda", "medium"), 
                  res = res_unshrunken, 
                  type = "normal")  

#Shrinkage of the log2 fold change
#res <- lfcShrink(dds2,
                 #coef = "Activation_lambda_vs_medium",
                 #type= "apeglm")

#plot the results
plotMA(res1_shrunken, ylim=c(-4, 4), xlab="mean expression", colSig = "red3", colLine = "#ff000080", cex= 0.5,  main = "DEG - Medium vs Lambda 12h")
```

Now keep only the gene wich are actual DEG:  

```{r}
summary(res1_shrunken)
```
**In total, 3129 DEG when comparing Lambda with medium.**  
**1695 of these DEG are upregulated.**  
**1434 of these DED are down reulated.**

Whiwh is way better compared to 6h:  

*In total, 819 DEG when comparing Lambda with medium.*   
*549 of these DEG are upregulated.*  
*270 of these DED are down reulated.*   

Gene Annotation: 
In order to understand which genesbare up/downregulated we should add the symbols and the gene's names to our table. 

```{r}
#Use the library org.Hs.eg.db: 
head(keytypes(org.Hs.eg.db))
uniKeys <- keys(org.Hs.eg.db, keytype="ENSEMBL")
cols <- c("SYMBOL", "GENENAME")
head(select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL"))
#okai it worked - I check on Ensemble website. 

annots <- select(org.Hs.eg.db, keys=rownames(res1_shrunken),
                 columns=c("SYMBOL","GENENAME"), keytype="ENSEMBL")
head(annots)
```

```{r}
#Results transformation into df:
res1_df <- res1_shrunken %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL")
head(res1_df)

# Join the res df with the gene symbols
res2 <- dplyr::inner_join(res1_df, annots, by = c("ENSEMBL"))
head(res2)
write_csv2(res2, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Lambda/GSEA_ML12.csv")
```

Recuperer les DEG only:  

```{r}
#Create a table retreving only DEG: 
ML_res2_sig <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 0.58)
nrow(ML_res2_sig)
```
**1765 DEG** found when comparing IFN-L to the Medium.   

Number of upregulated or downregulated DEG:
```{r}
ML_up <- subset(res2, padj < 0.05 & log2FoldChange > 0.58)
nrow(ML_up)
write_csv2(ML_up, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Lambda/DEG_ML_up.csv")

ML_down <- subset(res2, padj < 0.05 & log2FoldChange < (-0.58))
nrow(ML_down)
write_csv2(ML_down, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Lambda/DEG_ML_down.csv")
```

### 2. Volcano plot: 
To create a Volcano plot we first need to create a df readable for ggplot: 
```{r}
#Obtain logical vector where TRUE values denote padj value < 0.05 and log fold change > 0.5 in either direction  
res1_table_volcano <- res1_df %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 0.58)

#Reoder df with padj: 
res1_table_volcano <-  res1_table_volcano %>%
  arrange(padj)
head(res1_table_volcano)

#Add genes label: 
res1_table_volcano <- dplyr::inner_join(res1_table_volcano , annots, by = c("ENSEMBL"))

#Add a column to the df filled with "NA" to each row: 
res1_table_volcano$diffexpressed <- "NO"

#Now if the gene has log2FoldChange > 0.5 &  a pvalue < 0.05, the gene is significantly Upregulated: 
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange > 0.5 & res1_table_volcano$pvalue < 0.05] <- "UP"
#If the gene has log2FoldChange > -0.5 &  a pvalue < 0.05, the gene is significantly Downregulated: 
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange < -0.5 & res1_table_volcano$pvalue < 0.05] <- "DOWN"
```


We would like to label some genes that are DEG in the volcano plot: 
```{r}
#Add a column with "NA"
res1_table_volcano$delabel <- NA
#Add the labels only if the gene are NOT non-DEG: 
res1_table_volcano$delabel[res1_table_volcano$diffexpressed == "UP" | res1_table_volcano$diffexpressed == "DOWN"] <- res1_table_volcano$SYMBOL[res1_table_volcano$diffexpressed == "UP" | res1_table_volcano$diffexpressed == "DOWN"]

#plot the volcano plot:
ML <- ggplot(data= res1_table_volcano, aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed, label=delabel))+
  geom_point()+
  ggtitle("Medium vs Lambda - 12h")+
  scale_colour_manual(values = c("blue", "grey60","red")) +
  geom_text_repel() +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))
ML
```

```{r}
pdf("Volcano_ML.pdf", width = 10, height = 6)
ML
dev.off()
```

### 3. Creation of DEG tables 
#### a) All DEG
```{r}
#Extract the significative gene detected with the res2 from the normalized counts!
norm_counts <- read.delim("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/1_Processing/normalized_counts12.txt")
colnames(norm_counts)[colnames(norm_counts) == "Gene"] <- "ENSEMBL"
norm_counts_sig_ml <- dplyr::inner_join(ML_res2_sig, norm_counts, by = "ENSEMBL")
DEG_nc_IFNci <- norm_counts_sig_ml[, c(1,3,7,8, 16,23,30,37,44, 10,17,24,31,38, 11,18,25,32,39, 14,21,28,35,42, 15,22,29,36,43)]

norm_counts_sig_ml <- dplyr::inner_join(ML_res2_sig, norm_counts, by = "ENSEMBL")

#Create a table and merge it with annotate.
DEG_nc_IFNci <- DEG_nc_IFNci %>% 
  as.data.frame()
#Save the table in excel file: 
write_csv2(DEG_nc_IFNci, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_ML_12_nc_IFNci.csv")

#Selected only the column Medium and Lambda
norm_counts_sig_ml <- norm_counts_sig_ml[,c(1,8,16,23,30,37,44,15,22,29,36,43)]
look(norm_counts_sig_ml)

#Create a table and merge it with annotate.
norm_counts_sig_ml <- norm_counts_sig_ml %>% 
  as.data.frame()
head(norm_counts_sig_ml)

#Save the table in excel file: 
write_csv2(norm_counts_sig_ml, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_ML_nc_justML.csv")
```

#### b) All DEG - Z score
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_norm_counts_sig_ml <- norm_counts_sig_ml[,3:12]
MvsL_Z_score <- t(scale(t(Zscore_norm_counts_sig_ml)))
look(MvsL_Z_score)

MvsL_Z_score <- MvsL_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Z_score$SYMBOL <- norm_counts_sig_ml[,2]
look(MvsL_Z_score)

#Write and export the excel
write_csv2(MvsL_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsL_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  


#### c) Top 50 DEG
```{r}
#Select the top 50 DEG: 
res2_sig_Top50 <- ML_res2_sig %>%
  arrange(padj)
head(res2_sig_Top50)
res2_sig_Top50 <- res2_sig_Top50[1:50,]
res2_sig_Top50 <- res2_sig_Top50 %>%
  arrange(log2FoldChange)
look(res2_sig_Top50)

#Revtreive norm counts: 
res2_sig_Top50 <- dplyr::inner_join(res2_sig_Top50, norm_counts, by = "ENSEMBL")
look(res2_sig_Top50)

#Selected only the column Medium and Lambda
MvcL_res2_sig_Top50 <- res2_sig_Top50[,c(1,8,16,23,30,37,44,15,22,29,36,43)]
look(norm_counts_sig_ml)

#Create a table and merge it with annotate.
MvcL_res2_sig_Top50 <- MvcL_res2_sig_Top50 %>% 
  as.data.frame()
look(MvcL_res2_sig_Top50)

#Save the table: 
write_csv2(MvcL_res2_sig_Top50, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvcL_norm_counts_Top50.csv")
```

#### d) Z score - Top 50 DEG
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_Top_50_norm_counts <- MvcL_res2_sig_Top50[,3:12]
MvsL_Top_50_Z_score <- t(scale(t(Zscore_Top_50_norm_counts)))
look(MvsL_Top_50_Z_score)


MvsL_Top_50_Z_score <- MvsL_Top_50_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Top_50_Z_score$SYMBOL <- MvcL_res2_sig_Top50[,2]
look(MvsL_Top_50_Z_score)

#Write and export the excel
write_csv2(MvsL_Top_50_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsL_Top_50_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  

###4.Summary: 

**Nb DEG:**    
**Total: 1765**  
**UP: 755**  
**DOWN: 1010**  


Plot the Volcano plot Lambda vs Medium: 
```{r}
ML
```

Norm counts DEG ML IFNs condition, file called: "DEG_ML_12_nc_IFNci.csv"  
Norm counts DEG ML JUST Medium vs Lambda, file called : "DEG_ML_nc_justML.csv"  
Norm counts DEG ML JUST Medium vs Lambda *Top 50*, file called: "MvcL_norm_counts_Top50.csv"

Zscore ALL DEG ML JUST Medium vs Lambda, file called:"MvsL_Z_score.csv"
Zscore 50 DEG ML JUST Medium vs Lambda, file called : "MvsL_Top_50_Z_score.csv" 

Top 50 DEG Lambda vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Lambda/ML_50DEG.png")
```

ALL DEG Lambda vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Lambda/ML_DEG.png")
```

## B) Alpha vs Medium
### 1. Process the data

Extraction of the result table : Alpha vs Medium  

```{r}
#Compare "IFN-a" with "medium" as the base level 
res_unshrunken <- results(dds12_2, 
                         contrast = c("Stim", "alpha", "medium"), 
                         alpha = 0.05)
look(res_unshrunken)

plotMA(res_unshrunken, ylim=c(-8, 8), xlab="mean expression", colSig = "red3", colLine = "#ff000080", main = "DEG - Medium vs Alpha")
```
We can see that there is way more DEG when pDCs are activated with IFN-a compared with IFN-L.  


```{r}
#We should technically use this formula below, but I don't know why it does not work. So we have to use another formula. 
res1_shrunken <- lfcShrink(dds12_2,
                  contrast = c("Stim", "alpha", "medium"), 
                  res = res_unshrunken, 
                  type = "normal")  

#plot the results
plotMA(res1_shrunken, ylim=c(-4, 4), xlab="mean expression", colSig = "red3", colLine = "#ff000080", cex= 0.5,  main = "DEG - Medium vs alpha - 12h")
```

Now keep only the gene wich are actual DEG:  

```{r}
summary(res1_shrunken)
```


```{r}
#Results transformation into df:
res1_df <- res1_shrunken %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL")
head(res1_df)

# Join the res df with the gene symbols
res2 <- dplyr::inner_join(res1_df, annots, by = c("ENSEMBL"))
head(res2) 
```

```{r}
#Create a table retreving only DEG: 
ML_res2_sig <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 0.58)
nrow(ML_res2_sig)
```
**5549 genes** have a differential expression.  

Number of upregulated or downregulated DEG:
```{r}
MA_up <- subset(res2, padj < 0.05 & log2FoldChange > 0.58)
nrow(MA_up)
write_csv2(MA_up, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Alpha/DEG_MA_up.csv")

MA_down <- subset(res2, padj < 0.05 & log2FoldChange < (-0.58))
nrow(MA_down)
write_csv2(MA_down, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Alpha/DEG_MA_down.csv")
```

### 2. Volcano plot: 
Visualization ar volcano plot, but first create a df readable by ggplot:  

check this link to help:  
https://hbctraining.github.io/Training-modules/Visualization_in_R/lessons/03_advanced_visualizations.html  

```{r}
#Reoder df with padj: 
res1_table_volcano <-  res2 %>%
  arrange(padj)

#Create a column to indicate which genes to label
res1_table_volcano$diffexpressed <- "NO"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange > 0.5 & res1_table_volcano$pvalue < 0.05] <- "UP"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange < -0.5 & res1_table_volcano$pvalue < 0.05] <- "DOWN"

res1_table_volcano$delabel <- NA
res1_table_volcano$delabel[res1_table_volcano$diffexpressed != "NO"] <- res1_table_volcano$SYMBOL[res1_table_volcano$diffexpressed != "NO"]
look(res1_table_volcano)

#plot the volcano plot:
MA <- ggplot(data= res1_table_volcano, aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed, label=delabel))+
  geom_point()+
  ggtitle("Medium vs Alpha - 12h")+
  scale_colour_manual(values = c("blue", "grey60","red")) +
  geom_text_repel(max.overlaps = 8) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))
MA
```

```{r}
pdf("Volcano_MA.pdf", width = 10, height = 6)
MA
dev.off()
```

### 3. Creation of DEG tables 
#### a) All DEG

```{r}
#Extract the significative gene detected with the res2 from the normalized counts!
norm_counts_sig_ma <- dplyr::inner_join(ML_res2_sig, norm_counts, by = "ENSEMBL")
DEG_nc_IFNci <- norm_counts_sig_ma[, c(1,3,7,8, 16,23,30,37,44, 10,17,24,31,38, 11,18,25,32,39, 14,21,28,35,42, 15,22,29,36,43)]
look(DEG_nc_IFNci)
#Create a table and merge it with annotate.
DEG_nc_IFNci <- DEG_nc_IFNci %>% 
  as.data.frame()
#Save the table in excel file: 
write_csv2(DEG_nc_IFNci, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_MA_12_nc_IFNci.csv")


#Selected only the column Medium and Lambda
norm_counts_sig_ma <- norm_counts_sig_ma[,c(1,8, 16,23,30,37,44, 10,17,24,31,38)]
look(norm_counts_sig_ma)

#Create a table and merge it with annotate.
norm_counts_sig_ma <- norm_counts_sig_ma %>% 
  as.data.frame()
look(norm_counts_sig_ma)

#Save the table in excel file: 
write_csv2(norm_counts_sig_ma, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_MA_nc_justMA.csv")
```

#### b) All DEG - Z score
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_norm_counts_sig_ma <- norm_counts_sig_ma[,3:12]
MvsA_Z_score <- t(scale(t(Zscore_norm_counts_sig_ma)))
look(MvsA_Z_score)

MvsA_Z_score <- MvsA_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsA_Z_score$SYMBOL <- norm_counts_sig_ma[,2]
look(MvsA_Z_score)

#Write and export the excel
write_csv2(MvsA_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsA_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  


#### c) Top 50 DEG
```{r}
#Select the top 50 DEG: 
res2_sig_Top50 <- ML_res2_sig %>%
  arrange(padj)
head(res2_sig_Top50)
res2_sig_Top50 <- res2_sig_Top50[1:50,]
res2_sig_Top50 <- res2_sig_Top50 %>%
  arrange(log2FoldChange)
look(res2_sig_Top50)

#Revtreive norm counts: 
res2_sig_Top50 <- dplyr::inner_join(res2_sig_Top50, norm_counts, by = "ENSEMBL")
look(res2_sig_Top50)

#Selected only the column Medium and Lambda
MvcA_res2_sig_Top50 <- res2_sig_Top50[,c(1,8, 16,23,30,37,44, 10,17,24,31,38)]
look(MvcA_res2_sig_Top50)

#Create a table and merge it with annotate.
MvcA_res2_sig_Top50 <- MvcA_res2_sig_Top50 %>% 
  as.data.frame()
look(MvcA_res2_sig_Top50)

#Save the table: 
write_csv2(MvcA_res2_sig_Top50, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvcA_norm_counts_Top50.csv")
```

#### d) Z score - Top 50 DEG
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_Top_50_norm_counts <- MvcA_res2_sig_Top50[,3:12]
MvsL_Top_50_Z_score <- t(scale(t(Zscore_Top_50_norm_counts)))
look(MvsL_Top_50_Z_score)


MvsL_Top_50_Z_score <- MvsL_Top_50_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Top_50_Z_score$SYMBOL <- MvcA_res2_sig_Top50[,2]
look(MvsL_Top_50_Z_score)

#Write and export the excel
write_csv2(MvsL_Top_50_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsA_Top_50_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  

###4.Summary: 

**Nb DEG:**    
**Total: 5549**  
**UP: 2087**  
**DOWN: 3462**  


Plot the Volcano plot Alpha vs Medium: 
```{r}
MA
```

Norm counts DEG MA IFNs condition, file called: "DEG_MA_12_nc_IFNci.csv"  
Norm counts DEG MA JUST Medium vs alpha, file called : "DEG_MA_nc_justMA.csv"
Norm counts DEG MA JUST Medium vs Lambda *Top 50*, file called: "MvcA_norm_counts_Top50.csv"

Zscore ALL DEG ML JUST Medium vs Lambda, file called:"MvsA_Z_score.csv"
Zscore 50 DEG ML JUST Medium vs Lambda, file called : "MvsA_Top_50_Z_score.csv"  

Top 50 DEG Alpha vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Alpha/MA_50DEG.png")
```

ALL DEG Alpha vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Alpha/MA_DEG.png")
```

## C) Beta vs Medium
### 1. Process the data

Extraction of the result table : Beta vs Medium

```{r}
#Compare "IFN-a" with "medium" as the base level 
res_unshrunken <- results(dds12_2, 
                         contrast = c("Stim", "beta", "medium"), 
                         alpha = 0.05)
look(res_unshrunken)

plotMA(res_unshrunken, ylim=c(-8, 8), xlab="mean expression", colSig = "red3", colLine = "#ff000080", main = "DEG - Medium vs Alpha")
```
We can see that there is way more DEG when pDCs are activated with IFN-a compared with IFN-L.  


```{r}
#We should technically use this formula below, but I don't know why it does not work. So we have to use another formula. 
res1_shrunken <- lfcShrink(dds12_2,
                  contrast = c("Stim", "beta", "medium"), 
                  res = res_unshrunken, 
                  type = "normal")  

#plot the results
plotMA(res1_shrunken, ylim=c(-4, 4), xlab="mean expression", colSig = "red3", colLine = "#ff000080", cex= 0.5,  main = "DEG - Medium vs beta - 12h")
```

Now keep only the gene wich are actual DEG:  

```{r}
summary(res1_shrunken)
```

```{r}
#Results transformation into df:
res1_df <- res1_shrunken %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL")
head(res1_df)

# Join the res df with the gene symbols
res2 <- dplyr::inner_join(res1_df, annots, by = c("ENSEMBL"))
head(res2) 
```

```{r}
#Create a table retreving only DEG: 
ML_res2_sig <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 0.58)
nrow(ML_res2_sig)
```
**5769 genes** have a differential expression.  

Number of upregulated or downregulated DEG:
```{r}
MB_up <- subset(res2, padj < 0.05 & log2FoldChange > 0.58)
nrow(MB_up)
write_csv2(MB_up, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Beta/DEG_MB_up.csv")


MB_down <- subset(res2, padj < 0.05 & log2FoldChange < (-0.58))
nrow(MB_down)
write_csv2(MB_down, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Beta/DEG_MB_down.csv")
```

### 2. Volcano plot: 
Visualization ar volcano plot, but first create a df readable by ggplot:  

check this link to help:  
https://hbctraining.github.io/Training-modules/Visualization_in_R/lessons/03_advanced_visualizations.html  

```{r}
#Reoder df with padj: 
res1_table_volcano <-  res2 %>%
  arrange(padj)

#Create a column to indicate which genes to label
res1_table_volcano$diffexpressed <- "NO"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange > 0.5 & res1_table_volcano$pvalue < 0.05] <- "UP"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange < -0.5 & res1_table_volcano$pvalue < 0.05] <- "DOWN"

res1_table_volcano$delabel <- NA
res1_table_volcano$delabel[res1_table_volcano$diffexpressed != "NO"] <- res1_table_volcano$SYMBOL[res1_table_volcano$diffexpressed != "NO"]
look(res1_table_volcano)

#plot the volcano plot:
MB <- ggplot(data= res1_table_volcano, aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed, label=delabel))+
  geom_point()+
  ggtitle("Medium vs Beta - 12h")+
  scale_colour_manual(values = c("blue", "grey60","red")) +
  geom_text_repel(max.overlaps = 8) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))

MB
```

```{r}
pdf("Volcano_MB.pdf", width = 10, height = 6)
MB
dev.off()
```

### 3. Creation of DEG tables 
#### a) All DEG

```{r}
#Extract the significative gene detected with the res2 from the normalized counts!
norm_counts_sig_ma <- dplyr::inner_join(ML_res2_sig, norm_counts, by = "ENSEMBL")
DEG_nc_IFNci <- norm_counts_sig_ma[, c(1,3,7,8, 16,23,30,37,44, 10,17,24,31,38, 11,18,25,32,39, 14,21,28,35,42, 15,22,29,36,43)]
look(DEG_nc_IFNci)
#Create a table and merge it with annotate.
DEG_nc_IFNci <- DEG_nc_IFNci %>% 
  as.data.frame()
#Save the table in excel file: 
write_csv2(DEG_nc_IFNci, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_MB_12_nc_IFNci.csv")


#Selected only the column Medium and Lambda
norm_counts_sig_ma <- norm_counts_sig_ma[,c(1,8, 16,23,30,37,44, 11,18,25,32,39)]
look(norm_counts_sig_ma)

#Create a table and merge it with annotate.
norm_counts_sig_ma <- norm_counts_sig_ma %>% 
  as.data.frame()
look(norm_counts_sig_ma)

#Save the table in excel file: 
write_csv2(norm_counts_sig_ma, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_MB_nc_justMB.csv")
```

#### b) All DEG - Z score
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_norm_counts_sig_ma <- norm_counts_sig_ma[,3:12]
MvsA_Z_score <- t(scale(t(Zscore_norm_counts_sig_ma)))
look(MvsA_Z_score)

MvsA_Z_score <- MvsA_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsA_Z_score$SYMBOL <- norm_counts_sig_ma[,2]
look(MvsA_Z_score)

#Write and export the excel
write_csv2(MvsA_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsB_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  


#### c) Top 50 DEG
```{r}
#Select the top 50 DEG: 
res2_sig_Top50 <- ML_res2_sig %>%
  arrange(padj)
look(res2_sig_Top50)
res2_sig_Top50 <- res2_sig_Top50[1:50,]
res2_sig_Top50 <- res2_sig_Top50 %>%
  arrange(log2FoldChange)
look(res2_sig_Top50)

#Revtreive norm counts: 
res2_sig_Top50 <- dplyr::inner_join(res2_sig_Top50, norm_counts, by = "ENSEMBL")
look(res2_sig_Top50)

#Selected only the column Medium and Lambda
MvcA_res2_sig_Top50 <- res2_sig_Top50[,c(1,8, 16,23,30,37,44, 11,18,25,32,39)]
look(MvcA_res2_sig_Top50)

#Create a table and merge it with annotate.
MvcA_res2_sig_Top50 <- MvcA_res2_sig_Top50 %>% 
  as.data.frame()
look(MvcA_res2_sig_Top50)

#Save the table: 
write_csv2(MvcA_res2_sig_Top50, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvcB_norm_counts_Top50.csv")
```

#### d) Z score - Top 50 DEG
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_Top_50_norm_counts <- MvcA_res2_sig_Top50[,3:12]
MvsL_Top_50_Z_score <- t(scale(t(Zscore_Top_50_norm_counts)))
look(MvsL_Top_50_Z_score)


MvsL_Top_50_Z_score <- MvsL_Top_50_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Top_50_Z_score$SYMBOL <- MvcA_res2_sig_Top50[,2]
look(MvsL_Top_50_Z_score)

#Write and export the excel
write_csv2(MvsL_Top_50_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsB_Top_50_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  

###4.Summary: 

**Nb DEG:**    
**Total: 5769**  
**UP: 2185**  
**DOWN: 3584**  


Plot the Volcano plot Beta vs Medium: 
```{r}
MB
```

Norm counts DEG MA IFNs condition, file called: "DEG_MB_12_nc_IFNci.csv"  
Norm counts DEG MA JUST Medium vs alpha, file called : "DEG_MB_nc_justMB.csv"
Norm counts DEG MA JUST Medium vs Lambda *Top 50*, file called: "MvcB_norm_counts_Top50.csv"

Zscore ALL DEG ML JUST Medium vs Lambda, file called:"MvsB_Z_score.csv"
Zscore 50 DEG ML JUST Medium vs Lambda, file called : "MvsB_Top_50_Z_score.csv"  

Top 50 DEG Alpha vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Beta/MB_50DEG.png")
```

ALL DEG Alpha vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Beta/MB_DEG.png")
```

## D) Gamma vs Medium
### 1. Process the data

Extraction of the result table : Gamma vs Medium  

```{r}
#Compare "IFN-a" with "medium" as the base level 
res_unshrunken <- results(dds12_2, 
                         contrast = c("Stim", "gamma", "medium"), 
                         alpha = 0.05)
look(res_unshrunken)

plotMA(res_unshrunken, ylim=c(-8, 8), xlab="mean expression", colSig = "red3", colLine = "#ff000080", main = "DEG - Medium vs Gamma")
```
We can see that there is way more DEG when pDCs are activated with IFN-a compared with IFN-L.  


```{r}
#We should technically use this formula below, but I don't know why it does not work. So we have to use another formula. 
res1_shrunken <- lfcShrink(dds12_2,
                  contrast = c("Stim", "gamma", "medium"), 
                  res = res_unshrunken, 
                  type = "normal")  

#plot the results
plotMA(res1_shrunken, ylim=c(-4, 4), xlab="mean expression", colSig = "red3", colLine = "#ff000080", cex= 0.5,  main = "DEG - Medium vs Gamma - 12h")
```

Now keep only the gene wich are actual DEG:  

```{r}
summary(res1_shrunken)
```

```{r}
#Results transformation into df:
res1_df <- res1_shrunken %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL")
head(res1_df)

# Join the res df with the gene symbols
res2 <- dplyr::inner_join(res1_df, annots, by = c("ENSEMBL"))
head(res2) 
```

```{r}
#Create a table retreving only DEG: 
ML_res2_sig <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 0.58)
nrow(ML_res2_sig)
```
**182 genes** have a differential expression.  

Number of upregulated or downregulated DEG:
```{r}
MG_up <- subset(res2, padj < 0.05 & log2FoldChange > 0.58)
nrow(MG_up)
write_csv2(MG_up, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Gamma/DEG_MG_up.csv")

MG_down <- subset(res2, padj < 0.05 & log2FoldChange < (-0.58))
nrow(MG_down)
write_csv2(MG_down, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Gamma/DEG_MG_down.csv")
```

### 2. Volcano plot: 
Visualization ar volcano plot, but first create a df readable by ggplot:  

check this link to help:  
https://hbctraining.github.io/Training-modules/Visualization_in_R/lessons/03_advanced_visualizations.html  

```{r}
#Reoder df with padj: 
res1_table_volcano <-  res2 %>%
  arrange(padj)

#Create a column to indicate which genes to label
res1_table_volcano$diffexpressed <- "NO"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange > 0.5 & res1_table_volcano$pvalue < 0.05] <- "UP"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange < -0.5 & res1_table_volcano$pvalue < 0.05] <- "DOWN"

res1_table_volcano$delabel <- NA
res1_table_volcano$delabel[res1_table_volcano$diffexpressed != "NO"] <- res1_table_volcano$SYMBOL[res1_table_volcano$diffexpressed != "NO"]
look(res1_table_volcano)

#plot the volcano plot:
MG <- ggplot(data= res1_table_volcano, aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed, label=delabel))+
  geom_point()+
  ggtitle("Medium vs Gamma - 12h")+
  scale_colour_manual(values = c("blue", "grey60","red")) +
  geom_text_repel(max.overlaps = 8) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))
MG
```

```{r}
pdf("Volcano_MG.pdf", width = 10, height = 6)
MG
dev.off()
```

### 3. Creation of DEG tables 
#### a) All DEG

```{r}
#Extract the significative gene detected with the res2 from the normalized counts!
norm_counts_sig_ma <- dplyr::inner_join(ML_res2_sig, norm_counts, by = "ENSEMBL")
DEG_nc_IFNci <- norm_counts_sig_ma[, c(1,3,7,8, 16,23,30,37,44, 10,17,24,31,38, 11,18,25,32,39, 14,21,28,35,42, 15,22,29,36,43)]
look(DEG_nc_IFNci)
#Create a table and merge it with annotate.
DEG_nc_IFNci <- DEG_nc_IFNci %>% 
  as.data.frame()
#Save the table in excel file: 
write_csv2(DEG_nc_IFNci, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_MG_12_nc_IFNci.csv")


#Selected only the column Medium and Lambda
norm_counts_sig_ma <- norm_counts_sig_ma[,c(1,8, 16,23,30,37,44, 14,21,28,35,42)]
look(norm_counts_sig_ma)

#Create a table and merge it with annotate.
norm_counts_sig_ma <- norm_counts_sig_ma %>% 
  as.data.frame()
look(norm_counts_sig_ma)

#Save the table in excel file: 
write_csv2(norm_counts_sig_ma, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_MG_nc_justMG.csv")
```

#### b) All DEG - Z score
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_norm_counts_sig_ma <- norm_counts_sig_ma[,3:12]
MvsA_Z_score <- t(scale(t(Zscore_norm_counts_sig_ma)))
look(MvsA_Z_score)

MvsA_Z_score <- MvsA_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsA_Z_score$SYMBOL <- norm_counts_sig_ma[,2]
look(MvsA_Z_score)

#Write and export the excel
write_csv2(MvsA_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsG_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  


#### c) Top 50 DEG
```{r}
#Select the top 50 DEG: 
res2_sig_Top50 <- ML_res2_sig %>%
  arrange(padj)
look(res2_sig_Top50)
res2_sig_Top50 <- res2_sig_Top50[1:50,]
res2_sig_Top50 <- res2_sig_Top50 %>%
  arrange(log2FoldChange)
look(res2_sig_Top50)

#Revtreive norm counts: 
res2_sig_Top50 <- dplyr::inner_join(res2_sig_Top50, norm_counts, by = "ENSEMBL")
look(res2_sig_Top50)

#Selected only the column Medium and Lambda
MvcA_res2_sig_Top50 <- res2_sig_Top50[,c(1,8, 16,23,30,37,44, 14,21,28,35,42)]
look(MvcA_res2_sig_Top50)

#Create a table and merge it with annotate.
MvcA_res2_sig_Top50 <- MvcA_res2_sig_Top50 %>% 
  as.data.frame()
look(MvcA_res2_sig_Top50)

#Save the table: 
write_csv2(MvcA_res2_sig_Top50, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvcG_norm_counts_Top50.csv")
```

#### d) Z score - Top 50 DEG
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_Top_50_norm_counts <- MvcA_res2_sig_Top50[,3:12]
MvsL_Top_50_Z_score <- t(scale(t(Zscore_Top_50_norm_counts)))
look(MvsL_Top_50_Z_score)


MvsL_Top_50_Z_score <- MvsL_Top_50_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Top_50_Z_score$SYMBOL <- MvcA_res2_sig_Top50[,2]
look(MvsL_Top_50_Z_score)

#Write and export the excel
write_csv2(MvsL_Top_50_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/MvsG_Top_50_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  

###4.Summary: 

**Nb DEG:**    
**Total: 182**  
**UP: 96**  
**DOWN: 86**  


Plot the Volcano plot Gamma vs Medium: 
```{r}
MG
```

Norm counts DEG MA IFNs condition, file called: "DEG_MG_12_nc_IFNci.csv"  
Norm counts DEG MA JUST Medium vs alpha, file called : "DEG_MG_nc_justMG.csv"
Norm counts DEG MA JUST Medium vs Lambda *Top 50*, file called: "MvcG_norm_counts_Top50.csv"

Zscore ALL DEG ML JUST Medium vs Lambda, file called:"MvsG_Z_score.csv"
Zscore 50 DEG ML JUST Medium vs Lambda, file called : "MvsG_Top_50_Z_score.csv"  

Top 50 DEG Alpha vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Gamma/MG_50DEG.png")
```

ALL DEG Alpha vs Medium: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Med_Gamma/MG_DEG.png")
```



#II. DEG Lambda vs Alpha
### 1. Process the data

Extraction of the result table : Lambda vs Alpha  

```{r}
#Compare "IFN-L" with "IFN-a"  as the base level 
res_unshrunken <- results(dds12_2, 
                         contrast = c("Stim", "lambda", "alpha"), 
                         alpha = 0.05)
look(res_unshrunken)

plotMA(res_unshrunken, ylim=c(-8, 8), xlab="mean expression", colSig = "red3", colLine = "#ff000080", main = "DEG - Alpha vs Lambda - 12h")
```
We can see that there is way more DEG when pDCs are activated with IFN-a compared with IFN-L.  


```{r}
#We should technically use this formula below, but I don't know why it does not work. So we have to use another formula. 
res1_shrunken <- lfcShrink(dds12_2,
                  contrast = c("Stim", "lambda", "alpha"), 
                  res = res_unshrunken, 
                  type = "normal")  

#plot the results
plotMA(res1_shrunken, ylim=c(-4, 4), xlab="mean expression", colSig = "red3", colLine = "#ff000080", cex= 0.5,  main = "DEG - Alfa vs Lambda - 12h")
```

Now keep only the gene wich are actual DEG:  

```{r}
summary(res1_shrunken)
```


```{r}
#Results transformation into df:
res1_df <- res1_shrunken %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL")
look(res1_df)

# Join the res df with the gene symbols
res2 <- dplyr::inner_join(res1_df, annots, by = c("ENSEMBL"))
look(res2) 
```

```{r}
#Create a table retreving only DEG: 
ML_res2_sig <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 0.58)
nrow(ML_res2_sig)
```
**2291 genes** have a differential expression.  

Number of upregulated or downregulated DEG:
```{r}
MG_up <- subset(res2, padj < 0.05 & log2FoldChange > 0.58)
nrow(MG_up)
write_csv2(MG_up, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Alpha_Lambda/DEG_AL_up.csv")

MG_down <- subset(res2, padj < 0.05 & log2FoldChange < (-0.58))
nrow(MG_down)
write_csv2(MG_down, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/12H/Alpha_Lambda/DEG_AL_down.csv")
```

### 2. Volcano plot: 
Visualization ar volcano plot, but first create a df readable by ggplot:  

```{r}
#Reoder df with padj: 
res1_table_volcano <-  res2 %>%
  arrange(padj)

#Create a column to indicate which genes to label
res1_table_volcano$diffexpressed <- "NO"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange > 0.5 & res1_table_volcano$pvalue < 0.05] <- "UP"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange < -0.5 & res1_table_volcano$pvalue < 0.05] <- "DOWN"

res1_table_volcano$delabel <- NA
res1_table_volcano$delabel[res1_table_volcano$diffexpressed != "NO"] <- res1_table_volcano$SYMBOL[res1_table_volcano$diffexpressed != "NO"]
head(res1_table_volcano)

#plot the volcano plot:
AL <- ggplot(data= res1_table_volcano, aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed, label=delabel))+
  geom_point()+
  ggtitle("Alpha vs Lambda - 12h")+
  scale_colour_manual(values = c("blue", "grey60","red")) +
  geom_text_repel(max.overlaps = 8) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))
AL
```

```{r}
pdf("Volcano_AL.pdf", width = 10, height = 6)
AL
dev.off()
```

### 3. Creation of DEG tables 
#### a) All DEG

```{r}
#Extract the significative gene detected with the res2 from the normalized counts!
norm_counts_sig_ma <- dplyr::inner_join(ML_res2_sig, norm_counts, by = "ENSEMBL")
look(norm_counts_sig_ma)

#Selected only the column Medium and Lambda
norm_counts_sig_ma <- norm_counts_sig_ma[,c(1,8,10,17,24,31,38,15,22,29,36,43)]
look(norm_counts_sig_ma)

#Create a table and merge it with annotate.
norm_counts_sig_ma <- norm_counts_sig_ma %>% 
  as.data.frame()
look(norm_counts_sig_ma)

#Save the table in excel file: 
write_csv2(norm_counts_sig_ma, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/DEG_AL_nc_justAL.csv")
```

#### b) All DEG - Z score
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_norm_counts_sig_ma <- norm_counts_sig_ma[,3:12]
MvsA_Z_score <- t(scale(t(Zscore_norm_counts_sig_ma)))
look(MvsA_Z_score)

MvsA_Z_score <- MvsA_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsA_Z_score$SYMBOL <- norm_counts_sig_ma[,2]
look(MvsA_Z_score)

#Write and export the excel
write_csv2(MvsA_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/AvsL_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  


#### c) Top 50 DEG
```{r}
#Select the top 50 DEG: 
res2_sig_Top50 <- ML_res2_sig %>%
  arrange(padj)
res2_sig_Top50 <- res2_sig_Top50[1:50,]
res2_sig_Top50 <- res2_sig_Top50 %>%
  arrange(log2FoldChange)
look(res2_sig_Top50)

#Revtreive norm counts: 
res2_sig_Top50 <- dplyr::inner_join(res2_sig_Top50, norm_counts, by = "ENSEMBL")
look(res2_sig_Top50)

#Selected only the column Medium and Lambda
MvcA_res2_sig_Top50 <- res2_sig_Top50[,c(1,8,10,17,24,31,38,15,22,29,36,43)]
look(norm_counts_sig_ml)

#Create a table and merge it with annotate.
MvcA_res2_sig_Top50 <- MvcA_res2_sig_Top50 %>% 
  as.data.frame()
look(MvcA_res2_sig_Top50)

#Save the table: 
write_csv2(MvcA_res2_sig_Top50, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/AvsL_norm_counts_Top50.csv")
```

#### d) Z score - Top 50 DEG
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_Top_50_norm_counts <- MvcA_res2_sig_Top50[,3:12]
MvsL_Top_50_Z_score <- t(scale(t(Zscore_Top_50_norm_counts)))
look(MvsL_Top_50_Z_score)


MvsL_Top_50_Z_score <- MvsL_Top_50_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Top_50_Z_score$SYMBOL <- MvcA_res2_sig_Top50[,2]
look(MvsL_Top_50_Z_score)

#Write and export the excel
write_csv2(MvsL_Top_50_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/AvsL_Top_50_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  

###4.Summary: 

**Nb DEG:**    
**Total: 2291**  
**UP: 1391**  
**DOWN: 900**  


Plot the Volcano plot Gamma vs Medium: 
```{r}
AL
```

Norm counts DEG AL JUST Lambda vs Alpha, file called : "DEG_AL_nc_justAL.csv"
Norm counts DEG AL JUST Lambda vs Alpha *Top 50*, file called: "AvsL_norm_counts_Top50.csv"

Zscore ALL DEG AL JUST Medium vs Lambda, file called:"AvsL_Z_score.csv"
Zscore 50 DEG AL JUST Lambda vs Alpha, file called : "AvsL_Top_50_Z_score.csv"  

Top 50 DEG Lambda vs Alpha: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/Alpha_Lambda/AL_50DEG.png")
```

ALL DEG Lambda vs Alpha: 
```{r}
knitr::include_graphics("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/Alpha_Lambda/AL_DEG.png")
```



#III. DEG CpG-L vs CpG
### 1. Process the data

Extraction of the result table : CpG-L vs CpG

```{r}
#For now I would like to compare "IFN-L" with "medium" as the base level 
res_unshrunken <- results(dds12_2, 
                         contrast = c("Stim", "CpG_L", "CpG"), 
                         alpha = 0.05)
look(res_unshrunken)

plotMA(res_unshrunken, ylim=c(-8, 8), xlab="mean expression", colSig = "red3", colLine = "#ff000080", main = "DEG - CpG-L vs CpG 12h")
```
 

```{r}
#We should technically use this formula below, but I don't know why it does not work. So we have to use another formula. 
res1_shrunken <- lfcShrink(dds12_2,
                  contrast = c("Stim", "CpG_L", "CpG"), 
                  res = res_unshrunken, 
                  type = "normal")  

#plot the results
plotMA(res1_shrunken, ylim=c(-4, 4), xlab="mean expression", colSig = "red3", colLine = "#ff000080", cex= 0.5,  main = "DEG - CpG-L vs CpG 12h")
```

Now keep only the gene wich are actual DEG:  

```{r}
summary(res1_shrunken)
```


```{r}
#Results transformation into df:
res1_df <- res1_shrunken %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "ENSEMBL")
look(res1_df)

# Join the res df with the gene symbols
res2 <- dplyr::inner_join(res1_df, annots, by = c("ENSEMBL"))
look(res2) 
```


```{r}
#Create a table retreving only DEG: 
ML_res2_sig <- subset(res2, padj < 0.05 & abs(log2FoldChange) > 0.58)
nrow(ML_res2_sig)
summary(ML_res2_sig)
```
**1328 DEG** found when comparing IFN-L to the Medium.   

```{r}
ML_up <- subset(res2, padj < 0.05 & log2FoldChange > 0.58)
nrow(ML_up)

ML_down <- subset(res2, padj < 0.05 & log2FoldChange < (-0.58))
nrow(ML_down)
```

### 2. Volcano plot: 
Visualization ar volcano plot, but first create a df readable by ggplot:  

```{r}
#Reoder df with padj: 
res1_table_volcano <-  res2 %>%
  arrange(padj)

#Create a column to indicate which genes to label
res1_table_volcano$diffexpressed <- "NO"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange > 0.5 & res1_table_volcano$pvalue < 0.05] <- "UP"
res1_table_volcano$diffexpressed[res1_table_volcano$log2FoldChange < -0.5 & res1_table_volcano$pvalue < 0.05] <- "DOWN"

res1_table_volcano$delabel <- NA
res1_table_volcano$delabel[res1_table_volcano$diffexpressed != "NO"] <- res1_table_volcano$SYMBOL[res1_table_volcano$diffexpressed != "NO"]
head(res1_table_volcano)

#plot the volcano plot:
CCL <- ggplot(data= res1_table_volcano, aes(x=log2FoldChange, y=-log10(padj), colour=diffexpressed, label=delabel))+
  geom_point()+
  ggtitle("CpG-L vs CpG 12h")+
  scale_colour_manual(values = c("blue", "grey60","red")) +
  geom_text_repel(max.overlaps = 8) +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))
CCL
ggsave(filename = "CCL_Volcano.png", path="C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG", last_plot(), width = 10, height = 6, units = "in", dpi=300)
```


```{r}
pdf("Volcano_CCL.pdf", width = 10, height = 6)
CCL
dev.off()
```


### 3. Creation of DEG tables 
#### a) All DEG
```{r}
#Extract the significative gene detected with the res2 from the normalized counts!
norm_counts <- read.delim("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/1_Processing/normalized_counts12.txt")
colnames(norm_counts)[colnames(norm_counts) == "Gene"] <- "ENSEMBL"
norm_counts_sig_ml <- dplyr::inner_join(ML_res2_sig, norm_counts, by = "ENSEMBL")
norm_counts_sig_ml <- norm_counts_sig_ml[,c(1,8, 13,20,27,34,41, 12,19,26,33,40)]
look(norm_counts_sig_ml)

#Create a table and merge it with annotate.
norm_counts_sig_ml <- norm_counts_sig_ml %>% 
  as.data.frame()
head(norm_counts_sig_ml)

#Save the table in excel file: 
write_csv2(norm_counts_sig_ml, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/CCL_normcounts.csv")
```

#### b) All DEG - Z score
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_norm_counts_sig_ml <- norm_counts_sig_ml[,3:12]
MvsL_Z_score <- t(scale(t(Zscore_norm_counts_sig_ml)))
look(MvsL_Z_score)

MvsL_Z_score <- MvsL_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Z_score$SYMBOL <- norm_counts_sig_ml[,2]
look(MvsL_Z_score)

#Write and export the excel
write_csv2(MvsL_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/CCL_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  


#### c) Top 50 DEG
```{r}
#Select the top 50 DEG: 
res2_sig_Top50 <- ML_res2_sig %>%
  arrange(padj)
head(res2_sig_Top50)
res2_sig_Top50 <- res2_sig_Top50[1:50,]
res2_sig_Top50 <- res2_sig_Top50 %>%
  arrange(log2FoldChange)
look(res2_sig_Top50)

#Revtreive norm counts: 
res2_sig_Top50 <- dplyr::inner_join(res2_sig_Top50, norm_counts, by = "ENSEMBL")
look(res2_sig_Top50)

#Selected only the column Medium and Lambda
MvcL_res2_sig_Top50 <- res2_sig_Top50[,c(1,8, 13,20,27,34,41, 12,19,26,33,40)]
look(norm_counts_sig_ml)

#Create a table and merge it with annotate.
MvcL_res2_sig_Top50 <- MvcL_res2_sig_Top50 %>% 
  as.data.frame()
look(MvcL_res2_sig_Top50)

#Save the table: 
write_csv2(MvcL_res2_sig_Top50, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/CCL_norm_counts_Top50.csv")
```

#### d) Z score - Top 50 DEG
In order to plot our DEG, we should use Z-score:  
Z-score = (Sample value - mean)/SD

```{r}
#Perform  the transformation into Z-score: 
#1. Withdraw the column containing chaacters...: 
Zscore_Top_50_norm_counts <- MvcL_res2_sig_Top50[,3:12]
MvsL_Top_50_Z_score <- t(scale(t(Zscore_Top_50_norm_counts)))
look(MvsL_Top_50_Z_score)


MvsL_Top_50_Z_score <- MvsL_Top_50_Z_score %>%
  as.data.frame()

#Re add the SYMBOL Column: 
MvsL_Top_50_Z_score$SYMBOL <- MvcL_res2_sig_Top50[,2]
look(MvsL_Top_50_Z_score)

#Write and export the excel
write_csv2(MvsL_Top_50_Z_score, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/2_DEG/CCL_Top_50_Z_score.csv")
```

Quick checking of what we just did on Morpheus:  
https://software.broadinstitute.org/morpheus/  


