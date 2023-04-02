---
title: "pDCs Microbulk RNA seq 12h de stim"
author: "Candice SAKREF"
date: "21/10/2021"
output: 
  html_document:
    toc: true
    toc_float: true
    toc_depth: 5
---
We received the  15/12/2020 the project 134 gene expression table created by  Maude ARDIN.  
This project gathered 35 RNA samples extracted from pDCs microbulks stimulated in different conditions. 

```{r}
#Definition de l'environnement de travail 
setwd("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/1_Processing")

#Overview of the first elements of the df
look <- function(x){
  head(x)[1:5, 1:5]
}
```

```{r message = FALSE}
library(dplyr)
library(tidyverse)
library(ggplot2)
library(vsn)
library(DESeq2)
library(ggfortify)
```

# I. Data description  
### a) Samples 
**pDCs** were purified with a miltenyi purification kit and 100 000 pDCs were stimulated 6h or 12h in this conditions:  
1: Medium  
2: IFN-alfa  
3: IFN-beta  
4: IFN-gamma  
5: IFN-lambda
6: CpG B
7: CpG B + IFN-lambda

Thus, we have 5 x 7 **2** conditions.

## b) Alignment and QCs
Performed by Maude ARDIN

# II. Import of expression and meta data
### a) Expression data

Fisrt test the dataframe given by Maude: 
This "rds" document contains 3 table *(~Numpy array)*:  
1. counts  
2. lenghts  
3. Abundance

We can now import of the expression data: 
```{r}
txi12 <- read_rds("salmon_txi_Run_PGC204.rds")
txi6 <- read_rds("salmon_txi_projet134.rds")

#Transform the txi "R object" into a dataframe that we can process.
# Abundance = TPM  
TPM12 <- as.data.frame(txi12$abundance)
TPM6 <- as.data.frame(txi6$abundance)

#head(TPM)
#here we used the function "head()" instead of "look" because I wanted to see the names of all the columns! In order to create a correct meta table. See.Below. 
look(TPM12)
look(TPM6)
```
*Why do we look at abundance and not counts? Doesn't matter for now I think*     

The colmun names give you the order of your sample.  
I had to create a meta-data table that will  have the **SAME ORDER** as the column order. 

**CAREFUL** 
It will also be interesting to fuse the data from pDCs stimulated 6h and 12h hours togather. We have to combine the expression matrix for the 2 runs: 
```{r}
#First with the abundance table
#tranforme the matrix in df to merge the by row: 
df_tpm_run1 <- txi12$abundance %>% as.data.frame() %>% rownames_to_column("ensgene")
df_tpm_run2 <- txi6$abundance %>% as.data.frame() %>% rownames_to_column("ensgene")
abundance <- dplyr::full_join(df_tpm_run1, df_tpm_run2, by = "ensgene")
nrow(abundance)

#(!) as Deseq does not accept "NA" will have to replace them by 0: 
abundance[is.na(abundance)] = 0

#process the table so that we can have a correct matrix
row.names(abundance) <- abundance$ensgene
abundance <- abundance[,2:71]
abundance <- as.matrix(abundance)
look(abundance)


#Second with the counts: 
df_counts_run1 <- txi12$counts %>% as.data.frame() %>% rownames_to_column("ensgene")
df_counts_run2 <- txi6$counts %>% as.data.frame() %>% rownames_to_column("ensgene")
counts <- dplyr::full_join(df_counts_run1, df_counts_run2, by = "ensgene")
counts[is.na(counts)] = 0
row.names(counts) <- counts$ensgene
counts <- counts[,2:71]
counts <- as.matrix(counts)
look(counts)

#Third with the length: 
df_length_run1 <- txi12$length %>% as.data.frame() %>% rownames_to_column("ensgene")
df_length_run2 <- txi6$length %>% as.data.frame() %>% rownames_to_column("ensgene")
length <- dplyr::full_join(df_length_run1, df_length_run2, by = "ensgene")
length[is.na(length)] = 1 #(!) DESeq only accept length > 0
row.names(length) <- length$ensgene
length <- length[,2:71]
length <- as.matrix(length)
look(length)

#Fourth, DO NOT forget this elements:
countsFromAbundance <- "no"

#Once we hav all this files we should merge them in a list of element => txi: 
txiall <- list(abundance, counts, length, countsFromAbundance)
names(txiall) <- c("abundance", "counts", "length", "countsFromAbundance")

#Create the TPM df: 
TPM_all <- as.data.frame(txiall$abundance)
```



### b) Meta data table
For the following steps it is very important to create a Meta_data  table containing all the "Biological" information that could be of great importance on the analysis. 

```{r}
#Register the meta data inforation in R. 
meta12 <- read.table("Meta_data_Run2.txt", sep = "\t", header = TRUE)
str(meta12)

meta6 <- read.table("Meta_data_Run1.txt", sep = "\t", header = TRUE)
str(meta6)

metaall <- read.table("Meta_Data_All.txt", sep = "\t", header = TRUE)
str(metaall)
```
After the meta_data table has been create we can now fuse it with the txi document to create an **R object**.   
But first we have to check that our Meta_data actually matches with the expression expression data table to avoid any incorrect association.  

Reoder meta:  
```{r}
#"medium", "alpha", "beta", "gamma", "lambda", "CpG", "CpG_L"
meta12$Stim<- factor(meta12$Stim, levels = c("medium", "alpha", "beta", "gamma", "lambda", "CpG", "CpG_L"))
meta6$Stim<- factor(meta6$Stim, levels = c("medium", "alpha", "beta", "gamma", "lambda", "CpG", "CpG_L"))
metaall$tim<- factor(metaall$Stim, levels = c("medium", "alpha", "beta", "gamma", "lambda", "CpG", "CpG_L"))
str(meta6$Stim)
```


```{r}
# Are the colnames(TPM) the same AND in the same order as IDs given by Maude ? 
all(colnames(TPM12) %in% meta12$ID_Maude)
all(colnames(TPM6) %in% meta6$ID_Maude)
all(colnames(TPM_all) %in% metaall$ID_Maude)
```
The colnames of the expression data do match with the ID_Maude colmun in the meta_data table AND are in the same order. Perfect. Keep proceeding.  

In order to have a expression df clearer, it would be useful to rename its columns: 
```{r}
#12h
colnames(txi12$abundance) <- meta12$Sample_ID
colnames(txi12$counts) <- meta12$Sample_ID
colnames(txi12$length) <- meta12$Sample_ID

#6h
colnames(txi6$abundance) <- meta6$Sample_ID
colnames(txi6$counts) <- meta6$Sample_ID
colnames(txi6$length) <- meta6$Sample_ID

#Both 
colnames(txiall$abundance) <- metaall$Sample_ID
colnames(txiall$counts) <- metaall$Sample_ID
colnames(txiall$length) <- metaall$Sample_ID
```

#II. DESeq-2: 
### a) Creation of the "R object"

For the analysis of our RNA transxript, we are going to use the **DESeq2 package** which works with object called a "DESeqDataSet" from the expression datat table and the meta_data table containing the sample informations. 


A `DESeqDataSet` object must have an associated *design formula* where you indicated the sources of variation (added in the Meta_Table) you have to control and the factor of interest to test during differential expression testing.   
The design formula should have all of the factors in the meta df that account for major sources of (unwanted) variation, here the "Donor" information. The last factor entered in the formula shoud be the condition of interest, here the "Stim" information.   

```{r}
dds12 <- DESeqDataSetFromTximport(txi12,
                                   colData = meta12,
                                   design =  ~ Donor + Stim)

dds6 <- DESeqDataSetFromTximport(txi6,
                                   colData = meta6,
                                   design =  ~ Donor + Stim)

ddsall <- DESeqDataSetFromTximport(txiall,
                                   colData = metaall,
                                   design =  ~ Flow_cell_ID + Stim)

```


In this formula the batch effect hasn't be taken into account yet.  

*The design formula can be changed later, however then all differential analysis steps should be repeated, as the design formula is used to estimate the dispersions and to estimate the log2 fold changes of the model*.  
  
```{r}
dds12
```
**For pDCs stimulated 12h => 57 084 genes were retreive** in 35 samples  

```{r}
dds6
```

**For pDCs stimulated 12h => 21 386 genes were retreive** in 35 samples  
*As a reference Anais with her NK => 60 199 genes were retreive*  

*Here Margaux retrieve the normalized counts before the pre-filtering. I don't really understand why we are doing now, and it is a step that is not mentioned in the different documentation that I have. So I am not going to do this step. I will come back to it if necessary.*   

```{r}
ddsall
```
**For pDCs stimulated 12h => 58864 genes were retreive** in 70 samples 
Which means that only **1 780** were expressed in 6h hours and not in 12h stim.   

## b) Pre-filtering the dataset

```{r}
#Here it is recommended in to only filtered the gene that have counts = 0.  
keep12 <- rowSums(counts(dds12)) >= 10
dds12_2 <- dds12[keep12, ]
dds12_2
```
**22 542 genes left!** x 35 samples.  
So we removed (57084 - 22542 = ) **34 542 genes** which have counts =< 10.

```{r}
#Here it is recommended in to only filtered the gene that have counts = 0.  
keep6 <- rowSums(counts(dds6)) >= 10
dds6_2 <- dds6[keep6, ]
dds6_2
```
**16 922 genes left!** x 35 samples.  
So we removed (21386 - 16922 = ) **4 464 genes** which have counts =< 10.

```{r}
#Here it is recommended in to only filtered the gene that have counts = 0.  
keepall <- rowSums(counts(ddsall)) >= 10
ddsall_2 <- ddsall[keepall, ]
ddsall_2
```
**24 837 genes left!** x 70 samples.  
So we removed (58864 - 24837 = ) **34 027 genes** which have counts =< 10.


*It would be useful to save the dds2 df*  
```{r}
write_rds(dds6_2, "dds6_2.rds")
write_rds(dds12_2, "dds12_2.rds")
write_rds(ddsall_2, "ddsall_2.rds")
```


## c) Retreive the normalized counts:  
For further analysis it is recommanded to retrieve the normalized counts. 
```{r}
#12h
dds12_2 <- estimateSizeFactors(dds12_2)
norm_counts12 <- counts(dds12_2, normalize =TRUE)
look(norm_counts12)

norm_counts_table12 <- as.data.frame(counts(dds12_2, normalize =TRUE))
look(norm_counts_table12)

# Save normalized count df  
norm_counts_table12 %>%
  rownames_to_column(var = "Gene") %>%
  write.table("normalized_counts12.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
```


```{r}
#6h
dds6_2 <- estimateSizeFactors(dds6_2)
norm_counts6 <- counts(dds6_2, normalize =TRUE)
look(norm_counts6)

norm_counts_table6 <- as.data.frame(counts(dds6_2, normalize =TRUE))
look(norm_counts_table6)

# Save normalized count df  
norm_counts_table6 %>%
  rownames_to_column(var = "Gene") %>%
  write.table("normalized_counts6.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
```

```{r}
#Both
ddsall_2 <- estimateSizeFactors(ddsall_2)
norm_countsall <- counts(ddsall_2, normalize =TRUE)
look(norm_countsall)

norm_counts_tableall <- as.data.frame(counts(ddsall_2, normalize =TRUE))
look(norm_counts_tableall)

# Save normalized count df  
norm_counts_tableall %>%
  rownames_to_column(var = "Gene") %>%
  write.table("normalized_countsall_Flowcell.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
```

## d) The variance stabilizing transformamtion (vst) and the **rlog**  
If it is plan to perform a PCA in the following step it is mendatory to pass throught a 'vsn' or a 'rlog' transformation. Why?  
If a PCA is perform with the raw datas, the genes bringing the highest variation will be the genes that have the highest counts. Because they will have the highest variation between samples.   

That is why it is recommanded to pass by a "vsn" or a "rlog" transformation. 

**The rlog tends to work well on small datasets (n<30, we are at n=5), potenyially outperforming the "vst" when there is a wide range of sequencing depth across samples (it is not our case with this conditions, but good to know).**

```{r}
rld12 <- rlog(dds12_2, blind = TRUE)
rld6 <- rlog(dds6_2, blind = TRUE)
rldall <- rlog(ddsall_2, blind = TRUE)
#Blind = TRUE, for a complete unsupervise transformation.
```

Assess if th tranformation have work correctly:  
1. Plot the untransfromed data:   
```{r}
meanSdPlot(assay(dds12))
meanSdPlot(assay(dds6))
meanSdPlot(assay(ddsall))
```

2. Plot the **transfromed** data: 
```{r}
meanSdPlot(assay(rld12))
meanSdPlot(assay(rld6))
meanSdPlot(assay(rldall))
```

The transformation of the data with the "rlog" function has worked. We can see that the standard deviation is quite homogenous between all the genes (low or highly expressed).  

## e) PCA plot
Plot the PCA with the top 500 most variant genes: Colors = Stim 
```{r}
PCA_act6 <- plotPCA(rld6, 
                     intgroup = "Stim", 
                     ntop = 500)

PCA_act6
```
```{r}
pdf("PCA_act6.pdf", width = 10, height = 6)
PCA_act6
dev.off()
```


**We can clearly see here 2 effects:**  
**1. The PC1 looks like it relays on the stimulation with IFN and CpG activation.**  
**2. The PC2  looks like it relays on the donnor variation or maybe on the Bactch effect.**   

**Hypothesis: 1. The "Donnor" variation hasn't been removed correclty**  
            **2. The "Batch-extraction" is really important and have to be included in the formula (right at the begining).**  

Moreover, the PCA1 and PCA2 togather only explain 50% of the variation between samples. which is not a lot... 


However, Whe we lookat what is happening at **12h**  

```{r}
PCA_act12 <- plotPCA(rld12, 
                     intgroup = "Stim", 
                     ntop = 500)

PCA_act12
```
**We can clearly see here 2 effects:**  
**1. The PC1 looks like it relays on the stimulation with IFN and CpG activation.**  
**2. The PC2  looks like it relays on someting else.**   

**Hypothesis PC2 could be due to: 1. Donor** 
$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$**2. Gender**  
$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$**3. Blood type**
$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$**4. Age**
$~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~$**5. D antigen**

```{r}
PCA_actall <- plotPCA(rldall, 
                     intgroup = "Stim", 
                     ntop = 500)

PCA_actall
```
```{r}
pdf("PCA_actall.pdf", width = 10, height = 6)
PCA_actall
dev.off()
```

Test the run: 
```{r}
PCA_Run <- plotPCA(rldall, 
                     intgroup = "Flow_cell_ID", 
                     ntop = 500)

PCA_Run
```
```{r}
pdf("PCA_Run.pdf", width = 10, height = 6)
PCA_Run
dev.off()
```
IT is clear now that the:  
**- PC1 => Stim with or with CpG**   
**- PC2 => the Run**  

However, we here have a confunding factor, because the Run 1 is the run havin all the 6h stim condition, and the run 2 have all the 12h stim conditions. Thus, we cannot know of the PC2 is due by the Run or by the time of stimulation. 

We could look now at the PC3: 
```{r}
ntop <- 500
rv <- rowVars(assay(rldall))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
df_plot <- t(assay(rldall)[select, ])
PC <- prcomp(df_plot)
Sex <- metaall$Sex
Flow_cell_ID <- metaall$Flow_cell_ID

PCi <- data.frame(PC$x, Sex, Flow_cell_ID)

# PCA with PC3 & 2 
ggplot(PCi, aes(x = PC3, y = PC2, col = Sex, shape = Flow_cell_ID)) +
  geom_point(size = 3) + 
  labs(x = "PC3", y = "PC2: 21% variance")
```
```{r}
pdf("PCA_3andSex.pdf", width = 10, height = 6)
ggplot(PCi, aes(x = PC3, y = PC2, col = Sex, shape = Flow_cell_ID)) +
  geom_point(size = 3) + 
  labs(x = "PC3", y = "PC2: 21% variance")
dev.off()
```

```{r}
ntop <- 500
rv <- rowVars(assay(rldall))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
df_plot <- t(assay(rldall)[select, ])
PC <- prcomp(df_plot)
Sex <- metaall$Sex
Age <- metaall$Age
Flow_cell_ID <- metaall$Flow_cell_ID

PCi <- data.frame(PC$x, Age, Flow_cell_ID)

# PCA with PC3 & 2 
ggplot(PCi, aes(x = PC3, y = PC2, col = Age, shape = Flow_cell_ID)) +
  geom_point(size = 3) + 
  labs(x = "PC3", y = "PC2: 21% variance")
```

```{r}
pdf("PCA_3andAge.pdf", width = 10, height = 6)
ggplot(PCi, aes(x = PC3, y = PC2, col = Age, shape = Flow_cell_ID)) +
  geom_point(size = 3) + 
  labs(x = "PC3", y = "PC2: 21% variance")
dev.off()
```
**- PC1 => Stim with or with CpG**   
**- PC2 => the Run** 
**- PC3 => Sex (confunding with Age too)**  

What about PC4? 
Could it be due to the donor?  
```{r}

```

