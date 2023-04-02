---
title: "GSEA"
author: "Candice SAKREF"
date: "18/03/2022"
output: 
  html_document:
    toc: true
    toc_float: true
    toc_depth: 5
---

After we obtain all the DEG lists we wanted "2_DEG", Gene Set Enrichment Analysis GSEA aims to understand which pathways or gene networks the differentially expressed genes are implicated in.

#I. How to Run GSEA

You can find help for the whole procedure here:  
https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html   

Before hand, you have to downlaod the GSEA software (see.link above)  

##A) Laod the data

Once you installed GSEA software you have to generate 3 file correclty edited. 
$~~~~~$1. Norm count file with your conditions of interest **.gct** 
$~~~~~$2. Phenotype Labels **.cls**
$~~~~~$3. Chip annotation *toujours le même fichier transmit par Margaux HUBERT*

###1. ".gct" file
This file is a Norm count file with your conditions of interest, and containing the total number of genes. It has a specific extension: **.gct**  

```{r}
#Call Libraries: 
library(fgsea)
library(dplyr)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(nlme)
```

```{r}
#Overview of the first elements of the df
look <- function(x){
  head(x)[1:5, 1:5]
}
```

Prepare the Norm count table obtained from
```{r}
#Import the norm counts
norm_counts <- read.delim("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/1_Processing/normalized_counts12.txt")

norm_counts <- norm_counts %>% dplyr::rename(ENSEMBL=Gene)

#Download the annotations: 
uniKeys <- keys(org.Hs.eg.db, keytype="ENSEMBL")
cols <- c("SYMBOL", "GENENAME")
annots <- select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL")

#Merge Norm counts with the anotation
norm_counts <- dplyr::inner_join(norm_counts, annots, by = c("ENSEMBL"))

#To perfor the GSEA you will need the number of total gene: 
nrow(norm_counts)
```

####a) Lambda vs Medium: 

```{r}
colnames(norm_counts)
normcL <- norm_counts[, c(37, 8,15,22,29,36, 7,14,21,28,35)]
colnames(normcL)
normcL <- normcL %>% 
  dplyr::rename(NAME=SYMBOL) %>% 
  mutate(Description=NAME, .after=NAME)
write.table(normcL, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/GCT/Lambda2.gct.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

Once you have this file you have to: 
1. Open it with **Sublime Text**
2. Make sure there is no added row at the botton of the file
3. Add a line and write **"#1.2"** *Its always this*
4. Add a second line: **"Number of row"** TABULATION **Number of samples**
5. Ctrl+A => click on **Tab Size:4** and select **Convert Indentation to Tabs**
6. Save the file

####a) Lambda vs Alpha: 

```{r}
colnames(norm_counts)
normcL <- norm_counts[, c(37, 2,9,16,23,30, 7,14,21,28,35)]
colnames(normcL)
normcL <- normcL %>% 
  dplyr::rename(NAME=SYMBOL) %>% 
  mutate(Description=NAME, .after=NAME)
write.table(normcL, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/GCT/AL.gct.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

###2. Chip Annotation:

Margaux gave me her Chip annotation file.  
This file should work with all of my future GSEA run. 

###3. Phenotype labels: 

To create tis file you should follow this following steps: 
1. Open a file with Sublime text
2. 1st line: **"Nb od sample"** TABS **"2"** *(Always)* TABS **1** *(Always)*
3. Write the phenotype of each of your sample seperated by tabulation
4. Save the file by adding CLS at the en of the name. 



## B) Run GSEA

Open GSEA Software  

Go to "Load data"  
Swipe your 3 files: .gct, Phenotype, and chip in the dropbox  
Click "Load these files"  
Check that the message "Files loaded successfully: 3/3", then OK   

Go to "Run GSEA"  

In **Required fields**  
- Expression dataset: *Load the data of interest*  
- Gene sets database:  *For this project select: **"h.all.v7.5.symbols.gmt [Hallmarks]"** *  
- Number of permutation: *Let it at 1000*  
- Phenotype labels: *Select the comparaison you want*   
$~~~~~~~~~$*We want to compare LAMBDA to the ctrl MEDIUM => Lambda vs Medium*  
- Collapse/Remap to gene symbols : **"No collapse"**  
- Permutation type: **Penotype**  
- Chip platform: *=> Chips (local .chip) => select your chip*  

In **Basic Fields**  
- Min size: exclude smaller sets : **10**  
- Save results in this folder: *Chose your directory*  

In **Advanced fields**:   
Nothing to change  

THEN CLICK **RUN**

## C) Results
###1. Processing the file: 
####a. Lambda vs Medium

```{r}
GSEAL <- read.table("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/Results/ML/_gsea_report_for_Lambda_1648218618087.tsv" , sep="\t", header = TRUE)

GSEAL <- GSEAL[, c(1,2,4,5,6,8,10)]
summary(GSEAL)

#Withdraw all Hallmark with a FDR q-val > 0,25
GSEAL <- GSEAL %>% filter(FDR.q.val<0.25)
GSEAL$NAME <- as.factor(GSEAL$NAME)
nrow(GSEAL)
```

####a. Lambda vs Alpha

```{r}
GSEAL <- read.table("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/Results/AL/_gsea_report_for_Alpha_1649087444476.tsv" , sep="\t", header = TRUE)

GSEAL <- GSEAL[, c(1,4,5,6,8,10)]
summary(GSEAL)

#Withdraw all Hallmark with a FDR q-val > 0,25
GSEAL <- GSEAL %>% filter(FDR.q.val<0.25)
GSEAL$NAME <- as.factor(GSEAL$NAME)
nrow(GSEAL)
```

###2. Plot 
```{r}
ML <- ggplot(data=GSEAL, aes(x = reorder(factor(NAME), NES), y= NES , fill= FDR.q.val)) +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red")+
  geom_bar(width=0.7, stat = "identity")+
  theme(axis.title = element_blank())

ML
```

```{r}
pdf("Lambda_vs_Medium.pdf", width = 10, height = 8)
ML
dev.off()
```



```{r}
AL <- ggplot(data=GSEAL, aes(x = reorder(factor(NAME), NES), y= NES , fill= FDR.q.val)) +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red")+
  geom_bar(width=0.7, stat = "identity")+
  theme(axis.title = element_blank())

AL
```

```{r}
pdf("Lambda_vs_Alpha.pdf", width = 10, height = 8)
AL
dev.off()
```


```{r}
ggdotchart(GSEAL, x = NAME, y = NES,
           color = FDR.FDR.q.val,                                # Color by groups
           palette = c("blue", "red"), # Custom color palette
           sorting = NES,                       # Sort value in descending order
           rotate = TRUE,                                # Rotate vertically
           dot.size = 2,                                 # Large dot size
           y.text.col = TRUE,                            # Color y text by groups
           ggtheme = theme_pubr()                        # ggplot2 theme
           )
```

https://rpkgs.datanovia.com/ggpubr/

```{r}
ggdotchart(GSEAL, x = reorder(factor(NAME), NES), y = NES,
           color = FDR.q.val,                                # Color by groups
           scale_fill_gradient(low = "blue", high = "red"),  # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 6,                                 # Large dot size
           label = round(GSEAL$Enriched_genes,1),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")

```



#II. Difference in nom counts: 

Withdraw all the genes from the "MTORC1_SIGNALING" enriched in Lambda and in Alpha

```{r}
#For Lambda
mtorL <- read.table("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/Results/ML/HALLMARK_MTORC1_SIGNALING.tsv", sep="\t", header = TRUE )
summary(mtorL)
mtorL$CORE.ENRICHMENT <- as.factor(mtorL$CORE.ENRICHMENT)

#Filter "Yes" gene in "CORE ENRICHMENT"
mtorL <- mtorL %>% filter(CORE.ENRICHMENT=="Yes")
nrow(mtorL)

#For alpha
mtorA <- read.table("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/Results/AL/HALLMARK_MTORC1_SIGNALING.tsv", sep="\t", header = TRUE )
summary(mtorA)
mtorL$CORE.ENRICHMENT <- as.factor(mtorL$CORE.ENRICHMENT)

#Filter "Yes" gene in "CORE ENRICHMENT"
mtorA <- mtorA %>% filter(CORE.ENRICHMENT=="Yes")
nrow(mtorA)
```

Compare genes between Lambda and Alpha: 
```{r}
mtorAL <- dplyr::inner_join(mtorL, mtorA, by = c("NAME"))
nrow(mtorAL)
mtorAL$SYMBOL.x
```

Plot these 4 genes norm counts: 

Import the norm counts of 12h: 
```{r}
norm_counts12 <- read.delim("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/1_Processing/normalized_counts12.txt")

norm_counts12 <- norm_counts12[,c(1,8,15,22,29,36,2,9,16,23,30,7,14,21,28,35)]
look(norm_counts12)
```

```{r}
Fourg <- norm_counts12 %>% filter(Gene  == "ENSG00000044574" | Gene  == "ENSG00000183735" | Gene  == "ENSG00000120053"|Gene  == "ENSG00000196262", preserve=TRUE)
look(Fourg)

write_csv2(Fourg, "C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/4genes.csv")
```

Modification du fichier à la main: 


```{r}
Fourg <- read.csv2("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/4genes_1.csv")
Fourg$Stim <- as.factor(Fourg$Stim) 
Fourg$Donor <- as.factor(Fourg$Donor)
Fourg$Gene <- as.factor(Fourg$Gene)
summary(Fourg)
```
```{r}
```



```{r}
ggplot(Fourg, aes(Stim, Counts, color = Gene)) +
  geom_jitter(position = position_jitter(0.2)) + 
  geom_line(aes(group = Gene),data = Fourg) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "pink", "Blue")) +
  scale_y_continuous(trans = 'log10') +
  theme(legend.position = "top")
```


```{r}
Theme_GreyBackground <- theme( axis.title = element_text(size=14,face="bold"), 
                               panel.grid = element_line(size = 0.5, linetype = 'solid', colour = "grey"),
                               panel.border = element_rect(fill = "NA", color="grey", size=1),
                               legend.title = element_text(size = 14),
                               legend.text = element_text(size = 14), 
                               axis.text.x = element_text(face="bold", size=14, colour = "black"), 
                               axis.text.y = element_text(face="bold", size=14))
```



4
```{r}
ggplot(data=Fourg, aes(x = Stim, y= Counts))+
  geom_boxplot(aes(color = Gene), width = 0.5, position = position_dodge(0.8)) + 
  geom_dotplot(aes(fill = Gene, color = Gene), binaxis='y', stackdir='center', dotsize = 0.8,position = position_dodge(0.8))+
  scale_y_continuous(trans = 'log10') +
  scale_fill_manual(values = c("#3288BD", "#F46D43", "pink", "blue"))+
  scale_color_manual(values = c("#3288BD", "#F46D43", "pink", "blue"))+
  Theme_GreyBackground

#Save it - OMG so pretty 
ggsave(filename = "4genes.png", last_plot(), width = 15, height = 8, units = "in", dpi=300)
```




2: 

```{r}
Fourg <- read.csv2("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/2genes_1.csv")
Fourg$Stim <- as.factor(Fourg$Stim) 
Fourg$Donor <- as.factor(Fourg$Donor)
Fourg$Gene <- as.factor(Fourg$Gene)
summary(Fourg)
```

```{r}
ggplot(data=Fourg, aes(x = Stim, y= Counts))+
  geom_boxplot(aes(color = Gene), width = 0.5, position = position_dodge(0.8)) + 
  geom_dotplot(aes(fill = Gene, color = Gene), binaxis='y', stackdir='center', dotsize = 0.8,position = position_dodge(0.8))+
  scale_fill_manual(values = c("#3288BD", "#F46D43", "pink", "blue"))+
  scale_color_manual(values = c("#3288BD", "#F46D43", "pink", "blue"))+
  Theme_GreyBackground

#Save it - OMG so pretty 
ggsave(filename = "2genes1.png", last_plot(), width = 15, height = 8, units = "in", dpi=300)
```



```{r}
Fourg <- read.csv2("C:/Users/CAUX/Documents/Candice/Bioinfo/Microbulk/B_Stim12h_pDCs/3_GSEA/2genes_2.csv")
Fourg$Stim <- as.factor(Fourg$Stim) 
Fourg$Donor <- as.factor(Fourg$Donor)
Fourg$Gene <- as.factor(Fourg$Gene)
summary(Fourg)
```

```{r}
ggplot(data=Fourg, aes(x = Stim, y= Counts))+
  geom_boxplot(aes(color = Gene), width = 0.5, position = position_dodge(0.8)) + 
  geom_dotplot(aes(fill = Gene, color = Gene), binaxis='y', stackdir='center', dotsize = 0.8,position = position_dodge(0.8))+
  scale_fill_manual(values = c("#3288BD", "#F46D43", "pink", "blue"))+
  scale_color_manual(values = c("#3288BD", "#F46D43", "pink", "blue"))+
  Theme_GreyBackground

#Save it - OMG so pretty 
ggsave(filename = "2genes2.png", last_plot(), width = 15, height = 8, units = "in", dpi=300)
```



