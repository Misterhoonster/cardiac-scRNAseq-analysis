# Description
A collection of the scRNA sequencing analyses I ran for the Choudhury lab in summer 2022.

## Research goals
The broader goal of my research was to uncover novel relationships between factors such as heart disease, aging, and smoking
on gene and transposable element (TE) expression in the human heart. While my initial goal was to focus on aging, my research
expanded to examine the effects of smoking and ploidy on gene expression.

![image](https://user-images.githubusercontent.com/11368469/183768859-d2987423-450b-41ac-a74f-66472da51f64.png)

## Smoking analysis
My first project examined the link between smoking and gene expression in the heart. Of the samples collected by the lab,
there were only middle-aged and old female samples for smoking, so I controlled for this by drawing control samples from
the same demographic. However, since there were only two control samples that met these criteria, I sourced additional data
from a Nature paper (Koenig et al., 2022) examining cell-type-diversification in the human heart. Using the Seurat and DESeq2
packages in R, I generated a list of differentially expressed genes (DEGs) in the smoking population and found a significant
difference in the distribution of cell types between the smoking and control samples. A concern I had while running the differential
expression analyses was that differences in gene expression could be due to variations in preparation between the Nature and Choudhury
samples rather than genuine biological factors. This was supported by the fact that the samples were separated by source
(Choudhury lab vs. from the Nature paper) rather than age or smoking status on the PCA plot. Because of this, I re-ran the analysis using
solely the data from my lab. I then compiled a list of DEGs found in the Choudhury-only analysis and the combined analysis. This list
was composed of just two genes: RGCC, a cell-cycle regulator that was down-regulated in smoking samples, and FABP4, a fatty-acid binding
protein which was up-regulated in smoking samples. Previous studies have suggested that down-regulation of RGCC could play a role in the
pathogenesis of non-small-cell lung cancer (Kim et al., 2011) while elevated levels of FABP4 have been linked to cardiovascular disease
and death (Furuhashi et al., 2015; Saito et al., 2021).

![image](https://user-images.githubusercontent.com/11368469/183768958-672bc3eb-781b-4941-8081-066e9edaa83c.png)

## Ploidy analysis
My second project explored the relationship between ploidy (diploid v. hexaploid) and gene expression with the goal of finding potential
fusion markers for cardiomyocytes. Out of our samples, there was only one individual that we had collected diploid and hexaploid samples from,
so I began by running differential expression analysis on that data. Next, I pooled together all the diploid and hexaploid samples
(not necessarily from the same individuals) and ran another round of analysis. I was worried that confounding factors from this pooled analysis
such as age and smoking status would interfere with the results, so I accounted for this by selecting only the differentially expressed genes
found in both the single individual and pooled runs. As a result of this analysis, I compiled a list of DEGs which I provided to a post-doc in
my lab as possible candidate markers.
![image](https://user-images.githubusercontent.com/11368469/183769099-f2200c17-66e9-4c96-b092-ecd2c0fd4d20.png)

## TE analysis
My last project examined the effects of smoking and aging on the expression of transposable elements (TEs). Through this analysis, I found that
the distribution of TE types (LINE, SINE, DNA, Satellite…etc) does not appear to be significantly impacted by age or smoking. In addition, I
compiled a list of DEGs that were up and down-regulated in smokers and up-regulated in each age group (young, middle, old).
![image](https://user-images.githubusercontent.com/11368469/183769170-510bd658-3654-4647-a2f8-4a011503ebd8.png)

