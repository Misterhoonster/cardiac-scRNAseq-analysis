# 🫀 Cardiac scRNA Analysis
A collection of scRNA sequencing analyses conducted for the Choudhury lab in summer 2022.

## 📝 Publications
High-quality nuclei isolation from postmortem human heart muscle tissues for single-cell studies
[https://www.sciencedirect.com/science/article/pii/S002228282300055X](url)

## 🔬 Research goals
The broader goal of my research was to uncover novel relationships between factors such as heart disease, aging, and smoking
on gene and transposable element (TE) expression in the human heart. While my initial goal was to focus on aging, my research
expanded to examine the effects of smoking and ploidy on gene expression.

## 🚬 Smoking analysis
My first project examined the link between smoking and gene expression in the heart. Of the samples collected by the lab,
there were only middle-aged and old female samples for smoking, so I controlled for this by drawing control samples from
the same demographic. However, since there were only two control samples that met these criteria, I sourced additional data
from a Nature paper (Koenig et al., 2022) examining cell-type-diversification in the human heart. Using the Seurat and DESeq2
packages in R, I generated a list of differentially expressed genes (DEGs) in the smoking population and found a significant
difference in the distribution of cell types between the smoking and control samples.<br></br>A concern I had while running the differential
expression analyses was that differences in gene expression could be due to variations in preparation between the Nature and Choudhury
samples rather than genuine biological factors. This was supported by the fact that the samples were separated by source
(Choudhury lab vs. from the Nature paper) rather than age or smoking status on the PCA plot. Because of this, I re-ran the analysis using
solely the data from my lab. I then compiled a list of DEGs found in the Choudhury-only analysis and the combined analysis.<br></br>This list
was composed of just two genes: RGCC, a cell-cycle regulator that was down-regulated in smoking samples, and FABP4, a fatty-acid binding
protein which was up-regulated in smoking samples. Previous studies have suggested that down-regulation of RGCC could play a role in the
pathogenesis of non-small-cell lung cancer (Kim et al., 2011) while elevated levels of FABP4 have been linked to cardiovascular disease
and death (Furuhashi et al., 2015; Saito et al., 2021).
<br></br>

<img src="https://raw.githubusercontent.com/Misterhoonster/cardiac-scRNAseq-analysis/main/images/smoking/rgcc_fabp4_vln_smoking_cum.png" width="400px"></img>
<br />
*RGCC and FABP4 were markers found to be differentially expressed in smoking individuals*
<br></br>

<img src="https://github.com/Misterhoonster/cardiac-scRNAseq-analysis/blob/main/images/smoking/rgcc_fabp4_feature.png?raw=true" width="400px"></img>
<br />
*Expression patterns of RGCC and FABP4*

## 🧫 Ploidy analysis
My second project explored the relationship between ploidy (diploid v. hexaploid) and gene expression with the goal of finding potential
fusion markers for cardiomyocytes. Out of our samples, there was only one individual that we had collected diploid and hexaploid samples from,
so I began by running differential expression analysis on that data. Next, I pooled together all the diploid and hexaploid samples
(not necessarily from the same individuals) and ran another round of analysis. I was worried that confounding factors from this pooled analysis
such as age and smoking status would interfere with the results, so I accounted for this by selecting only the differentially expressed genes
found in both the single individual and pooled runs. As a result of this analysis, I compiled a list of DEGs which I provided to a post-doc in
my lab as possible candidate markers.
<br></br>

<img src="https://github.com/Misterhoonster/cardiac-scRNAseq-analysis/blob/main/images/ploidy/2n_all_enriched_pathways.png?raw=true" width="400px"></img>
<br />
*Pathways enriched in diploid (2n) cells*
<br></br>

<img src="https://github.com/Misterhoonster/cardiac-scRNAseq-analysis/blob/main/images/ploidy/6n_all_enriched_pathways.png?raw=true" width="400px"></img>
<br />
*Pathways enriched in hexaploid (6n) cells*
<br></br>

## 🧬 TE analysis
My last project examined the effects of smoking and aging on the expression of transposable elements (TEs). Through this analysis, I found that
the distribution of TE types (LINE, SINE, DNA, Satellite…etc) does not appear to be significantly impacted by age or smoking. In addition, I
compiled a list of DEGs that were up and down-regulated in smokers and up-regulated in each age group (young, middle, old).
<br></br>

<img src="https://user-images.githubusercontent.com/11368469/183781035-0db6eec6-0add-43d9-a6f8-b240dc04a055.png" width="400px"></img>
<br />
*TEs up-regulated in young individuals*
<br></br>

<img src="https://user-images.githubusercontent.com/11368469/183781069-170cc46d-982d-4c99-9f0c-4ecbbe75446a.png" width="400px"></img>
<br />
*TEs up-regulated in middle-aged individuals*
<br></br>

<img src="https://user-images.githubusercontent.com/11368469/183781095-28435893-baf5-4830-b768-76420ed5c2df.png" width="400px"></img>
<br />
*TEs up-regulated in old individuals*
<br></br>

## Citations
Furuhashi, M., Saitoh, S., Shimamoto, K., & Miura, T. (2015). Fatty Acid-Binding Protein 4 (FABP4): Pathophysiological Insights and Potent Clinical     Biomarker of Metabolic and Cardiovascular Diseases. Clinical Medicine Insights. Cardiology, 8(Suppl 3), 23–33. https://doi.org/10.4137/CMC.S17067

Kim, D. S., Lee, J. Y., Lee, S. M., Choi, J. E., Cho, S., & Park, J. Y. (2011). Promoter methylation of the RGC32 gene in nonsmall cell lung cancer. Cancer, 117(3), 590–596. https://doi.org/10.1002/cncr.25451

Koenig, A. L., Shchukina, I., Amrute, J., Andhey, P. S., Zaitsev, K., Lai, L., Bajpai, G., Bredemeyer, A., Smith, G., Jones, C., Terrebonne, E., Rentschler, S. L., Artyomov, M. N., & Lavine, K. J. (2022). Single-cell transcriptomics reveals cell-type-specific diversification in human heart failure. Nature Cardiovascular Research, 1(3), 263–280. https://doi.org/10.1038/s44161-022-00028-6

Saito, N., Furuhashi, M., Koyama, M., Higashiura, Y., Akasaka, H., Tanaka, M., Moniwa, N., Ohnishi, H., Saitoh, S., Ura, N., Shimamoto, K., & Miura, T. (2021). Elevated circulating FABP4 concentration predicts cardiovascular death in a general population: A 12-year prospective study. Scientific Reports, 11(1), 4008. https://doi.org/10.1038/s41598-021-83494-5
