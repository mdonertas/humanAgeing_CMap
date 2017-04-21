# humanAgeing_CMap

## GTEx Ageing Data

### Data Preparation

Data is preprocessed following the steps explained in my [GTEx](https://github.com/mdonertas/GTEx) repository.

Among all tissues, only the ones having at least 20 subjects are considered. We also excluded 'Cells-Transformedfibroblasts' category. As a result 35 tissues<sup>\*1</sup> (17 major tissue type) are used for the downstream analysis.

Files are saved as .rds data objects under './data/processed/GTEx/' folder.

Exploratory data analysis results are saved under './results/GTEx/' folder.

### CMap analysis

#### All Tissues:

There was only 1 gene that has a consistent gene expression change with age: [ENSG00000269834 (ZNF528 antisense RNA 1)](http://mar2017.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000269834;r=19:52388842-52397766)

Since there is also not enough number of genes of which expression values are significantly<sup>\*2</sup> changing with age, we resorted to another approach:
* Genes that are not expressed in all tissues are excluded from the downstream analysis (19064 genes left).
* For each gene, we calculated the number of datasets having a positive correlation with age, irrespective of the effect size.
* Using that distribution we considered a gene as a consistently changing gene if it is in the lower 0.5% (Decreasing expression according to 31 datasets) or upper 0.5% (Increasing expression according to 32 datasets).
* As a result we have 112 up and 104 down regulated genes.

Next, we asked 'How many of these up and down regulated genes are the ones we determined using Brain ageing microarray data?'

||GTEx-up|GTEx-down|Brain_micro-up|Brain_micro-down|
|-|-|-|-|-|
|GTEx-up|X|0|4%|0|
|GTEx-down|0|X|0|3.4%|
|Brain_micro-up|3.6%|0|X|0|
|Brain_micro-down|0|3.8%|0|X|

Since the resulting numbers are quite low, I further checked the expression patterns of the genes determined by brain microarray studies. The resulting heatmap is './results/GTEx/Brain_microarrayUpsandDowns_inGTEx.pdf'. It seems like GTEx brain data shows similar pattern as the microarray but the other tissues differ.

##### CMap Query:

* The list of up and down regulated gene id s are converted to *affymetrix hg u133a* ids, using *getBM* function in **biomaRt**<sup>1</sup> R<sup>2</sup> package. This step is required by the CMap web service.
* If there are some probeset ids shared between up and down regulated genes, these are discarded.
* The list of up and down regulated genes are uploaded to the query system in [CMap website](https://portals.broadinstitute.org/cmap/) with the name 'gtex'.
* The result is downloaded as '.xlsx' file and converted to '.csv' files.
* The '.csv' files are read into R in order to correct the p values for multiple testing.

|rank|cmap.name|mean|n|enrichment|p|specificity|percent.non.null|padj|
|----|---------|----|-|----------|-|-----------|----------------|----|
|1|sirolimus |  0.412 | 44   |  0.361       | 0    | 0.1928            |  72 | 0.000000
|2|LY-294002 |  0.372 | 61   |  0.355       | 0    | 0.2282            |  65 | 0.000000
|3|isoxicam | -0.633 |  5   | -0.848 | 0.00024         | 0            | 100 | 0.039088
|4|quinostatin |  0.816 |  2   |  0.987 | 0.00028    | 0.0111            | 100 | 0.039088
|5|cephaeline |  0.725 |  5   |  0.833 | 0.00028    | 0.0959            | 100 | 0.039088

#### All Tissues except Brain:
*22 tissues*

There were only 4 genes that has a consistent gene expression change with age:
* [ENSG00000269834 (ZNF528 antisense RNA 1)](http://mar2017.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000269834;r=19:52388842-52397766)
* [ENSG00000134574 (damage specific DNA binding protein 2)](http://mar2017.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000134574;r=11:47214465-47239240;t=ENST00000256996)
* [ENSG00000170160 (coiled-coil domain containing 144A)](http://mar2017.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000170160;r=17:16689537-16777881)
* [ENSG00000162695 (solute carrier family 30 member 7)](http://mar2017.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000162695;r=1:100896076-100981753)

Since there is also not enough number of genes of which expression values are significantly<sup>\*</sup> changing with age, we resorted to another approach:
* Genes that are not expressed in all tissues are excluded from the downstream analysis (19456 genes left).
* For each gene, we calculated the number of datasets having a positive correlation with age, irrespective of the effect size.
* Using that distribution we considered a gene as a consistently changing gene if it is in the lower 0.5% (Decreasing expression according to 19 datasets) or upper 0.5% (Increasing expression according to 20 datasets).
* As a result we have 205 up and 254 down regulated genes.

Next, we asked 'How many of these up and down regulated genes are the ones we determined using Brain ageing microarray data?'

||GTEx-up|GTEx-down|Brain_micro-up|Brain_micro-down|
|-|-|-|-|-|
|GTEx-up|X|0|2%|0.9%|
|GTEx-down|0|X|1%|3.4%|
|Brain_micro-up|1%|0.4%|X|0|
|Brain_micro-down|0.5%|1.6%|0|X|

##### CMap Query:

* The list of up and down regulated gene id s are converted to *affymetrix hg u133a* ids, using *getBM* function in **biomaRt**<sup>1</sup> R<sup>2</sup> package. This step is required by the CMap web service.
* If there are some probeset ids shared between up and down regulated genes, these are discarded.
* The list of up and down regulated genes are uploaded to the query system in [CMap website](https://portals.broadinstitute.org/cmap/) with the name 'gtex'.
* The result is downloaded as '.xlsx' file and converted to '.csv' files.
* The '.csv' files are read into R in order to correct the p values for multiple testing.

|rank|cmap.name|mean|n|enrichment|p|specificity|percent.non.null|padj|
|----|---------|----|-|----------|-|-----------|----------------|----|
|  1          | anisomycin |  0.771 |  4   |  0.975       | 0    | 0.0155            | 100 | 0.000000000
|  2           | puromycin |  0.701 |  4   |  0.947       | 0    | 0.0449            | 100 | 0.000000000
|  3          | cephaeline |  0.806 |  5   |  0.944       | 0     | 0.048            | 100 | 0.000000000
|  4        | thioridazine |  0.531 | 20   |  0.637       | 0    | 0.0776            |  80 | 0.000000000
|  5        | tanespimycin |  0.241 | 62   |  0.289       | 0    | 0.3938            |  53 | 0.000000000
|  6           | LY-294002 |  0.262 | 61   |  0.285 | 0.00002    | 0.3826            |  52 | 0.001326667
|  7         | terfenadine |  0.738 |  3   |  0.979 | 0.00004    | 0.0049            | 100 | 0.001326667
|  8             | emetine |  0.746 |  4   |  0.921 | 0.00004    | 0.0211            | 100 | 0.001326667
|  9         | niclosamide |  0.567 |  5   |  0.921 | 0.00004    | 0.0105            | 100 | 0.001326667
| 10       | cicloheximide |  0.714 |  4   |  0.918 | 0.00004    | 0.0339            | 100 | 0.001326667
| 11    | prochlorperazine |  0.431 | 16   |  0.572 | 0.00004    | 0.0631            |  75 | 0.001326667
| 12     | trifluoperazine |  0.414 | 16   |  0.569 | 0.00004     | 0.125            |  75 | 0.001326667
| 13          | loperamide |  0.497 |  6   |  0.791 | 0.00018    | 0.0101            | 100 | 0.005510769
| 14          | indoprofen | -0.554 |  4   | -0.897 | 0.00022         | 0            | 100 | 0.006254286
| 15           | sirolimus |  0.266 | 44   |  0.309 | 0.00028    | 0.3072            |  54 | 0.007429333
| 16        | alvespimycin |  0.382 | 12   |  0.572 | 0.00032    | 0.0402            |  75 | 0.007960000
| 17        | lanatoside C |  0.503 |  6   |  0.743 | 0.00077    | 0.1009            |  83 | 0.018027059
| 18        | quinisocaine |  0.541 |  4   |  0.845 | 0.00086         | 0            | 100 | 0.019015556
| 19            | thiamine | -0.537 |  3   | -0.916 | 0.00104         | 0            | 100 | 0.021293000
| 20       | digitoxigenin |  0.535 |  4   |  0.838 | 0.00107    | 0.0429            | 100 | 0.021293000
| 21           | rottlerin |  0.566 |  3   |  0.917 | 0.00118    | 0.0518            | 100 | 0.021890000
| 22    | phenoxybenzamine |  0.483 |  4   |  0.832 | 0.00121    | 0.2525            | 100 | 0.021890000
| 23         | PNU-0251126 | -0.413 |  6   | -0.706 | 0.00157    | 0.0133            |  83 | 0.027167826
| 24      | nicotinic acid | -0.419 |  4   | -0.827 | 0.00165         | 0            |  75 | 0.027362500
| 25             | 5224221 |  0.619 |  2   |  0.961 | 0.00264    | 0.1341            | 100 | 0.039800000
| 26       | calmidazolium |  0.618 |  2   |  0.959 | 0.00276    | 0.0474            | 100 | 0.039800000
| 27         | ellipticine | -0.312 |  4   | -0.804 | 0.00282    | 0.0382            |  50 | 0.039800000
| 28                 | H-7 | -0.271 |  4   | -0.804 | 0.00284    | 0.2273            |  50 | 0.039800000
| 29        | DL-thiorphan |  0.596 |  2   |  0.959 | 0.00296         | 0            | 100 | 0.039800000
| 30             | 5255229 |  0.590 |  2   |  0.958 | 0.00308         | 0            | 100 | 0.039800000
| 31 | tetrahydroalstonine | -0.323 |  4   | -0.799 | 0.00318         | 0            |  75 | 0.039800000
| 32        | blebbistatin |  0.637 |  2   |  0.957  | 0.0032    | 0.0122            | 100 | 0.039800000
| 33           | cefotetan | -0.435 |  3   | -0.878 | 0.00357    | 0.0063            | 100 | 0.043056364

###### Footnotes:

<sup>\*1</sup> Tissues analysed in this analysis

[1] "Adipose-Subcutaneous"               
[2] "Adipose-Visceral-Omentum"           
[3] "Artery-Aorta"                       
[4] "Artery-Tibial"                      
[5] "Brain-Amygdala"                     
[6] "Brain-Anteriorcingulatecortex-BA24"
[7] "Brain-Caudate-basalganglia"         
[8] "Brain-CerebellarHemisphere"         
[9] "Brain-Cerebellum"                   
[10] "Brain-Cortex"                       
[11] "Brain-FrontalCortex-BA9"            
[12] "Brain-Hippocampus"                  
[13] "Brain-Hypothalamus"                 
[14] "Brain-Nucleusaccumbens-basalganglia"
[15] "Brain-Putamen-basalganglia"         
[16] "Brain-Spinalcord-cervicalc-1"       
[17] "Brain-Substantianigra"              
[18] "Breast-MammaryTissue_male"          
[19] "Colon-Sigmoid"                      
[20] "Esophagus-GastroesophagealJunction"
[21] "Esophagus-Mucosa"                   
[22] "Esophagus-Muscularis"               
[23] "Heart-AtrialAppendage"              
[24] "Heart-LeftVentricle"                
[25] "Liver"                              
[26] "Lung"                               
[27] "Muscle-Skeletal"                    
[28] "Nerve-Tibial"                       
[29] "Pituitary"                          
[30] "Prostate"                           
[31] "Skin-NotSunExposed-Suprapubic"      
[32] "Skin-SunExposed-Lowerleg"           
[33] "Testis"                             
[34] "Thyroid"                            
[35] "WholeBlood"

<sup>\*2</sup> FDR adjusted p value < 0.05.

| Tissue        | # of Genes|
| ------------- |:-------------:|
| Adipose-Subcutaneous |                   0 |
| Adipose-Visceral-Omentum |               0 |
| Artery-Aorta |                          60 |
| Artery-Tibial |                        271 |
| Brain-Amygdala |                        24 |
| Brain-Anteriorcingulatecortex-BA24 |     6 |
| Brain-Caudate-basalganglia |             0 |
| Brain-CerebellarHemisphere |           608 |
| Brain-Cerebellum |                       0 |
| Brain-Cortex |                        1260 |
| Brain-FrontalCortex-BA9 |               81 |
| Brain-Hippocampus |                   3892 |
| Brain-Hypothalamus |                   478 |
| Brain-Nucleusaccumbens-basalganglia | 1414 |
| Brain-Putamen-basalganglia |             1 |
| Brain-Spinalcord-cervicalc-1 |           1 |
| Brain-Substantianigra |                  0 |
| Breast-MammaryTissue_male |              0 |
| Colon-Sigmoid |                          0 |
| Esophagus-GastroesophagealJunction |     0 |
| Esophagus-Mucosa |                       0 |
| Esophagus-Muscularis |                   0 |
| Heart-AtrialAppendage |                  0 |
| Heart-LeftVentricle |                    1 |
| Liver |                                  0 |
| Lung |                                   0 |
| Muscle-Skeletal |                        1 |
| Nerve-Tibial |                           1 |
| Pituitary |                              1 |
| Prostate |                               0 |
| Skin-NotSunExposed-Suprapubic |          0 |
| Skin-SunExposed-Lowerleg |               0 |
| Testis |                                 0 |
| Thyroid |                                0 |
| WholeBlood |                             0 |

###### References:
<sup>1</sup> Durinck, S., Moreau, Y., Kasprzyk, A., Davis, S., De Moor, B., Brazma, A., & Huber, W. (2005). BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis. Bioinformatics (Oxford, England), 21(16), 3439–3440. https://doi.org/10.1093/bioinformatics/bti525

<sup>2</sup> R Core Team (2016). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL
  https://www.R-project.org/.
