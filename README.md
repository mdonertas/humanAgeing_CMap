# humanAgeing_CMap

## GTEx Ageing Data

### Data Preparation

Data is preprocessed following the steps explained in my [GTEx](https://github.com/mdonertas/GTEx) repository.

Among all tissues, only the ones having at least 20 subjects are considered. We also excluded 'Cells-Transformedfibroblasts' category. As a result 35 tissues (17 major tissue type) are used for the downstream analysis:

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

Files are saved as .rds data objects under './data/processed/GTEx/' folder.

Exploratory data analysis results are saved under './results/GTEx/' folder.

### CMap analysis

#### All Tissues:

There was only 1 gene that has a consistent gene expression change with age: [ENSG00000269834 (ZNF528 antisense RNA 1)](http://mar2017.archive.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=ENSG00000269834;r=19:52388842-52397766)

Since there is also not enough number of genes of which expression values are significantly<sup>\*</sup> changing with age, we resorted to another approach:
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

Since the resulting numbers are quite low, I further checked the expression patterns of the genes determined by brain microarray studies. The resulting heatmap is './results/GTEx/Brain_microarrayUpsandDowns_inGTEx.pdf'

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

##### Footnotes:
<sup>\*</sup> FDR adjusted p value < 0.05.

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
