# TaphonomeAnalyst 2.0

## Introduction

### 1. TaphonomeAnalyst 1.0 [_Initial version_]

Lacustrine taphocoenosis network analysis software.  

Related Article: A new method for examining the co-occurrence network of fossil assemblages  
Publication: https://doi.org/10.1038/s42003-023-05417-6  
Contract: 13820113071@163.com (Shilong Guo) / wma19952022@163.com (Wang Ma)

### 2. TaphonomeAnalyst 2.0 [_*current version_]

Integrated analysis software for lacustrine taphocoenosis network and geochemical data.

Article:   
Publication:   
Contract: 13820113071@163.com (Shilong Guo) / wma19952022@163.com (Wang Ma)

### 3. Version Comparison

| Module                                                                                         | Function                                                                                                                                                                                   | Command                    | 1.0              | 2.0              |
|------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------|------------------|------------------|
| Assessment of sampling effort and estimation of theoretical maximum biodiversity **(ModuleⅠ)** | Sobs<br>Chao1<br>ACE                                                                                                                                                                       | m1sobs<br>m1chao<br>m1ace | √<br>√<br>√      | √<br>√<br>√      |
| Relative abundance of OTU analysis **(Module Ⅱ)**                                              |                                                                                                                                                                                            | m2abundplots                 | √                | √                |
| Proportion of taphonomic preservational grade of species analysis **(Module Ⅲ)**               |                                                                                                                                                                                            | m3otus<br>m3plots          | √                | √                |
| Taphonomic environment analysis **(Module Ⅳ)**                                                 | Including the creation of Veen diagrams that compare the diversity found across sedimentary environments or outcrops.                                                                      | m4cluster<br>m4venn      | √                | √                |
| Visualization of geochemical data **(Module Ⅴ)**                                               |                                                                                                                                                                                            | m5cluster                 | ×                | √                |
| Assembling dissimilarity-environmental distance test **(Module Ⅵ)**                            |                                                                                                                                                                                            | m6dissenvtest                | ×                | √                |
| Mantel Test between species abundance and ecological environmental variables **(Module Ⅶ)**    |                                                                                                                                                                                            | m7mantel                     | ×                | √                |
| Species correlation semi-matrix graphics **(Module Ⅷ)**                                        |                                                                                                                                                                                            | m8corrotus                   | √                | √                |
| Correlational Network Visualization **(Module Ⅸ)**                                             | SparCC<br>Pearson<br>Spearman<br>Kendall                                                                                                                                                   | m9cooccurnet                 | ×<br>√<br>√<br>√ | √<br>√<br>√<br>√ |
| Comparison of networks **(Module Ⅹ)**                                                          | The capability to compare networks under different groups of plots. Visualization of total nodes, total linked nodes, total edges, density, modularity, complexity, degree and robustness. | m10netVC                      | ×                | √                |

## Development Environments

This project requires the following software versions:

- **Python**: 3.9.21
- **R**: 4.3.3

### Dependencies

#### 1. Python Packages

This project uses the following Python packages:

- `rpy2`: 3.5.16
- `community`: 1.0.0b1
- `dask`: 2023.5.0
- `h5py`: 3.6.0
- `matplotlib`: 3.4.1
- `networkx`: 2.5.1
- `numba`: 0.58.1
- `numpy`: 1.23.0
- `pandas`: 1.5.0
- `python_igraph`: 0.11.6
- `python_louvain`: 0.15
- `scikit_bio`: 0.6.2
- `scipy`: 1.9.0
- `seaborn`: 0.11.1
- `venn`: 0.1.3

#### 2. R Packages

This project uses the following R packages:

- `dplyr`: 1.1.4
- `ggplot2`: 3.4.4
- `linkET`: 0.0.7.4

## Installation

There are two ways to install TaphonomeAnalyst 2.0.

### 1. Use Conda (Recommend)

#### * Note : How to download and install conda? [Documentation](https://docs.conda.io/projects/conda/en/stable/).

Command:

    conda create -n taphonomeAnalyst2 python=3.9 r-base=4.3.3
    conda activate taphonomeAnalyst2
    conda install r-dplyr r-ggplot2 r-devtools r-vegan openpyxl
    pip install -r ./requirements.txt
    Rscript ./install_packages.R


### 2. Only use pip and Rscript

Command:

    pip install -r ./requirements.txt
    Rscript ./install_packages.R

#### * Note : In this way, you need to install Python 3.8.8 and R 4.3.2 in advance and add them to the environment variables.

## Documentation

Using `python ./TaphonomeAnalyst2.py -h` and `python ./TaphonomeAnalyst2.py [command] -h` to read the brief.  

Or check the detail in [User Guide](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Supplementary%20material1%20The%20User%20Guide%20of%20TaphonomeAnalyst%202.0.pdf)

## Information and Data

### 1. User Guide

[Supplementary material1 The User Guide of TaphonomeAnalyst 2.0.pdf](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Supplementary%20material1%20The%20User%20Guide%20of%20TaphonomeAnalyst%202.0.pdf) 

### 2. OTUs Data

[Supplementary material2.xlsx](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Supplementary%20material2.xlsx)  

### 3. Geochemical Data

[Supplementary material3.xlsx](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Supplementary%20material3.xlsx)  

## Usage and Result

### 1. Hierarchical clustering-sedimentary environment. **(Module Ⅳ)** [clustermap]

Command:
        
    python ./TaphonomeAnalyst2.py m4cluster --input "./Supplementary material2.xlsx" --aquatic "Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l)" --level species

Result:

![clusterenv_aquatic.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/clusterenv_aquatic.png)

![clusterenv_tree.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/clusterenv_tree.png)

### 2. Visualization of geochemical data. **(Module Ⅴ)** [clustermap]

Command:
        
    python ./TaphonomeAnalyst2.py m5cluster --input "./Supplementary material2.xlsx" --aquatic "Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l)" --geochem "./Supplementary material3.xlsx" --level species

Result:

![clusterenv_geochem.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/clusterenv_geochem.png)

### 3. Venn diagram-sampling locations or environments. **(Module Ⅳ)** [venn]

Command:
        
    python ./TaphonomeAnalyst2.py m4venn --input "./Supplementary material2.xlsx" --groups plot3-3/plot3-1/plot3-2,plot1-1/plot2-2/plot2-1/plot2-3/plot1-2/plot1-3 --level family

Result:

![divvenn.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/divvenn.png)

### 4. Taphonomic grades-taxa. **(Module Ⅲ)** [barh]

Command:

    python ./TaphonomeAnalyst2.py m3otus --input "./Supplementary material2.xlsx" --level order

Result:

![TGotus.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/TGotus.png)

### 5. Taphonomic grades-sampling plots (in customized order). **(Module Ⅲ)** [barh]

Command:

    python ./TaphonomeAnalyst2.py m3plots --input "./Supplementary material2.xlsx"

    python ./TaphonomeAnalyst2.py m3plots --input "./Supplementary material2.xlsx" --groups plot3-3/plot3-1/plot3-2,plot1-1/plot2-2/plot2-1/plot2-3/plot1-2/plot1-3

Result:

![TGplots.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/TGplots.png)

![TGplots_groups.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/TGplots_groups.png)

### 6. Abundance-sampling plots. **(Module Ⅱ)** [barh]

Command:

    python ./TaphonomeAnalyst2.py m2abundplots --input "./Supplementary material2.xlsx" --level order

Result:

![abundplots.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/abundplots.png)

### 7. Co-occurrence networks. **(Module Ⅸ)** [network/venn]

Command:

    python ./TaphonomeAnalyst2.py m9cooccurnet --input "./Supplementary material2.xlsx" --level family --corr pearson --corr_coef 0.7 --p_value 0.1

Result:

![cooccurnet.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/cooccurnet.png)

![cooccurnet_venn.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/cooccurnet_venn.png)

### 8. Sampling coverage curve. **(Module Ⅰ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py m1sobs --input "./Supplementary material2.xlsx" --level family --groups plot1:plot1-1/plot1-2/plot1-3,plot2:plot2-1/plot2-2/plot2-3,plot3:plot3-1/plot3-2/plot3-3

Result:

![samplecurve.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/samplecurve.png)

### 9. Chao1 potential diversity curve. **(Module Ⅰ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py m1chao --input "./Supplementary material2.xlsx" --level family --groups plot1:plot1-1/plot1-2/plot1-3,plot2:plot2-1/plot2-2/plot2-3,plot3:plot3-1/plot3-2/plot3-3

Result:

![chao.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/chao.png)

### 10. ACE potential diversity curve. **(Module Ⅰ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py m1ace --input "./Supplementary material2.xlsx" --level family --groups plot1:plot1-1/plot1-2/plot1-3,plot2:plot2-1/plot2-2/plot2-3,plot3:plot3-1/plot3-2/plot3-3 --rare 10

Result:

![ace.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/ace.png)

### 11. Heatmap-OTUs correlation analysis. **(Module Ⅷ)** [heatmap]

Command:

    python ./TaphonomeAnalyst2.py m8corrotus --input "./Supplementary material2.xlsx" --level family --corr pearson --p_value 0.1

Result:

![corrotus.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/corrotus.png)

### 12. Mantel Test between species abundance and ecological environmental variables. **(Module Ⅶ)** [multiplot]

Command:

    python ./TaphonomeAnalyst2.py m7mantel --input "./Supplementary material2.xlsx" --rhome "C:\Program Files\R\R-4.3.3" --geochem "./Supplementary material3.xlsx" --aquatic Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l) --level_aquatic species --level_terrestrial family --corr pearson

Result:

![mantel.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/mantel.png)

### 13. Assembling dissimilarity- environmental distance test. **(Module Ⅵ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py m6dissenvtest --input "./Supplementary material2.xlsx" --aquatic Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l) --level_aquatic species --level_terrestrial family --geochem "./Supplementary material3.xlsx"

Result:

![dissenvtest.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/dissenvtest.png)

### 14. Generate a unified layout network for comparison. (Module Ⅹ)\t[network/boxplot/barplot]

Command:

    python ./TaphonomeAnalyst2.py m10netVC --input "./Supplementary material2.xlsx" --aquatic Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l) --level_aquatic species --level_terrestrial family --groups plot3-3/plot3-1/plot3-2,plot1-1/plot2-2/plot2-1/plot2-3/plot1-2/plot1-3 --corr spearman --corr_coef 0.7 --p_value 0.01

Result:

![netVC_Aquatic_A.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Aquatic_A.png)

![netVC_Aquatic_B.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Aquatic_B.png)

![netVC_Aquatic_Robustness.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Aquatic_Robustness.png)

![netVC_Terrestrial_A.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Terrestrial_A.png)

![netVC_Terrestrial_B.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Terrestrial_B.png)

![netVC_Terrestrial_Robustness.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Terrestrial_Robustness.png)

![netVC_Total_nodes.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Total_nodes.png)

![netVC_Total_links.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Total_links.png)

![netVC_Total_edges.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Total_edges.png)

![netVC_Degree.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Degree.png)

![netVC_Density.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Density.png)

![netVC_Complexity.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Complexity.png)

![netVC_Modularity.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/netVC_Modularity.png)


[//]: # (## Citation)
