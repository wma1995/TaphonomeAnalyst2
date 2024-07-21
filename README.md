# TaphonomeAnalyst 2.0

## Introduction

### 1. TaphonomeAnalyst 1.0 [_Initial version_]

Lacustrine taphocoenosis network analysis software.  

Related Article: A new method for examining the co-occurrence network of fossil assemblages  
Publication: https://doi.org/10.1038/s42003-023-05417-6  
Contract: 13820113071@163.com (shilong Guo) / wma19952022@163.com (wang Ma)

### 2. TaphonomeAnalyst 2.0 [_*current version_]

Integrated analysis software for lacustrine taphocoenosis network and geochemical data.

Article:   
Publication:   
Contract: 13820113071@163.com (shilong Guo) / wma19952022@163.com (wang Ma)

### 3. Version Comparison

| Module                                                                                         | Function                                                                                                                                                                                   | Command                    | 1.0              | 2.0              |
|------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------|------------------|------------------|
| Assessment of sampling effort and estimation of theoretical maximum biodiversity **(ModuleⅠ)** | Sobs<br>Chao1<br>ACE                                                                                                                                                                       | samplecurve<br>chao<br>ace | √<br>√<br>√      | √<br>√<br>√      |
| Relative abundance of OTU analysis **(Module Ⅱ)**                                              |                                                                                                                                                                                            | abundplots                 | √                | √                |
| Proportion of taphonomic preservational grade of species analysis **(Module Ⅲ)**               |                                                                                                                                                                                            | TGotus / TGplots           | √                | √                |
| Taphonomic environment analysis **(Module Ⅳ)**                                                 | Including the creation of Veen diagrams that compare the diversity found across sedimentary environments or outcrops.                                                                      | clusterenv / divvenn       | √                | √                |
| Visualization of geochemical data **(Module Ⅴ)**                                               |                                                                                                                                                                                            | clusterenv                 | ×                | √                |
| Assembling dissimilarity-environmental distance test **(Module Ⅵ)**                            |                                                                                                                                                                                            | dissenvtest                | ×                | √                |
| Mantel Test between species abundance and ecological environmental variables **(Module Ⅶ)**    |                                                                                                                                                                                            | mantel                     | ×                | √                |
| Species correlation semi-matrix graphics **(Module Ⅷ)**                                        |                                                                                                                                                                                            | corrotus                   | √                | √                |
| Correlational Network Visualization **(Module Ⅸ)**                                             | SparCC<br>Pearson<br>Spearman<br>Kendall                                                                                                                                                   | cooccurnet                 | ×<br>√<br>√<br>√ | √<br>√<br>√<br>√ |
| Comparison of networks **(Module Ⅹ)**                                                          | The capability to compare networks under different groups of plots. Visualization of total nodes, total linked nodes, total edges, density, modularity, complexity, degree and robustness. | netVC                      | ×                | √                |

## Development Environments

This project requires the following software versions:

- **Python**: 3.8.8
- **R**: 4.3.2

### Dependencies

#### 1. Python Packages

This project uses the following Python packages:

- `community`: 1.0.0b1
- `dask`: 2023.5.0
- `h5py`: 3.6.0
- `matplotlib`: 3.4.1
- `networkx`: 2.5.1
- `numba`: 0.58.1
- `numpy`: 1.23.0
- `pandas`: 1.2.4
- `python_igraph`: 0.11.6
- `python_louvain`: 0.15
- `rpy2`: 3.5.16
- `scikit_bio`: 0.6.2
- `scipy`: 1.6.2
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

Command:

    conda create -n taphonomeAnalyst2 python=3.8.8 r-base=4.3.2
    conda activate taphonomeAnalyst2
    pip install -r ./requirements.txt
    Rscript ./install_packages.R

#### * Note : How to download and install conda? [Documentation](https://docs.conda.io/projects/conda/en/stable/).

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

### 1. Hierarchical clustering-sedimentary environment. **(Module Ⅳ and Ⅴ)** [clustermap]

Command:
        
    python ./TaphonomeAnalyst2.py clusterenv --input ./Supplementary material2 --aquatic 'Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l)' --geochem './Supplementary material3.xlsx' --level 'species'

Result:

![clusterenv_aquatic.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/clusterenv_aquatic.png)

![clusterenv_tree.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/clusterenv_tree.png)

![clusterenv_geochem.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/clusterenv_geochem.png)

### 2. Venn diagram-sampling locations or environments. **(Module Ⅳ)** [venn]

Command:
        
    python ./TaphonomeAnalyst2.py divvenn --input ./Supplementary material2.xlsx --groups 'plot3-3/plot3-1/plot3-2,plot1-1/plot2-2/plot2-1/plot2-3/plot1-2/plot1-3' --level 'family

Result:

![divvenn.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/divvenn.png)

### 3. Taphonomic grades-taxa. **(Module Ⅲ)** [barh]

Command:

    python ./TaphonomeAnalyst2.py TGotus --input ./Supplementary material2.xlsx --level 'order'

Result:

![TGotus.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/TGotus.png)

### 4. Taphonomic grades-sampling plots (in customized order). **(Module Ⅲ)** [barh]

Command:

    python ./TaphonomeAnalyst2.py TGplots --input ./Supplementary material2.xlsx

    python ./TaphonomeAnalyst2.py TGplots --input ./Supplementary material2.xlsx --groups 'plot3-3/plot3-1/plot3-2,plot1-1/plot2-2/plot2-1/plot2-3/plot1-2/plot1-3'

Result:

![TGplots.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/TGplots.png)

![TGplots_groups.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/TGplots_groups.png)

### 5. Abundance-sampling plots. **(Module Ⅱ)** [barh]

Command:

    python ./TaphonomeAnalyst2.py abundplots --input ./Supplementary material2.xlsx --level 'order'

Result:

![abundplots.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/abundplots.png)

### 6. Co-occurrence networks. **(Module Ⅸ)** [network/venn]

Command:

    python ./TaphonomeAnalyst2.py cooccurnet --input ./Supplementary material2.xlsx --level 'family' --corr 'pearson' --corr_coef 0.7 --p_value 0.1

Result:

![cooccurnet.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/cooccurnet.png)

![cooccurnet_venn.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/cooccurnet_venn.png)

### 7. Sampling coverage curve. **(Module Ⅰ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py samplecurve --input ./Supplementary material2.xlsx --level 'family' --groups 'plot1:plot1-1/plot1-2/plot1-3,plot2:plot2-1/plot2-2/plot2-3,plot3:plot3-1/plot3-2/plot3-3'

Result:

![samplecurve.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/samplecurve.png)

### 8. Chao1 potential diversity curve. **(Module Ⅰ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py chao --input ./Supplementary material2.xlsx --level 'family' --groups 'plot1:plot1-1/plot1-2/plot1-3,plot2:plot2-1/plot2-2/plot2-3,plot3:plot3-1/plot3-2/plot3-3'

Result:

![chao.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/chao.png)

### 9. ACE potential diversity curve. **(Module Ⅰ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py ace --input ./Supplementary material2.xlsx --level 'family' --groups 'plot1:plot1-1/plot1-2/plot1-3,plot2:plot2-1/plot2-2/plot2-3,plot3:plot3-1/plot3-2/plot3-3' --rare 10

Result:

![ace.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/ace.png)

### 10. Heatmap-OTUs correlation analysis. **(Module Ⅷ)** [heatmap]

Command:

    python ./TaphonomeAnalyst2.py corrotus --input ./Supplementary material2.xlsx --level 'family' --corr 'pearson' --p_value 0.1

Result:

![corrotus.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/corrotus.png)

### 11. Mantel Test between species abundance and ecological environmental variables. **(Module Ⅶ)** [multiplot]

Command:

    python ./TaphonomeAnalyst2.py mantel --input ./Supplementary material2.xlsx --rhome 'C:\Program Files\R\R-4.3.2' --geochem './Supplementary material3.xlsx' --aquatic 'Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l)' --level_aquatic 'species' --level_terrestrial 'family' --corr 'pearson'

Result:

![mantel.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/mantel.png)

### 12. Assembling dissimilarity- environmental distance test. **(Module Ⅵ)** [regplot]

Command:

    python ./TaphonomeAnalyst2.py dissenvtest --input ./Supplementary material2.xlsx --aquatic 'Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l)' --level_aquatic 'species' --level_terrestrial 'family' --geochem './Supplementary material3.xlsx'

Result:

![dissenvtest.png](https://github.com/wma1995/TaphonomeAnalyst2/blob/main/Result/dissenvtest.png)

### 13. Generate a unified layout network for comparison. (Module Ⅹ)\t[network/boxplot/barplot]

Command:

    python ./TaphonomeAnalyst2.py netVC --input ./Supplementary material2.xlsx --aquatic 'Daohugounectes primitinus(l),Triglypta haifanggouensis,Triglypta haifanggouensis,Yanliaocorixa chinensis,Karataviella popovi,Samarura gigantea(l),Anisoptera fam. gen. sp1.(l),Platyperla platypoda(l),Ferganoconcha sibirica,Qiyia jurassica(l),Mesomyzon sp1.,Triops sp1.,Chirocephalidae gen. sp1.,Eurythoracalis mirabilis(l),Shantous lacustris(l),Foliomimus latus(l),Furvoneta viriosus(l),Furvoneta raucus(l),Mesobaetis sibirica(l),Clavineta eximia(l)' --level_aquatic 'species' --level_terrestrial 'family' --groups 'plot3-3/plot3-1/plot3-2,plot1-1/plot2-2/plot2-1/plot2-3/plot1-2/plot1-3' --corr 'spearman' --corr_coef 0.7 --p_value 0.01

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
