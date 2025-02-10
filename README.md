# High dimensional Immune profiling following Autologous Hematopoietic Stem Cell Transplantation in Relapsing-Remitting Multiple Sclerosis   

We present the analysis of ~100 million single cells profiled with mass cytometry from 25 RRMS patients enrolled in one arm of the RAM-MS randomised clinical trial. The samples were profiled using a 41-antibody panel and different time points were screened inlcuding before aHSCT and 100 days, 6 months, 12 months and 24 months after aHSCT. 
Phenotyping was performed in a hierarchical scheme: First manual gating was performed to detect CD45+CD66- , CD45-CD66b+ and CD45-CD66b- cells. Then the CD45+CD66b- gate was undergone unsuperived clustering based on FlowSOM and ConsensusClustPlus algorithms. After two rounds of FlowSOM/ConsensusClustPlus we detected 17 different cell populations and we described their reconstitution dynamics. Diffcyt analysis was performed to quantify changes in the abundance of these cell populations before and after aHSCT. Next to fully utilise all antibodies in the panel, the detected 17 cell populations were further analysed using functional markers to detect distinct functional states within these populations. We present the discovery of 107 different functional states using a computational framework that combines unsupervised clustering with supervised learning based on Linear Discriminant Analysis (LDA). Finally, we present the development of a computational framework for multidimensional analysis of cell populations using RadViz/FreeViz and SARA scores.    

In this page we provide codes and datasets required to reproduce the analysis presented in our relevant publication

## Publication

Title: High dimensional Immune profiling following Autologous Hematopoietic Stem Cell Transplantation in Relapsing-Remitting Multiple Sclerosis   

Journal: The paper is submitted to Nature Communications jourmal 

Published: pre-print available at bioRxiv 

## Code Description

The framework is presented using different R Markdowns that implement parts of the analysis as follows:

1. Part1.Rmd : immune reconsitution trajectories of CD45+CD66b- cells before and after aHSCT (17 cell populations)

2. Part2.Rmd : analysis of functional states (107 cell populations)

3. Part3.Rmd : multidimensional analysis of cell populations using RadViz/FreeViz and SARA scores

Important note for users

Please modify all paths found in the R markdowns and change them to your computer's file system. Since we are not allowed to share patient personal information, we provide anonymised and aggregated single-cell data. For some parts of the analysis we also provide already processed data that have been made non-identifiable. 
For more specific information please contact the authors.

## Data Sources

Below we provide the sources of the datasets used in the study:

Annotated mass cytometry datasets are available at: XXXXXXX

In the same repo we provide trained LDA models for the cellular states as well as other pre-processed data that have are non-identifiable.


#### Releases

15-Feb-2025 : Beta version 1

## Contact

Comments and bug reports are welcome, please email: Dimitrios Kleftogiannis (dimitrios.kleftogiannis@uib.no)

We are also interested to know about how you have used our source code, including any improvements that you have implemented.
 
You are free to modify, extend or distribute our source code, as long as our copyright notice remains unchanged and included in its entirety. 

## License

This project is licensed under the MIT License.

Copyright 2024 Department of Clinical Sciences, University of Bergen (UiB) and Neuro-SysMed center for clinical treatment research, Norway

You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENSE.md".
