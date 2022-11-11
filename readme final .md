#Final Project Outline
####Title Differential Gene Expression in TCGA within Stage 1 kidney Adenocarcinomas comparing heavy and light smokers using DeSEQ2
#####Author Siqi Tao

##Overview of project 
I will identify differentially expressed genes between kidney Cancer Adenocarcomas for alcohol drinkers and vs. non drinkers. This analysis will utilize the package DeSEQ2 and follow the specific vignette: (http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) For this analysis, I'll use the TCGA cohort and have identified 66 ht-seq counts files for tumors that fit within my cohort with 18 light smokers（less than 1 cigarettes） and 48 heavy smokers. Within the analysis, I will control for race and ethnicity.

##Data

I will use the data from https://portal.gdc.cancer.gov/repository. Examining clinical data, there are 66 tumor samples, and 18 are defined by me as light smokers and 48 are identified as heavy smokers. The specific files are available are here [data](https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.demographic.race%22%2C%22value%22%3A%5B%22white%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22kidney%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TCGA%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.analysis.workflow_type%22%2C%22value%22%3A%5B%22STAR%20-%20Counts%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%5D%7D&searchTableTab=files).

##Milestone 1

You should define your milestone based on your initial view of the data. Due Date: Thursday November 11th

** Data fully loaded into vignette through HT-SEQ steps ** I will complete an entire first draft of analysis analyzed through the vignette.

##Milestone 2

You should define your milestone 2 which should be an initial running of the entire vignette for feedback. An initial completion of vignette. I will complete an entire first draft of analysis analyzed through the vignette.Data loaded into vignette (through htseq), for seeking feedback. Not all sections in the writing will be completed, but will be final project.

Deliverable

Due Date: December 3rd

A complete repository with clear documentation and description of your analysis and results.