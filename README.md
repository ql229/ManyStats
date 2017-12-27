# ManyStats
Power functions for calculating statistics/graphics for complex metabolomics projects with several experimental groups

Several metabolomics projects can have a complex design with several experimental groups. Performing statistical analyses for those projects using online statistical tools such as metaboanalyst is tedious. 

ManyStats has some power functions (mostly wrappers) for doing statistics for those complex projects. 

Step 1 : Prepare the inputs

Metabolomics data table need to be split into three files -

a. Data matrix
b. Sample meta-data
c. Data dictionary

Step 2: Prepare the parameters file :

For each statistical test, an user need to create a parameter file in excel which will be used by Many Stats functions. The file will be in CSV format and can also be generated automatically. 

Step 3: Use the function :

Call the functions using the files created in step 1-2. 


# Installation and Usage Instructions

if (!require("devtools"))

install.packages('devtools', repos="http://cran.rstudio.com/")

if (!require("opencpu"))

install.packages('opencpu', repos="http://cran.rstudio.com/")

if (!require("RCurl"))

install.packages('RCurl', repos="http://cran.rstudio.com/")

library(devtools)

library(RCurl)

source('https://bioconductor.org/biocLite.R')

install_github('barupal/ManyStats')

library(ManyStats)

pacman::p_load(gridExtra,ggplot2,officer,magrittr,rvg,flextable,ggplot2,plotly,ggpubr)

createCSVFiles("A metabolomics dataset")










