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

# R version must be 3.4.3 or latest.

Run below commands in the R-console. It will installed the required packages and files. 

```
if (!require("devtools"))
install.packages('devtools', repos="http://cran.rstudio.com/")

if (!require("RCurl"))
install.packages('RCurl', repos="http://cran.rstudio.com/")

if (!require("opencpu"))
install.packages('opencpu', repos="http://cran.rstudio.com/")

if (!require("pacman"))
install.packages('pacman', repos="http://cran.rstudio.com/")

library(devtools)
library(RCurl)
library(pacman)

pacman::p_load(htmlwidgets,DT,RCurl,RJSONIO,ape,devEMF,dynamicTreeCut,extrafont,flextable,ggplot2,ggpubr,ggrepel,grid,gridExtra,installr,magrittr,officer,openxlsx,phytools,plotly,plotrix,rvg,mixOmics)

installr::install.Rtools()

install_github('barupal/ManyStats')
library(ManyStats)

```
# Step 1. Setup a new RStudio project.

For better organization and efficiency, create a new RStudio project. To create a new project, click on the File--> New project and then give it a name eg Study1 or whatever you like. By default R-Projects are save in your home directory which in "My Documents" for Windows. This will create a folder eg Study1 in your home directory. 

To know how to setup a new RStudio project - watch this youtube video. 
How to setup a new RStudio project[www.Youtube.com]

# Step 2. Move the metabolomcis dataset to the rproject folder. 

Copy and paste the metabolomics dataset (Excel sheet) that you have received from WCMC, Metabolon or other service centers into this folder (eg Study1). 

# Step 3. Create a new R-script file.
Now create a new R-script file and save it as "ManyStats.R". You will type all the codes for a project in this script. Watch this video on how to setup and run scripts in R. 

# Step 4. Load the ManyStats package.

In the script file type this command - 

```
pacman::p_load(htmlwidgets,DT,RCurl,RJSONIO,ape,devEMF,dynamicTreeCut,extrafont,flextable,ggplot2,ggpubr,ggrepel,grid,gridExtra,installr,magrittr,officer,openxlsx,phytools,plotly,plotrix,rvg,mixOmics) 
library(ManyStats)
```
It will load all the necessary packages that you will need to run the ManyStats functions. 

# Step 5. Generating the CSV files. 
Three CSV files representing data matrix, data dictionary and sample metadata are needed to perform statistical analyses using ManyStats functions. To generate these files run this command in R console. 

```
createCSVFiles(input="mx 69088_HepG2 cells_Hirahatake & Meissen_high fructose_summer course_08-2015_submitDATA.xlsx")
```
This will create three files 

*data_matrix.csv
*sample_metadata.csv
*data_dictionary.csv

in the project directory. You will use all these files as input for ManyStats functions. 

# Step 6. Calculating Statistics

# Calculate t-tests

To calculate t-tests for your study, you need to prepare a t-test parameter file ttest_param.csv. 
Watch this video on how to prepare this file for your project. 

```
calculate.ttests(
   numericData = "data_matrix.csv",
   sampleInfo ="sample_metadata.csv",
   cpdInfo="data_dictionary.csv",
   ttestgroups="ttest_param.csv"
 )

```











