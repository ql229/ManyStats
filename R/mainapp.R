### ManyStats analysis
## Why a package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY

## Import libararies
require("pacman")
pacman::p_load(gridExtra,ggplot2,officer,magrittr,rvg,flextable,ggplot2,plotly,ggpubr,openxlsx,installr,ggrepel)
pacman::p_load(officer,openxlsx,grid,RJSONIO,RCurl,dynamicTreeCut,ape,ggplot2,ggrepel,phytools,plotrix,plotly, htmlwidgets,DT,extrafont,devEMF,rvg,magrittr)# better idea will be to load the packages as needed.

####################
## Import the data #
####################

createCSVFiles <- function(input="An excel file") {
  datatable <- openxlsx::read.xlsx(input, colNames=F)
  xind <- 1
  for(i in 1:ncol(datatable)) {
    if(is.na(datatable[1,i])) {
      xind <- c(xind,i)
    } else {
      break
    }
  }
  xind <- xind[-1]
  yind <- 1
  for(i in 1:nrow(datatable)) {
    if(is.na(datatable[i,1])) {
      yind <- c(yind,i)
    } else {
      break
    }
  }
  yind <- yind[-1]
  if(datatable[ (tail(yind,1)+1),1]!="CPDID") { stop( "First column label does not have CPDID annotation in it")}
  if(datatable[1, (tail(xind,1)+1)]!="file id") { stop( "First row  does not have 'file id' annotation in it")}

  ## Data dictionary export
  datadict <- datatable[(tail(yind,1)+1):nrow(datatable),xind]
  write.table(datadict, file="data_dictionary.csv",row.names=FALSE, col.names=FALSE, sep=",")

  ## Sample meta data export
  sdf <- datatable[yind,(tail(xind,1)+1):ncol(datatable)]
  write.table(t(sdf), file="sample_metadata.csv",row.names=FALSE, col.names=FALSE, sep=",")

  ## data matrix export
  ndf <- datatable[(tail(yind,1)+2):nrow(datatable),(tail(xind,1)+2):ncol(datatable)]
  write.table(ndf, file="data_matrix.csv",row.names=FALSE, col.names=FALSE, sep=",")
}
#createCSVFiles(input="mx 69088_HepG2 cells_Hirahatake & Meissen_high fructose_summer course_08-2015_submitDATA.xlsx")








































