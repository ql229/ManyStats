### ManyStats analysis
## Why a new stat package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY

## Import libararies
require("pacman")
pacman::p_load(gridExtra,ggplot2,officer,magrittr,rvg,flextable,ggplot2,plotly,ggpubr)

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

# example createCSVFiles(input="mx 69088_HepG2 cells_Hirahatake & Meissen_high fructose_summer course_08-2015_submitDATA.xlsx")

## Student TTest for multiple groups.

calculate.ttests <- function(
  numericData = "",
  sampleInfo ="",
  cpdInfo="",
  ttestgroups="" ) {
  ndf <- read.csv(file = numericData, stringsAsFactors = F,header = F)
  cdf <- read.csv(file = cpdInfo, stringsAsFactors = F,check.names = F)
  sdf <- read.csv(file = sampleInfo, stringsAsFactors = F,check.names = F)
  tgroups <- read.csv(file = ttestgroups, stringsAsFactors = F, header = F)

  if(nrow(sdf)!=ncol(ndf)) { stop("Sample metadata file or data matrix is not complete.")   }

  ## Missing value computation
  row.min <- sapply(1:nrow(ndf), function (x) {
    vec <- ndf[x,]
    vec <- vec[!is.na(vec)]
    vec <- vec[which(vec!=0)]
    min(vec)
  })

  xrows <- which(row.min!=Inf)  ### select the rows without inf as min. Get rid of rows with Inf as minimum value.

  ndf <- ndf[xrows,] # selecting only rows and all columns
  cdf <- cdf[xrows,] # subset the compounds.

  row.min  <- row.min[xrows]
  # misdf <- misdf[xrows,]

  for (i in 1:nrow(ndf)) {
    ndf[i,][is.na(ndf[i,])] <- row.min[i] ## if it is NA, we replace it with minimum
    ndf[i,][which(ndf[i,]==0)] <- row.min[i] ## if is is 0, we replace it with min.
  }

  print("Missing value computation finished")

  ### T-Test calculations (pvalues, FDR and fold-change)
  ttest_pvals <- lapply(1:nrow(tgroups) , function(y) { # ttest pvalues
    ttestCol <- tgroups[y,1]
    pvec <- sapply(1:nrow(ndf), function (x) {
      res <- list()
      res$p.value <- 1
      tryCatch(res <- t.test(ndf[x,which(sdf[,ttestCol]==tgroups[y,3])], ndf[x,which(sdf[,ttestCol]==tgroups[y,2])]),
               error=function(e) {})
      res$p.value
    })

    fcvec <- sapply(1:nrow(ndf), function (x) {
      res <- 1
      tryCatch(res <- median(as.numeric(  ndf[x,which(sdf[,ttestCol]==tgroups[y,3])]  ))/median(as.numeric( ndf[x,which(sdf[,ttestCol]==tgroups[y,2])]     )),
               error=function(e) {})
      res
    })
    xdf <- data.frame(pvec,p.adjust(pvec,method = "fdr"),fcvec,stringsAsFactors = F)
    names(xdf) <- c(paste0(ttestCol,"_",tgroups[y,3],"_vs_",tgroups[y,2],"_pval"),paste0(ttestCol,"_",tgroups[y,3],"_vs_",tgroups[y,2],"_pval.fdr"), paste0(ttestCol,"_",tgroups[y,3],"_vs_",tgroups[y,2],"_foldchange"))
    xdf
  })

  print("TTest calculation finished")

  ttest.res.df <- do.call(cbind,ttest_pvals)
  ttest.res.df <- cbind(cdf,ttest.res.df)

  write.table(ttest.res.df,"student_ttest_results.txt",col.names = T,sep="\t",row.names = F)

  ### Export the excel output

  ttest.results <- read.delim("student_ttest_results.txt", header = T, stringsAsFactors = F)
  l <- list(TTESTResults = ttest.results)
  openxlsx::write.xlsx(l, file = "Student_t_test_results.xlsx", asTable = TRUE)
  print("TTest calculations were successfull, check out the results")

}

### Usage For the TTest.

calculate.ttests(
  numericData = "data_matrix.csv",
  sampleInfo ="sample_metadata.csv",
  cpdInfo="data_dictionary.csv",
  ttestgroups="ttest_param.csv"
)
