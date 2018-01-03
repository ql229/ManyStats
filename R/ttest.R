### ManyStats analysis
## Why a package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY

## Student TTest for multiple groups.

## roxygen parameters.
#' Calcualte multiple t-tests
#'
#' @param numericData an metabolomics dataset
#' @param sampleInfo sample metadata file
#' @param cpdInfo compound data dictionary
#' @param ttestgroups ttest parameter file
#' @return student t-test result in an excel sheet - Student_t_test_results.xlsx
#' @examples
#' calculate.ttests(numericData = system.file("/inst/Examples", "data_matrix.csv", package="ManyStats"), sampleInfo =system.file("/inst/Examples", "sample_metadata.csv", package="ManyStats"), cpdInfo=system.file("/inst/Examples", "data_dictionary.csv", package="ManyStats"), ttestgroups=system.file("/inst/Examples", "ttest_param.csv", package="ManyStats"))
#'
#'
#'
#'
#'

calculate.ttests <- function(
  numericData = "",
  sampleInfo ="",
  cpdInfo="",
  ttestgroups="" ) {

  pacman::p_load(gridExtra,ggplot2,officer,magrittr,rvg,flextable,ggplot2,plotly,ggpubr,openxlsx,installr)

  ndf <- read.csv(file = numericData, stringsAsFactors = F,header = F)
  cdf <- read.csv(file = cpdInfo, stringsAsFactors = F,check.names = F)
  sdf <- read.csv(file = sampleInfo, stringsAsFactors = F,check.names = F)
  tgroups <- read.csv(file = ttestgroups, stringsAsFactors = F, header = T)

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

# ### Usage For the TTest.
#
# calculate.ttests(
#   numericData = "data_matrix.csv",
#   sampleInfo ="sample_metadata.csv",
#   cpdInfo="data_dictionary.csv",
#   ttestgroups="ttest_param.csv"
# )
calculate.ttests(numericData = system.file("Examples", "data_matrix.csv", package="ManyStats"), sampleInfo =system.file("Examples", "sample_metadata.csv", package="ManyStats"), cpdInfo=system.file("Examples", "data_dictionary.csv", package="ManyStats"), ttestgroups=system.file("Examples", "ttest_param.csv", package="ManyStats"))


