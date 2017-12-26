### ManyStats analysis
## Why a package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
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
#createCSVFiles(input="mx 69088_HepG2 cells_Hirahatake & Meissen_high fructose_summer course_08-2015_submitDATA.xlsx")

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

# ### Usage For the TTest.
#
# calculate.ttests(
#   numericData = "data_matrix.csv",
#   sampleInfo ="sample_metadata.csv",
#   cpdInfo="data_dictionary.csv",
#   ttestgroups="ttest_param.csv"
# )

###################################
## Create box and whisker plots ###
###################################

generateBoxPlots <- function(
  numericData = "data_matrix.csv",
  sampleInfo ="sample_metadata.csv",
  cpdInfo="data_dictionary.csv",
  bw.param="boxplot_param.csv"
) {
  ndf <- read.csv(file = numericData, stringsAsFactors = F,header = F)
  cdf <- read.csv(file = cpdInfo, stringsAsFactors = F,check.names = F)
  sdf <- read.csv(file = sampleInfo, stringsAsFactors = F,check.names = F)
  bwGroups <- read.csv(file = bw.param, stringsAsFactors = F,check.names = F)

  for(i in 1:nrow(bwGroups)) {
    ndf1 <- ndf[cdf$CPDID%in%bwGroups$metabolite[i],sdf[bwGroups$groupingVar[i]][,1]%in%strsplit(bwGroups$groupsInclude[i],";")[[1]]]
    sdf1 <- sdf [sdf[bwGroups$groupingVar[i]][,1]%in%strsplit(bwGroups$groupsInclude[i],";")[[1]],]
    cdf1 <- cdf[cdf$CPDID%in%bwGroups$metabolite[i],]
    df1 <- data.frame(sdf1,metabolite=as.numeric(ndf1))

    gg <- ggboxplot(df1, bwGroups$groupingVar[i], "metabolite", color = bwGroups$colorVar[i], palette =strsplit(bwGroups$colors[1],";")[[1]],
                    add = "jitter",ylab = paste0("metabolite levels ", "(",bwGroups$label[i],")"),order=strsplit(bwGroups$groupsInclude[i],";")[[1]])
    read_pptx() %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(gg), type = "body", width = 10,height = 7, offx = 2, offy = .3) %>%
      print(target = paste0(bwGroups$metabolite[i],"_bw_plot.pptx")) %>%
      invisible()
  }
}

# generateBoxPlots(
#   numericData = "data_matrix.csv",
#   sampleInfo ="sample_metadata.csv",
#   cpdInfo="data_dictionary.csv",
#   bw.param="boxplot_param.csv"
# )

#####################
#### PCA analysis ###
#####################

calculatePCAs <- function (
  numericData = "data_matrix.csv",
  sampleInfo ="sample_metadata.csv",
  cpdInfo="data_dictionary.csv",
  pcaParam ="pca_param.csv" # PCA parameters file.

){

  ndf <- read.csv(file = numericData, stringsAsFactors = F,header = F)
  cdf <- read.csv(file = cpdInfo, stringsAsFactors = F,check.names = F)
  sdf <- read.csv(file = sampleInfo, stringsAsFactors = F,check.names = F)
  pca.param <- read.csv(file = pcaParam, stringsAsFactors = F,header = T, check.names = F)

  ##################
  ## Calculate PCA #
  ##################

  for (i in 1:nrow(pca.param)) {

    ndf1 <- ndf[,sdf[pca.param$groupingVar[i]][,1]%in%strsplit(pca.param$groupsInclude[i],";")[[1]]]
    sdf1 <- sdf [sdf[pca.param$groupingVar[i]][,1]%in%strsplit(pca.param$groupsInclude[i],";")[[1]],]

    lvar <- ""
    tryCatch(lvar <- sdf1[pca.param$labelVar[i]][,1], error=function(e) {
      stop("Sample label not provided in the PCA param file")
    })

    cvar <- "black"      # color
    tryCatch(cvar <- sdf1[pca.param$groupingVar[i]][,1], error=function(e) {
    })

    if(length(unique(cvar))==1) {
      cvar = "black"
    } else {
      cvar <- as.factor(cvar)
    }

    shvar <- 16 # defualt shape is filled black circles.

    evar <- "" # elipse variable, we want to export them and user can edit them in the pptx
    tryCatch(evar <- sdf1[pca.param$groupingVar[i]][,1], error=function(e) {
      stop("Elipse draw variable is not provided in the PCA param file")
    })

    pcamat <- ndf1
    pcamat <- do.call("cbind",lapply(pcamat,as.numeric))
    pcamat[which(is.na(pcamat)==TRUE)] <- min(pcamat[which(is.na(pcamat)==FALSE)])
    pcamat[which(is.na(pcamat)==TRUE)] <- 100
    pca.res <- mixOmics::pca(pcamat, ncomp=2,max.iter=100,center = T, scale = F)
    pc_scores <- pca.res$loadings[[1]]

    data_bw <- data.frame(snames = lvar, PC1 = pc_scores[,1], PC2=pc_scores[,2], stringsAsFactors = F)

    data_bw$PC1 <- as.numeric(pc_scores[,1])
    data_bw$PC2 <- as.numeric(pc_scores[,2])

    f2 <- ggplot(data_bw, aes(PC2,PC1, label=snames)) +
      scale_y_continuous(paste("PC2 - variance explained : ",signif(pca.res$explained_variance[2],digit=4)," ",sep="")) +
      scale_x_continuous(paste("PC1 - variance explained : ",signif(pca.res$explained_variance[1],digit=4)," ",sep="")) +
      theme_bw() +
      labs(title = "Principal component analysis (PCA)") +
      theme(
        plot.title = element_text(face="bold", size=20),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20, angle=90),
        panel.grid.major = element_blank(), # switch off major gridlines
        panel.grid.minor = element_blank(), # switch off minor gridlines
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "bottom", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
        legend.title = element_blank(), # switch off the legend title
        legend.text = element_text(size=20),
        legend.key.size = unit(1.5, "lines"),
        legend.key = element_blank(), # switch off the rectangle around symbols in the legend
        legend.spacing = unit(.05, "cm"),
        axis.text.x = element_text(size=0,angle = 45, hjust = 1),
        axis.text.y = element_text(size=0,angle = 45, hjust = 1)
      )
    f3 <- f2 + geom_point(aes(size=5,colour=cvar)) + geom_text(aes(label=lvar),hjust=0, vjust=0)
    f4 <- f3 +  stat_ellipse( aes(colour=evar ), linetype = 1 )

    ## Export PCA slides

    read_pptx() %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(f4), type = "body", width = 7,height = 7, offx = 2, offy = .3) %>%
      print(target = paste0(pca.param[i,1],".pptx")) %>%
      invisible()
  }
}

###############
## PCA usage. #
###############

# calculatePCAs(
#   numericData = "data_matrix.csv",
#   sampleInfo ="sample_metadata.csv",
#   cpdInfo="data_dictionary.csv",
#   pcaParam ="pca_param.csv" # PCA parameters file.
# )

##########################
## Create volcano plots  #
##########################



##################
##### PLS-DA #####
##################






############################
### Box and Whisker plots ##
############################






##########################
### ChemRICH Analysis ####
##########################






##########################
### Pathway enrichment  ##
##########################


