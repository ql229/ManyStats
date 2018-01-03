### ManyStats analysis
## Why a package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY

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

    if(pca.param$labelVar[i]!=""){
      tryCatch(lvar <- sdf1[pca.param$labelVar[i]][,1], error=function(e) {
        stop("Sample label not provided in the PCA param file")
      })
    }

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

    library(ggplot2)

    f2 <- ggplot(data_bw, aes(PC2,PC1, label=snames)) +
      scale_y_continuous(paste("PC2 - variance explained : ",signif(pca.res$explained_variance[2],digit=4)," ",sep="")) +
      scale_x_continuous(paste("PC1 - variance explained : ",signif(pca.res$explained_variance[1],digit=4)," ",sep="")) +
      theme_bw() +
      guides(size=FALSE)+
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

    library(officer)
    library(magrittr)
    library(rvg)

    read_pptx() %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(f4), type = "body", width = 7,height = 7, offx = 2, offy = .3) %>%
      print(target = paste0(pca.param[i,1],".pptx")) %>%
      invisible()
    print(paste0(pca.param[i,1],".pptx created"))
  }
}

###############
## PCA usage. #
###############

calculatePCAs(
  numericData = "data_matrix.csv",
  sampleInfo ="sample_metadata.csv",
  cpdInfo="data_dictionary.csv",
  pcaParam ="pca_param.csv" # PCA parameters file.
)
