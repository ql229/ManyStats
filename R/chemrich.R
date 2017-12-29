### ManyStats analysis
## Why a package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY

##################################################
### ChemRICH Analysis for Logistic Regression ####
##################################################

calculateChemRICHlogit <- function (
  statres = "adni_logit_result.csv",
  classinfo = "adni_lipid_clusters.csv",
) {

  stat_res <- read.csv(file = statres, stringsAsFactors = F,header = T)
  classinfo <- read.csv(file = classinfo, stringsAsFactors = F,check.names = F)

  if(length(grep("N/A",classinfo$SMILES))>0) {
    stat_res <- stat_res[-grep("N/A",classinfo$SMILES),]
    classinfo <- classinfo[-grep("N/A",classinfo$SMILES),]
  }

  if(nrow(stat_res)!=length(which((stat_res$FLDNAME==classinfo$FLDNAME)==TRUE)))  { stop("FLDNAME field has different values in the class and the statistical result file.They shall be exactly same and in the same order.") }
  if( (length(which(is.na(classinfo$SMILES)==TRUE)))) { stop("Missing SMILES codes. Please check the input.") }
  if( length(which(table(classinfo$Compound.Name)>1))>0 ) { stop("Please remove the duplicates compound names") }

  classinfo$xlogp <- as.numeric(sapply(classinfo$SMILES, function(x)  { rcdk::get.xlogp(rcdk::parse.smiles(x)[[1]]) }))

  ## get the value locations.

  pval.ind <- grep("pvalue.$",names(stat_res))
  fc.ind <- grep("beta.$",names(stat_res))

  for (i in 1:length(pval.ind)){
    df1 <- data.frame(foldchange=stat_res[,fc.ind[i]],pvalue=stat_res[,pval.ind[i]],SMILES=classinfo$SMILES,class=classinfo$Cluster, stringsAsFactors = F)
    df1$xlogp <- classinfo$xlogp

    ## Code break rules.
    if( length(table(is.na(df1$pvalue)))==2 ) { stop("Missing pvalue") }
    if( length(table(is.na(df1$foldchange)))==2 ) { stop("Missing foldchange or effect size") }
    pacman::p_load(officer,openxlsx,grid,RJSONIO,RCurl,dynamicTreeCut,ape,ggplot2,ggrepel,phytools,plotrix,plotly, htmlwidgets,DT,extrafont,devEMF,rvg,magrittr)# better idea will be to load the packages as needed.

    #############################################
    ### Calculation of Enrichment Statistics#####
    #############################################

    clusterids <- names(which(table(df1$class)>2))
    clusterids <- clusterids[which(clusterids!="")] ## get rid of the empty cells.

    cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
      cl.member <- which(df1$class==x)
      if( length(which(df1$pvalue[cl.member]<.05)) >1 ){
        pval.cl.member <- df1$pvalue[cl.member]
        p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
        p.test.results$p.value
      } else {
        1
      }
    })


    cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
    clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)

    clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
      dfx <- df1[which(df1$class==x),]
      dfx$Compound.Name[which.min(dfx$pvalue)]
    })

    altrat <- sapply(clusterdf$name, function (k) {
      length(which(df1$class==k & df1$pvalue<0.05))/length(which(df1$class==k))
    })

    uprat <-sapply(clusterdf$name, function (k) {
      length(which(df1$class==k & df1$pvalue<0.05 & df1$foldchange > 0))/length(which(df1$class==k))
    })

    clust_s_vec <- sapply(clusterdf$name, function (k) {
      length(which(df1$class==k))
    })

    clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(df1$class==k & df1$pvalue<0.05))})
    clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(df1$class==k & df1$pvalue<0.05 & df1$foldchange > 0))})
    clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(df1$class==k & df1$pvalue<0.05 & df1$foldchange < 0))})
    clusterdf$upratio <- uprat
    clusterdf$altratio <- altrat
    clusterdf$csize <- clust_s_vec
    clusterdf <- clusterdf[which(clusterdf$csize>2),]
    clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")

    clusterdf$xlogp <- as.numeric(sapply(clusterdf$name, function(x) {  median(df1$xlogp[which(df1$class==x)]) }))

    clustdf <- clusterdf[which(cluster.pvalues!=1),]

    #################################################
    ########## Impact Visualization Graph ###########
    #################################################

    clustdf.alt.impact <- clustdf[which(clustdf$pvalues<0.05 & clustdf$csize>1 & clustdf$alteredMetabolites>1) ,]
    #clustdf.alt.impact <- clustdf
    clustdf.alt.impact <- clustdf.alt.impact[order(clustdf.alt.impact$xlogp),]
    clustdf.alt.impact$order <- 1:nrow(clustdf.alt.impact) ### Order is decided by the hclust algorithm.
    clustdf.alt.impact$logPval <- -log(clustdf.alt.impact$pvalues)

    p2 <- ggplot(clustdf.alt.impact,aes(x=xlogp,y=-log(pvalues)))
    p2 <- p2 + geom_point(aes(size=csize, color=upratio)) +
      #labs(subtitle = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
      scale_color_gradient(low = "blue", high = "red")+
      scale_size(range = c(5, 30)) +
      scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+4  )) +
      scale_x_continuous(" median XlogP of clusters ") +
      theme_bw() +
      #labs(title = "ChemRICH cluster impact plot") +
      geom_label_repel(aes(label = name), color = "gray20",family="Arial",data=subset(clustdf.alt.impact, csize>2),force = 5)+
      theme(text=element_text(family="Arial Black"))+
      theme(
        plot.title = element_text(face="bold", size=30,hjust = 0.5),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20, angle=90),
        panel.grid.major = element_blank(), # switch off major gridlines
        panel.grid.minor = element_blank(), # switch off minor gridlines
        legend.position = "none", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
        legend.title = element_blank(), # switch off the legend title
        legend.text = element_text(size=12),
        legend.key.size = unit(1.5, "lines"),
        legend.key = element_blank(), # switch off the rectangle around symbols in the legend
        legend.spacing = unit(.05, "cm"),
        axis.text.x = element_text(size=10,angle = 0, hjust = 1),
        axis.text.y = element_text(size=15,angle = 0, hjust = 1)
      )

    read_pptx() %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(p2), type = "body", width=10, height=8, offx =0.0 , offy = 0.0) %>%
      print(target = paste0("chemrich_impact ",i,"_",names(stat_res)[pval.ind[i]],".pptx")) %>%
      invisible()
  }
}


##################################################
### ChemRICH Analysis for Logistic Regression ####
##################################################

calculateChemRICHttest <- function (
  statres = "ttest_results",
  classinfo = "compound_clusterinfo.csv",
) {

  stat_res <- read.csv(file = statres, stringsAsFactors = F,header = T)
  classinfo <- read.csv(file = classinfo, stringsAsFactors = F,check.names = F)

  if(length(grep("N/A",classinfo$SMILES))>0) {
    stat_res <- stat_res[-grep("N/A",classinfo$SMILES),]
    classinfo <- classinfo[-grep("N/A",classinfo$SMILES),]
  }

  if(nrow(stat_res)!=length(which((stat_res$FLDNAME==classinfo$FLDNAME)==TRUE)))  { stop("FLDNAME field has different values in the class and the statistical result file.They shall be exactly same and in the same order.") }
  if( (length(which(is.na(classinfo$SMILES)==TRUE)))) { stop("Missing SMILES codes. Please check the input.") }
  if( length(which(table(classinfo$Compound.Name)>1))>0 ) { stop("Please remove the duplicates compound names") }

  classinfo$xlogp <- as.numeric(sapply(classinfo$SMILES, function(x)  { rcdk::get.xlogp(rcdk::parse.smiles(x)[[1]]) }))

  ## get the value locations.

  pval.ind <- grep("pvalue.$",names(stat_res))
  fc.ind <- grep("beta.$",names(stat_res))

  for (i in 1:length(pval.ind)){
    df1 <- data.frame(foldchange=stat_res[,fc.ind[i]],pvalue=stat_res[,pval.ind[i]],SMILES=classinfo$SMILES,class=classinfo$Cluster, stringsAsFactors = F)
    df1$xlogp <- classinfo$xlogp

    ## Code break rules.
    if( length(table(is.na(df1$pvalue)))==2 ) { stop("Missing pvalue") }
    if( length(table(is.na(df1$foldchange)))==2 ) { stop("Missing foldchange or effect size") }
    pacman::p_load(officer,openxlsx,grid,RJSONIO,RCurl,dynamicTreeCut,ape,ggplot2,ggrepel,phytools,plotrix,plotly, htmlwidgets,DT,extrafont,devEMF,rvg,magrittr)# better idea will be to load the packages as needed.

    #############################################
    ### Calculation of Enrichment Statistics#####
    #############################################

    clusterids <- names(which(table(df1$class)>2))
    clusterids <- clusterids[which(clusterids!="")] ## get rid of the empty cells.

    cluster.pvalues <- sapply(clusterids, function(x) { # pvalues were calculated if the set has at least 2 metabolites with less than 0.10 pvalue.
      cl.member <- which(df1$class==x)
      if( length(which(df1$pvalue[cl.member]<.05)) >1 ){
        pval.cl.member <- df1$pvalue[cl.member]
        p.test.results <- ks.test(pval.cl.member,"punif",alternative="greater")
        p.test.results$p.value
      } else {
        1
      }
    })


    cluster.pvalues[which(cluster.pvalues==0)] <- 2.2e-20 ### All the zero are rounded to the double.eps pvalues.\
    clusterdf <- data.frame(name=clusterids,pvalues=cluster.pvalues, stringsAsFactors = F)

    clusterdf$keycpdname <- sapply(clusterdf$name, function(x) {
      dfx <- df1[which(df1$class==x),]
      dfx$Compound.Name[which.min(dfx$pvalue)]
    })

    altrat <- sapply(clusterdf$name, function (k) {
      length(which(df1$class==k & df1$pvalue<0.05))/length(which(df1$class==k))
    })

    uprat <-sapply(clusterdf$name, function (k) {
      length(which(df1$class==k & df1$pvalue<0.05 & df1$foldchange > 0))/length(which(df1$class==k))
    })

    clust_s_vec <- sapply(clusterdf$name, function (k) {
      length(which(df1$class==k))
    })

    clusterdf$alteredMetabolites <- sapply(clusterdf$name, function (k) {length(which(df1$class==k & df1$pvalue<0.05))})
    clusterdf$upcount <- sapply(clusterdf$name, function (k) {length(which(df1$class==k & df1$pvalue<0.05 & df1$foldchange > 0))})
    clusterdf$downcount <- sapply(clusterdf$name, function (k) {length(which(df1$class==k & df1$pvalue<0.05 & df1$foldchange < 0))})
    clusterdf$upratio <- uprat
    clusterdf$altratio <- altrat
    clusterdf$csize <- clust_s_vec
    clusterdf <- clusterdf[which(clusterdf$csize>2),]
    clusterdf$adjustedpvalue <- p.adjust(clusterdf$pvalues, method = "fdr")

    clusterdf$xlogp <- as.numeric(sapply(clusterdf$name, function(x) {  median(df1$xlogp[which(df1$class==x)]) }))

    clustdf <- clusterdf[which(cluster.pvalues!=1),]

    #################################################
    ########## Impact Visualization Graph ###########
    #################################################

    clustdf.alt.impact <- clustdf[which(clustdf$pvalues<0.05 & clustdf$csize>1 & clustdf$alteredMetabolites>1) ,]
    #clustdf.alt.impact <- clustdf
    clustdf.alt.impact <- clustdf.alt.impact[order(clustdf.alt.impact$xlogp),]
    clustdf.alt.impact$order <- 1:nrow(clustdf.alt.impact) ### Order is decided by the hclust algorithm.
    clustdf.alt.impact$logPval <- -log(clustdf.alt.impact$pvalues)

    p2 <- ggplot(clustdf.alt.impact,aes(x=xlogp,y=-log(pvalues)))
    p2 <- p2 + geom_point(aes(size=csize, color=upratio)) +
      #labs(subtitle = "Figure Legend : Point size corresponds to the count of metabolites in the group. Point color shows that proportion of the increased metabolites where red means high and blue means low number of upregulated compounds.")+
      scale_color_gradient(low = "blue", high = "red")+
      scale_size(range = c(5, 30)) +
      scale_y_continuous("-log (pvalue)",limits = c(0, max(-log(clustdf.alt.impact$pvalues))+4  )) +
      scale_x_continuous(" median XlogP of clusters ") +
      theme_bw() +
      #labs(title = "ChemRICH cluster impact plot") +
      geom_label_repel(aes(label = name), color = "gray20",family="Arial",data=subset(clustdf.alt.impact, csize>2),force = 5)+
      theme(text=element_text(family="Arial Black"))+
      theme(
        plot.title = element_text(face="bold", size=30,hjust = 0.5),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=20, angle=90),
        panel.grid.major = element_blank(), # switch off major gridlines
        panel.grid.minor = element_blank(), # switch off minor gridlines
        legend.position = "none", # manually position the legend (numbers being from 0,0 at bottom left of whole plot to 1,1 at top right)
        legend.title = element_blank(), # switch off the legend title
        legend.text = element_text(size=12),
        legend.key.size = unit(1.5, "lines"),
        legend.key = element_blank(), # switch off the rectangle around symbols in the legend
        legend.spacing = unit(.05, "cm"),
        axis.text.x = element_text(size=10,angle = 0, hjust = 1),
        axis.text.y = element_text(size=15,angle = 0, hjust = 1)
      )

    read_pptx() %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(p2), type = "body", width=10, height=8, offx =0.0 , offy = 0.0) %>%
      print(target = paste0("chemrich_impact ",i,"_",names(stat_res)[pval.ind[i]],".pptx")) %>%
      invisible()
  }
}




