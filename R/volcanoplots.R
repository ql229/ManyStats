### ManyStats analysis
## Why a package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY

################################################
## Create volcano plots for a logistic model  ##
################################################


createVolcanoLogit(
  inputData = "ADNI_Diag_lipidomics_logistic_model_final.txt",
  pvalueCut = 0.05, ## This will draw a horizontal line on the vlocano plot
  fdrCutff = 0.20, # This will draw another line
  betacutoff = 0.2, # vertical lines for the eta
  labelOnly = 5, # only the top 5% significany compounds will be labeled.
  labalVar = "Compound.Name",
  interactive =TRUE # if true an interactive plot using ggplotly will be exported.
){
  stat_res <- read.delim(inputData, header = T, stringsAsFactors = F)

  ## get the value locations.

  pval.ind <- grep("pvalue.$",names(stat_res))
  pvalfdr.ind <- grep("pval.fdr.$",names(stat_res))
  fc.ind <- grep("beta.$",names(stat_res))



  for (i in 1:length(pval.ind)){
    df2 <- data.frame(clabel=stat_res[labalVar][,1],pval=stat_res[,pval.ind[i]],pval.fdr=  p.adjust(stat_res[,pval.ind[i]],method="fdr"),fc=stat_res[,fc.ind[i]], stringsAsFactors = F)

    df2$Changed <- "No Change"
    df2$Changed[which(df2$pval<0.05 & df2$fc>0)] <- "UP"
    df2$Changed[which(df2$pval<0.05 & df2$fc<0)] <- "DOWN"
    df2$Changed <- as.factor(df2$Changed)
    #df2$fc <- round(sapply(df2$fc, function(x) { if(x>1) {x} else {1/x} }), digits = 1)
    #df2$fc[ df2$fc>10] <- 10

    p2 <-   ggplot(df2, aes(label=clabel,x=fc, y=-log(pval,base = 10),colour = Changed)) +
      #geom_line(position=pd, size=2)+
      #geom_errorbar(aes(ymin = V2-V3 , ymax=V2+V3), width=.3,size=2,position=pd) +
      geom_point(size=5) + # 21 is filled circle
      #geom_bar(stat="identity", size=.1,position=position_dodge()) +
      scale_y_continuous("pvalue (-log10)") +
      scale_x_continuous("Logistic regression beta coefficient") +
      scale_color_manual("Beta direction",values=c("blue", "yellow", "red","white")) +
      #scale_fill_manual("",values=c("white", "yellow", "red","white")) +
      #scale_shape_manual("Pathway found",values=c(1,16))+
      #scale_shape(solid = FALSE) +
      theme_bw() +
      labs(title = "Metabolite volcano plot") +
      theme(
        plot.title = element_text(face="bold", size=30,hjust = 0.5),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=30, angle=90),
        panel.grid.major = element_blank(), # switch off major gridlines
        panel.grid.minor = element_blank(), # switch off minor gridlines
        #legend.justification=c(1,0),
        #legend.position=c(1,.6),
        #legend.position = "none",
        #legend.title = element_blank(), # switch off the legend title
        legend.text = element_text(size=25),
        #legend.key.size = unit(1.5, "lines"),
        #legend.key = element_blank(), # switch off the rectangle around symbols in the legend
        #legend.spacing = unit(.05, "cm"),
        #axis.text.x = element_text(size=15,angle = 45, hjust = 1.0),
        axis.text.x= element_text(size=15,angle = 0, hjust = 0.5),
        axis.text.y = element_text(size=15,angle = 0, hjust = 0.5)
      )

    p3 <- p2 + geom_hline(yintercept = -(log(pvalueCut,base = 10)) ) +
      geom_label_repel(aes(label = clabel), color = "gray20", family = "Arial", data = subset(df2,pval.fdr< fdrCutff), force = 5) +
      #geom_text(data= subset(df2,pval.fdr< fdrCutff),aes(label=clabel),hjust=0.5, vjust=-1) +  # label only the compound below FDR cutoff.
      annotate("text", max(df2$fc), -(log(pvalueCut,base = 10)), vjust = -1, label = paste0("pvalue ",pvalueCut)) +
      geom_hline(yintercept = -(log(max(df2$pval[which(df2$pval.fdr<fdrCutff)]),base = 10)) ) +
      annotate("text", max(df2$fc), -(log(max(df2$pval[which(df2$pval.fdr<fdrCutff)]),base = 10)), vjust = -1, label = paste0("FDR ",fdrCutff))
    plot(p3)
    read_pptx() %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(p3), type = "body", width = 10,height = 7, offx = .3, offy = .3) %>%
      print(target = gsub("_pval$","_volcano_plot.pptx",names(stat_res)[pval.ind[i]])) %>%
      invisible()

    p4 <- ggplotly(p3, width = 1600, height = 1000)
    htmlwidgets::saveWidget(p4, gsub("_pval$","_volcano_plot.html",names(stat_res)[pval.ind[i]]), selfcontained = T)
  }
}


################################################
## Create volcano plots for a Student ttest.  ##
################################################


createVolcano.ttest(
  inputData = "ADNI_Diag_lipidomics_logistic_model_final.txt",
  pvalueCut = 0.05, ## This will draw a horizontal line on the vlocano plot
  fdrCutff = 0.20, # This will draw another line
  betacutoff = 0.2, # vertical lines for the eta
  labelOnly = 5, # only the top 5% significany compounds will be labeled.
  labalVar = "Compound.Name",
  interactive =TRUE # if true an interactive plot using ggplotly will be exported.
){
  stat_res <- read.delim(inputData, header = T, stringsAsFactors = F)

  ## get the value locations.

  pval.ind <- grep("pvalue.$",names(stat_res))
  pvalfdr.ind <- grep("pval.fdr.$",names(stat_res))
  fc.ind <- grep("beta.$",names(stat_res))



  for (i in 1:length(pval.ind)){
    df2 <- data.frame(clabel=stat_res[labalVar][,1],pval=stat_res[,pval.ind[i]],pval.fdr=  p.adjust(stat_res[,pval.ind[i]],method="fdr"),fc=stat_res[,fc.ind[i]], stringsAsFactors = F)

    df2$Changed <- "No Change"
    df2$Changed[which(df2$pval<0.05 & df2$fc>0)] <- "UP"
    df2$Changed[which(df2$pval<0.05 & df2$fc<0)] <- "DOWN"
    df2$Changed <- as.factor(df2$Changed)
    #df2$fc <- round(sapply(df2$fc, function(x) { if(x>1) {x} else {1/x} }), digits = 1)
    #df2$fc[ df2$fc>10] <- 10

    p2 <-   ggplot(df2, aes(label=clabel,x=fc, y=-log(pval,base = 10),colour = Changed)) +
      #geom_line(position=pd, size=2)+
      #geom_errorbar(aes(ymin = V2-V3 , ymax=V2+V3), width=.3,size=2,position=pd) +
      geom_point(size=5) + # 21 is filled circle
      #geom_bar(stat="identity", size=.1,position=position_dodge()) +
      scale_y_continuous("pvalue (-log10)") +
      scale_x_continuous("Logistic regression beta coefficient") +
      scale_color_manual("Beta direction",values=c("blue", "yellow", "red","white")) +
      #scale_fill_manual("",values=c("white", "yellow", "red","white")) +
      #scale_shape_manual("Pathway found",values=c(1,16))+
      #scale_shape(solid = FALSE) +
      theme_bw() +
      labs(title = "Metabolite volcano plot") +
      theme(
        plot.title = element_text(face="bold", size=30,hjust = 0.5),
        axis.title.x = element_text(face="bold", size=20),
        axis.title.y = element_text(face="bold", size=30, angle=90),
        panel.grid.major = element_blank(), # switch off major gridlines
        panel.grid.minor = element_blank(), # switch off minor gridlines
        #legend.justification=c(1,0),
        #legend.position=c(1,.6),
        #legend.position = "none",
        #legend.title = element_blank(), # switch off the legend title
        legend.text = element_text(size=25),
        #legend.key.size = unit(1.5, "lines"),
        #legend.key = element_blank(), # switch off the rectangle around symbols in the legend
        #legend.spacing = unit(.05, "cm"),
        #axis.text.x = element_text(size=15,angle = 45, hjust = 1.0),
        axis.text.x= element_text(size=15,angle = 0, hjust = 0.5),
        axis.text.y = element_text(size=15,angle = 0, hjust = 0.5)
      )

    p3 <- p2 + geom_hline(yintercept = -(log(pvalueCut,base = 10)) ) +
      geom_label_repel(aes(label = clabel), color = "gray20", family = "Arial", data = subset(df2,pval.fdr< fdrCutff), force = 5) +
      #geom_text(data= subset(df2,pval.fdr< fdrCutff),aes(label=clabel),hjust=0.5, vjust=-1) +  # label only the compound below FDR cutoff.
      annotate("text", max(df2$fc), -(log(pvalueCut,base = 10)), vjust = -1, label = paste0("pvalue ",pvalueCut)) +
      geom_hline(yintercept = -(log(max(df2$pval[which(df2$pval.fdr<fdrCutff)]),base = 10)) ) +
      annotate("text", max(df2$fc), -(log(max(df2$pval[which(df2$pval.fdr<fdrCutff)]),base = 10)), vjust = -1, label = paste0("FDR ",fdrCutff))
    plot(p3)
    read_pptx() %>% add_slide(layout = "Title and Content", master = "Office Theme") %>%
      ph_with_vg(code = print(p3), type = "body", width = 10,height = 7, offx = .3, offy = .3) %>%
      print(target = gsub("_pval$","_volcano_plot.pptx",names(stat_res)[pval.ind[i]])) %>%
      invisible()

    p4 <- ggplotly(p3, width = 1600, height = 1000)
    htmlwidgets::saveWidget(p4, gsub("_pval$","_volcano_plot.html",names(stat_res)[pval.ind[i]]), selfcontained = T)
  }
}
