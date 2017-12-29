### ManyStats analysis
## Why a package : For running statistical analysis of larger projects where multiple groups needs to be compared and large number of plots and graphics needs to be created and exported.
## Author : Dinesh Kumar Barupal dinkumar@ucdavis.edu
## Lincense : CC-BY

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
