#'Excel File of Gene Expression Differences 
#'
#'Loops geDiffs() over top n mutated genes making sheets of top N Limma-called differentially expressed genes.
#'
#'@param d  Clinical dataframe made by muts() after tidyClin() (set number of most mutated genes in muts call
#'@param v     Variant dataframe v made by mkBeatAML
#'@param cpm     Counts per million dataframe cpm made by mkBeatAML
#'@param f     Name of Excel file to write DGE.
#'@param N     Top N most differentially expressed genes in each sheet
#'@return  Returns Limma's topTable(X, coef=2, adjust="BH",number=n)
#'@author Tom Radivoyevitch
#'@examples
#' library(AMLbeatR)
#' load("~/data/BeatAML/BeatAML.RData") # 672  specimens (rows in clin) 
#' geDiffs("TET2",v,cpm)
#'@name geXL
#'@export
#'@import  dplyr openxlsx  
#'@importFrom limma voom lmFit eBayes topTable 
#'@importFrom Biobase ExpressionSet pData<- 
#'@importFrom stats model.matrix median 
#'@importFrom tidyr nest unnest unite
#'@importFrom stringr str_c str_length
#'@importFrom purrr map map_chr map_dbl

geXL<-function(d,v,cpm,f="~/Results/AML/topGE.xlsx",N=1000)  {
  # labId=lid=vaf=t_vaf=symbol=data=muts=vafs=ref=alt=mut1=cnts=insLen=NULL
  labId=lid=vaf=t_vaf=symbol=data=muts=vafs=NULL
  
  # library(tidyverse)
  # library(AMLbeatR)
  # library(openxlsx)
  # load("~/data/BeatAML/BeatAML.RData")
  # (d=tidyClin(clin))
  # n=30
  # d=muts(d,v,av,n)
  # N=1000
  # f="~/Results/AML/topGE.xlsx"
  
  hs1=createStyle(fgFill="#DDDDDD",halign="CENTER",textDecoration="bold")
  # hs2=createStyle(fgFill="#FFFFFF",halign="CENTER",textDecoration="bold")
  unlink(f)
  wb <- createWorkbook() 
  addWorksheet(wb,"clinical")
  writeData(wb,"clinical",data.frame("Beat AML Clinical Data"), startRow=1,startCol=1,colNames=F)
  writeData(wb,"clinical",data.frame(paste("Top",n,"genes")), startRow=1,startCol=17,colNames=F)
  writeData(wb,"clinical",data.frame("Point muts in all variants in Table S16"), startRow=1,startCol=19,colNames=F)
  writeData(wb,"clinical",d, startRow=2,startCol=1,headerStyle = hs1)
  freezePane(wb, "clinical" ,firstActiveRow = 3,  firstActiveCol = 8)
  # setColWidths(wb, "clinical", cols=1:24, widths ="auto")
  setColWidths(wb, "clinical", cols=1:16, widths = 6)
  setColWidths(wb, "clinical", cols=11, widths = 8)
  setColWidths(wb, "clinical", cols=17, widths = 26)
  setColWidths(wb, "clinical", cols=18, widths = 16)
  setColWidths(wb, "clinical", cols=19:24, widths = 4)
  
  attributes(d)
  (loop=names(attr(d,"topp")))
  for (i in 1:length(loop)) {
    # i="TET2"
    g=loop[i]
    cat("Working on ",i," out of top ",n," mutated genes:",g,"\n")
    tb=geDiffs(g,v,cpm,N=N)
    # head(tb)
    addWorksheet(wb,g)
    # writeData(wb,i,tb, startRow=1,startCol=1,headerStyle = hs1,rowNames=TRUE)
    writeData(wb,g,cbind(gene=row.names(tb),tb), startRow=1,startCol=1,headerStyle = hs1)
    freezePane(wb,g,firstActiveRow = 2,  firstActiveCol = 2)
    setColWidths(wb, g, cols=1, widths = 16)
  }
  
  saveWorkbook(wb,file=f,overwrite = TRUE)
  
  
  
}
