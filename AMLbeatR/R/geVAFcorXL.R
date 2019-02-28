#'Excel File of Gene Expression Correlations with VAFs 
#'
#'Loops geVAFcor() over top n mutated genes making sheets of top N Limma-called VAF-correlated genes.
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
#' load("~/data/BeatAML/BeatAML.RData") 
#' (d=tidyClin(clin)) #672 rows in clin, one for each measurement; 562 rows in d, one for each patient
#' (d=muts(d,v,av,n=10))# adds muts and vafs of top n genes and point mut counts from all variants
#' geVAFcorXL(d,v,cpm,f="~/Results/AML/topGEvafCor.xlsx")#1st sht = clinical, rest = expr vs vafs
#'@name geVAFcorXL
#'@export
#'@import  dplyr openxlsx  
#'@importFrom limma voom lmFit eBayes topTable 
#'@importFrom Biobase ExpressionSet pData<- 
#'@importFrom stats model.matrix median 
#'@importFrom tidyr nest unnest unite
#'@importFrom stringr str_c str_length
#'@importFrom purrr map map_chr map_dbl

geVAFcorXL<-function(d,v,cpm,f="~/Results/AML/topGEvafCor.xlsx",N=1000)  {
  # labId=lid=vaf=t_vaf=symbol=data=muts=vafs=ref=alt=mut1=cnts=insLen=NULL
  labId=lid=vaf=t_vaf=symbol=data=muts=vafs=NULL
  
  # library(tidyverse)
  # library(AMLbeatR)
  # library(openxlsx)
  # load("~/data/BeatAML/BeatAML.RData")
  # (d=tidyClin(clin))
  # n=5
  # d=muts(d,v,av,n)
  # N=1000
  # f="~/Results/AML/topGEvafCor.xlsx"
  
  hs1=createStyle(fgFill="#DDDDDD",halign="CENTER",textDecoration="bold")
  # hs2=createStyle(fgFill="#FFFFFF",halign="CENTER",textDecoration="bold")
  unlink(f)
  wb <- createWorkbook() 
  addWorksheet(wb,"clinical")
  writeData(wb,"clinical",data.frame("Beat AML Clinical Data"), startRow=1,startCol=1,colNames=F)
  writeData(wb,"clinical",data.frame(paste("Top",N,"genes")), startRow=1,startCol=17,colNames=F)
  writeData(wb,"clinical",data.frame("Point muts in all variants in Table S16"), startRow=1,startCol=19,colNames=F)
  writeData(wb,"clinical",d, startRow=2,startCol=1,headerStyle = hs1)
  freezePane(wb, "clinical" ,firstActiveRow = 3,  firstActiveCol = 8)
  # setColWidths(wb, "clinical", cols=1:24, widths ="auto")
  setColWidths(wb, "clinical", cols=1:16, widths = 6)
  setColWidths(wb, "clinical", cols=11, widths = 8)
  setColWidths(wb, "clinical", cols=17, widths = 26)
  setColWidths(wb, "clinical", cols=18, widths = 16)
  setColWidths(wb, "clinical", cols=19:24, widths = 4)
  
  # attributes(d)
  (loop=names(attr(d,"topv")))
  n=length(loop)
  L=NULL
  L[["clin"]]=d
  for (i in 1:n) {
    # i=1
    g=loop[i]
    cat("Working on ",i," out of top ",n," mutated genes:",g,"\n")
    tb=geVAFcor(g,v,cpm,N=N)
    L[[g]]=tb
    nms=names(attr(tb,"counts"))
    cnts=as.numeric(attr(tb,"counts"))
    hd=paste("Number of GE Samples: ",nms[1]," = ",cnts[1],"  ",nms[2]," = ",cnts[2],sep="",collapse="")
    # head(tb)
    addWorksheet(wb,g)
    writeData(wb,g,data.frame(hd), startRow=1,startCol=1,colNames=F)
    # writeData(wb,i,tb, startRow=1,startCol=1,headerStyle = hs1,rowNames=TRUE)
    writeData(wb,g,cbind(gene=row.names(tb),tb), startRow=2,startCol=1,headerStyle = hs1)
    freezePane(wb,g,firstActiveRow = 3)
    setColWidths(wb, g, cols=1, widths = 16)
  }
  
  saveWorkbook(wb,file=f,overwrite = TRUE)
  invisible(L)
  L
}
