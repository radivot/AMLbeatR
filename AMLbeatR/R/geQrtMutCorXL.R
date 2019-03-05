#'Excel File of Gene Expression 1st vs 4th quartile correlations with mutations 
#'
#'Loops over gene list to find mutations that correlate with low levels of gene expression. Odd Ratios are postive if
#'mutations correlate with low levels of expression of the gene. 
#'
#'@param d  Clinical dataframe made by muts() after tidyClin() (set number of most mutated genes in muts call
#'@param v     Variant dataframe v made by mkBeatAML
#'@param cpm     Counts per million dataframe cpm made by mkBeatAML
#'@param f     Name of Excel file to write Fisher.test() OR estimates and P values
#'@param genes     Loops through 1st vs 4th quartiles of expression of genes included in this character vector
#'@return  Invisibly, a list of dataframes, one for each sheet in the Excel file. 
#'@author Tom Radivoyevitch
#'@examples
#' library(AMLbeatR)
#' load("~/data/BeatAML/BeatAML.RData") 
#' d=tidyClin(clin) 
#' (d=muts(d,v,av,n=14))# crashes with n>14 
#' geQrtMutCorXL(d,v,cpm,genes=genesE)
#'@name geQrtMutCorXL
#'@export
#'@import  dplyr openxlsx  
#'@importFrom tidyr nest unnest unite
#'@importFrom stats quantile fisher.test
#'@importFrom stringr str_c str_detect
#'@importFrom purrr map map_chr map_dbl

geQrtMutCorXL<-function(d,v,cpm,f="~/Results/AML/geQrtMutCorXL.xlsx",genes=names(attr(d,"topv")))  {
  # labId=lid=vaf=t_vaf=symbol=data=muts=vafs=ref=alt=mut1=cnts=insLen=NULL
  labId=lid=vaf=t_vaf=symbol=data=muts=vafs=NULL
  
  # rm(list=ls()) 
  # library(tidyverse)
  # library(openxlsx)
  # library(AMLbeatR)
  # load("~/data/BeatAML/BeatAML.RData")
  # d=tidyClin(clin)
  # n=14   # top n most mutated genes in AML
  # d=muts(d,v,av,n) #make d with muts col and attributes
  # attributes(d)
  # (genes=names(attr(d,"topv")))
  # genesE=c("DDX41",sort(genes)) #one sheet for each of these gene expression quartiles
  # 
  # 
  # ###### code below is now in geXL
  # f="~/Results/AML/qrts.xlsx"
  

  (D0=tibble(genes))
  (d0=tibble(genes))
  
  qrts1=function(gene1,v,cpm) {
    Symbol=Gene=geneExpr=M=f=P=OR=hit=NULL
    # gene1="DDX41"
    
    x=cpm%>%filter(Symbol==gene1)%>%select(-Symbol,-Gene)
    x=unlist(x)
    x=sort(x)
    qs=quantile(x)
    (low=x[x<qs[2]])
    (hi=x[x>qs[4]])
    (lowS=names(low))
    (hiS=names(hi))
    (d=d%>%select(lid,muts)%>%filter(lid%in%c(lowS,hiS))%>%mutate(geneExpr=ifelse(lid%in%lowS,"low","hi")))
    (d=d%>%mutate(muts=ifelse(is.na(muts),"none",muts))) #missing = no mutations
    
    getD=function(gene2) {
      D=d%>%mutate(hit=ifelse(str_detect(muts,gene2),"yes","no"))
      D%>%group_by(geneExpr,hit)%>%summarize(n=n())
    }
    d0=D0
    (d0=d0%>%mutate(D=map(genes,getD)))
    d0$D
    
    getM=function(D) {
      matrix(D$n,ncol=2,dimnames=list(hit=c("no","yes"),expr=c("hi","low")))
    }
    d0=d0%>%mutate(M=map(D,getM))
    d0=d0%>%mutate(f=map(M,fisher.test))
    d0=d0%>%mutate(P=map_dbl(f,function(x) x$p.value))
    d0=d0%>%mutate(OR=map_dbl(f,function(x) x$estimate))
    d0=d0%>%mutate(ORci=map_chr(f,function(x) 
      str_c(round(x$estimate,2)," (",round(x$conf.int[1],2),", ",round(x$conf.int[2],2),")") ))
    d0=d0%>%select(-D,-M,-f)%>%arrange(P)
    d0
  }
  
  
  hs1=createStyle(fgFill="#DDDDDD",halign="CENTER",textDecoration="bold")
  # hs2=createStyle(fgFill="#FFFFFF",halign="CENTER",textDecoration="bold")
  unlink(f)
  wb <- createWorkbook() 
  addWorksheet(wb,"clinical")
  writeData(wb,"clinical",data.frame("Beat AML Clinical Data"), startRow=1,startCol=1,colNames=F)
  # writeData(wb,"clinical",data.frame(paste("Top",d$n,"genes")), startRow=1,startCol=17,colNames=F)
  # writeData(wb,"clinical",data.frame("Point muts in all variants in Table S16"), startRow=1,startCol=19,colNames=F)
  writeData(wb,"clinical",d, startRow=2,startCol=1,headerStyle = hs1)
  freezePane(wb, "clinical" ,firstActiveRow = 3,  firstActiveCol = 8)
  setColWidths(wb, "clinical", cols=1:16, widths = 6)
  setColWidths(wb, "clinical", cols=11, widths = 8)
  setColWidths(wb, "clinical", cols=17, widths = 26)
  setColWidths(wb, "clinical", cols=18, widths = 16)
  setColWidths(wb, "clinical", cols=19:24, widths = 4)
  N=length(genes)
  L=NULL
  L[["clin"]]=d
  for (i in 1:N) {
    # i=1
    g=genes[i]
    cat("Working on ",i," out of top ",N," mutated genes:",g,"\n")
    D=qrts1(g,v,cpm)
    L[[g]]=D
    addWorksheet(wb,g)
    # writeData(wb,i,tb, startRow=1,startCol=1,headerStyle = hs1,rowNames=TRUE)
    writeData(wb,g,D, startRow=1,startCol=1,headerStyle = hs1)
    freezePane(wb,g,firstActiveRow = 2,  firstActiveCol = 2)
    setColWidths(wb, g, cols=1, widths = 16)
  }
  saveWorkbook(wb,file=f,overwrite = TRUE)
  invisible(L)
  L
}
  
  
  
  
