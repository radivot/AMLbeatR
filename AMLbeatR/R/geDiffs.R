#'Gene Expression Differences 
#'
#'Uses Limma to find genes differential expressed when a gene is mutated vs. not (reference). 
#'
#'@param gene  Gene commonly mutated in AML    
#@param d  Clinical dataframe made by tidyClin
#'@param v     Variant dataframe v made by mkBeatAML
#'@param cpm     Counts per million dataframe cpm made by mkBeatAML
#'@param N     Top N most differentially expressed genes
#'@param cutOff     Counts per million reads below which a gene is ignored (median across samples)
#'@return  Returns Limma's topTable(X, coef=2, adjust="BH",number=n)
#'@author Tom Radivoyevitch
#'@examples
#' library(AMLbeatR)
#' load("~/data/BeatAML/BeatAML.RData") # 672  specimens (rows in clin) 
#' geDiffs("TET2",v,cpm)
#'@name geDiffs
#'@export
#'@import  dplyr   
#'@importFrom limma voom lmFit eBayes topTable 
#'@importFrom Biobase ExpressionSet pData<- 
#'@importFrom stats model.matrix median 
#'@importFrom tidyr nest unnest unite
#'@importFrom stringr str_c str_length
#'@importFrom purrr map map_chr map_dbl

geDiffs<-function(gene,v,cpm,N=2000,cutOff=1)  {
  # labId=lid=vaf=t_vaf=symbol=data=muts=vafs=ref=alt=mut1=cnts=insLen=NULL
  labId=lid=vaf=t_vaf=symbol=data=muts=vafs=NULL
  
  # library(tidyverse)
  # library(AMLbeatR)
  # library(Biobase)
  # library(limma)
  # load("~/data/BeatAML/BeatAML.RData")
  # # (d=tidyClin(clin)) # collected from 562 patients (rows in d)
  # n=2000
  # cutOff=1
  # gene="TET2"
  
  names(v)
  # v%>%filter(symbol=="TET2")%>%select(lid=labId,t_vaf,symbol)
  getPD=function(v,gene) {
    sv=v%>%mutate(sym=ifelse(symbol==gene,"Mutant","WT"))%>%select(lid=labId,t_vaf,sym)
    sv=sv%>%group_by(lid,sym)%>%summarize(vaf=mean(t_vaf,na.rm=T))
    sv=sv%>%group_by(lid)%>%nest()
    getM=function(x) if("Mutant"%in%x$sym) "Mutant" else "WT" 
    sv=sv%>%mutate(State=map_chr(data,getM))%>%select(-data)%>%arrange(lid)
    pD=data.frame(sv)
    rownames(pD)=sv$lid
    pD
  }
  (pD=getPD(v,gene))
  eids=names(cpm)[-c(1:2)] #451 RNAseq measurements
  int=intersect(eids,pD$lid) #399 
  pD=pD%>%filter(lid%in%int)
  (meta=cpm[names(cpm)[1:2]])
  expr=as.matrix(cpm[int])
  rownames(expr)=meta[,2]
  eset=ExpressionSet(assayData=expr)
  I=apply(expr,1,median)
  I=I>cutOff
  eset=eset[I,]
  # eset
  rownames(pD)=int
  pD$State=factor(pD$State,c("WT","Mutant"))
  pData(eset)=pD
  (design=model.matrix(~pD$State))
  v <- voom(eset,design,plot=TRUE,normalize.method="quantile")
  fit <- lmFit(v,design)
  efitM <- eBayes(fit)
  topTable(efitM, coef=2, adjust.method="BH",number=N)
}
