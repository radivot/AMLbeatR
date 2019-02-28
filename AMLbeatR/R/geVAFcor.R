#'Gene Expression Correlations with VAFs 
#'
#'Uses Limma to find genes whosw expression correlates with VAF of a specific mutation. 
#'
#'@param gene  Gene commonly mutated in AML    
#@param d  Clinical dataframe made by tidyClin
#'@param v     Variant dataframe v made by mkBeatAML
#'@param cpm     Counts per million dataframe cpm made by mkBeatAML
#'@param N     Top N genes most correlated with VAF 
#'@param cutOff     Counts per million reads below which a gene is ignored (median across samples)
#'@return  Returns Limma's topTable(X, coef=2, adjust="BH",number=n) where coef=2 is the linear slope parameter
#'@author Tom Radivoyevitch
#'@examples
#' library(AMLbeatR)
#' load("~/data/BeatAML/BeatAML.RData") # 672  specimens (rows in clin) 
#' geVAFcor("TET2",v,cpm)
#'@name geVAFcor
#'@export
#'@import  dplyr   
#'@importFrom limma voom lmFit eBayes topTable 
#'@importFrom Biobase ExpressionSet pData<- 
#'@importFrom stats model.matrix median 
#'@importFrom tidyr nest unnest unite
#'@importFrom stringr str_c str_length
#'@importFrom purrr map map_chr map_dbl

geVAFcor<-function(gene,v,cpm,N=2000,cutOff=1)  {
  # labId=lid=vaf=t_vaf=symbol=data=muts=vafs=ref=alt=mut1=cnts=insLen=NULL
  labId=lid=vaf=t_vaf=symbol=data=muts=vafs=NULL
  
  # library(tidyverse)
  # library(AMLbeatR)
  # library(Biobase)
  # library(limma)
  # load("~/data/BeatAML/BeatAML.RData")
  # # (d=tidyClin(clin)) # collected from 562 patients (rows in d)
  # N=2000
  # cutOff=1
  # gene="TET2"
  # gene="NPM1"
  # gene="FLT3"
  
  getPDv=function(v,gene) {
    sv=v%>%filter(symbol==gene)%>%select(lid=labId,t_vaf)
    sv=sv%>%group_by(lid)%>%summarize(vaf=mean(t_vaf,na.rm=T))
    sv=sv%>%mutate(vaf=ifelse(vaf>0.5,0.5,vaf))
    sv=sv%>%mutate(svaf=factor(ifelse(vaf>0.25,"High","Low"),levels=c("Low","High")))
    pD=data.frame(sv)
    rownames(pD)=sv$lid
    pD
  }
  (pD=getPDv(v,gene))
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
  # pD$vaf 
  pData(eset)=pD
  (design=model.matrix(~pD$vaf))
  # (design=model.matrix(~pD$svaf))
  v <- voom(eset,design,plot=TRUE,normalize.method="quantile")
  fit <- lmFit(v,design)
  efitM <- eBayes(fit)
  tb=topTable(efitM, coef=2, adjust.method="BH",number=N)
  attr(tb,"counts")<-table(pD$svaf)
  # head(tb,10)
  tb
}

