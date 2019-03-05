#'Get Expression and Mutation Status for Two Genes 
#'
#'Returns a dataframe with expression values and mutation status for two genes 
#'
#'@param geneE  Gene Expression of interest    
#'@param geneM  Gene commonly mutated in AML    
#'@param v     Variant dataframe v made by mkBeatAML
#'@param cpm     Counts per million dataframe cpm made by mkBeatAML
#'@return  Returns dataframe with gene expression and mutation status for geneE and geneM
#'@author Tom Radivoyevitch
#'@examples
#' library(AMLbeatR)
#' load("~/data/BeatAML/BeatAML.RData") # 672  specimens (rows in clin) 
#' getExpr("FLT3","FLT3",v,cpm)
#' 
#'@name getExpr
#'@export
#'@import  dplyr   
#'@importFrom tidyr nest unnest unite
#'@importFrom stringr str_c str_length
#'@importFrom purrr map map_chr map_dbl

getExpr<-function(geneE,geneM,v,cpm)  {
  # labId=lid=vaf=t_vaf=symbol=data=muts=vafs=ref=alt=mut1=cnts=insLen=NULL
  labId=lid=vaf=t_vaf=symbol=data=muts=vafs=NULL
  
  # library(tidyverse)
  # library(AMLbeatR)
  # load("~/data/BeatAML/BeatAML.RData")
  # geneE=geneM="FLT3"
  
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
  (pD=getPD(v,geneM))
  eids=names(cpm)[-c(1:2)] #451 RNAseq measurements
  int=intersect(eids,pD$lid) #399 
  pD=pD%>%filter(lid%in%int)
  (meta=cpm[names(cpm)[1:2]])
  expr=as.matrix(cpm[int])
  rownames(expr)=meta[,2]
  pD$E=expr[geneE,]
  pD
}
