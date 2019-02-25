#'Add Mutations to Clinical Data 
#'
#'Adds mutation information from the variance dataframe v to the genetics list column of clin 
#'
#'@param d  Clinical dataframe made by tidyClin
#'@param v     Variant dataframe made by mkBeatAML
#'@param n     Top n most frequent mutations are encodes in the muts field
#'@return  d with additional columns of muts and vafs separated by / and point mutantion counts
#'@author Tom Radivoyevitch
#'@examples
#' library(AMLbeatR)
#' load("~/data/BeatAML/BeatAML.RData") # 672  specimens (rows in clin) 
#' (d=tidyClin(clin)) # collected from 562 patients (rows in d)
#' (d=muts(d,v)) 
#' (D=d%>%unnest()) 
#'@name muts
#'@export
#'@import  dplyr   
#'@importFrom tidyr nest unnest unite
#'@importFrom stringr str_c str_length
#'@importFrom purrr map map_chr map_dbl

muts<-function(d,v,n=10)  {
  labId=lid=vaf=t_vaf=symbol=data=muts=vafs=ref=alt=mut1=cnts=insLen=NULL
  
  # library(tidyverse)
  # library(AMLbeatR)
  # load("~/data/BeatAML/BeatAML.RData")
  # (d=tidyClin(clin)) # collected from 562 patients (rows in d)
  # n=10
  
  
  (t=sort(table(v$symbol),decreasing=T)[1:n])
  (top=names(t))
  sv=v%>%select(lid=labId,vaf=t_vaf,sym=symbol)%>%filter(sym%in%top)
  sv=sv%>%group_by(lid,sym)%>%summarize(vaf=max(vaf,na.rm=T))
  (tp=sort(table(sv$sym),decreasing=T)[1:n])
  sv=sv%>%group_by(lid)%>%nest()
  sv=sv%>%mutate(muts=map_chr(data,function(x) str_c(x$sym,collapse="/")),
                 vafs=map_chr(data,function(x) str_c(round(x$vaf,1),collapse="/"))
                 )%>%select(-data)
  # sort(table(sv$State))
  # sv
  D=d%>%unnest()
  D=left_join(D,sv)
  # x=c("d","sss","dsd")
  # x=c(NA,"sss","dsd")
  # x=c(NA,NA)
  getLong=function(x) {
   if (all(is.na(x))) return (x) else {
     k=which(str_length(x)==max(str_length(x),na.rm=T))
     k=k[length(k)]
     x[is.na(x)]=x[k]
     return(x)
   }
  }
  # getLong(x)
  # D
  D=D%>%nest(t:vafs)%>%mutate(data=map(data,function(x){x$muts=getLong(x$muts);x$vafs=getLong(x$vafs);return(x)}))
  D=D%>%unnest()
  D
  
  sv=v%>%select(ref,alt,lid=labId)%>%mutate(insLen=str_length(alt)-str_length(ref))%>%filter(insLen==0)
  sv=sv%>%unite(mut1,ref:alt,sep=">")%>%select(-insLen)
  map6=function(x) {
    # x=sv$mut1
    x[str_detect(x,"A>G")]="T>C"
    x[str_detect(x,"A>T")]="T>A"
    x[str_detect(x,"A>C")]="T>G"
    x[str_detect(x,"G>A")]="C>T"
    x[str_detect(x,"G>T")]="C>A"
    x[str_detect(x,"G>C")]="C>G"
    x
  }
  sv=sv%>%mutate(mut1=map6(mut1))
  # head(sv)
  mkTab=function(x) {
    # x=sv
    t=table(x$mut1)
    y=as.list(as.numeric(t))
    names(y)=names(t)
    as_tibble(y)
  }  
  
  ord=c("C>A","C>G","C>T","T>A","T>G","T>C")
  (tots=mkTab(sv))
  sv=sv%>%group_by(lid)%>%nest()%>%mutate(cnts=map(data,mkTab))%>%unnest(cnts)%>%select(-data)
  sv=sv[c("lid",ord)]
  # sv
  sv[is.na(sv)]=0
  D=left_join(D,sv)
  
  attr(D,"n")=n
  attr(D,"topv")=t
  attr(D,"topp")=tp
  # sort(t[names(tp)]/tp,decreasing=TRUE)
  attr(D,"vpRatio")=sort(t[names(tp)]/tp)
  attr(D,"tots")=tots
  # attributes(D)
  D
}
