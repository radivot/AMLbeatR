#'Add Mutations to Clinical Data 
#'
#'Adds mutation information from the variance dataframe v to the genetics list column of clin 
#'
#'@param d  Clinical dataframe made by tidyClin
#'@param v     Variant dataframe made by mkBeatAML
#'@param n     Top n most frequent mutations are encodes in the muts field
#'@return  d with additional columns od muts and vafs separated by / 
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
#'@importFrom tidyr nest unnest
#'@importFrom stringr str_c
#'@importFrom purrr map map_chr map_dbl

muts<-function(d,v,n=10)  {
  labId=lid=vaf=t_vaf=symbol=data=NULL
  
  # library(tidyverse)
  # library(AMLbeatR)
  # load("~/data/BeatAML/BeatAML.RData")
  # (d=tidyClin(clin)) # collected from 562 patients (rows in d)
  # n=10
  
  
  (t=sort(table(v$symbol),decreasing=T)[1:n])
  (top=names(t))
  sv=v%>%select(lid=labId,vaf=t_vaf,sym=symbol)%>%filter(sym%in%top)
  sv=sv%>%group_by(lid,sym)%>%summarize(vaf=max(vaf,na.rm=T))
  sv=sv%>%group_by(lid)%>%nest()
  sv=sv%>%mutate(muts=map_chr(data,function(x) str_c(x$sym,collapse="/")),
                 vafs=map_chr(data,function(x) str_c(round(x$vaf,1),collapse="/"))
                 )%>%select(-data)
  # sort(table(sv$State))
  # sv
  D=d%>%unnest()
  D=left_join(D,sv)
  D
}
