#'Tidy up Clinical Data
#'
#'Organizes clin as a tibble sorted by patient ids, with one row per patient and a list column data for multiple measurements
#'
#'@param clin  Clinical dataframe made by mkBeatAML
#'@return  Tibble with one patient per row, sorted by patient ids, with lab ids et al nested in a list column called data
#'@author Tom Radivoyevitch
#'@note  Inspired by Chapter 25 in R4DS, with countries as patients.
#'@examples
#' library(AMLbeatR)
#' load("~/data/BeatAML/BeatAML.RData") # 672  specimens (rows in clin) 
#' (d=mkListColTib(clin,v,cpm)) # collected from 562 patients (rows in d)
#' (D=d%>%unnest()) 
#'@name tidyClin
#'@export
#'@import  dplyr   
#'@importFrom stringr str_detect
#'@importFrom tidyr nest
#'@importFrom purrr map map_chr map_dbl

tidyClin<-function(clin=clin)  {
  sex=PatientId=id=LabId=consensus_sex=age=ageAtDiagnosis=inferred_ethnicity=race=NULL
  surv=overallSurvival=status=vitalStatus=t=timeOfSampleCollectionRelativeToInclusion=NULL
  LDH=PBB=BMB=FLT3ITD=FLT3c=NPM1=NPM1c=tis=specimenType=RNA=DNA=data=trt=mostRecentTreatmentType=NULL
  
  library(tidyverse)
  load("~/data/BeatAML/BeatAML.RData")
  
  d=clin%>%select(id=PatientId,lid=LabId,sex=consensus_sex,age=ageAtDiagnosis,race=inferred_ethnicity,
                  surv=overallSurvival,status=vitalStatus,
                  LDH,PBB="%.Blasts.in.PB",BMB="%.Blasts.in.BM",tis=specimenType,DNA,RNA,trt=mostRecentTreatmentType,t=timeOfSampleCollectionRelativeToInclusion,
                  FLT3c="FLT3-ITD",NPM1c=NPM1)%>%
                  mutate(status=ifelse(status=="Dead",1,ifelse(status=="Alive",0,NA)),
                        NPM1c=ifelse(NPM1c=="positive",1,ifelse(NPM1c=="negative",0,NA)),
                        FLT3c=ifelse(FLT3c=="positive",1,ifelse(FLT3c=="negative",0,NA))
                        )
  d$race[which(str_detect(d$race,"White"))]="White"
  d$race[which(str_detect(d$race,"Black"))]="Black"
  d$race[which(str_detect(d$race,"Asian"))]="Other"  #races as in SEERaBomb
  d$race[which(str_detect(d$race,"Hisp"))]="Other"
  d$tis[which(str_detect(d$tis,"Bone"))]="BM"
  d$tis[which(str_detect(d$tis,"Peripheral"))]="PB"
  d$tis[which(str_detect(d$tis,"Leukapheresis"))]="PB"
  d$trt[which(str_detect(d$trt,"Transplant"))]="BMT"
  d$trt[which(str_detect(d$trt,"DLI"))]="BMT"
  d$trt[which(str_detect(d$trt,"Kinase"))]="TKI"
  d$trt[which(str_detect(d$trt,"NONE"))]="None"
  d$trt[which(str_detect(d$trt,"Unknown"))]=NA
  d$trt[which(str_detect(d$trt,"Intrathecal"))]="Other"
  d$trt[which(str_detect(d$trt,"Standard"))]="Chemo"
  d$trt[which(str_detect(d$trt,"Targeted Therapy - Other"))]="TTO"
  d$trt[which(str_detect(d$trt,"Palliative"))]="SPC"
  d$LDH=as.numeric(d$LDH)
  d$PBB=as.numeric(d$PBB)
  d$BMB=as.numeric(d$BMB)
  d=d%>%group_by(id,sex,status)%>%nest()%>%arrange(id)
  d=d%>%mutate(race=map_chr(data,function(x) x$race[!is.na(x$race)][1] ), #grab first non-missing race in patient group
            age=map_dbl(data,function(x) mean(x$age,na.rm=T) ), #fix two age at diagnosis values for 1 pt
            surv=map_dbl(data,function(x) mean(x$surv,na.rm=T) ), #fix two surv times for same pt
            n=map_dbl(data,function(x) dim(x)[1]), #fix two surv times for same pt
            data=map(data,function(x) x%>%select(-race,-age,-surv) )
             )%>%select(-data,everything())
  D=d%>%unnest()
  # Nature Absract: 672 tumour specimens (rows in D) collected from 562 patients (rows in d)
  (D=D%>%select(id,sex,race,age,surv,status,n,lid,t,tis,DNA,RNA,trt,LDH,PBB,BMB,FLT3c,NPM1c))
  # dput(names(D))
  (D=D%>%group_by(id,sex,race,age,surv,status,n)%>%nest(lid:NPM1c))
  # system.time(
  # D%>%mutate(timeNsrc=map(data,function(x) x%>%select(lid:RNA) ) ,  #### this takes 2.5 secs!!!
  #            u=map(data,function(x) x%>%select(trt) ),
  #            x=map(data,function(x) x%>%select(FLT3c,NPM1c) ),
  #            y=map(data,function(x) x%>%select(LDH) ))%>%select(-data)
  # )
  SH=c("lid", "t", "tis", "DNA", "RNA")
  U=c("trt")
  Y=c("LDH","PBB","BMB")
  X=c("FLT3c", "NPM1c")
  system.time(
  d<-D%>%mutate(timeNsrc=map(data,function(x) x[SH] ) ,  #### this takes 0.9 secs!!!
             inputs=map(data,function(x) x[U] ),
             outputs=map(data,function(x) x[Y] ),
             genetics=map(data,function(x) x[X] ))%>%select(-data)
  )
  d
  # system.time({
  #   d<-D%>%mutate(timeNsrc=map(data,function(x) x[SH]))  #### this takes 0.9 secs!!!
  #     d<-d%>%mutate(inputs=map(data,function(x) x[U] ))
  #     d<-d%>%mutate(outputs=map(data,function(x) x[Y] ))
  #     d<-d%>%mutate(genetics=map(data,function(x) x[X] ))%>%select(-data) }
  # ) #same time
  
  
  # SH=c("lid", "t", "tis", "DNA", "RNA","trt")
  # Y=c("FLT3c", "NPM1c","LDH")
  # system.time(
  # D%>%mutate(timeNsrc=map(data,function(x) x[SH] ) ,  #### this takes 0.4 secs!!!
  #            y=map(data,function(x) x[Y] ))%>%select(-data)
  # )
  # d
  d%>%unnest()
  d%>%unnest(timeNsrc)
  d%>%unnest(inputs)
  d%>%unnest(outputs)
  d%>%unnest(genetics)
  d
}
