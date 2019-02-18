#'Make Beat AML R Binary
#'
#'Converts sheets 5,7, 9, 12 1nd 13 in the Supplementary Tables excel file 
#'of J.W. Tyner et al. (Nature 562 2018) into tibbles in an R binary file.
#'
#'@param clin  Clinical tibble made by mkBeatAML
#'@param v  Variant tibble made by mkBeatAML 
#'@param cpm  Counts per million reads tibble made by mkBeatAML
#'@return  Tibble with patients as rows and lab ids nested within a tibbles using list columns
#'@author Tom Radivoyevitch
#'@note  Inspired by Chapter 25 in R4DS, with countries as patients.
#'@examples
#' library(AMLbeatR)
#' mkBeatAML()
#'@name mkListColTib
#'@export
#'@import openxlsx dplyr tidyr stringr

mkListColTib<-function(clin=clin,v=v,cpm=cpm)  {
  sex=PatientId=id=LabId=inferred_sex=age=ageAtDiagnosis=inferred_ethnicity=race=NULL
  surv=overallSurvival=status=vitalStatus=t=timeOfSampleCollectionRelativeToInclusion=NULL
  LDH=FLT3ITD=FLT3c=NPM1=NPM1c=tis=specimenType=rnaSeq=exomeSeq=trt=mostRecentTreatmentType=NULL
  library(AMLbeatR)
  load("~/data/BeatAML/BeatAML.RData")
  d=clin%>%select(id=PatientId,lid=LabId,sex=inferred_sex,age=ageAtDiagnosis,race=inferred_ethnicity,
                  surv=overallSurvival,status=vitalStatus,
                  LDH,tis=specimenType,rnaSeq,exomeSeq,trt=mostRecentTreatmentType,t=timeOfSampleCollectionRelativeToInclusion,
                  FLT3c="FLT3-ITD",NPM1c=NPM1)%>%
                  mutate(status=ifelse(status=="Dead",1,ifelse(status=="Alive",0,NA)),
                        NPM1c=ifelse(NPM1c=="positive",1,ifelse(NPM1c=="negative",0,NA)),
                        FLT3c=ifelse(FLT3c=="positive",1,ifelse(FLT3c=="negative",0,NA)),
                         rna=ifelse(rnaSeq=="y",1,0),dna=ifelse(exomeSeq=="y",1,0))%>%select(-rnaSeq,-exomeSeq)
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
  # d
  # table(d$trt)
  d%>%group_by(id,sex,age,race,surv,status)%>%nest()
}
