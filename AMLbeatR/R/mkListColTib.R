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
#'@import openxlsx dplyr tidyr

mkListColTib<-function(clin=clin,v=v,cpm=cpm)  {
  # d=v=clin=cpm=
  sex=PatientId=id=LabId=inferred_sex=age=ageAtDiagnosis=inferred_ethnicity=race=NULL
  surv=overallSurvival=status=vitalStatus=NULL
  LDH=FLT3ITD=NPM1=NULL
  # library(AMLbeatR)
  # load("~/data/BeatAML/BeatAML.RData")
  d=clin%>%select(id=PatientId,lid=LabId,sex=inferred_sex,age=ageAtDiagnosis,race=inferred_ethnicity,
                  surv=overallSurvival,status=vitalStatus,
                  LDH,FLT3c="FLT3-ITD",NPM1c=NPM1)%>%
                  mutate(status=ifelse(status=="Dead",1,ifelse(status=="Alive",0,NA)))%>%
    group_by(id,sex,age,race,surv,status)%>%nest()
  d
}
