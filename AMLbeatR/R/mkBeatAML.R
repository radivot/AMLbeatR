#'Make Beat AML R Binary
#'
#'Converts sheets in the 24-Table 290-MB Supplementary Tables excel file found
#'in Supplementary Information of J.W. Tyner et al. (Nature 562 2018) into a list-column tibble.
#'@param beatHome  Location of 290 MB Beat AML excel file
#'@param beatFile  Name of 290 MB Beat AML excel file. Default is 41586_2018_623_MOESM3_ESM.xlsx
#'@param outFile  Name of binary output file, which will be placed in \code{beatHome}.
#'@return  Tibble (big data file) saved as d in \code{outFile}.
#'@note  Inspired by Chapter 25 in R4DS, with countries as patients.
#'@examples
#' library(AMLbeatR)
#' mkBeatTib()
#'@name mkBeatTib
#'@export
#'@import tibble
#'@import stringi
#'@import dplyr
#'@import labelled

mkBeatTib<-function(beatHome="~/data/BeatAML",beatFile="41586_2018_623_MOESM3_ESM.xlsx",outFile="beatAML")  {
  library(dplyr)
  outFile="beatAML"
  beatFile="beatAML"
  beatHome="~/data/BeatAML"
  beatFile="41586_2018_623_MOESM3_ESM.xlsx"

  race=sex=age=casenum=NULL

  outF=file.path(beatHome,paste0(outFile,".RData"))
  inF=file.path(beatHome,beatFile)

  wesIDs=read.xlsx(inF,sheet=12) #622
  rnaIDs=read.xlsx(inF,sheet=13) # 451
  cpm=read.xlsx(inF,sheet=9) #counts per million reads
  v=read.xlsx(inF,sheet=7)
  clin=read.xlsx(inF,sheet=5)
  save(d,clin,v,cpm,wesIDs,rnaIDs,file=outF)
  cat("Beat AML data has been written to: ",outF,"\n")
}
