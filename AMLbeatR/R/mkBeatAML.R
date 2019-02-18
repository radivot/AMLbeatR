#'Make Beat AML R Binary
#'
#'Converts sheets 5,7, 9, 12 1nd 13 in the Supplementary Tables excel file 
#'of J.W. Tyner et al. (Nature 562 2018) into tibbles in an R binary file.
#'
#'@param beatHome  Location of Beat AML excel file.
#'@param beatFile  Name of Beat AML file. Default is 41586_2018_623_MOESM3_ESM.xlsx.
#'@param outFile  Name of binary output file, which will be placed in \code{beatHome}.
#'@return  Nothing returned. Tibbles v (WES variants), cpm (RNAseq counts per million) and clin are in \code{outFile}.
#'@author Tom Radivoyevitch
#'@note  The goal here is simply to make reading the files into R faster 
#'@examples
#' library(AMLbeatR)
#' mkBeatAML()
#'@name mkBeatAML
#'@export
#'@import openxlsx dplyr

mkBeatAML<-function(beatHome="~/data/BeatAML",beatFile="41586_2018_623_MOESM3_ESM.xlsx",outFile="beatAML")  {
  # library(dplyr)
  # outFile="beatAML"
  # beatFile="beatAML"
  # beatHome="~/data/BeatAML"
  # beatFile="41586_2018_623_MOESM3_ESM.xlsx"


  outF=file.path(beatHome,paste0(outFile,".RData"))
  inF=file.path(beatHome,beatFile)

  # wesIDs=read.xlsx(inF,sheet=12) #622
  # rnaIDs=read.xlsx(inF,sheet=13) # 451
  cpm=read.xlsx(inF,sheet=9) #counts per million reads
  v=read.xlsx(inF,sheet=7)
  clin=read.xlsx(inF,sheet=5)
  save(clin,v,cpm,file=outF)
  cat("Beat AML data has been written to: ",outF,"\n")
}
