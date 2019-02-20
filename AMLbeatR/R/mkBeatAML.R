#'Make Beat AML R Binary
#'
#'Converts sheets 5, 7, and 9 in the Supplementary Tables excel file 
#'of J.W. Tyner et al. (Nature 562 2018) into dataframes clin, v and cpm in an R binary file.
#'
#'@param beatHome  Location of Beat AML excel file.
#'@param beatFile  Name of Beat AML file. Default is 41586_2018_623_MOESM3_ESM.xlsx.
#'@param outFile  Name of binary output file, which will be placed in \code{beatHome}.
#'@return  Nothing returned. The dataframes v (WES variants), cpm (RNAseq counts per million) and clin are in \code{outFile}.
#'@author Tom Radivoyevitch
#'@note  The goal here is simply to make reading the files into R faster via a binary file. Slapped on 
#'to clin are DNAseq and RNAseq TRUE/FALSE columns via sheets 12 and 13. This is to replace
#'exomeSeq and rnaSeq clin columns.  
#'@examples
#' library(AMLbeatR)
#' mkBeatAML()
#' load("~/data/BeatAML/BeatAML.RData")
#' (vids=sort(unique(v$labId))) #608 have at least one variant
#' (eids=names(cpm)[-c(1:2)]) #451 RNAseq measurements go 12-xx to 16-xx
#' (cvids=clin$LabId[clin$DNAseq]) #622 had DNAseq done on them =>14 AMLs with no variants
#' (int=intersect(eids,cvids)) #405 samples have RNA and DNA readouts (group of interest)
#' (int1=intersect(eids,vids)) #399 => 6 with both had no variants 

#'@name mkBeatAML
#'@export
#'@importFrom openxlsx read.xlsx

mkBeatAML<-function(beatHome="~/data/BeatAML",beatFile="41586_2018_623_MOESM3_ESM.xlsx",outFile="beatAML")  {
  # library(dplyr)
  # library(openxlsx)
  # outFile="beatAML"
  # beatFile="beatAML"
  # beatHome="~/data/BeatAML"
  # beatFile="41586_2018_623_MOESM3_ESM.xlsx"


  outF=file.path(beatHome,paste0(outFile,".RData"))
  inF=file.path(beatHome,beatFile)

  wesIDs=read.xlsx(inF,sheet=12) #622
  rnaIDs=read.xlsx(inF,sheet=13) # 451
  cpm=read.xlsx(inF,sheet=9) #counts per million reads
  v=read.xlsx(inF,sheet=7)
  clin=read.xlsx(inF,sheet=5)
  # u=unique(wesIDs$labId) #622=>all are unique, so skip SeqID
  # u=unique(clin$LabId) #672=>all clin lids are unique
  clin$DNAseq=FALSE
  clin$RNAseq=FALSE
  rownames(clin)=clin$LabId
  clin[wesIDs$labId,"DNAseq"]=TRUE
  clin[rnaIDs$labId,"RNAseq"]=TRUE
  # names(clin)
  save(clin,v,cpm,file=outF)
  cat("Beat AML data has been written to: ",outF,"\n")
}
