# AMLbeatR
The goal of this R package is to convert Beat AML data in 
Tyner et al, Functional genomic landscape of acute myeloid leukaemia,
*Nature* **562**, 526–531 (2018) into a list-column tibble 
as described in Chapter 25 of [R4DS](https://r4ds.had.co.nz/).  

To install AMLbeatR use:  
```
devtools::install_github("radivot/AMLbeatR",subdir="AMLbeatR")
```

## Beat AML Data
To use AMLbeatR you must first download a large excel file from  [Supplementary information](https://www.nature.com/articles/s41586-018-0623-z#Sec38) of the 2018 Nature paper cited above. Therein [Supplementary Information link](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM1_ESM.pdf)summarizes the sheets in the file and [Supplementary Table](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx) is the 290-MB excel data file of interest.

## Make R Binary Data
Place the large file (41586_2018_623_MOESM3_ESM.xlsx) in the folder ~/data/BeatAML, where ~ is your home directory. 

```
library(AMLbeatR)  #loads installed package AMLbeatR into memory 
mkBeatAML()        #makes data file ~/data/BeatAML/BeatAML.Rdata
```


Check the R binary using 
```
load("~/data/BeatAML/BeatAML.RData") #loads d, clinical, variant, and expression tibbles into memory 
head(d,3)   #shows top 3 rows of d. Rows are patients and list columns handle multiple samples per patient.
``` 







