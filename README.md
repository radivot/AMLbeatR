# AMLbeatR
AMLbeatR converts Beat AML data in 
[Tyner et al, Functional genomic landscape of acute myeloid leukaemia,
*Nature* **562**, 526–531 (2018)](https://www.nature.com/articles/s41586-018-0623-z) into a list-column tibble,
described in Chapter 25 of [R4DS](https://r4ds.had.co.nz/).  

To install AMLbeatR use:  
```
devtools::install_github("radivot/AMLbeatR",subdir="AMLbeatR")
```

## Beat AML Data
To use AMLbeatR you must first download a [290 MB excel file](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx) and a [pdf](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM1_ESM.pdf) that summarizes the sheets/tables in it. Sheets/Tables 5 (clinical), 7 (variant) and 9 (expression) are used by AMLbeatR. Place the 290 MB file 41586_2018_623_MOESM3_ESM.xlsx in the folder ~/data/BeatAML where ~ is your home directory. Now run the following  and explore the tibble in the upper left panel of R Studio.

```
library(AMLbeatR)  #loads installed package AMLbeatR into memory 
mkBeatAML() #makes ~/data/BeatAML/BeatAML.Rdata containing clinical, variant, and expression tibbles
load("~/data/BeatAML/BeatAML.RData") #loads  clinical (clin), variant (v), and expression (cpm) tibbles  
d=mkListColTib(clin,v,cpm)  #make list column tibble d 
View(d)   #Rows are patients and data is a list column of dataframes where rows are multiple samples per patient.
D=d%>%unnest() #reverse nesting to get back to 672 samples being rows
View(D)
eids=names(cpm)[-c(1:2)] #451 RNAseq measurements
vids=unique(v$labId) #608 variant samples
int=intersect(eids,vids) #399 samples have RNA and DNA readouts
De=D%>%filter(lid%in%eids)
De=De%>%group_by(id)%>%nest() #411 patients have 451 RNAseqs
Dv=D%>%filter(lid%in%vids)
Dv=Dv%>%group_by(id)%>%nest() #519 patients have 608 DNAseqs
intp=intersect(De$id,Dv$id) #376 patients have RNA and DNA reads 
``` 

## Example 1
It was recently [reported](https://www.ncbi.nlm.nih.gov/pubmed/29665898) 
that Lactate Dehydrogenase levels in the plasma are higher in AML patients with FLT3-ITD mutations. Using only data 
in the tibble clin the following code shows minimal differences
![](docs/LdhFLT3C.png)

```
D=d%>%unnest()%>%filter(!is.na(FLT3c))
table(D$FLT3c,useNA="always")
D$FLT3c=c("WT","Mutant")[D$FLT3c+1]
tc=function(sz) theme_classic(base_size=sz);
D%>%ggplot(aes(x=FLT3c,y=LDH))+scale_y_log10()+geom_boxplot(alpha=0)+geom_jitter(width=.15)+
  ylab("Lactate Dehydrogenase")+tc(14)
ggsave("~/Results/AML/LdhFLT3C.png",width=4,height=3)
wilcox.test(D$LDH~D$FLT3c) #P=0.06
t.test(log10(D$LDH)~D$FLT3c) #P=0.17
``` 







