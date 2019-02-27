# AMLbeatR
An R package for analyses of  Beat AML data in
[Tyner et al, *Nature* **562**, 526–531 (2018)](https://www.nature.com/articles/s41586-018-0623-z).  

To install it use:  
```
devtools::install_github("radivot/AMLbeatR",subdir="AMLbeatR")
```

## Beat AML Data
Place the 290 MB file  [41586_2018_623_MOESM3_ESM.xlsx ](https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0623-z/MediaObjects/41586_2018_623_MOESM3_ESM.xlsx) in ~/data/BeatAML where ~ is your home directory and run  

```
rm(list=ls()) 
library(tidyverse)
library(AMLbeatR)  #loads installed package AMLbeatR into memory 
mkBeatAML() #Do only once (this takes some time). Makes ~/data/BeatAML/BeatAML.Rdata 
load("~/data/BeatAML/BeatAML.RData") #loads  clinical (clin), variant (v, av), and expression (rpkm, cpm) data  
(d=tidyClin(clin)) #672 rows in clin, one for each measurement; 562 rows in d, one for each patient
(d=muts(d,v,av,n=10)) # adds muts and vafs of top n genes and point mutation counts from all variants
attributes(d) # muts() also added some summary attributes
```

## Differentially Expressed Genes
To create an excel file with 10 sheets of differentially expressed genes, one for each mutated gene, run
```
geXL(d,v,cpm,f="~/Results/AML/topGE.xlsx") #WT vs mut differential gene expression (ge), one per mutated gene
``` 

## Survival
An independent group recently [reported](https://www.nejm.org/doi/full/10.1056/NEJMoa1516192) 
that among the four combinations of DNMT3A (D) and NPM1 (N) mutations, FLT3 (F) mutations change survival 
the most when D and N are both mutated. This result is replicated nicely below.
![](docs/survALL.png)

```
###### Survival
library(survival);library(survminer)
sv=v%>%filter(symbol%in%c("FLT3","DNMT3A","NPM1"))%>%select(lid=labId,sym=symbol,t_vaf)
sv=sv%>%group_by(lid,sym)%>%summarize(vaf=max(t_vaf,na.rm=T))
(sv=sv%>%mutate(sym=str_sub(sym,1,1)))
sv=sv%>%group_by(lid)%>%nest()
sv=sv%>%mutate(State=map(data,function(x) str_c(x$sym,collapse="")))%>%select(-data)
sv=sv%>%unnest()
D=left_join(d,sv)
D=D%>%group_by(id,surv,status)%>%nest()
getLong=function(x) x$State[str_length(x$State)==max(str_length(x$State))]
getFirst=function(x) x[1]
D=D%>%mutate(State=map(data,getLong))%>%select(-data)
D=D%>%mutate(State1=map(State,getFirst))%>%select(-State)
D=D%>%unnest()
D=D%>%mutate(surv=surv/365.25) # make in Years
D$State1[is.na(D$State1)]="WT"
D1=D%>%filter(State1%in%c("WT","F"))
D2=D%>%filter(State1%in%c("D","DF"))
D3=D%>%filter(State1%in%c("N","FN"))
D4=D%>%filter(State1%in%c("DN","DFN"))
D1$grp="None"; D2$grp="D"; D3$grp="N"; D4$grp="DN";
D=bind_rows(D1,D2,D3,D4)
labs=c("Fwt","Fmut")
D$F=labs[str_detect(D$State1,"F")+1]
D=as.data.frame(D) #fixes survplot error 
D$grp=factor(D$grp,c("None","D","N","DN"))
D$F=factor(D$F,c("Fwt","Fmut"))
fit=survfit(Surv(surv,status)~F,data=D)

sbb=theme(strip.background=element_blank())
gy=ylab("Survival Probability")
gx=xlab("Years")
svts=scale_x_continuous(breaks=c(0,2,4,6,8,10))#surv times
lg=theme(legend.margin=margin(0,0,0,0),legend.position=c(0.9,0.85))#,legend.direction="horizontal")
ggsurvplot_facet(fit,D,legend.title="",facet.by="grp",nrow=1,short.panel.labs=T,
                 xlim=c(0,10),pval=T,pval.coord=c(3,0.45))+svts+sbb+gy+gx+lg
ggsave("~/Results/AML/survALL.pdf",width=10,height=3)
ggsave("~/Results/AML/survALL.png",width=7,height=3)
```

