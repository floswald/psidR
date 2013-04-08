
you should change into ./inst of this package
---------------------------------------------

fam.vars = data.frame(year=c(2001,2003),house.value=c("ER17044","ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all")	
using small examples manually specyfing paths and files
d <- build.panel(datadir=NULL,fam.vars,fam.files=c("short2001.csv","short2003.csv"),ind.file="shortIND.csv",design="all")
d <- build.panel(datadir=NULL,fam.vars,fam.files=c("short2001.csv","short2003.csv"),ind.file="shortIND.csv",design="balanced")
d <- build.panel(datadir=NULL,fam.vars,fam.files=c("short2001.csv","./inst/short2003.csv"),ind.file="./inst/shortIND.csv",design=1,heads.only=FALSE,core=FALSE)
d <- build.panel(datadir=NULL,fam.vars,fam.files=c("stata-data/short2001.dta","stata-data/short2003.dta"),ind.file="stata-data/shortIND.dta",design="all")
d <- build.panel(datadir=NULL,fam.vars,fam.files=c("stata-data/short2001.dta","stata-data/short2003.dta"),ind.file="stata-data/shortIND.dta",design="balanced")
d <- build.panel(datadir=NULL,fam.vars,fam.files=c("stata-data/short2001.dta","stata-data/short2003.dta"),ind.file="stata-data/shortIND.dta",design=1)
 
 you can specify if a variable is missing in some years
---------------------------------------------
 
fam.vars = data.frame(year=c(2001,2003),house.value=c(NA,"ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
d <- build.panel(datadir=NULL,fam.vars,fam.files=c("./inst/short2001.csv","./inst/short2003.csv"),ind.file="./inst/shortIND.csv",design="all")
