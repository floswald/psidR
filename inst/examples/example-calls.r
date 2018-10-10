#
# examples of fam.vars declarations
#
famvars1 = data.frame(year=c(2001,2003),
                       age=c("ER17013","ER21017"),
                       house.value=c("ER17044","ER21043"),
                       total.income=c("ER20456","ER24099"),
                       education=c("ER20457","ER24148"))
famvars2 = data.frame(year=c(1986,1987),
                      house.value=c("V12524","V13725"),
                      total.income=c("V13623","V14670"),
                      education=c("V13640","V14687"))
famvars3 <- data.frame(year=c(1985,1986),
                       faminc=c("V12371","V13623"),
                       house.value=c("V11125","V12524"),
                       educ=c("V12400","V13640"))
#
# examples of ind.vars declarations
#
indvars1 = data.frame(year=c(2001,2003),longitud.wgt=c("ER33637","ER33740"))
indvars2 = data.frame(year=c(1986,1987),ind.weight=c("ER30534","ER30569"))
indvars3 <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
#
# assume your data directory is ~/datasets/PSID/fam-files
#
d1 <- build.panel(datadir="~/datasets/PSID/fam-files",
                 fam.vars=famvars1,
                 ind.vars=indvars1,
                 design="all")  
d2 <- build.panel(datadir="~/datasets/PSID/fam-files",
                  fam.vars=famvars2,
                  ind.vars=indvars2,
                  design="balanced",
                  verbose=TRUE)	
d3 <- build.panel(datadir="~/datasets/PSID/fam-files",
                  fam.vars=famvars3,
                  ind.vars=indvars3,
                  heads.only=FALSE,
                  core=FALSE)	
d4 <- build.panel(fam.vars=famvars3,SAScii=TRUE,verbose=TRUE)	
# with SAScii=TRUE data.dir not required
# 
# you can specify if a variable is missing in some years
#
fam.vars = data.frame(year=c(2001,2003),
                     house.value=c(NA,"ER21043"),
                     total.income=c("ER20456","ER24099"),
                     education=c("ER20457","ER24148"))
d <- build.panel(datadir="~/datasets/PSID/fam-files",famvars1,design="all")
# 
# suppose datadir is empty and you want to download from PSID directly 
# (caution: takes a lot of time)
#
famvars4 = data.frame(year=c(2001,2003),
                      house.value=c(NA,"ER21043"),
                      total.income=c("ER20456","ER24099"),
                      education=c("ER20457","ER24148"))
d <- build.panel(datadir="~/datasets/PSID/from-SAS/",
                 fam.vars=famvars4,design="all")	
                 # will store .rda dataframes into datadir
} 
