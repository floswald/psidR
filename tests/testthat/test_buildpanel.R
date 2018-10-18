

context("test build.panel")

n=1000
attr = 7
my.dir = system.file(package = "psidR","testdata")
# test data was created like this:
# td = testPSID(N=n,N.attr=attr)
# fam1985 <- copy(td$famvars1985)
# fam1986 <- copy(td$famvars1986)
# IND2015ER <- copy(td$IND2009ER)
# save(fam1985,file=paste0("inst/testdata","/FAM1985ER.rda"))
# save(fam1986,file=paste0("inst/testdata","/FAM1986ER.RData"))
# save(IND2015ER,file=paste0("inst/testdata","/IND2015ER.RData"))
# 

load(file.path(my.dir,"FAM1985ER.rda"))
load(file.path(my.dir,"FAM1986ER.RData"))
load(file.path(my.dir,"IND2015ER.RData"))



test_that("check balanced sample design", {
  
    famvars <- data.frame(year=c(1985,1986),money=c("Money85","Money86"),age=c("age85","age86"))
    
    # and ind.vars 
    indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
    d <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample=NULL,heads.only=FALSE,design="all",loglevel=DEBUG)   

    dd = d

    # full sample design: keep all obs
    # for all attrited, there should only be a 1985 row
    attrited <- subset(IND2015ER,(ER30463!=0) & (ER30498 == 0) )
    attrited$pernum <- attrited$ER30001*1000 + attrited$ER30002
    expect_true(nrow(attrited) == attr)
    expect_true(dd[pid %in% attrited$pernum,unique(year)] == 1985)

    # check that age_t+1 = age_t + 1
    dd = dd[!(pid %in% attrited$pernum)]
    setkey(dd,pid,year)
    expect_true( all( dd[,list(dage=diff(age)),by=pid][,dage] == 1 ))

    # check that year_t+1 = year_t + 1
    setkey(dd,pid,year)
    expect_true( all( dd[,list(dage=diff(year)),by=pid][,dage] == 1 ))

    # balanced design: only keep people who are in both waves
    d <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample=NULL,heads.only=FALSE,design="balanced")   
    dd = d
    expect_true(any(dd[,pid] %in% attrited$pernum) == FALSE)
    
    # not subset to heads:
    expect_true( "relation.head" %in% names(dd) )
    # check sequence numbers
    expect_true( !all( dd[,sequence == 1]))
    
    # check relationship to head
    expect_true( !all( dd[,relation.head == 10]))
} )

test_that("check subsetting to head and wife sample", {
  
  famvars <- data.frame(year=c(1985,1986),money=c("Money85","Money86"),age=c("age85","age86"))
  
  # and ind.vars 
  indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
  core <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample=NULL,heads.only=TRUE,design="all")   
  cored = core
  
  # check sequence numbers
  expect_true( all( cored[,(sequence >0) & (sequence < 21) ]))
  
  # check relationship to head
  expect_true( all( cored[,relation.head == 10]))

})

test_that("check subsetting to current heads only", {
  
  famvars <- data.frame(year=c(1985,1986),money=c("Money85","Money86"),age=c("age85","age86"))
  
  # and ind.vars 
  indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
  core <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample=NULL,current.heads.only=TRUE,design="all")   
  cored = core
  
  # check sequence numbers
  expect_true( all( cored[,sequence==1]))
  
  # check relationship to head
  expect_true( all( cored[,relation.head == 10]))
  
})

test_that("check subsetting to core/immigrant/latino", {
  
  famvars <- data.frame(year=c(1985,1986),money=c("Money85","Money86"),age=c("age85","age86"))
  
  # and ind.vars 
  indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
  src <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample="SRC",heads.only=FALSE,design="all")   
  seo <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample="SEO",heads.only=FALSE,design="all")   
  lat <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample="latino",heads.only=FALSE,design="all")   
  imm <- build.panel(datadir=my.dir,fam.vars=famvars,ind.vars=indvars,sample="immigrant",heads.only=FALSE,design="all")   
  
  # check interview numbers
  expect_true( all( src[,ID1968 < 3000 ]))
  expect_true( all( seo[,ID1968 > 5000 & ID1968 < 7000 ]))
  expect_true( all( lat[,ID1968 > 7000 & ID1968 < 9308]))
  expect_true( all( imm[,ID1968 > 3000 & ID1968 < 7000  ]))
  
} )


