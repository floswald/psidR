



test_that("check testPSID function", {
  
    n=15
    attr = 2
    td = testPSID(N=n,N.attr=attr)
    
    expect_true(all(names(td) == c("famvars1985","famvars1986", "IND2009ER")))
    expect_that(td, is_a("list"))
    expect_true(nrow(td$famvars1985) == n)
    expect_true(nrow(td$famvars1986) == n)
    expect_true(all(unlist(lapply(td,function(x) !is.null(x)))))

    expect_true(nrow( subset(td$IND2009ER,ER30498!=0) ) == n-attr)
  
} )


test_that("check build.panel on fake data", {

    n=150
    attr = 10
    td = testPSID(N=n,N.attr=attr)
    fam1985 <- copy(td$famvars1985)
    fam1986 <- copy(td$famvars1986)
    IND2009ER <- copy(td$IND2009ER)

    # create a temporary datadir
    my.dir <- tempdir()
    # save those in the datadir
    save(fam1985,file=paste0(my.dir,"/FAM1985ER.rda"))
    save(fam1986,file=paste0(my.dir,"/FAM1986ER.RData"))
    save(IND2009ER,file=paste0(my.dir,"/IND2009ER.RData"))

    famvars <- data.frame(year=c(1985,1986),
                          money=c("Money85","Money86"),
                          age=c("age85","age86"))
    
    # and ind.vars 
    indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
    d <- build.panel(datadir=my.dir,fam.vars=famvars,
                     ind.vars=indvars,core=FALSE,
                     heads.only=FALSE,verbose=FALSE,design="all")   

    expect_true(all(names(d) == c("data","dict")))

    dd = d$data

    # full sample design: keep all obs
    # for all attrited, there should only be a 1985 row
    attrited <- subset(IND2009ER,(ER30463!=0) & (ER30498 == 0) )
    attrited$pernum <- attrited$ER30001*1000 + attrited$ER30002
    expect_true(nrow(attrited) == attr)
    expect_true(dd[pid %in% attrited$pernum,unique(year)] == 1985)

    # check that age_t+1 = age_t + 1
    setkey(dd,pid,year)
    expect_true( all( dd[!(pid %in% attrited$pernum),list(dage=diff(age)),by=pid][,dage] == rep(1,n-attr)) )

    # check that year_t+1 = year_t + 1
    setkey(dd,pid,year)
    expect_true( all( dd[!(pid %in% attrited$pernum),list(dage=diff(year)),by=pid][,dage] == rep(1,n-attr)) )

    # balanced design: only keep people who are in both waves
    d <- build.panel(datadir=my.dir,fam.vars=famvars,
                     ind.vars=indvars,core=FALSE,
                     heads.only=FALSE,verbose=FALSE,design="balanced")   
    dd = d$data
    expect_true(any(dd[,pid] %in% attrited$pernum) == FALSE)




} )



