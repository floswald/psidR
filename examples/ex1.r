#' 
#' ## Start psidR example
#' 
#' ## make reproducible family data sets for 2 years
#' ## variables are: family income (Money) and age 
#' 
#' # suppose there are N individuals in year 1 and year 2. zero attrition. 
#' 
#' N <- 10
#' 
#' fam <- data.frame(int85 = 1:N,int86=sample(1:N),Money85=rlnorm(n=N,10,1),age85=sample(20:80,size=N,replace=T))
#' fam$Money86 <- fam$Money85+rnorm(N,500,30)
#' fam$age86 <- fam$age85+1
#' fam
#' 
#' # separate into data.frames.
#' # you would download files like those two:
#' fam1985 <- subset(fam,select = c(int85,Money85,age85))
#' fam1986 <- subset(fam,select = c(int86,Money86,age86))
#' names(fam1985)[1] <- "V11102"	# assign correct PSID varname of "family interview 1985"
#' names(fam1986)[1] <- "V12502"
#' 
#' 
#' # Individual index file
#' # needs to have a unique person number (ER30001) and an indicator for whether from core etc, as well as the interview number for each year
#' # for sake of illustration, let's take 2N people (i.e. N are neither in year1 nor year2)
#' IND2009ER <- data.frame(ER30001=sample((2*N):(4*N),size=2*N),ER30002=sample(1:(2*N),size=2*N))
#' 
#' # if a person is observed, they have an interview number in both years. otherwise zero. randomly allocate persons to ER30001.
#' tmp <- rbind(fam[,1:2],data.frame(int85=rep(0,N),int86=rep(0,N)))
#' 
#' IND2009ER <- cbind(IND2009ER,tmp[sample(1:(2*N)),])
#' names(IND2009ER)[3:4] <- c("ER30463","ER30498")
#' 
#' # also need relationship to head in each year in the index
#' IND2009ER$ER30465 <- sample(c(10,20),prob=c(0.5,0.5),size=2*N,replace=T)	# 50% prob of being head in year1
#' IND2009ER$ER30500 <- sample(c(10,20),prob=c(0.9,0.1),size=2*N,replace=T)
#' IND2009ER
#' 
#' # create a temporary datadir
#' my.dir <- tempdir()
#' # save those in the datadir
#' save(fam1985,file=paste0(my.dir,"/FAM1985ER.rda"))
#' save(fam1986,file=paste0(my.dir,"/FAM1986ER.RData"))	# notice different R formats admissible
#' save(IND2009ER,file=paste0(my.dir,"/IND2009ER.RData"))
#' 
#' # now famvars
#' famvars <- data.frame(year=c(1985,1986),money=c("Money85","Money86"),age=c("age85","age86"))
#' 
#' # call the builder
#' # need to set core==FALSE because person numbering indicates that all ids<2931 are not core. 
#' # set heads to FALSE to have a clear count
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,core=FALSE,heads=FALSE,verbose=TRUE)	
#' print(d$data[order(pid)],nrow=Inf)	# check the age column
#' 
#' # see what happens if we drop non-heads
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,core=FALSE,heads=TRUE)	
#' print(d$data[order(pid)],nrow=Inf)
#' 
#' file.remove(paste0(my.dir,"/FAM1985ER.rda"),paste0(my.dir,"/FAM1986ER.RData"),paste0(my.dir,"/IND2009ER.RData"))
#' 
#' # END psidR example
#' 
