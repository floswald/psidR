\name{build.panel}
\alias{build.panel}
\title{build.panel: Build PSID panel data set}
\usage{
  build.panel(datadir, fam.vars, ind.vars = NULL,
    SAScii = FALSE, heads.only = TRUE, core = TRUE,
    design = "balanced", verbose = FALSE)
}
\arguments{
  \item{datadir}{directory containing family files
  ("FAMyyyy.dta") and individual file ("IND2009ER.dta") in
  Stata or other admissible formats (naming convention
  required for stata files)}

  \item{fam.vars}{data.frame of variable to retrieve from
  family files. see example for required format.}

  \item{ind.vars}{optional data.frame of non-default
  variables to get from individual file.}

  \item{SAScii}{logical TRUE if you want to directly
  download data into Rda format (no dependency on
  STATA/SAS/SPSS). may take a long time.}

  \item{heads.only}{logical TRUE if user wants household
  heads only. if FALSE, data contains a row with value of
  "relation to head" variable.}

  \item{core}{logical TRUE if user wants core sample only.
  if FALSE, data will oversample poverty sample.}

  \item{design}{either character "balanced" or "all" or
  integer. "Balanced" means only individuals who appear in
  each wave are considered. "All" means all are taken. An
  integer value stands for minimum consecutive years of
  participation, i.e. design=3 means at least 3 waves.}

  \item{verbose}{logical TRUE if you want verbose output.}
}
\value{
  \item{data}{resulting \code{data.table}. the variable
  \code{pid} is the unique person identifier, constructed
  from ID1968 and pernum.} \item{dict}{data dictionary if
  stata data was supplied, NULL else}
}
\description{
  Builds a panel data set in wide format with id variables
  \code{personID} and \code{period} from individual PSID
  family files.
}
\details{
  takes family files for specified years in folder
  \code{datadir} and merges using the id information in
  \code{ind.vars}, which must be in the same directory.
  There is an option to directly download the data from the
  PSID server to folder \code{datadir}. The user can change
  subsetting criteria as well as sample designs. If there
  are N individuals in each of T waves, the individual file
  contains NT rows. If an individual has non-response in a
  given wave, values in the family file are NA. the
  variables \code{interview number} in each family file map
  to the \code{interview number} variable of a given year
  in the individual file. Run \code{example(build.panel)}
  for a demonstration. Accepted input data are stata format
  .dta, .csv files or R data formats .rda and RData.
  Similar in usage to stata module \code{psiduse}.
}
\examples{
\dontrun{
#
# examples of fam.vars declarations
#
famvars1 = data.frame(year=c(2001,2003),age=c("ER17013","ER21017"),house.value=c("ER17044","ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
famvars2 = data.frame(year=c(1986,1987),house.value=c("V12524","V13725"),total.income=c("V13623","V14670"),education=c("V13640","V14687"))
famvars3 <- data.frame(year=c(1985,1986),faminc=c("V12371","V13623"),house.value=c("V11125","V12524"),educ=c("V12400","V13640"))
#
# assume your data directory is ~/datasets/PSID/fam-files
#
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all")
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars2,design="balanced",verbose=TRUE)
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars2,heads.only=FALSE,core=FALSE)
d <- build.panel(datadir="~/datasets/PSID/psidR-test/",famvars3)
#
# you can specify if a variable is missing in some years
#
fam.vars = data.frame(year=c(2001,2003),house.value=c(NA,"ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all")
#
# suppose datadir is empty and you want to download from PSID directly (caution: takes a lot of time)
#
fam.vars = data.frame(year=c(2001,2003),house.value=c(NA,"ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all",SAScii=TRUE)
}
# ###################################################
# Here is a reproducible example you can actually run
# ###################################################

## Start psidR example

## make reproducible family data sets for 2 years
## variables are: family income (Money) and age

# suppose there are N individuals in year 1 and year 2. zero attrition.

N <- 10

fam <- data.frame(int85 = 1:N,int86=sample(1:N),Money85=rlnorm(n=N,10,1),age85=sample(20:80,size=N,replace=TRUE))
fam$Money86 <- fam$Money85+rnorm(N,500,30)
fam$age86 <- fam$age85+1
fam

# separate into data.frames.
# you would download files like those two:
fam1985 <- subset(fam,select = c(int85,Money85,age85))
fam1986 <- subset(fam,select = c(int86,Money86,age86))
names(fam1985)[1] <- "V11102"	# assign correct PSID varname of "family interview 1985"
names(fam1986)[1] <- "V12502"


# Individual index file
# needs to have a unique person number (ER30001)
# and an indicator for whether from core etc, as well as the interview number for each year
#
# for sake of illustration, let's take 2N people (i.e. N are neither in year1 nor year2)
IND2009ER <- data.frame(ER30001=sample((2*N):(4*N),size=2*N),ER30002=sample(1:(2*N),size=2*N))

# if a person is observed, they have an interview number in both years. otherwise zero. randomly allocate persons to ER30001.
tmp <- rbind(fam[,1:2],data.frame(int85=rep(0,N),int86=rep(0,N)))

IND2009ER <- cbind(IND2009ER,tmp[sample(1:(2*N)),])
names(IND2009ER)[3:4] <- c("ER30463","ER30498")

# also need relationship to head in each year in the index
IND2009ER$ER30465 <- sample(c(10,20),prob=c(0.5,0.5),size=2*N,replace=TRUE)	# 50\% prob of being head in year1
IND2009ER$ER30500 <- sample(c(10,20),prob=c(0.9,0.1),size=2*N,replace=TRUE)
IND2009ER

# create a temporary datadir
my.dir <- tempdir()
# save those in the datadir
save(fam1985,file=paste0(my.dir,"/FAM1985ER.rda"))
save(fam1986,file=paste0(my.dir,"/FAM1986ER.RData"))	# notice different R formats admissible
save(IND2009ER,file=paste0(my.dir,"/IND2009ER.RData"))

# now famvars
famvars <- data.frame(year=c(1985,1986),money=c("Money85","Money86"),age=c("age85","age86"))

# call the builder
# need to set core==FALSE because person numbering indicates that all ids<2931 are not core.
# set heads to FALSE to have a clear count
d <- build.panel(datadir=my.dir,fam.vars=famvars,core=FALSE,heads=FALSE,verbose=TRUE)
print(d$data[order(pid)],nrow=Inf)	# check the age column

# see what happens if we drop non-heads
d <- build.panel(datadir=my.dir,fam.vars=famvars,core=FALSE,heads=TRUE)
print(d$data[order(pid)],nrow=Inf)

file.remove(paste0(my.dir,"/FAM1985ER.rda"),paste0(my.dir,"/FAM1986ER.RData"),paste0(my.dir,"/IND2009ER.RData"))

# END psidR example
}

