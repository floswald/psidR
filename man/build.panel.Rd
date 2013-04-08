\name{build.panel}
\alias{build.panel}
\title{Build PSID panel data set}
\usage{
  build.panel(datadir, fam.vars, ind.vars = NULL,
    fam.files = NULL, ind.file = NULL, heads.only = TRUE,
    core = TRUE, design = "balanced", verbose = FALSE)
}
\arguments{
  \item{datadir}{directory containing family files
  ("FAMyyyy.dta") and individual file ("IND2009ER.dta").
  fixed data format required.}

  \item{fam.vars}{data.frame of variable to retrieve from
  family files. see example for required format.}

  \item{ind.vars}{optional data.frame of non-default
  variables to get from individual file.}

  \item{fam.files}{optional character vector of file
  locations of family files}

  \item{ind.file}{optional character vector of file
  location of individual file}

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
  \itemize{ \item Default behaviour: takes family files for
  specified years in folder \code{datadir} and merges using
  the id information in \code{ind.vars}, which must be in
  the same directory. \item Alternative: location of files
  may be manually specified. \item in both: the user can
  change subsetting criteria as well as sample designs.
  \item Merge: If there are N individuals in each of T
  waves, the individual file contains NT rows. If an
  individual has non-response in a given wave, values in
  the family file are NA. the variables \code{interview
  number} in each family file map to the \code{interview
  number} variable of a given year in the individual file.
  \item Accepted input data are stata format .dta or .csv
  files. \item See examples for usage. Note that examples
  use fake data, i.e. this is not real PSID data. Data
  dictionaries do not work therefore. }
}
\examples{

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
}

