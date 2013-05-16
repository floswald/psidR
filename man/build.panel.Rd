\name{build.panel}
\alias{build.panel}
\title{Build PSID panel data set}
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
  in the individual file. Accepted input data are stata
  format .dta, .csv files or R data formats .rda and RData.
  Similar in usage to stata module \code{psiduse}.
}
\examples{
\dontrun{
fam.vars = data.frame(year=c(2001,2003),age=c("ER17013","ER21017"),house.value=c("ER17044","ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
fam.vars2 = data.frame(year=c(1986,1987),house.value=c("V12524","V13725"),total.income=c("V13623","V14670"),education=c("V13640","V14687"))
famvars3 <- data.frame(year=c(1985,1986),faminc=c("V12371","V13623"),house.value=c("V11125","V12524"),educ=c("V12400","V13640"))
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all")
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars2,design="balanced")
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars2,heads.only=FALSE,core=FALSE)
d <- build.panel(datadir="~/datasets/PSID/psidR-test/",famvars3)
 ##
 ## you can specify if a variable is missing in some years
 ##
fam.vars = data.frame(year=c(2001,2003),house.value=c(NA,"ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all")
}
}

