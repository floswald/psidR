

#' build.panel: Build PSID panel data set
#' 
#' @description Builds a panel data set in wide format with id variables \code{personID} and \code{period} from individual PSID family files.
#' @details 
#' takes desired variables from family files for specified years in folder \code{datadir} and merges using the id information in \code{IND2011ER.xyz}, which must be in the same directory. Note that only one IND file may be present in the directory (each PSID shipping comes with a new IND file).
#' There is an option to directly download the data from the PSID server to folder \code{datadir} or \code{tmpdir}.
#' The user can change subsetting criteria as well as sample designs. 
#' Merge: the variables \code{interview number} in each family file map to 
#' the \code{interview number} variable of a given year in the individual file. Run \code{example(build.panel)} for a demonstration.
#' Accepted input data are stata format .dta, .csv files or R data formats .rda and RData. Similar in usage to stata module \code{psiduse}.
#' @param datadir either \code{NULL}, in which case saves to tmpdir or path to directory containing family files ("FAMyyyy.xyz") and individual file ("IND2009ER.xyz") in admissible formats .xyz. Admissible are .dta, .csv, .RData, .rda. Please follow naming convention.
#' @param fam.vars data.frame of variable to retrieve from family files. see example for required format.
#' @param ind.vars data.frame of variables to get from individual file. In almost all cases this will be the type of survey weights you want to use. don't include id variables ER30001 and ER30002.
#' @param SAScii logical TRUE if you want to directly download data into Rda format (no dependency on STATA/SAS/SPSS). may take a long time.
#' @param heads.only logical TRUE if user wants household heads only. if FALSE, data contains a row with value of "relation to head" variable.
#' @param core logical TRUE if user wants core sample only. if FALSE, data will oversample poverty sample.
#' @param design either character "balanced" or "all" or integer. "Balanced" means only individuals who appear in each wave are considered. "All" means all are taken. An integer value stands for minimum consecutive years of participation, i.e. design=3 means present in at least 3 consecutive waves.
#' @param verbose logical TRUE if you want verbose output.
#' @import SAScii RCurl data.table
#' @return
#' \item{data}{resulting \code{data.table}. the variable \code{pid} is the unique person identifier, constructed from ID1968 and pernum.}
#' \item{dict}{data dictionary if stata data was supplied, NULL else}
#' @export
#' @examples 
#' \dontrun{
#' # specify variables from family files you want
#' 
#' myvars <- data.frame(year=c(2001,2003),
#'                      house.value=c("ER17044","ER21043"),
#'                      total.income=c("ER20456","ER24099"),
#'                      education=c("ER20457","ER24148"))
#'                      
#' # specify variables from individual index file       
#'                
#' indvars = data.frame(year=c(2001,2003),
#'                      longitud.wgt=c("ER33637","ER33740"))
#' 
#' # call builder
#' # mydir is a directory that contains FAM2001ER.dta, 
#' # FAM2003ER.dta and IND2011ER.dta
#' 
#' # default
#' d <- build.panel(datadir=mydir,
#'                 fam.vars=myvars,
#'                 ind.vars=indvars)
#'                 
#' # also non-heads               
#' d <- build.panel(datadir=mydir,
#'                 fam.vars=myvars,
#'                 ind.vars=indvars,
#'                 heads.only=FALSE)
#'                 
#' # non-balanced panel design               
#' d <- build.panel(datadir=mydir,
#'                 fam.vars=myvars,
#'                 ind.vars=indvars,
#'                 heads.only=FALSE,
#'                 design=2) # keep if stay 2+ periods
#' } 
#' 
#' # ######################################
#' # reproducible example on artifical data. 
#' # run this with example(build.panel).
#' # ######################################
#' 
#' ## make reproducible family data sets for 2 years
#' ## variables are: family income (Money) and age 
#' 
#' # suppose there are N individuals in year 1 and year 2. 
#' # zero attrition. 
#' 
#' N <- 10
#' 
#' fam <- data.frame(int85 = 1:N,int86=sample(1:N),
#'                  Money85=rlnorm(n=N,10,1),
#'                  age85=sample(20:80,size=N,replace=TRUE))
#' fam$Money86 <- fam$Money85+rnorm(N,500,30)
#' fam$age86 <- fam$age85+1
#' fam
#' 
#' # separate into data.frames.
#' # you would download files like those two:
#' fam1985 <- subset(fam,select = c(int85,Money85,age85))
#' fam1986 <- subset(fam,select = c(int86,Money86,age86))
#' 
#' # assign correct PSID varname of "family interview 1985"
#' names(fam1985)[1] <- "V11102"	
#' names(fam1986)[1] <- "V12502"
#' 
#' 
#' # construct an Individual index file: that would be IND2009ER
#' # needs to have a unique person number (ER30001) 
#' # and an indicator for whether from core etc, 
#' # as well as the interview number for each year
#' # 
#' # for sake of illustration, suppose the PSID has a total
#' # of 2N people (i.e. N are neither in year1 nor year2, 
#' # but in some other years)
#' IND2009ER <- data.frame(ER30001=sample((2*N):(4*N),size=2*N),
#'                         ER30002=sample(1:(2*N),size=2*N))
#' 
#' # if a person is observed, they have an interview number 
#' # in both years. if not observed, it's zero. 
#' # randomly allocate persons to ER30001.
#' tmp <- rbind(fam[,1:2],data.frame(int85=rep(0,N),int86=rep(0,N)))
#' 
#' IND2009ER <- cbind(IND2009ER,tmp[sample(1:(2*N)),])
#' names(IND2009ER)[3:4] <- c("ER30463","ER30498")
#' 
#' # also need relationship to head in each year in the index
#' # 50% prob of being head in year1
#' IND2009ER$ER30465 <- sample(c(10,20),prob=c(0.5,0.5),
#'                             size=2*N,replace=TRUE)	
#' IND2009ER$ER30500 <- sample(c(10,20),prob=c(0.9,0.1),
#'                            size=2*N,replace=TRUE)
#' # and a survey weight
#' IND2009ER$ER30497 <- runif(20)
#' IND2009ER$ER30534 <- runif(20)
#' IND2009ER
#' 
#' # setup the ind.vars data.frame
#' indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
#' 
#' # create a temporary datadir
#' my.dir <- tempdir()
#' # save those in the datadir
#' # notice different R formats admissible
#' save(fam1985,file=paste0(my.dir,"/FAM1985ER.rda"))
#' save(fam1986,file=paste0(my.dir,"/FAM1986ER.RData"))	
#' save(IND2009ER,file=paste0(my.dir,"/IND2009ER.RData"))
#' 
#' # now famvars
#' famvars <- data.frame(year=c(1985,1986),
#'                       money=c("Money85","Money86"),
#'                       age=c("age85","age86"))
#' 
#' # call the builder
#' # need to set core==FALSE because person numbering indicates
#' # that all ids<2931 are not core. 
#' # set heads to FALSE to have a clear count. 
#' # data will contain column "relation.head" holding the relationship code.
#' 
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,core=FALSE,
#'                  heads=FALSE,verbose=TRUE)	
#'
#' # notice: all 2*N individuals are present
#' print(d$data[order(pid)],nrow=Inf)	# check the age column
#' 
#' # see what happens if we drop non-heads
#' # only the ones who are heads in BOTH years 
#' # are present (since design='balanced' by default)
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,core=FALSE,
#'                  heads=TRUE,verbose=FALSE)	
#' print(d$data[order(pid)],nrow=Inf)
#' 
#' # change sample design to "all": 
#' # we'll keep individuals if they are head in one year,
#' # and drop in the other
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,core=FALSE,heads=TRUE,
#'                  verbose=FALSE,design="all")	
#' print(d$data[order(pid)],nrow=Inf)
#' 
#' file.remove(paste0(my.dir,"/FAM1985ER.rda"),
#'             paste0(my.dir,"/FAM1986ER.RData"),
#'             paste0(my.dir,"/IND2009ER.RData"))
#' 
#' # END psidR example
#' 
#' # #####################################################################
#' # Please go to https://github.com/floswald/psidR for more example usage
#' # #####################################################################
build.panel <- function(datadir=NULL,fam.vars,ind.vars=NULL,SAScii=FALSE,heads.only=TRUE,core=TRUE,design="balanced",verbose=FALSE){
	

	# locally bind all variables to be used in a data.table

	interview <- headyes <- .SD <- fam.interview <- ind.interview <- ind.head <- ER30001 <- ind.head.num <- pid <- ID1968 <- pernum <- isna <- present <- always <- enough <- NULL

	stopifnot(is.numeric(fam.vars$year))
	years <- fam.vars$year

	s <- .Platform$file.sep
	
	# sort out data storage: either to tmp or datadir
	if (is.null(datadir)){
		datadir <- tempdir()	
		datadir <- paste0(datadir,s)
	} else {
		if (substr(datadir,nchar(datadir),nchar(datadir))!=s) datadir <- paste0(datadir,s)
	}

	# data acquisition
	# ----------------

	if (SAScii){
    	confirm <- readline("This can take several hours/days to download.\n want to go ahead? give me 'yes' or 'no'.")
		if (confirm=="yes"){


			ftype <- "Rdata"
			
			user <- readline("please enter your PSID username: ")
			pass <- readline("please enter your PSID password: ")
			
			curl = getCurlHandle()
			curlSetOpt(cookiejar = 'cookies.txt', followlocation = TRUE, autoreferer = TRUE, curl = curl)

			html <- getURL('http://simba.isr.umich.edu/u/Login.aspx', curl = curl)

			viewstate <- as.character(sub('.*id="__VIEWSTATE" value="([0-9a-zA-Z+/=]*).*', '\\1', html))

			params <- list(
				'ctl00$ContentPlaceHolder3$Login1$UserName'    = paste(user),
				'ctl00$ContentPlaceHolder3$Login1$Password'    = paste(pass),
				'ctl00$ContentPlaceHolder3$Login1$LoginButton' = 'Log In',
				'__VIEWSTATE'                                  = viewstate
				)
				
			family    <- data.frame(year = c( 1968:1997 , seq( 1999 , 2011 , 2 ) ),file = c( 1056 , 1058:1082 , 1047:1051 , 1040 , 1052 , 1132 , 1139 , 1152  , 1156 ))
			psidFiles <- data.frame(year=c(family[family$year %in% years,]$year,"2011" ),file=c(family[family$year %in% years,]$file, 1053))

			for ( i in seq( nrow(psidFiles ) -1 )) get.psid( psidFiles[ i , 'file' ] ,name= paste0(datadir, "FAM" , psidFiles[ i , 'year' ], "ER") , params , curl )
			get.psid( psidFiles[ nrow(psidFiles ) , 'file' ] ,name= paste0(datadir, "IND" , psidFiles[ nrow(psidFiles ) , 'year' ], "ER") , params , curl )

			cat('finished downloading files to', datadir,'. continuing now to build the dataset.\n')
						
		} else if (confirm=="no") {
			break
		}
	}


	# work out file types
	# -------------------

	# figure out filestypes in datadir
	l <- list.files(datadir)
	if (length(l)==0) stop('there is something wrong with the data directory. please check path')
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "dta") ftype   <- "stata"
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "rda") ftype   <- "Rdata"
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "RData") ftype <- "Rdata"
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "csv") ftype   <- "csv"

	if (verbose) cat('psidR: loading data\n')
	if (ftype=="stata"){
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE))
        fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		tmp <- grep("IND",l,value=TRUE)
		if (length(tmp)>1) {
		    warning(cat("Warning: you have more than one IND file in your datadir.\nI take the last one:",tail(tmp,1),"\n"))
			ind.file <- paste0(datadir,tail(tmp,1))	# needs to be updated with next data delivery.
		} else {
			ind.file <- paste0(datadir,grep("IND",l,value=TRUE))	# needs to be updated with next data delivery.
		}
		ind      <- read.dta(file=ind.file)
		ind.dict <- data.frame(code=names(ind),label=attr(ind,"var.labels"))
		ind      <- data.table(ind)
	} else if (ftype=="Rdata") {
		# data downloaded directly into a dataframe
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE))
		fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		tmp <- grep("IND",l,value=TRUE)
		if (length(tmp)>1) {
		    warning(cat("Warning: you have more than one IND file in your datadir.\nI take the last one:",tail(tmp,1),"\n"))
			ind.file <- paste0(datadir,tail(tmp,1))	# needs to be updated with next data delivery.
		} else {
			ind.file <- paste0(datadir,grep("IND",l,value=TRUE))	# needs to be updated with next data delivery.
		}
		tmp.env  <- new.env()
		load(file=ind.file,envir=tmp.env)
		ind      <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
		ind.dict <- NULL
		ind      <- data.table(ind)
	} else if (ftype=="csv") {
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE))
		fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		tmp <- grep("IND",l,value=TRUE)
		if (length(tmp)>1) {
		    warning(cat("Warning: you have more than one IND file in your datadir.\nI take the last one:",tail(tmp,1),"\n"))
			ind.file <- paste0(datadir,tail(tmp,1))	# needs to be updated with next data delivery.
		} else {
			ind.file <- paste0(datadir,grep("IND",l,value=TRUE))	# needs to be updated with next data delivery.
		}
		ind      <- fread(input=ind.file)
		ind.dict <- NULL
	}
	
	if (verbose){
		cat('loaded individual file:',ind.file,'\n')
		cat('total memory load in MB:\n')
		vvs = ceiling(object.size(ind)/1024^2)
		print(as.numeric(vvs))
	}

	# data dictionaries
	fam.dicts <- vector("list",length(years))

	# output data.tables
	datas <- vector("list",length(years))

	# convert fam.vars to data.table
	stopifnot(is.data.frame(fam.vars))
	fam.vars <- data.table(fam.vars)
	fam.vars <- copy(fam.vars[,lapply(.SD,make.char)])
	setkey(fam.vars,year)
	
	# convert ind.vars to data.table if not null
	if (!is.null(ind.vars)){
		stopifnot(is.data.frame(ind.vars))
		ind.vars <- data.table(ind.vars)
		ind.vars <- copy(ind.vars[,lapply(.SD,make.char)])
		setkey(ind.vars,year)
	}


	# make an index of interview numbers for each year
	ids <- makeids()
	if (verbose){
		cat('here is the list of hardcoded PSID variables\n')
		cat('The merge is based on equal values in ind.interview and fam.interview\n')
		print(ids)
	}
	
	# which vars to keep from ind.files?
	if (!is.null(ind.vars))	stopifnot(is.list(ind.vars))

	# add compulsory vars to fam.vars
	fam.vars[,interview := ids[fam.vars][,fam.interview]]

	# loop over years
	for (iy in 1:length(years)){
		
		if (verbose) {
			cat('=============================================\n')
			cat('currently working on data for year',years[iy],'\n')
		}

		# keeping only relevant columns from individual file
		# subset for core sample and heads only if requested.
		curr <- ids[list(years[iy])]
		ind.subsetter <- as.character(curr[,list(ind.interview,ind.head)])	# keep from ind file
		def.subsetter <- c("ER30001","ER30002")	# must keep those in all years
		yind <- copy(ind[,c(def.subsetter,unique(c(ind.subsetter,as.character(ind.vars[list(years[iy]),which(names(ind.vars)!="year"),with=FALSE])))),with=FALSE])

		if (core) {
		   n    <- nrow(yind)
		   yind <- copy(yind[ER30001>2930])	# individuals 1-2930 are from poor sample
		   if (verbose){
			   cat('full',years[iy],'sample has',n,'obs\n')
			   cat('dropping non-core individuals leaves',nrow(yind),'obs\n')
		   }
		   if (nrow(yind)==0) stop('you dropped all observations by selecting core only.\n This means you supplied a family file only from the poor sample. unusual.')
		}

		if (heads.only) {
		   n    <- nrow(yind)
		   yind <- yind[,headyes := yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num]]
		   yind <- copy(yind[headyes==TRUE])
		   if (verbose){
			   cat('dropping non-heads leaves',nrow(yind),'obs\n')
		   }
		   yind[,c(curr[,ind.head],"headyes") := NULL]
			# set names on individual index
			setnames(yind,c("ID1968","pernum","interview",names(ind.vars)[-1]))
		} else {
			# set names on individual index
			setnames(yind,c("ID1968","pernum","interview","relation.head",names(ind.vars)[-1]))
		}
		yind[,pid := ID1968*1000 + pernum]	# unique person identifier
		setkey(yind,interview)

		# bring in family files, subset them
		# load data for current year, make data dictionary for subsets and save data as data.table
		if (ftype=="stata") {
			tmp             <- read.dta(file=fam.dat[iy])
			fam.dicts[[iy]] <- data.frame(code=names(tmp),label=attr(tmp,"var.labels"))
			tmp             <- data.table(tmp)
		} else if (ftype=="csv") {
			tmp <- fread(input=fam.dat[iy])
			fam.dicts[[iy]] <- NULL
		} else if (ftype=="Rdata") {
			rm(list=ls(envir=tmp.env),envir=tmp.env)
		   	load(file=fam.dat[iy],envir=tmp.env)
			tmp             <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
			fam.dicts[[iy]] <- NULL
			tmp             <- data.table(tmp)
		}

		if (verbose){
			cat('loaded family file:',fam.dat[iy],'\n')
			cat('current memory load in MB:\n')
			vs = ceiling(object.size(tmp))
			print(vs,units="Mb")
		}
	

		# add vars from from file that the user requested.
		curvars <- fam.vars[list(years[iy]),which(names(fam.vars)!="year"),with=FALSE]
		curnames <- names(curvars)
		# current set of variables
		# caution if there are specified NAs
		if (curvars[,any(is.na(.SD))]) {
			na      <- curvars[,which(is.na(.SD))]
			codes   <- as.character(curvars)
			nanames <- curnames[na]
			tmp     <- copy(tmp[,codes[-na],with=FALSE])
			tmp[,nanames := NA_real_,with=FALSE]
			setnames(tmp,c(curnames[-na],nanames))
			setkey(tmp,interview)
		} else {
			codes <- as.character(curvars)
			tmp   <- copy(tmp[,codes,with=FALSE])
			setnames(tmp,curnames)
			setkey(tmp,interview)
		}
		
		# merge
		m <- copy(tmp[yind])
		m[,year := years[iy] ]
	
		# note: a person who does not respond in wave x has an interview number in that wave, but NAs in the family file variables. remove does records.
		idx <- which(!is.na(unlist(fam.vars[list(years[iy])][,curnames,with=FALSE])))[1]	# index of first non NA variable
		m[,isna := is.na(m[,curnames[idx],with=FALSE])]
		m <- copy(m[isna == FALSE])
		m[,isna := NULL]
		# all remaining NAs are NAs which the user knows about and actually requested when specifying fam.vars
		if (iy>1)	setcolorder(m,names(datas[[1]]))
		datas[[iy]] <- copy(m)


	}
	data2 <- rbindlist(datas)	# glue together
	rm(datas)

	# design
	# keep all
	# keep only obs who show up in each wave: balanced
	# keep only obs who show up at least in g consecutive waves

	data2[,present := length(year), by=pid]
	if (design == "balanced"){
		n <- nrow(data2)
		data2[,always := max(present) == length(years),by=pid]
		data2 <- copy(data2[always==TRUE])
		data2[,always := NULL]
		if (verbose) cat("balanced design reduces sample from",n,"to",nrow(data2),"\n")
	} else if (is.numeric(design)){
		n <- nrow(data2)
		data2[,enough := max(present) >= design,by=pid]
		data2 <- copy(data2[enough==TRUE])
		data2[,enough := NULL]
		if (verbose) cat("design choice reduces sample from",n,"to",nrow(data2),"\n")
	} else if (design=="all"){
		# do nothing
	}
	data2[,present := NULL]

	#     setkey(datas,pid,year)
	#     rm(ind)
	out <- list(data=data2,dict=fam.dicts)
	if (verbose){
		cat('\n\n\nEnd of build.panel(verbose=TRUE)\n\n\n')
		cat('==================================\n')
	}

	return(out)
}








	  
	  
	  
	  








