

#' build.panel: Build PSID panel data set
#' 
#' @description Builds a panel data set in wide format with id variables \code{personID} and \code{period} from individual PSID family files.
#' @details 
#' takes family files for specified years in folder \code{datadir} and merges using the id information in \code{ind.vars}, which must be in the same directory. 
#' There is an option to directly download the data from the PSID server to folder \code{datadir}.
#' The user can change subsetting criteria as well as sample designs. If there are N individuals in each of T waves, the individual file contains NT rows. 
#' If an individual has non-response in a given wave, values in the family file are NA. the variables \code{interview number} in each family file map to 
#' the \code{interview number} variable of a given year in the individual file. Run \code{example(build.panel)} for a demonstration.
#' Accepted input data are stata format .dta, .csv files or R data formats .rda and RData. Similar in usage to stata module \code{psiduse}.
#' @param datadir directory containing family files ("FAMyyyy.dta") and individual file ("IND2009ER.dta") in Stata or other admissible formats (naming convention required for stata files)
#' @param fam.vars data.frame of variable to retrieve from family files. see example for required format.
#' @param ind.vars optional data.frame of non-default variables to get from individual file.
#' @param SAScii logical TRUE if you want to directly download data into Rda format (no dependency on STATA/SAS/SPSS). may take a long time.
#' @param heads.only logical TRUE if user wants household heads only. if FALSE, data contains a row with value of "relation to head" variable.
#' @param core logical TRUE if user wants core sample only. if FALSE, data will oversample poverty sample.
#' @param design either character "balanced" or "all" or integer. "Balanced" means only individuals who appear in each wave are considered. "All" means all are taken. An integer value stands for minimum consecutive years of participation, i.e. design=3 means at least 3 waves.
#' @param verbose logical TRUE if you want verbose output.
#' @import SAScii RCurl data.table
#' @return
#' \item{data}{resulting \code{data.table}. the variable \code{pid} is the unique person identifier, constructed from ID1968 and pernum.}
#' \item{dict}{data dictionary if stata data was supplied, NULL else}
#' @export
#' @examples \dontrun{
#' #
#' # examples of fam.vars declarations
#' #
#' famvars1 = data.frame(year=c(2001,2003),age=c("ER17013","ER21017"),house.value=c("ER17044","ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
#' famvars2 = data.frame(year=c(1986,1987),house.value=c("V12524","V13725"),total.income=c("V13623","V14670"),education=c("V13640","V14687"))
#' famvars3 <- data.frame(year=c(1985,1986),faminc=c("V12371","V13623"),house.value=c("V11125","V12524"),educ=c("V12400","V13640"))
#' #
#' # assume your data directory is ~/datasets/PSID/fam-files
#' #
#' d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all")	
#' d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars2,design="balanced",verbose=TRUE)	
#' d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars2,heads.only=FALSE,core=FALSE)	
#' d <- build.panel(datadir="~/datasets/PSID/psidR-test/",famvars3)	
#' # 
#' # you can specify if a variable is missing in some years
#' #
#' fam.vars = data.frame(year=c(2001,2003),house.value=c(NA,"ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
#' d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all")
#' # 
#' # suppose datadir is empty and you want to download from PSID directly (caution: takes a lot of time)
#' #
#' fam.vars = data.frame(year=c(2001,2003),house.value=c(NA,"ER21043"),total.income=c("ER20456","ER24099"),education=c("ER20457","ER24148"))
#' d <- build.panel(datadir="~/datasets/PSID/fam-files",fam.vars,design="all",SAScii=TRUE)
#' } 
#' # ###################################################
#' # Here is a reproducible example you can actually run
#' # ###################################################
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
#' fam <- data.frame(int85 = 1:N,int86=sample(1:N),Money85=rlnorm(n=N,10,1),age85=sample(20:80,size=N,replace=TRUE))
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
#' # needs to have a unique person number (ER30001) 
#' # and an indicator for whether from core etc, as well as the interview number for each year
#' # 
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
#' IND2009ER$ER30465 <- sample(c(10,20),prob=c(0.5,0.5),size=2*N,replace=TRUE)	# 50% prob of being head in year1
#' IND2009ER$ER30500 <- sample(c(10,20),prob=c(0.9,0.1),size=2*N,replace=TRUE)
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
build.panel <- function(datadir,fam.vars,ind.vars=NULL,SAScii=FALSE,heads.only=TRUE,core=TRUE,design="balanced",verbose=FALSE){
	
	years <- fam.vars$year

	# data acquisition
	# ----------------

	if (SAScii){
    	confirm <- readline("This can take several hours/days to download.\n want to go ahead? give me 'yes' or 'no'.")
		if (confirm=="yes"){

			if (is.null(datadir)) stop("please supply datadir argument. will store datasets there.")
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
				
			family    <- data.frame(year = c( 1968:1997 , seq( 1999 , 2009 , 2 ) ),file = c( 1056 , 1058:1082 , 1047:1051 , 1040 , 1052 , 1132 , 1139 , 1152 ))
			psidFiles <- data.frame(year=c(family[family$year %in% years,]$year,"2009" ),file=c(family[family$year %in% years,]$file, 1053))

			for ( i in seq( nrow(psidFiles ) -1 )) get.psid( psidFiles[ i , 'file' ] ,name= paste0(datadir, "/FAM" , psidFiles[ i , 'year' ], "ER") , params , curl )
			get.psid( psidFiles[ nrow(psidFiles ) , 'file' ] ,name= paste0(datadir, "/IND" , psidFiles[ nrow(psidFiles ) , 'year' ], "ER") , params , curl )

			cat('finished downloading files to', datadir,'. continuing now to build the dataset.\n')
						
		} else if (confirm=="no") {
			break
		}
	}


	# work out file types
	# -------------------

	# if last char of datadir is not "/", add it
	if (substr(datadir,nchar(datadir),nchar(datadir))!="/") datadir <- paste0(datadir,"/")

	# figure out filestypes in datadir
	l <- list.files(datadir)
	if (length(l)==0) stop('there is something wrong with the data directory. please check path')
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "dta") ftype   <- "stata"
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "rda") ftype   <- "Rdata"
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "RData") ftype <- "Rdata"
	if (tail(strsplit(l[1],"\\.")[[1]],1) == "csv") ftype   <- "csv"

	if (ftype=="stata"){
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE))
        fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		ind.file <- paste0(datadir,grep("IND",l,value=TRUE))	# needs to be updated with next data delivery.
		ind      <- read.dta(file=ind.file)
		ind.dict <- data.frame(code=names(ind),label=attr(ind,"var.labels"))
		ind      <- data.table(ind)
	} else if (ftype=="Rdata") {
		# data downloaded directly into a dataframe
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE))
		fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		ind.file <- paste0(datadir,grep("IND",l,value=TRUE))	# needs to be updated with next data delivery.
		tmp.env  <- new.env()
		load(file=ind.file,envir=tmp.env)
		ind      <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
		ind.dict <- NULL
		ind      <- data.table(ind)
	} else if (ftype=="csv") {
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE))
		fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		ind.file <- paste0(datadir,grep("IND",l,value=TRUE))	# needs to be updated with next data delivery.
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
	if (!is.data.table(fam.vars)) {
		fam.vars <- data.table(fam.vars)
		fam.vars <- copy(fam.vars[,lapply(.SD,make.char)])
		setkey(fam.vars,year)
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
		curr <- ids[.(years[iy])]
		ind.subsetter <- as.character(curr[,list(ind.interview,ind.head)])	# keep from ind file
		def.subsetter <- c("ER30001","ER30002")	# must keep those in all years
		yind <- copy(ind[,c(def.subsetter,unique(c(ind.subsetter,ind.vars[[iy]]))),with=FALSE])

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
			   cat(years[iy],'sample with heads has',n,'obs\n')
			   cat('dropping non-heads leaves',nrow(yind),'obs\n')
		   }
		   yind[,c(curr[,ind.head],"headyes") := NULL]
			# set names on individual index
			setnames(yind,c("ID1968","pernum","interview"))
		} else {
			# set names on individual index
			setnames(yind,c("ID1968","pernum","interview","relation.head"))
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
			rm(list=ls(all=T,envir=tmp.env),envir=tmp.env)
		   	load(file=fam.dat[iy],envir=tmp.env)
			tmp             <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
			fam.dicts[[iy]] <- NULL
			tmp             <- data.table(tmp)
		}

		if (verbose){
			cat('loaded family file:',fam.dat,'\n')
			cat('current memory load in MB:\n')
			vs = ceiling(object.size(tmp))
			print(vs,units="Mb")
		}
	

		# add vars from from file that the user requested.
		curvars <- fam.vars[.(years[iy]),which(names(fam.vars)!="year"),with=FALSE]
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
		idx <- which(!is.na(unlist(fam.vars[.(years[iy])][,curnames,with=FALSE])))[1]	# index of first non NA variable
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
		data2[,always := max(present) == length(years),by=pid]
		data2 <- copy(data2[always==TRUE])
		data2[,always := NULL]
	} else if (is.numeric(design)){
		data2[,enough := max(present) >= design,by=pid]
		data2 <- copy(data2[enough==TRUE])
		data2[,enough := NULL]
	} else if (design=="all"){
		# do nothing
	}
	data2[,present := NULL]

	#     setkey(datas,pid,year)
	#     rm(ind)
	out <- list(data=data2,dict=fam.dicts)
	return(out)
}








	  
	  
	  
	  








