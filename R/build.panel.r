

#' build.panel: Build PSID panel data set
#' 
#' @description Builds a panel data set in wide format with id variables \code{pid} (unique person identifier) and \code{year} from individual PSID family files.
#' @details 
#' takes desired variables from family files for specified years in folder \code{datadir} and merges using the id information in \code{IND2013ER.xyz}, which must be in the same directory. The raw data can be supplied in stata .dta format or it can be directly downloaded from the PSID server to folders \code{datadir} or \code{tmpdir}. Notice that currently only stata format <= 12 is supported (so do \code{saveold} in stata). The user can change subsetting criteria as well as sample designs. The package allows the missing variables in certain waves to be accounted for automatically, i.e. the variables are inserted in the missing year as \code{NA}.
#' Merge: the variables \code{interview number} in each family file map to 
#' the \code{interview number} variable of a given year in the individual file. Run \code{example(build.panel)} for a demonstration.
#' Accepted input data are stata format .dta, .csv files or R data formats .rda and RData. Similar in usage to stata module \code{psiduse}.
#' @param datadir either \code{NULL}, in which case saves to tmpdir or path to directory containing family files ("FAMyyyy.xyz") and individual file ("IND2009ER.xyz") in admissible formats .xyz. Admissible are .dta, .csv, .RData, .rda. Please follow naming convention. Only .dta version <= 12 supported. Recommended usage is to specify \code{datadir}.
#' @param fam.vars data.frame of variable to retrieve from family files. Can contain see example for required format.
#' @param ind.vars data.frame of variables to get from individual file. In almost all cases this will be the type of survey weights you want to use. don't include id variables ER30001 and ER30002.
#' @param SAScii logical TRUE if you want to directly download data into Rda format (no dependency on STATA/SAS/SPSS). may take a long time, but downloads only once if you specify \code{datadir}.
#' @param heads.only logical TRUE if user wants current household heads only. 
#' @param sample string indicating which sample to select: "SRC" (survey research center), "SEO" (survey for economic opportunity), "immigrant" (immigrant sample), "latino" (Latino family sample). Defaults to NULL, so no subsetting takes place.
#' @param design either character \emph{balanced} or \emph{all} or integer. \emph{balanced} means only individuals who appear in each wave are considered. \emph{All} means all are taken. An integer value stands for minimum consecutive years of participation, i.e. design=3 means present in at least 3 consecutive waves.
#' @param verbose logical TRUE if you want verbose output.
#' @import SAScii RCurl data.table
#' @return
#' \item{data}{resulting \code{data.table}. the variable \code{pid} is the unique person identifier, constructed from ID1968 and pernum.}
#' \item{dict}{data dictionary if stata data was supplied, NULL else}
#' @export
#' @examples 
#' \dontrun{
#' # ################################################
#' # Real-world example: not run because takes long.
#' # Build panel with income, wage, age and education
#' # ################################################
#' 
#' r = system.file(package="psidR")
#' f = fread(file.path(r,"psid-lists","famvars.txt"))
#' i = fread(file.path(r,"psid-lists","indvars.txt"))
#' 
#' f[1:38,vgroup := "wage"]
#' f[39:76,vgroup := "earnings"]
#' setkey(f,vgroup)
#' 
#' i[1:38, vgroup := "age"]
#' i[39:76, vgroup := "educ"]  # caution about 2 first years: no educ data
#' i[77:114, vgroup := "weight"]
#' setkey(i,vgroup)
#' 
#' ind = cbind(i[J("age"),list(year,age=variable)],
#' 			   i[J("educ"),list(educ=variable)],
#' 			   i[J("weight"),list(weight=variable)])
#' fam = cbind(f[J("wage"),list(year,wage=variable)],
#' 			   f[J("earnings"),list(earnings=variable)])
#' 
#' # caution: this step will take many hours
#' d = build.panel(datadir="~/data",
#'                 fam.vars=fam,
#' 				   ind.vars=ind,
#'                 SAScii = TRUE, 
#'                 heads.only = TRUE,
#'                 sample="SRC",
#'                 design=2)
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
#' ## Data acquisition step: you download data or
#' ## run build.panel with sascii=TRUE
#' 
#' # testPSID creates artifical PSID data
#' td <- testPSID(N=12,N.attr=0)
#' fam1985 <- copy(td$famvars1985)
#' fam1986 <- copy(td$famvars1986)
#' IND2009ER <- copy(td$IND2009ER)
#' 
#' # create a temporary datadir
#' my.dir <- tempdir()
#' #save those in the datadir
#' # notice different R formats admissible
#' save(fam1985,file=paste0(my.dir,"/FAM1985ER.rda"))
#' save(fam1986,file=paste0(my.dir,"/FAM1986ER.RData"))
#' save(IND2009ER,file=paste0(my.dir,"/IND2009ER.RData"))
#' 
#' ## end Data acquisition step.
#' 
#' # now define which famvars
#' famvars <- data.frame(year=c(1985,1986),
#'                       money=c("Money85","Money86"),
#'                       age=c("age85","age86"))
#' 
#' # create ind.vars
#' indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
#' 
#' # call the builder
#' # data will contain column "relation.head" holding the relationship code.
#' 
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,
#'                  heads.only=FALSE,verbose=TRUE)	
#' 
#' # see what happens if we drop non-heads
#' # only the ones who are heads in BOTH years 
#' # are present (since design='balanced' by default)
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,
#'                  heads.only=TRUE,verbose=FALSE)	
#' print(d$data[order(pid)],nrow=Inf)
#' 
#' # change sample design to "all": 
#' # we'll keep individuals if they are head in one year,
#' # and drop in the other
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,heads.only=TRUE,
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
build.panel <- function(datadir=NULL,fam.vars,ind.vars=NULL,SAScii=FALSE,heads.only=FALSE,sample=NULL,design="balanced",verbose=FALSE){
	

	# locally bind all variables to be used in a data.table

	interview <- headyes <- .SD <- fam.interview <- ind.interview <- ind.head <- ER30001 <- ind.head.num <- pid <- ID1968 <- pernum <- isna <- present <- always <- enough <- ind.seq <- NULL

	stopifnot(is.numeric(fam.vars$year))
	years <- fam.vars$year

	# figure out platform
	s <- .Platform$file.sep
	if ( .Platform$OS.type != 'windows' ) {
		warning("I'm setting your encoding to windows now")
		options( encoding = "windows-1252" )		# # only macintosh and *nix users need this line
	}
	
	
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

		lf = list.files(datadir)

		# all psid family files
		family    <- data.frame(year = c( 1968:1997 , seq( 1999 , 2013 , 2 ) ),file = c( 1056 , 1058:1082 , 1047:1051 , 1040 , 1052 , 1132 , 1139 , 1152  , 1156, 1164 ))

		#subset to the years we want
		family <- family[family$year %in% years, ]

		families.down <- rep(FALSE,nrow(family))
		ind.down <- FALSE
		for ( i in 1:nrow( family )) {
			if ((paste0("FAM" , family[ i , 'year' ], "ER.rda") %in% lf)) {
				families.down[i] <- TRUE
				cat('found ',paste0("FAM" , family[ i , 'year' ], "ER.rda"), 'already downloaded \n')
			}
		}

		# check if datadir contains individual index already
		if (("IND2013ER.rda" %in% lf)) {
			#download latest individual index
			ind.down = TRUE
		}
		if (all(all(families.down),ind.down)) {
			cat("everything already downloaded. Build dataset now.\n")
		} else {
			cat("Will download missing datasets now.\n")
	    	confirm <- readline("This can take several hours/days to download.\n want to go ahead? give me 'yes' or 'no'.")
			if (confirm=="yes"){

				ftype <- "Rdata"
				
				user <- readline("please enter your PSID username: ")
				pass <- readline("please enter your PSID password: ")
				
				curl = getCurlHandle()
				curlSetOpt(cookiejar = 'cookies.txt', followlocation = TRUE, autoreferer = TRUE, curl = curl)

				html <- getURL('http://simba.isr.umich.edu/u/Login.aspx', curl = curl)

				viewstate <- as.character(sub('.*id="__VIEWSTATE" value="([0-9a-zA-Z+/=]*).*', '\\1', html))

				# extract the `eventvalidation` string
				eventvalidation <- 
					as.character(
						sub(
							'.*id="__EVENTVALIDATION" value="([0-9a-zA-Z+/=]*).*' , 
							'\\1' , 
							html
						)
					)

				# construct a list full of parameters to pass to the umich website
				params <- 
					list(
						'ctl00$ContentPlaceHolder1$Login1$UserName'    = user ,
						'ctl00$ContentPlaceHolder1$Login1$Password'    = pass ,
						'ctl00$ContentPlaceHolder1$Login1$LoginButton' = 'Log In' ,
						'__VIEWSTATE'                                  = viewstate ,
						'__EVENTVALIDATION'                            = eventvalidation
				    )

				#file number 1053 is always the individual cross year index file
				#it must always be the last file in this list.
				#you always want to download that.

				for ( i in 1:nrow( family )) {
					if (!(paste0("FAM" , family[ i , 'year' ], "ER.rda") %in% lf)) {
						get.psid( family[ i , 'file' ] ,name= paste0(datadir, "FAM" , family[ i , 'year' ], "ER") , params , curl )
					}
				}

				# check if datadir contains individual index already
				if (!("IND2013ER.rda" %in% lf)) {
					#download latest individual index
					get.psid( 1053 ,name= paste0(datadir, "IND2013ER") , params , curl )
				}

				cat('finished downloading files to', datadir,'. continuing now to build the dataset.\n')
							
			} else if (confirm=="no") {
				break
			}
		}  # end download data
	}


	# work out file types
	# -------------------

	# figure out filestypes in datadir
	l <- list.files(datadir)
	if (length(l)==0) stop('there is something wrong with the data directory. please check path')
	for (i in 1:length(l)) {
  		if (tail(strsplit(l[i],"\\.")[[1]],1) == "dta") {ftype   <- "stata"; break}
  		if (tail(strsplit(l[i],"\\.")[[1]],1) == "rda") {ftype   <- "Rdata"; break}
  		if (tail(strsplit(l[i],"\\.")[[1]],1) == "RData") {ftype <- "Rdata"; break}
  		if (tail(strsplit(l[i],"\\.")[[1]],1) == "csv") {ftype   <- "csv"; break}
	}
	if (!(exists("ftype"))) stop("No .dta, .rda, .RData or .csv files found in directory")

	if (verbose) cat('psidR: loading data\n')
	if (ftype=="stata"){
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE,ignore.case=TRUE))
        fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		tmp <- grep("IND",l,value=TRUE,ignore.case=TRUE)
		if (length(tmp)>1) {
		    warning(cat("Warning: you have more than one IND file in your datadir.\nI take the last one:",tail(tmp,1),"\n"))
			ind.file <- paste0(datadir,tail(tmp,1))	# needs to be updated with next data delivery.
		} else {
			ind.file <- paste0(datadir,grep("IND",l,value=TRUE,ignore.case=TRUE))	# needs to be updated with next data delivery.
		}
		# TODO
		# maybe should use 
		# ind      <- read.dta13(file=ind.file)
		# if version 13?
		ind      <- read.dta(file=ind.file)
		ind.dict <- data.frame(code=names(ind),label=attr(ind,"var.labels"))
		ind      <- data.table(ind)
	} else if (ftype=="Rdata") {
		# data downloaded directly into a dataframe
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE,ignore.case=TRUE))
		fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		tmp <- grep("IND",l,value=TRUE,ignore.case=TRUE)
		if (length(tmp)>1) {
		    warning(cat("Warning: you have more than one IND file in your datadir.\nI take the last one:",tail(tmp,1),"\n"))
			ind.file <- paste0(datadir,tail(tmp,1))	# needs to be updated with next data delivery.
		} else {
			ind.file <- paste0(datadir,grep("IND",l,value=TRUE,ignore.case=TRUE))	# needs to be updated with next data delivery.
		}
		tmp.env  <- new.env()
		load(file=ind.file,envir=tmp.env)
		ind      <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
		setnames(ind,names(ind), sapply(names(ind), toupper))	## convert all column names to uppercase
		ind.dict <- NULL
		ind      <- data.table(ind)
	} else if (ftype=="csv") {
		fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE,ignore.case=TRUE))
		fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)
		tmp <- grep("IND",l,value=TRUE,ignore.case=TRUE)
		if (length(tmp)>1) {
		    warning(cat("Warning: you have more than one IND file in your datadir.\nI take the last one:",tail(tmp,1),"\n"))
			ind.file <- paste0(datadir,tail(tmp,1))	# needs to be updated with next data delivery.
		} else {
			ind.file <- paste0(datadir,grep("IND",l,value=TRUE,ignore.case=TRUE))	# needs to be updated with next data delivery.
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


	# make an index of interview numbers for each year
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
		# subset only if requested.
		curr <- ids[list(years[iy])]
		if (years[iy] == 1968){
		  # there is no sequence variable
		  ind.subsetter <- as.character(curr[,list(ind.interview,ind.head)])	# keep from ind file
		} else {
		  ind.subsetter <- as.character(curr[,list(ind.interview,ind.seq,ind.head)])	# keep from ind file
		}
		def.subsetter <- c("ER30001","ER30002")	# must keep those in all years


		# select required variables from ind 
		# ----------------------------------
		# the default column order is: 
		# "ER30001","ER30002", "current year interview var", "current year sequence number", "current year head indicator"
		# get a character vector from ind.vars with variable names for this year
		ind.vars.yr <- ind.vars[list(years[iy]),which(names(ind.vars)!="year"),with=FALSE]

		# issue https://github.com/floswald/psidR/issues/4
		# ------------------------------------------------
		# check for NA in ind.vars: these are years when a certain variable isn not available in the individual index file.
		# adjust for first year (1968) when `sequence` was not available
		ind.nas <- subset(ind.vars.yr,select=names(ind.vars.yr)[is.na(ind.vars.yr)])
		yind    <- copy(ind[,c(def.subsetter,c(ind.subsetter,as.character(subset(ind.vars.yr,select=names(ind.vars.yr)[!is.na(ind.vars.yr)])))),with=FALSE])	
		# add NA columns
		if (length(ind.nas)>0){
		    yind[,(names(ind.nas)) := NA]
		}
		if (years[iy]==1968){
		    yind[,sequence := NA]
		    setcolorder(yind,c(1:3,ncol(yind),4:(ncol(yind)-1)))
		}


    
    	# sample selection
    	# ----------------

    	# based on: https://psidonline.isr.umich.edu/Guide/FAQ.aspx?Type=ALL#250

    	if (!is.null(sample)){

    		if (sample == "SRC"){
			   n    <- nrow(yind)
			   yind <- copy(yind[ER30001<3000])	# individuals 1-2999 are from SRC sample
			   if (verbose){
				   cat('full',years[iy],'sample has',n,'obs\n')
				   cat(sprintf('you selected %d obs belonging to %s \n',nrow(yind),sample))
			   }
    		} else if (sample == "SEO"){
			   n    <- nrow(yind)
			   yind <- copy(yind[ER30001<7000 & ER30001>5000])
			   if (verbose){
				   cat('full',years[iy],'sample has',n,'obs\n')
				   cat(sprintf('you selected %d obs belonging to %s \n',nrow(yind),sample))
			   }
    		} else if (sample == "immigrant"){
			   n    <- nrow(yind)
			   yind <- copy(yind[ER30001<5000 & ER30001>3000])	# individuals 1-2999 are from SRC sample
			   if (verbose){
				   cat('full',years[iy],'sample has',n,'obs\n')
				   cat(sprintf('you selected %d obs belonging to %s \n',nrow(yind),sample))
			   }
    		} else if (sample == "latino"){
			   n    <- nrow(yind)
			   yind <- copy(yind[ER30001<9309 & ER30001>7000])	# individuals 1-2999 are from SRC sample
			   if (verbose){
				   cat('full',years[iy],'sample has',n,'obs\n')
				   cat(sprintf('you selected %d obs belonging to %s \n',nrow(yind),sample))
			   }
    		}
    	}

    	# heads only selection
    	# --------------------

    	# https://github.com/floswald/psidR/issues/2
		# for current heads only need to subset "relationship to head" AS WELL AS "sequence number" == 1 
		# otherwise a head who died between last and this wave is still head, so there would be two heads in that family.
		# https://psidonline.isr.umich.edu/Guide/FAQ.aspx?Type=ALL#150
		if (heads.only) {
		   n    <- nrow(yind)
		   if (years[iy]==1968){
		     yind <- yind[,headyes := (yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num])]
		   } else {
		     yind <- yind[,headyes := (yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num]) & (yind[,curr[,ind.seq],with=FALSE]== 1)]
		   }
		   yind <- copy(yind[headyes==TRUE])
		   if (verbose){
			   cat('dropping non-current-heads leaves',nrow(yind),'obs\n')
		   }
		   yind[,headyes := NULL]
		}
		

		if (length(ind.nas)>0){
		    setnames(yind,c("ID1968","pernum","interview","sequence","relation.head",
		                  (names(ind.vars)[-1])[-which(names(ind.vars)[-1] %in% names(ind.nas))],
		                  names(ind.nas)))  
		} else {
		    setnames(yind,c("ID1968","pernum","interview","sequence","relation.head",
		                  (names(ind.vars)[-1])))
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

		# convert all variable names to lower case in both fam.vars and data file
		curvars <- fam.vars[list(years[iy]),which(names(fam.vars)!="year"),with=FALSE]
		tmpnms = tolower(as.character(curvars))
		for (i in 1:length(tmpnms)){
			curvars[[i]] <- tmpnms[i]
		}
		setnames(tmp,tolower(names(tmp)))
	
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



#' Build full example PSID
#' 
#' @description Builds a panel from the full PSID dataset (i.e. all years)
#' @export
#' @return a data.table with panel data
build.psid <- function(){
	r = system.file(package="psidR")
	f = fread(file.path(r,"psid-lists","famvars.txt"))
	i = fread(file.path(r,"psid-lists","indvars.txt"))

	# this is just to make R CMD CHECK happy
	J <- vgroup <- variable <- NULL
	f[1:38,vgroup := "wage"]
	f[39:76,vgroup := "earnings"]
	setkey(f,vgroup)

	i[1:38, vgroup := "age"]
	i[39:76, vgroup := "educ"]  # caution about 2 first years: no educ data
	i[77:114, vgroup := "weight"]
	setkey(i,vgroup)

	ind = cbind(i[J("age"),list(year,age=variable)],i[J("educ"),list(educ=variable)],i[J("weight"),list(weight=variable)])
	fam = cbind(f[J("wage"),list(year,wage=variable)],f[J("earnings"),list(earnings=variable)])

	d = build.panel(datadir="~/git/bk/data/",fam.vars=fam,ind.vars=ind,SAScii = TRUE, heads.only = TRUE,sample="SRC",design=2,verbose=TRUE)
	save(d,file="~/psid.RData")
	return(d$data)
}






	  
	  
	  
	  








