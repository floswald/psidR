

#' build.panel: Build PSID panel data set
#' 
#' @description Builds a panel data set with id variables \code{pid} (unique person identifier) and \code{year} from individual PSID family files and supplemental wealth files.
#' @details 
#' There are several supported approches. Approach one downloads stata data, uses stata to build each wave, then puts it together with `psidR`. The second (recommended) approach downloads all data directly from the psid servers (no Stata needed. For this approach you need to supply the precise names of psid variables - those variable names vary by year. E.g. \emph{total family income} will have different names in different waves. The function \code{\link{getNamesPSID}} greatly helps collecting names for all waves.
#' Merge: the variables \code{interview number} in each family file map to 
#' the \code{interview number} variable of a given year in the individual file. Run \code{example(build.panel)} for a demonstration.
#' For approach one, accepted input data are stata format .dta, .csv files or R data formats .rda and RData. This usage is similar to stata module \code{psiduse}. 
#' Approach two follows the strategy introduced at \url{http://asdfree.com}. In fact, both approaches use the same function \code{save.psid} to download the data, \code{psidR} automates the merge and subsetting proceedure for you. 
#' @param datadir either \code{NULL}, in which case saves to tmpdir or path to directory containing family files ("FAMyyyy.RData") and individual file ("IND2009ER.RData").
#' @param fam.vars data.frame of variable to retrieve from family files. Can contain see example for required format.
#' @param ind.vars data.frame of variables to get from individual file. In almost all cases this will be the type of survey weights you want to use. don't include id variables ER30001 and ER30002.
#' @param wealth.vars data.frame of variables to get from the wealth supplement files.
#' @param current.heads.only logical TRUE if user wants current household heads only. Distinguishes mover outs heads.
#' @param heads.only logical TRUE if user wants household heads only. Household heads in sample year.
#' @param sample string indicating which sample to select: "SRC" (survey research center), "SEO" (survey for economic opportunity), "immigrant" (immigrant sample), "latino" (Latino family sample). Defaults to NULL, so no subsetting takes place.
#' @param design either character \emph{balanced} or \emph{all} or integer. \emph{balanced} means only individuals who appear in each wave are considered. \emph{All} means all are taken. An integer value stands for minimum consecutive years of participation, i.e. design=3 means present in at least 3 consecutive waves.
#' @param loglevel one of INFO, WARN and DEBUG. INFO by default.
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
#' # optionally: add wealth supplements!
#' # ################################################
#' 
#' # The package is installed with a list of variables
#' # Alternatively, search for names with \code{\link{getNamesPSID}}
#' # This is the body of function build.psid()
#' # (so why not call build.psid() and see what happens!)
#' r = system.file(package="psidR")
#' if (small){
#'   f = fread(file.path(r,"psid-lists","famvars-small.txt"))
#'   i = fread(file.path(r,"psid-lists","indvars-small.txt"))
#'   if (wealth){
#'     w = fread(file.path(r,"psid-lists","wealthvars-small.txt"))
#'   }
#' } else {
#'   f = fread(file.path(r,"psid-lists","famvars.txt"))
#'   i = fread(file.path(r,"psid-lists","indvars.txt"))
#'   if (wealth){
#'     w = fread(file.path(r,"psid-lists","wealthvars.txt"))
#'   }
#' }
#' setkey(i,"name")
#' setkey(f,"name")
#' if (wealth) setkey(w,"name")
#' i = dcast(i[,list(year,name,variable)],year~name)
#' f = dcast(f[,list(year,name,variable)],year~name)
#' if (wealth) {
#'   w = dcast(w[,list(year,name,variable)],year~name)
#'   d = build.panel(datadir="~/datasets/psid/",fam.vars=f,
#'                  ind.vars=i,wealth.vars=w, 
#'                  heads.only =TRUE,sample="SRC",
#'                  design="all")
#'   save(d,file="~/psid_wealth.RData")
#' } else {
#'   d = build.panel(datadir="~/datasets/psid/",fam.vars=f,
#'                  ind.vars=i, 
#'                  heads.only =TRUE,sample="SRC",
#'                  design="all")
#'   save(d,file="~/psid_no_wealth.RData")
#' }
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
#' ## Data acquisition step:
#' ## run build.panel with sascii=TRUE
#' 
#' # testPSID creates artifical PSID data
#' td <- testPSID(N=12,N.attr=0)
#' fam1985 <- data.table::copy(td$famvars1985)
#' fam1986 <- data.table::copy(td$famvars1986)
#' IND2009ER <- data.table::copy(td$IND2009ER)
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
#'                  heads.only=FALSE)	
#' 
#' # see what happens if we drop non-heads
#' # only the ones who are heads in BOTH years 
#' # are present (since design='balanced' by default)
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,
#'                  heads.only=TRUE)	
#' print(d[order(pid)],nrow=Inf)
#' 
#' # change sample design to "all": 
#' # we'll keep individuals if they are head in one year,
#' # and drop in the other
#' d <- build.panel(datadir=my.dir,fam.vars=famvars,
#'                  ind.vars=indvars,heads.only=TRUE,
#'                  design="all")	
#' print(d[order(pid)],nrow=Inf)
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
build.panel <- function(datadir=NULL,fam.vars,ind.vars=NULL,wealth.vars=NULL,heads.only=FALSE,current.heads.only=FALSE,sample=NULL,design="balanced",loglevel=INFO){
	
  flog.threshold(loglevel)

	# locally bind all variables to be used in a data.table
	# or R CMD CHECK complains.

	interview <- headyes <- .SD <- fam.interview <- ind.interview <- ind.head <- ER30001 <- ind.head.num <- pid <- ID1968 <- pernum <- isna <- present <- always <- enough <- ind.seq <- name <- variable <- NULL

	stopifnot(is.numeric(fam.vars$year))
	years <- fam.vars$year
	flog.debug("years:",years,capture=TRUE)

	# figure out platform
	s <- .Platform$file.sep
	if ( .Platform$OS.type != 'windows' ) {
		# warning("I'm setting your encoding to windows now")
		options( encoding = "windows-1252" )		# # only macintosh and *nix users need this line
	}
	
	
	# sort out data storage: either to tmp or datadir
	if (is.null(datadir)){
		datadir <- tempdir()	
		datadir <- paste0(datadir,s)
	} else {
		if (substr(datadir,nchar(datadir),nchar(datadir))!=s) datadir <- paste0(datadir,s)
	}
	flog.debug("datadir: %s",datadir)

	# are we processing wealth supplements?
	any.wealth = is.data.frame(wealth.vars)
	flog.debug("any.wealth? %d",any.wealth)
	

	# data acquisition
	# ----------------
	lf = list.files(datadir)
	flog.debug("datadir contains: ",lf,capture=TRUE)

	wlth.down <- TRUE  # initiate to something

	# all psid family files
	family    <- data.frame(year = c( 1968:1997 , seq( 1999 , 2015 , 2 ) ),file = c( 1056 , 1058:1082 , 1047:1051 , 1040 , 1052 , 1132 , 1139 , 1152  , 1156, 1164 , 1183 ))
	# family    <- data.frame(year = c( 1968:1997 , seq( 1999 , 2013 , 2 ) ),file = c( 1056 , 1058:1082 , 1047:1051 , 1040 , 1052 , 1132 , 1139 , 1152  , 1156, 1164  ))

	#subset to the years we want
	family <- family[family$year %in% years, ]
	flog.debug("family df:",family,capture=TRUE)

	families.down <- rep(FALSE,nrow(family))
	for ( i in 1:nrow( family )) {
		if ((paste0("FAM" , family[ i , 'year' ], "ER.rda") %in% lf) | (paste0("FAM" , family[ i , 'year' ], "ER.RData") %in% lf)) {
			families.down[i] <- TRUE
			flog.info('found %s already downloaded',paste0("FAM" , family[ i , 'year' ], "ER.rda"))
		}
	}

	if (any.wealth){
		# if any of 1984, 1989, 1994, 1999, 2001, 2003, 2005, 2007 in years, also download the associated wealth supplement
		wlth = data.frame(year=c(1984, 1989, 1994, 1999, 2001, 2003, 2005, 2007),file=c(1147,1148,1149,1150,1130,1131,1133,1140))
		flog.debug("working on wealth")
		wlth = wlth[wealth.vars$year %in% years, ]
		flog.debug("wlth: ",wlth, capture=TRUE)

		if ( nrow(wlth) > 0 ){
		  wlth.down <- rep(FALSE,nrow(wlth))
		    for ( i in 1:nrow( wlth )) {
		      if ((paste0("WEALTH" , wlth[ i , 'year' ], "ER.rda") %in% lf) | (paste0("WEALTH" , wlth[ i , 'year' ], "ER.RData") %in% lf)) {
		        wlth.down[i] <- TRUE
		        flog.info('found %s already downloaded',paste0("WEALTH" , wlth[ i , 'year' ]))
		      } else {
		        flog.info('will download as %s',paste0("WEALTH" , wlth[ i , 'year' ], "ER.rda"))
		      }
		    }
			} else {
			  wlth.down = TRUE
			  flog.info('All Wealth files already downloaded')
			  
			  # saying we already downloaded even if we didn't look for wealth vars.
			}
		
	}

	ind.down <- FALSE
	# check if datadir contains individual index already
	if (("IND2015ER.rda" %in% lf) | ("IND2015ER.RData" %in% lf)) {
		#download latest individual index
		ind.down = TRUE
	}
	if (all(all(families.down),ind.down,all(wlth.down))) {
		flog.info("everything already downloaded. Build dataset now")
	} else {
		flog.info("Will download missing datasets now")
		if (!all(families.down)) {
			flog.info("will download family files: %s",paste(family[!families.down,"year"],collapse=", "))
		} 
		if (!ind.down) {
			flog.info("will download: IND2015ER")
		} 
		if (!wlth.down) {
			flog.info("will download missing wealth files.")
		}

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
      if (any.wealth){
        if (nrow(wlth)>0){
          for ( i in 1:nrow(wlth)){
            if (!(wlth.down[i])){
              get.psid( wlth[ i , 'file' ] ,name= paste0(datadir, "WEALTH" , wlth[ i , 'year' ], "ER") , params , curl )
            }
          }
        }
      }
			
			# check if datadir contains individual index already
			if (!("IND2015ER.rda" %in% lf)) {
				#download latest individual index
				get.psid( 1053 ,name= paste0(datadir, "IND2015ER") , params , curl )
			}

			flog.info('finished downloading files to %s', datadir)
			flog.info('continuing now to build the dataset')
						
		} else if (confirm=="no") {
			return(0)
		}
	}  # end download data

	l <- list.files(datadir)
	if (length(l)==0) stop('there is something wrong with the data directory. please check path')
	for (i in 1:length(l)) {
  		if (tail(strsplit(l[i],"\\.")[[1]],1) == "rda") {ftype   <- "Rdata"; break}
  		if (tail(strsplit(l[i],"\\.")[[1]],1) == "RData") {ftype <- "Rdata"; break}
	}
	if (!(exists("ftype"))) stop("No .rda or .RData files found in directory\n you need to download data first.")

	flog.info('psidR: Loading Family data from .rda files')
	# familiy data downloaded directly into a dataframe
	fam.dat  <- paste0(datadir,grep("FAM",l,value=TRUE,ignore.case=TRUE))
	fam.dat  <- grep(paste(years,collapse="|"),fam.dat,value=TRUE)

	# wealth data downloaded directly into a dataframe
	if (any.wealth){
		wlth.dat  <- paste0(datadir,grep("WEALTH",l,value=TRUE,ignore.case=TRUE))
		wlth.dat  <- grep(paste(years,collapse="|"),wlth.dat,value=TRUE)
	}

	# individual index
	tmp <- grep("IND",l,value=TRUE,ignore.case=TRUE)
	if (length(tmp)>1) {
	    flog.warn("Warning: you have more than one IND file in your datadir. I take the last one:",tail(tmp,1))
		ind.file <- paste0(datadir,tail(tmp,1))	# needs to be updated with next data delivery.
	} else {
		ind.file <- paste0(datadir,grep("IND",l,value=TRUE,ignore.case=TRUE))	# needs to be updated with next data delivery.
	}
	tmp.env  <- new.env()
	load(file=ind.file,envir=tmp.env)
	ind      <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
	setnames(ind,names(ind), sapply(names(ind), toupper))	## convert all column names to uppercase
	ind      <- data.table(ind)
	

	flog.info('psidR: loaded individual file: %s',ind.file)
	vvs = ceiling(object.size(ind)/1024^2)
	flog.info("psidR: total memory load in MB: %d",as.numeric(vvs))

	# output data.tables
	datas <- vector("list",length(years))

	# make an index of interview numbers for each year
	ids <- makeids()
	w_ids <- makeids.wealth()

		flog.debug('psidR: here is the list of hardcoded PSID variables')
		flog.debug('psidR: The merge is based on equal values in ind.interview and fam.interview')
		flog.debug("ids:",ids,capture=TRUE)
	
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
	# convert ind.vars to data.table if not null
	if (any.wealth){
		stopifnot(is.data.frame(wealth.vars))
		wealth.vars <- data.table(wealth.vars)
		wealth.vars <- copy(wealth.vars[,lapply(.SD,make.char)])
		setkey(wealth.vars,year)
		wealth.vars[,interview := w_ids[wealth.vars][,interview]]
	}

	# which vars to keep from ind.files?
	if (!is.null(ind.vars))	stopifnot(is.list(ind.vars))

	# add compulsory vars to fam.vars
	fam.vars[,interview := ids[fam.vars][,fam.interview]]

	# loop over years
	for (iy in 1:length(years)){
		
			flog.info('psidR: currently working on data for year %d',years[iy])
		
 
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
		ind.nas <- NULL
		if (!is.null(ind.vars)){
			ind.vars.yr <- ind.vars[list(years[iy]),which(names(ind.vars)!="year"),with=FALSE]
			flog.debug("ind.vars.yr:",ind.vars.yr,capture=TRUE)

			# issue https://github.com/floswald/psidR/issues/4
			# ------------------------------------------------
			# check for NA in ind.vars: these are years when a certain variable isn not available in the individual index file.
			# adjust for first year (1968) when `sequence` was not available
			ind.notnas <- NULL
			if (any(is.na(ind.vars.yr))){
				ind.nas <- ind.vars.yr[,which(is.na(ind.vars.yr)),with=FALSE]
				flog.debug("ind.nas:",ind.nas,capture=TRUE)
			}
			if (any(!is.na(ind.vars.yr))){
				ind.notnas <- ind.vars.yr[,which(!is.na(ind.vars.yr)),with=FALSE]
				flog.debug("ind.notnas:",ind.notnas,capture=TRUE)
				yind    <- copy(ind[,c(def.subsetter,c(ind.subsetter,as.character(ind.notnas))),with=FALSE])	
			} else {
				yind    <- copy(ind[,c(def.subsetter,c(ind.subsetter)),with=FALSE])	
			}
			# add NA columns
			if (!is.null(ind.nas)){
			    yind[,(as.character(ind.nas)) := NA]
			}
		} else {
			yind <- copy(ind[,c(def.subsetter,c(ind.subsetter)),with=FALSE])	
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
			   
			   flog.info('full %d sample has %d obs',years[iy],n)
			   flog.info('you selected %d obs belonging to %s',nrow(yind),sample)
			   
    		} else if (sample == "SEO"){
			   n    <- nrow(yind)
			   yind <- copy(yind[ER30001<7000 & ER30001>5000])
			   
			   flog.info('full %d sample has %d obs',years[iy],n)
			   flog.info('you selected %d obs belonging to %s',nrow(yind),sample)
			   
    		} else if (sample == "immigrant"){
			   n    <- nrow(yind)
			   yind <- copy(yind[ER30001<5000 & ER30001>3000])	# individuals 1-2999 are from SRC sample
			   flog.info('full %d sample has %d obs',years[iy],n)
			   flog.info('you selected %d obs belonging to %s',nrow(yind),sample)
    		} else if (sample == "latino"){
			   n    <- nrow(yind)
			   yind <- copy(yind[ER30001<9309 & ER30001>7000])	# individuals 1-2999 are from SRC sample
			   flog.info('full %d sample has %d obs',years[iy],n)
			   flog.info('you selected %d obs belonging to %s',nrow(yind),sample)
    		}
    	}

    	# current heads only selection
    	# --------------------

    	# https://github.com/floswald/psidR/issues/2
		# for current heads only need to subset "relationship to head" AS WELL AS "sequence number" == 1 
		# otherwise a head who died between last and this wave is still head, so there would be two heads in that family.
		# https://psidonline.isr.umich.edu/Guide/FAQ.aspx?Type=ALL#150
		if (current.heads.only) {
			n    <- nrow(yind)
			if (years[iy]==1968){
				yind <- yind[,headyes := (yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num])]
			} else {
				yind <- yind[,headyes := (yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num]) & (yind[,curr[,ind.seq],with=FALSE]== 1)]
			}
			yind <- copy(yind[headyes==TRUE])
			flog.info('dropping non-current-heads leaves %d obs',nrow(yind))
			yind[,headyes := NULL]
		} else if (heads.only){
			# https://psidonline.isr.umich.edu/Guide/FAQ.aspx?Type=ALL#250
			# To create a single year Head/Wife file: 
			# Select individuals with Relationship to Head of "Head" (a code value of 1 for 1968-1982; code 10 from 1983 onward) 
			# and with values for Sequence Number in the range 1-20. 
		   n    <- nrow(yind)
		   if (years[iy]==1968){
		     yind <- yind[,headyes := (yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num])]
		   } else {
		     yind <- yind[,headyes := (yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num]) & ((yind[,curr[,ind.seq],with=FALSE]> 0) & (yind[,curr[,ind.seq],with=FALSE]< 21))]
		   }
		   yind <- copy(yind[headyes==TRUE])
		   
		   flog.info('dropping non-heads leaves %d obs',nrow(yind))
		   
		   yind[,headyes := NULL]

		}

		

		if (!is.null(ind.nas)){
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
		# ==================================

		# load data for current year, make data dictionary for subsets and save data as data.table
		rm(list=ls(envir=tmp.env),envir=tmp.env)
	   	load(file=fam.dat[iy],envir=tmp.env)
		tmp             <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
		tmp             <- data.table(tmp)
		
		vs = ceiling(object.size(tmp)/1024^2)
		flog.debug('loaded family file: ',fam.dat[iy])
		flog.debug('current memory load in MB: %d',vs)



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
			tmp[,(nanames) := NA_real_,with=FALSE]
			setnames(tmp,c(curnames[-na],nanames))
			setkey(tmp,interview)
		} else {
			codes <- as.character(curvars)
			tmp   <- copy(tmp[,codes,with=FALSE])
			setnames(tmp,curnames)
			setkey(tmp,interview)
		}
		
		# merge family and yind
		m <- copy(tmp[yind])
		m[,year := years[iy] ]
		setkey(m,interview)

		# bring in wealth files, subset them
		# ==================================
		# merge m and wealth
		if (any.wealth){
			# check if there is a wealth file for this year
			iw = grep(years[iy],wlth.dat,value=TRUE)
			if (length(iw)>0){
				rm(list=ls(envir=tmp.env),envir=tmp.env)
			   	load(file=iw,envir=tmp.env)
				tmp             <- get(ls(tmp.env),tmp.env)	# assign loaded dataset a new name
				tmp             <- data.table(tmp)
				flog.debug("wealth tmp: ",head(tmp),capture=TRUE)
			
				# convert all variable names to lower case in both fam.vars and data file
				curvars <- wealth.vars[list(years[iy]),which(names(wealth.vars)!="year"),with=FALSE]
				flog.debug("wealth curvars: ",curvars,capture=TRUE)
				curvars[,name := tolower(name)]
				curvars[,variable := tolower(variable)]
				curvars[,interview := tolower(interview)]

				setnames(tmp,tolower(names(tmp)))
				flog.debug("wealth tmp: ",head(tmp),capture=TRUE)
			
				# current set of variables
				codes <- c(curvars[,variable],curvars[,interview][[1]])
				flog.debug("wealth codes: ",codes,capture=TRUE)
				tmp   <- copy(tmp[,codes,with=FALSE])
				setnames(tmp,c(curvars[,name],"interview"))
				flog.debug("wealth tmp: ",head(tmp),capture=TRUE)
				setkey(tmp,interview)
				
				# merge m and wealthfile
				m <- merge(m,tmp,all.x=TRUE)

			}  # end wealth files
		}
	
		# note: a person who does not respond in wave x has an interview number in that wave, but NAs in the family file variables. remove those records.
		idx <- which(!is.na(unlist(fam.vars[list(years[iy])][,curnames,with=FALSE])))[1]	# index of first non NA variable
		m[,isna := is.na(m[,curnames[idx],with=FALSE])]
		m <- copy(m[isna == FALSE])
		m[,isna := NULL]
		# all remaining NAs are NAs which the user knows about and actually requested when specifying fam.vars
		# if (iy>1)	setcolorder(m,names(datas[[1]]))
		datas[[iy]] <- copy(m)


	}  # end year
	
	data2 <- rbindlist(datas,use.names=TRUE,fill=TRUE)	# glue together
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
		flog.info("balanced design reduces sample from %d to %d",n,nrow(data2))
	} else if (is.numeric(design)){
		n <- nrow(data2)
		data2[,enough := max(present) >= design,by=pid]
		data2 <- copy(data2[enough==TRUE])
		data2[,enough := NULL]
		flog.info("design choice reduces sample from %d to %d",n,nrow(data2))
	} else if (design=="all"){
		# do nothing
	}
	data2[,present := NULL]

	#     setkey(datas,pid,year)
	#     rm(ind)
	flog.info("End of build.panel")
	return(data2)
}

#' one year test, no ind file
#' 
#' @param dd Data Dictionary location. If NULL, 
#' use temp dir and force download
#' @export
small.test.noind <- function(dd=NULL){
 	cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
 	head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003))
 	famvars = data.frame(year=c(2003),age=head_age_var_name)
 	build.panel(fam.vars=famvars,datadir=dd)
}

#' one year test, ind file
#' 
#' @param dd Data Dictionary location. If NULL, 
#' use temp dir and force download
#' @export
small.test.ind <- function(dd=NULL){
 	cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
 	head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003))
 	educ = getNamesPSID("ER30323",cwf,years=2003)
 	famvars = data.frame(year=c(2003),age=head_age_var_name)
 	indvars = data.frame(year=c(2003),educ=educ)
 	build.panel(fam.vars=famvars,ind.vars=indvars,datadir=dd)
}

#' three year test, ind file
#' 
#' @param dd Data Dictionary location. If NULL, 
#' use temp dir and force download
#' @export
medium.test.ind <- function(dd=NULL){
	cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
	head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003,2005,2007))
 	educ = getNamesPSID("ER30323",cwf,years=c(2003,2005,2007))
	famvars = data.frame(year=c(2003,2005,2007),age=head_age_var_name)
	indvars = data.frame(year=c(2003,2005,2007),educ=educ)
	build.panel(fam.vars=famvars,ind.vars=indvars,datadir=dd)
}

#' three year test, no ind file
#' 
#' @param dd Data Dictionary location
#' @export
medium.test.noind <- function(dd=NULL){
  cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
  head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003,2005,2007))
  famvars = data.frame(year=c(2003,2005,2007),age=head_age_var_name)
  build.panel(fam.vars=famvars,datadir=dd)
}

#' three year test, ind file and one NA variable
#' 
#' @param dd Data Dictionary location. If NULL, 
#' use temp dir and force download
#' @export
medium.test.ind.NA <- function(dd=NULL){
	cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
	head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2003,2005,2007))
 	educ = getNamesPSID("ER30323",cwf,years=c(2003,2005,2007))
 	educ[2] = NA
	famvars = data.frame(year=c(2003,2005,2007),age=head_age_var_name)
	indvars = data.frame(year=c(2003,2005,2007),educ=educ)
	build.panel(fam.vars=famvars,ind.vars=indvars,datadir=dd,loglevel = DEBUG)
}

#' three year test, ind file and one NA variable and wealth
#' 
#' @param dd Data Dictionary location. If NULL, 
#' use temp dir and force download
#' @export
medium.test.ind.NA.wealth <- function(dd=NULL){
	cwf = openxlsx::read.xlsx(system.file(package="psidR","psid-lists","psid.xlsx"))
	head_age_var_name <- getNamesPSID("ER17013", cwf, years=c(2005,2007))
 	educ = getNamesPSID("ER30323",cwf,years=c(2005,2007))
 	educ[2] = NA
	r = system.file(package="psidR")
    w = fread(file.path(r,"psid-lists","wealthvars-small.txt"))
	famvars = data.frame(year=c(2005,2007),age=head_age_var_name)
	indvars = data.frame(year=c(2005,2007),educ=educ)
	build.panel(fam.vars=famvars,ind.vars=indvars,wealth.vars=w,datadir=dd,loglevel = DEBUG)
}


#' Build example PSID
#' 
#' @description Builds a panel from the full PSID dataset
#' @export
#' @param datadr string of the data directory
#' @param small logical TRUE if only use years 2013 and 2015.
#' @param wealth logical TRUE if want to use wealth supplements
#' @return a data.table with panel data
build.psid <- function(datadr="~/datasets/psid/",small=TRUE,wealth=FALSE){
  variable <- name <- NULL
  r = system.file(package="psidR")
  if (small){
    f = fread(file.path(r,"psid-lists","famvars-small.txt"))
    i = fread(file.path(r,"psid-lists","indvars-small.txt"))
    if (wealth){
      w = fread(file.path(r,"psid-lists","wealthvars-small.txt"))
    }
    
  } else {
    f = fread(file.path(r,"psid-lists","famvars.txt"))
    i = fread(file.path(r,"psid-lists","indvars.txt"))
    if (wealth){
      w = fread(file.path(r,"psid-lists","wealthvars.txt"))
    }
    
  }
  setkey(i,"name")
  setkey(f,"name")
  if (wealth) setkey(w,"name")
  
  i = dcast(i[,list(year,name,variable)],year~name)
  f = dcast(f[,list(year,name,variable)],year~name)
  if (wealth) {
    w = dcast(w[,list(year,name,variable)],year~name)
    d = build.panel(datadir=datadr,fam.vars=f,ind.vars=i,wealth.vars=w, heads.only = TRUE,sample="SRC",design="all")
    save(d,file="~/psid_wealth.RData")
  } else {
    d = build.panel(datadir=datadr,fam.vars=f,ind.vars=i, heads.only = TRUE,sample="SRC",design="all")
    save(d,file="~/psid_no_wealth.RData")
  }
	return(d)
}



#' GetPSID variables names from various years
#'
#' The user can specify one variable name from any year. This function
#' will find that variable's correct name in any of the years
#' specified by the user. If user does not specify the \code{years}
#' variable, return will represent all years in which variable was
#' present.
#'
#' This uses the psid.xlsx crosswalk file from UMich, which is
#' available at http://psidonline.isr.umich.edu/help/xyr/psid.xlsx. In the 
#' example, the package openxlsx's read.xlsx is used to import the crosswalk
#' file.
#'
#' Ask for one variable at a time.
#' @param aname A variable name in any of the PSID years
#' @param cwf A data.frame representation of the cross-walk file,
#'     (the psid.xlsx file).
#' @param years A vector of years. If NULL, all years in which that
#'     variable existed are returned
#' @return A vector of names, one for each year.
#' @author Paul Johnson <pauljohn@@ku.edu>
#' @export
#' @examples
#' # read UMich crosswalk from installed file
#' r = system.file(package="psidR")
#' cwf = openxlsx::read.xlsx(file.path(r,"psid-lists","psid.xlsx"))
#' 
#' # or download directly
#' # cwf <- read.xlsx("http://psidonline.isr.umich.edu/help/xyr/psid.xlsx")
#' 
#' # then get names with
#' getNamesPSID("ER17013", cwf, years = 2001)
#' getNamesPSID("ER17013", cwf, years = 2003)
#' getNamesPSID("ER17013", cwf, years = NULL)
#' getNamesPSID("ER17013", cwf, years = c(2005, 2007, 2009))
getNamesPSID <- function(aname, cwf, years = NULL){
    myvar <- which(cwf == aname, arr.ind=TRUE)
    ## variables that begin with Y
    ynames.all <- grep("^Y", colnames(cwf))

    if (is.null(years)){
        yearkeep <- ynames.all
    } else {
        yearkeep <- paste0("Y", years)
        yearkeep <- yearkeep[yearkeep %in% colnames(cwf)]
    }
    ovalue <- transpose(cwf[myvar[1], yearkeep, drop = FALSE])
    ovalue$V1
}

