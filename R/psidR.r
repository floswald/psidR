

#' Build PSID panel data set
#' 
#' @description Builds a panel data set in wide format with id variables \code{personID} and \code{period} from individual PSID family files.
#' @details 
#' \itemize{
#' \item Default behaviour: takes family files for specified years in folder \code{datadir} and merges using the id information in \code{ind.vars}, which must be in the same directory. 
#' \item Alternative: location of files may be manually specified. 
#' \item in both: the user can change subsetting criteria as well as sample designs. 
#' \item Merge: If there are N individuals in each of T waves, the individual file contains NT rows. If an individual has non-response in a given wave, values in the family file are NA. the variables \code{interview number} in each family file map to the \code{interview number} variable of a given year in the individual file. 
#' \item Accepted input data are stata format .dta or .csv files.
#' \item See examples for usage. Note that examples use fake data, i.e. this is not real PSID data. Data dictionaries do not work therefore.
#' \item Similar in usage to stata module \code{psiduse}. 
#' }
#' @param datadir directory containing family files ("FAMyyyy.dta") and individual file ("IND2009ER.dta"). fixed data format required.
#' @param fam.vars data.frame of variable to retrieve from family files. see example for required format.
#' @param ind.vars optional data.frame of non-default variables to get from individual file.
#' @param fam.files optional character vector of file locations of family files
#' @param ind.file optional character vector of file location of individual file
#' @param heads.only logical TRUE if user wants household heads only. if FALSE, data contains a row with value of "relation to head" variable.
#' @param core logical TRUE if user wants core sample only. if FALSE, data will oversample poverty sample.
#' @param design either character "balanced" or "all" or integer. "Balanced" means only individuals who appear in each wave are considered. "All" means all are taken. An integer value stands for minimum consecutive years of participation, i.e. design=3 means at least 3 waves.
#' @param verbose logical TRUE if you want verbose output.
#' @return
#' \item{data}{resulting \code{data.table}. the variable \code{pid} is the unique person identifier, constructed from ID1968 and pernum.}
#' \item{dict}{data dictionary if stata data was supplied, NULL else}
#' @export
#' @example inst/examples.r
build.panel <- function(datadir,fam.vars,ind.vars=NULL,fam.files=NULL,ind.file=NULL,heads.only=TRUE,core=TRUE,design="balanced",verbose=FALSE){
	
	stata   <- FALSE
	csv     <- FALSE
	default <- FALSE

	years <- fam.vars$year

	# work out file types
	# -------------------

	if (is.null(datadir)){
		# user specifies file locations
		stopifnot(!is.null(fam.files) | !is.null(ind.file))
		stopifnot(is.character(fam.files) | is.character(ind.file))

		# get file type: .dta or .csv
		# load ind file. load fam file only at the point we need it to keep memory small.
		if (tail(strsplit(fam.files[1],"\\.")[[1]],1) == "dta"){
			stata    <- TRUE
			ind      <- read.dta(file=ind.file)
			ind.dict <- data.frame(code=names(ind),label=attr(ind,"var.labels"))
			ind      <- data.table(ind)
		} else if (tail(strsplit(fam.files[1],"\\.")[[1]],1) == "csv") {
			csv      <- TRUE
			ind      <- fread(input=ind.file)
			warning('no data dictionary from csv file')
		}
	} else {
		# default to stata data in datadir
		default <- TRUE
		fam.dta <- paste(datadir,"/FAM",years,".dta",sep="")
		ind.dta <- paste(datadir,"/IND2009ER.dta",sep="")	# needs to be updated with next data delivery.
		ind      <- read.dta(file=ind.dta)
		ind.dict <- data.frame(code=names(ind),label=attr(ind,"var.labels"))
		ind      <- data.table(ind)
	}
	
	if (verbose){
		if (!is.null(datadir)) {
			cat('loaded individual file:',ind.dta,'\n')
		} else {
			cat('loaded individual file:',ind.file,'\n')
		}
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
		if (stata) {
			tmp             <- read.dta(file=fam.files[iy])
			fam.dicts[[iy]] <- data.frame(code=names(tmp),label=attr(tmp,"var.labels"))
			tmp             <- data.table(tmp)
		} else if (csv) {
			tmp <- fread(input=fam.files[iy])
			fam.dicts[[iy]] <- NULL
		} else if (default) {
		   	tmp             <- read.dta(file=fam.dta[iy])
			fam.dicts[[iy]] <- data.frame(code=names(tmp),label=attr(tmp,"var.labels"))
			tmp             <- data.table(tmp)
		}

		if (verbose){
			if (default) {
				cat('loaded family file:',fam.dta,'\n')
			} else {
				cat('loaded family file:',fam.files[iy],'\n')
			}
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


#' Convert factor to character
#'
#' @description helper function to convert factor to character in a data.table
#' @export
make.char <- function(x){
	if (is.factor(x)){
		return(as.character(x))
	} else {
		return(x)
	}
}






#' ID list for mergeing PSID
#'
#' @description this list is taken from http://ideas.repec.org/c/boc/bocode/s457040.html
#' @export
#' @details this function hardcodes the PSID variable names of "interview number" from both family and individual file for each wave, as well as "sequence number", "relation to head" and numeric value x of that variable such that "relation to head" == x means the individual is the head. Varies over time.
makeids <- function(){

	id.list <- data.table(year=c(1968:1997,seq(1999,2009,by=2)))
	id.list$ind.interview <- c("ER30001","ER30020","ER30043","ER30067",
									"ER30091","ER30117","ER30138","ER30160", 
									"ER30188","ER30217","ER30246","ER30283",
									"ER30313","ER30343","ER30373","ER30399",
									"ER30429","ER30463","ER30498","ER30535",
									"ER30570","ER30606","ER30642","ER30689",
									"ER30733","ER30806","ER33101","ER33201",
									"ER33301","ER33401","ER33501","ER33601",
									"ER33701","ER33801","ER33901","ER34001")

	id.list$ind.seq <- c(NA,"ER30021","ER30044","ER30068","ER30092","ER30118","ER30139",
						 "ER30161","ER30189","ER30218","ER30247","ER30284","ER30314", 
						 "ER30344","ER30374","ER30400","ER30430","ER30464","ER30499", 
						 "ER30536","ER30571","ER30607","ER30643","ER30690","ER30734", 
						 "ER30807","ER33102","ER33202","ER33302","ER33402","ER33502", 
						 "ER33602","ER33702","ER33802","ER33902","ER34002")
	
	# name of variable "relationship to head"
	id.list$ind.head <- c("ER30003",
						  "ER30022",
						  "ER30045",
						  "ER30069",
						  "ER30093",
						  "ER30119",
						  "ER30140",
						  "ER30162",
						  "ER30190",
						  "ER30219",
						  "ER30248",
						  "ER30285",
						  "ER30315",
						  "ER30345",
						  "ER30375",
						  "ER30401",
						  "ER30431",
						  "ER30465",
						  "ER30500",
						  "ER30537",
						  "ER30572",
						  "ER30608",
						  "ER30644",
						  "ER30691",
						  "ER30735",
						  "ER30808",
						  "ER33103",
						  "ER33203",
						  "ER33303",
						  "ER33403",
						  "ER33503",
						  "ER33603",
						  "ER33703",
						  "ER33803",
						  "ER33903",
						  "ER34003")
						  
	# numeric code for "i am the head"
	id.list$ind.head.num <- c(rep(1,15),rep(10,21))

	id.list$fam.interview <- c("V3"      , "V442"    , "V1102"   , "V1802"   , "V2402"    , "V3002" ,
							   "V3402"   , "V3802"   , "V4302"   , "V5202"   , "V5702"    ,
							   "V6302"   , "V6902"   , "V7502"   , "V8202"   , "V8802"    ,
							   "V10002"  , "V11102"  , "V12502"  , "V13702"  , "V14802"   ,
							   "V16302"  , "V17702"  , "V19002"  , "V20302"  , "V21602"   ,
							   "ER2002"  , "ER5002"  , "ER7002"  , "ER10002" , "ER13002"  ,
							   "ER17002" , "ER21002" , "ER25002" , "ER36002" , "ER42002")
	setkey(id.list,year)
	return(id.list)
}


	  
	  
	  
	  
	  
	  








