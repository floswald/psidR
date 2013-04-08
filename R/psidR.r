

#' psidR: package to built panel data sets form PSID raw data

library(foreign)
library(data.table)

#' main function
#' @desc takes family files for specified years in folder datadir and merges using the id information in ind.vars, which must be in the same directory.
#' @param datadir directory containing family files ("FAMyyyy.dta") and individual file ("IND2009ER.dta")
#' @param fam.vars named list of variable names to retrieve from family files
#' @param ind.vars named list of variable names to retrieve form individual file
#' @param years numeric vector of years to consider
build.panel <- function(datadir,fam.vars,ind.vars=NULL,fam.files=NULL,ind.file=NULL,heads.only,core,design,missing.vars,verbose){

	
fam.vars = data.frame(year=c(1999,2001),cleaning=c("ER13027","ER17030"),oth.services=c("ER13028","ER17031"))
	
	years <- fam.vars$year

	stopifnot(is.data.frame(fam.vars))
	if (!is.data.table(fam.vars)) {
		fam.vars <- data.table(fam.vars)
		fam.vars <- copy(fam.vars[,lapply(.SD,make.char)])
		setkey(fam.vars,year)
	}

	#     if (!is.null(fam.files) & !is.null(ind.file)){
	#         if (strsplit(fam.files[1],"\\.")[[1]][2] == "dta"){
	#             stata    <- TRUE
	#             ind      <- read.dta(file=ind.file)
	#             ind.dict <- data.frame(code=names(ind),label=attr(ind,"var.labels"))
	#             ind      <- data.table(ind)
	#         
	#         } else if (strsplit(fam.files[1],"\\.")[[1]][2] == "csv") {
	#             csv      <- TRUE
	#             ind      <- fread(file=ind.file)
	#             warning('no data dictionary from csv file')
	# 
	#         }
	#     } else {
		# default to stata data in datadir
	#         default <- TRUE
		fam.dta <- paste(datadir,"/FAM",years,".dta",sep="")
		ind.dta <- paste(datadir,"/IND2009ER.dta",sep="")	# needs to be updated with next data delivery.
		ind      <- read.dta(file=ind.dta)
		ind.dict <- data.frame(code=names(ind),label=attr(ind,"var.labels"))
		ind      <- data.table(ind)
		#     }

	# data dictionaries
	fam.dicts <- vector("list",length(years))

	if (verbose){
		cat('loaded individual file:',ind.dta,'\n')
		cat('total memory load in MB:\n')
		vvs = ceiling(object.size(ind)/1024^2)
		print(as.numeric(vvs))
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

	# loop over years
	for (iy in 1:length(years)){
		
		if (verbose) {
			cat('=============================================\n')
			cat('currently working on data for year',years[iy],'\n')
		}

		# keeping only relevant columns from individual file
		# subset for core sample and heads only
		curr <- ids[.(years[iy])]
		ind.subsetter <- as.character(curr[,list(ind.interview,ind.head)])	# keep from ind file
		def.subsetter <- c("ER30001","ER30002")	# must keep those in all years

		yind            <- copy(ind[,c(def.subsetter,unique(c(ind.subsetter,ind.vars[[iy]]))),with=FALSE])
		if (core) {
		   n <- nrow(yind)
		   yind <- copy(yind[ER30001>2930])	# individuals 1-2930 are from poor sample
		   if (verbose){
			   cat('full sample has',n,'obs\n')
			   cat('dropping non-core individuals leaves',nrow(yind),'obs\n')
		   }
		}
		if (heads) {
		   n <- nrow(yind)
		   yind <- yind[,headyes := yind[,curr[,ind.head],with=FALSE]==curr[,ind.head.num]]
		   yind <- copy(yind[headyes==TRUE])
		   if (verbose){
			   cat('sample with heads has',n,'obs\n')
			   cat('dropping non-heads leaves',nrow(yind),'obs\n')
		   }
		   yind[,c(curr[,ind.head],"headyes") := NULL]
		}
		# set names on individual index
		setnames(yind,c("ID1968","pernum","interview"))
		yind[,pid := ID1968*1000 + pernum]	# unique person identifier
		setkey(yind,interview)

		# bring in family files, subset them
		# load data for current year, make data dictionary for subsets and save data as data.table
		#         if (stata) {
		#                tmp            <- read.dta(file=fam.files[iy])
		#             fam.dicts[[f]] <- data.frame(code=names(tmp),label=attr(tmp,"var.labels"))
		#             tmp            <- data.table(tmp)
		#         } else if (csv) {
		#             tmp <- fread(file=fam.files[iy])
		#         } else if (default) {
		   	tmp            <- read.dta(file=fam.dta[iy])
			fam.dicts[[iy]] <- data.frame(code=names(tmp),label=attr(tmp,"var.labels"))
			tmp            <- data.table(tmp)
			#         }

		if (verbose){
			cat('loaded family files:',fam.dta,'\n')
			cat('current memory load in MB:\n')
			vs = sapply(tmp, function(x) ceiling(object.size(x)/1024^2))
			print(vs)
		}
	
		# add compulsory vars to fam.vars
		fam.vars[,interview := ids[fam.vars][,fam.interview]]


		tmp <- copy(tmp[,as.character(fam.vars[.(years[iy]),which(names(fam.vars)!="year"),with=FALSE]),with=FALSE])
		tmpnames <- names(fam.vars)[-which(names(fam.vars)=="year")]
		setnames(tmp,tmpnames)	# call interview "interview"
		#         tmp[,year := years[iy]]
		setkey(tmp,interview)

		m <- copy(tmp[yind])

		#TODO name columns of tmp with names in fam.vars before merge to avoid confusion.

		setnames(m,tmpnames,paste(tmpnames,years[1],sep="_"))







	}
	





	}

	}


	# subset ind.vars for heads only

	# merge is done on fam.interview ind.interview


	# pid is 
	# gen pid = ER30001 * 1000 + ER30002

	# subset according to fam.vars and ind.vars
	# ind.vars will hold a hardcoded list of variables that is needed to make a panel
	# if the user adds more vars, get them as well.




}



make.char <- function(x){
	if (is.factor(x)){
		return(as.character(x))
	} else {
		return(x)
	}
}






#' makes id list for merges
#' this list is taken from http://ideas.repec.org/c/boc/bocode/s457040.html
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


	  
	  
	  
	  
	  
	  








