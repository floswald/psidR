
#' ID list for mergeing PSID
#'
#' @description this list is taken from http://ideas.repec.org/c/boc/bocode/s457040.html
#' @details this function hardcodes the PSID variable names of "interview number" from both family and individual file for each wave, as well as "sequence number", "relation to head" and numeric value x of that variable such that "relation to head" == x means the individual is the head. Varies over time.
makeids <- function(){

	id.list <- data.table(year=c(1968:1997,seq(1999,2011,by=2)))
	id.list$ind.interview <- c("ER30001","ER30020","ER30043","ER30067",
									"ER30091","ER30117","ER30138","ER30160", 
									"ER30188","ER30217","ER30246","ER30283",
									"ER30313","ER30343","ER30373","ER30399",
									"ER30429","ER30463","ER30498","ER30535",
									"ER30570","ER30606","ER30642","ER30689",
									"ER30733","ER30806","ER33101","ER33201",
									"ER33301","ER33401","ER33501","ER33601",
									"ER33701","ER33801","ER33901","ER34001","ER34101")

	id.list$ind.seq <- c(NA,"ER30021","ER30044","ER30068","ER30092","ER30118","ER30139",
						 "ER30161","ER30189","ER30218","ER30247","ER30284","ER30314", 
						 "ER30344","ER30374","ER30400","ER30430","ER30464","ER30499", 
						 "ER30536","ER30571","ER30607","ER30643","ER30690","ER30734", 
						 "ER30807","ER33102","ER33202","ER33302","ER33402","ER33502", 
						 "ER33602","ER33702","ER33802","ER33902","ER34002","ER34102")
	
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
						  "ER34003",
						  "ER34103")
						  
	# numeric code for "i am the head"
	id.list$ind.head.num <- c(rep(1,15),rep(10,22))

	id.list$fam.interview <- c("V3"      , "V442"    , "V1102"   , "V1802"   , "V2402"    , "V3002" ,
							   "V3402"   , "V3802"   , "V4302"   , "V5202"   , "V5702"    ,
							   "V6302"   , "V6902"   , "V7502"   , "V8202"   , "V8802"    ,
							   "V10002"  , "V11102"  , "V12502"  , "V13702"  , "V14802"   ,
							   "V16302"  , "V17702"  , "V19002"  , "V20302"  , "V21602"   ,
							   "ER2002"  , "ER5002"  , "ER7002"  , "ER10002" , "ER13002"  ,
							   "ER17002" , "ER21002" , "ER25002" , "ER36002" , "ER42002"  , "ER47302")
	setkey(id.list,year)
	return(id.list)
}


	  
	  
#' get.psid connects to PSID database and downloads into Rda
#'
#' see \url{http://www.asdfree.com/} for other usage and \url{http://stackoverflow.com/questions/15853204/how-to-login-and-then-download-a-file-from-aspx-web-pages-with-r}
#' @author Anthony Damico <ajdamico@@gmail.com>
#' @param file string psid file number
#' @param name string of filename on disc
#' @param params postForm{RCurl} parameters
#' @param curl postForm{RCurl} curl handle
get.psid <- function( file , name , params , curl ){

		html = postForm('http://simba.isr.umich.edu/u/Login.aspx', .params = params, curl = curl)
		
		if ( !grepl('Logout', html) ) stop( 'no longer logged in' )

	
		tf <- tempfile() ; td <- tempdir()
		
		file <- getBinaryURL( paste0( "http://simba.isr.umich.edu/Zips/GetFile.aspx?file=" , file ) , curl = curl )
		writeBin( file , tf )
		z <- unzip( tf , exdir = td )
		fn <- z[ grepl( ".txt" , tolower( z ) , fixed = TRUE ) & ! grepl( "_vdm|readme|doc|errata" , tolower( z ) ) ]
		sas_ri <- z[ grepl( '.sas' , z , fixed = TRUE ) ]

		cat('now reading and processing SAS file',name,'into R\n')
		x <- read.SAScii( fn , sas_ri )

		save( x , file = paste0( name , '.rda' ) )
	
		file.remove( tf , z )
	
		rm( x )
		
		gc()

		TRUE
	}



#' Convert factor to character
#'
#' @param x a \code{factor}
#' @return a character
#' @description helper function to convert factor to character in a data.table
make.char <- function(x){
	if (is.factor(x)){
		return(as.character(x))
	} else {
		return(x)
	}
}



#' Create a test PSID dataset
#'
#' makes artifical PSID data with variables \code{age} and \code{income}
#' for two consecutive years 1985 and 1986.
#' @param N how many people per wave
#' @param N.attr number of people lost to attrition
#' @return list with (fake) individual index file IND2009ER and
#' family files for 1985 and 1986
#' @export
testPSID <-function(N=10,N.attr = 0){

	interview85 <- NULL  # R CHECK

  # you would download files like those two data.frames:
	fam1985 <- data.frame(interview85 = 1:N,
	                 Money85=rlnorm(n=N,10,1),
	                 age85=sample(20:80,size=N,replace=TRUE))
	fam1986 <- data.frame(interview86 = sample(1:N,size=N,replace=FALSE))
	fam1986$Money86 <- subset(fam1985,interview85 %in% fam1986$interview86)$Money85+rnorm(nrow(fam1986),500,30)
	fam1986$age86 <- subset(fam1985,interview85 %in% fam1986$interview86)$age85+1

	# assign correct PSID varname of "family interview 1985/86"
	names(fam1985)[1] <- "V11102"	
	names(fam1986)[1] <- "V12502"
	
	# construct an Individual index file: that would be IND2009ER
	# needs to have a 1968 interview and person number (ER30001, ER30002) 
	# and an indicator for whether from core etc, 
	# as well as the interview number for each year
	# 
	# for sake of illustration, suppose the PSID has a total
	# of 2N people (i.e. N are neither in year1 nor year2, 
	# but in some other years)
	IND2009ER <- data.frame(ER30001=sample((2*N):(4*N),size=2*N),
	                        ER30002=sample(1:(2*N),size=2*N))
	
	# if a person is observed, they have an interview number 
	# in both years. if not observed, its zero. 
    tmp = data.frame(interview85=fam1985[,1], interview86=fam1986[,1])
    tmp = rbind(tmp,data.frame(interview85=rep(0,N), interview86=rep(0,N)))
	IND2009ER <- cbind(IND2009ER,tmp[sample(1:(2*N)),])
    if (N.attr>0){
        out = sample(which(IND2009ER[["interview85"]] != 0),size=N.attr,replace=FALSE)
        IND2009ER[out,"interview86"] <- 0
    }
	
	names(IND2009ER)[3:4] <- c("ER30463","ER30498")
	
	# also need relationship to head in each year in the index
	# 50% prob of being head in year1
	IND2009ER$ER30465 <- sample(c(10,20),prob=c(0.5,0.5),
	                            size=2*N,replace=TRUE)	
	IND2009ER$ER30500 <- sample(c(10,20),prob=c(0.9,0.1),
	                           size=2*N,replace=TRUE)
	
	# as well as the sequence number: 1 for current heads, > 50 for movers
	# 90% prob of being current head
	IND2009ER$ER30464 <- sample(c(1,20),prob=c(0.95,0.05),
	                            size=2*N,replace=TRUE)	
	IND2009ER$ER30499 <- sample(c(1,20),prob=c(0.95,0.05),
	                           size=2*N,replace=TRUE)
	# and a survey weight
	IND2009ER$ER30497 <- runif(2*N)
	IND2009ER$ER30534 <- runif(2*N)
	IND2009ER
	
	# setup the ind.vars data.frame
	indvars <- data.frame(year=c(1985,1986),ind.weight=c("ER30497","ER30534"))
  
  return(list(famvars1985=fam1985,famvars1986=fam1986,IND2009ER=IND2009ER))
}
