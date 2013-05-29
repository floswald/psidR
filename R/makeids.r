
#' ID list for mergeing PSID
#'
#' @description this list is taken from http://ideas.repec.org/c/boc/bocode/s457040.html
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
