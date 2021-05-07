
#' ID list for mergeing PSID
#'
#' @description this list is taken from http://ideas.repec.org/c/boc/bocode/s457040.html
#' @details this function hardcodes the PSID variable names of "interview number" from both family and individual file for each wave, as well as "sequence number", "relation to head" and numeric value x of that variable such that "relation to head" == x means the individual is the head. Varies over time.
makeids <- function(){

	id.list <- data.table(year=c(1968:1997,seq(1999,2019,by=2)))
	id.list$ind.interview <- c("ER30001","ER30020","ER30043","ER30067",
									"ER30091","ER30117","ER30138","ER30160", 
									"ER30188","ER30217","ER30246","ER30283",
									"ER30313","ER30343","ER30373","ER30399",
									"ER30429","ER30463","ER30498","ER30535",
									"ER30570","ER30606","ER30642","ER30689",
									"ER30733","ER30806","ER33101","ER33201",
									"ER33301","ER33401","ER33501","ER33601",
									"ER33701","ER33801","ER33901","ER34001",
									"ER34101","ER34201","ER34301","ER34501","ER34701")

	id.list$ind.seq <- c(NA,"ER30021","ER30044","ER30068","ER30092","ER30118","ER30139",
						 "ER30161","ER30189","ER30218","ER30247","ER30284","ER30314", 
						 "ER30344","ER30374","ER30400","ER30430","ER30464","ER30499", 
						 "ER30536","ER30571","ER30607","ER30643","ER30690","ER30734", 
						 "ER30807","ER33102","ER33202","ER33302","ER33402","ER33502", 
						 "ER33602","ER33702","ER33802","ER33902","ER34002","ER34102",
						 "ER34202","ER34302","ER34502","ER34702")
	
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
						  "ER34103",
						  "ER34203",
						  "ER34303",
						  "ER34503",
			     			  "ER34703")
						  
	# numeric code for "i am the head"
	id.list$ind.head.num <- c(rep(1,15),rep(10,26))

	id.list$fam.interview <- c("V3"      , "V442"    , "V1102"   , "V1802"   , "V2402"    , "V3002" ,
							   "V3402"   , "V3802"   , "V4302"   , "V5202"   , "V5702"    ,
							   "V6302"   , "V6902"   , "V7502"   , "V8202"   , "V8802"    ,
							   "V10002"  , "V11102"  , "V12502"  , "V13702"  , "V14802"   ,
							   "V16302"  , "V17702"  , "V19002"  , "V20302"  , "V21602"   ,
							   "ER2002"  , "ER5002"  , "ER7002"  , "ER10002" , "ER13002"  ,
							   "ER17002" , "ER21002" , "ER25002" , "ER36002" , "ER42002"  , 
							   "ER47302" , "ER53002" , "ER60002" , "ER66002" , "ER72002")
	id.list$stratum <- rep("ER31996",nrow(id.list))
	setkey(id.list,year)
	return(id.list)
}

makeids.wealth <- function(){

	id = data.table(year=c(1984, 1989, 1994, 1999, 2001, 2003, 2005, 2007),interview = c("S101","S201","S301","S401","S501","S601","S701","S801"))
	setkey(id,year)
	return(id)

}


	  
	  
#' get.psid connects to PSID database and downloads into Rda
#'
#' see \url{http://asdfree.com/} for other usage and \url{https://stackoverflow.com/questions/15853204/how-to-login-and-then-download-a-file-from-aspx-web-pages-with-r}
#' @author Anthony Damico <ajdamico@@gmail.com>
#' @param file string psid file number
#' @param name string of filename on disc
#' @param params postForm{RCurl} parameters
#' @param curl postForm{RCurl} curl handle
get.psid <- function( file , name , params , curl ){

		html = postForm('http://simba.isr.umich.edu/u/Login.aspx', .params = params, curl = curl)
		
		if ( !grepl('Logout', html) ) stop( 'no longer logged in' )

	
		tf <- tempfile() ; td <- tempdir()
		
		flog.info('downloading file %s',name)
		file <- getBinaryURL( paste0( "http://simba.isr.umich.edu/Zips/GetFile.aspx?file=" , file ) , curl = curl )
		writeBin( file , tf )
		z <- unzip( tf , exdir = td )
		fn <- z[ grepl( ".txt" , tolower( z ) , fixed = TRUE ) & ! grepl( "_vdm|readme|doc|errata" , tolower( z ) ) ]
		sas_ri <- z[ grepl( '.sas' , z , fixed = TRUE ) ]

		flog.info('now reading and processing SAS file %s into R',name)

		#SAScii version check SAScii_fork
		# if (!exists("SAScii_fork",mode="function")){
		# 	warning("you may run into trouble now. There was a change of file format on some PSID family files. \n If you get the an error \n toupper(SASinput) \n
		# 		then you need to re-install the SAScii package from my github fork at \n
		# 		https://github.com/floswald/SAScii \n
		# 		an easy way to do this is to use the devtools package. then:
		# 		install_github('floswald/SAScii')")
		# }


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
#' @param N number of people in each wave
#' @param N.attr number of people lost to attrition
#' @return list with (fake) individual index file IND2009ER and
#' family files for 1985 and 1986
#' @export
testPSID <-function(N=100,N.attr = 0){

	# need to bind some vars for R CHECK:
	ER30001 <- ER30002 <-intnum86 <- intnum85 <- Money85 <- Money86 <- age85 <- age86 <- smpls <- NULL

	stopifnot(N %% 4 == 0)
	# for sake of illustration, suppose the PSID has a total
	# of 2N people (i.e. N are neither in year1 nor year2, 
	# but in some other years)
	# draw interview numbers from realistic range: [1,9308]

	# in each wave, a quarter of observations is from
	# 1) core sample: interview number < 3000
	# 2) immigrant sample: interview number in [3000,5000)
	# 3) poor sample: interview number in [5000,7000)
	# 4) latino sample: interview number in [7000,9308)

	smpl <- ceiling(1:N / (N/4))
	IND2009ER <- data.table(smpls=c(smpl,smpl))  # get 2*N inds
    IND2009ER[, ER30001 := 0L]
	IND2009ER[smpls==1, ER30001 := sample(1:2999,size=sum(smpls==1))]
	IND2009ER[smpls==2, ER30001 := sample(3001:4999,size=sum(smpls==2))]
	IND2009ER[smpls==3, ER30001 := sample(5001:6999,size=sum(smpls==3))]
	IND2009ER[smpls==4, ER30001 := sample(7001:9308,size=sum(smpls==4))]

    IND2009ER[,c("intnum85","intnum86") := lapply(1:2,function(x) sample(1:(2*N)))]
	# IND2009ER[smpls==1,c("intnum85","intnum86") := lapply(1:2,function(x) sample(1:2999,size=sum(smpls==1)))]
	# IND2009ER[smpls==2,c("intnum85","intnum86") := lapply(1:2,function(x) sample(3001:4999,sum(smpls==2)))]
	# IND2009ER[smpls==3,c("intnum85","intnum86") := lapply(1:2,function(x) sample(5001:6999,sum(smpls==3)))]
	# IND2009ER[smpls==4,c("intnum85","intnum86") := lapply(1:2,function(x) sample(7001:9308,sum(smpls==4)))]
  
    # only N invidividuals show up in 1985 though.
    IND2009ER[sample(1:(2*N),size=N),c("intnum85") := 0]
  
    # and there is potential attrition from 1985 to 1986
	if (N.attr>0){
	    out = sample(which(IND2009ER[["intnum85"]] != 0),size=N.attr,replace=FALSE)
	    IND2009ER[out,"intnum86" := 0]
	}

    # add 1968 person id
	IND2009ER[,ER30002 := sample(1:(2*N))]

    # also need relationship to head in each year in the index
	# 50% prob of being head in year1
	IND2009ER$ER30465 <- sample(rep(c(10,20),c(N,N)),size=2*N,replace=TRUE)	
	IND2009ER$ER30500 <- sample(rep(c(10,20),c(N,N)),size=2*N,replace=TRUE)
	
	# as well as the sequence number: 1 for current heads, > 50 for movers
	# 75% prob of being current head
	IND2009ER$ER30464 <- sample(rep(c(1,20),c(ceiling(0.75*2*N),ceiling(0.5*N))),size=2*N,replace=TRUE)	
	IND2009ER$ER30499 <- sample(rep(c(1,20),c(ceiling(0.75*2*N),ceiling(0.5*N))), size=2*N,replace=TRUE)
	# and a survey weight
	IND2009ER$ER30497 <- runif(2*N)
	IND2009ER$ER30534 <- runif(2*N)
  
    # you would download files like those two data.frames:
    fam85 <- data.table(intnum85 = IND2009ER[intnum85>0,sample(intnum85,size=N)],Money85=rlnorm(n=N,10,1),age85=sample(20:80,size=N,replace=TRUE))
	fam86 <- data.table(intnum86 = IND2009ER[intnum86>0,sample(intnum86,size=sum(intnum86>0))],Money86 = rlnorm(n=IND2009ER[,sum(intnum86>0)],10,1),age86 = sample(20L:80L,size=IND2009ER[,sum(intnum86>0)],replace=TRUE))
  
    continuers85 = IND2009ER[intnum85>0 & intnum86>0][["intnum85"]]
	continuers86 = IND2009ER[intnum85>0 & intnum86>0][["intnum86"]]
  
    for (i in 1:nrow(fam86)){
        if (fam86[i,intnum86] %in% continuers86){
        int86 = fam86[i,intnum86]
        fam86[i,Money86 := fam85[intnum85==IND2009ER[intnum86==int86,intnum85], Money85 + rnorm(1,500,30)]]
        fam86[i,age86 := fam85[intnum85==IND2009ER[intnum86==int86,intnum85], age85 + 1]]
        }
    }
	
	# assign correct PSID varname of "family interview 1985/86"
	setnames(fam85,"intnum85", "V11102")
	setnames(fam86,"intnum86","V12502")
	
    # same on index file
	setnames(IND2009ER,c("intnum85","intnum86"), c("ER30463","ER30498"))

    return(list(famvars1985=fam85,famvars1986=fam86,IND2019ER=IND2009ER))
}
