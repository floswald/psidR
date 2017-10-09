

#' This function allows to build a data table containing the original PSID variable names (of type AB12345) for each desired year and checks
#' on the psidonline web site whether a variable is available or not. Information is obtained on https://simba.isr.umich.edu/VS/s.aspx
#' through the package XML.
#'
#' Instead of building the data table by hand and look up each variable on the online dictionary to understand in which year it is available, 
#' it is sufficient to look at the dictionary for a given wave, select some variables of interest and type in their original name for a given 
#' year y, as well as a list of years that one wants to be included in the panel. If the variable is not available for a given year, the 
#' function automatically types NA in the corresponding cell.
#'
#' Note that, in the list of variables provided, one could include a variable from year y and a variable from year y'. This is allowed and 
#' irrelevant to the functioning of the function. This could be useful if one wants to recover some variables that are only present in some
#' specific waves. Only remark: if one types in the same variable for both year y and year y', it is going to show up twice in the final data 
#' table.
#'
#' The following options are also available:
#' - if scrape_name equals TRUE, then the function recovers the variable names from the online dictionary and edits them to avoid the presence
#' of special characters such as "/", "#" or "-"; the dictionary names suggested by psidonline might be long and cryptic and are NOT always
#' consistent over years, even if the coding of the variable is the same;
#' - if the vector varnames is provided and scrape_name equals FALSE, then the function names the variables after the user-provided labels;
#' varnames must be of the same length as input vector vars.
#'
#' Here is shown how to obtain the ind and fam data tables used in the main example function. In this example, the function is required to
#' scrape names from psidonline.
#'
#' years=c(1968:1997, seq(from=1999, to=2011, by=2))
#' i <- c("ER30046", "ER30052", "ER30066") # Variable original names for the 1970 wave
#' f <- c("V1567", "V1514") # Variable original names for the 1970 wave
#' ind <- get.variable.names(i, years, scrape_name=T)
#' fam <- get.variable.names(f, years, scrape_name=T)
#'
#' > head(fam)
#'      year heads_avg_hrly_ern tot_fu_mon_inc_ov414
#' 1968 1968               V337                  V81
#' 1969 1969               V871                 V529
#' 1970 1970              V1567                V1514
#' 1971 1971              V2279                V2226
#' 1972 1972              V2906                V2852
#' 1973 1973              V3275                V3256
#' 
#' > head(ind)
#'      year age_of_individual grade_finished individual_weight
#' 1968 1968           ER30004        ER30010           ER30019
#' 1969 1969           ER30023             NA           ER30042
#' 1970 1970           ER30046        ER30052           ER30066
#' 1971 1971           ER30070        ER30076           ER30090
#' 1972 1972           ER30094        ER30100           ER30116
#' 1973 1973           ER30120        ER30126           ER30137
#' 

get.variable.names <- function(vars, years, scrape_name=F, varnames = NULL){
  
  # Check if the function is required to scrape variable names from the dictionary or not
  using_names = !(is.null(varnames))
  
  if (using_names & scrape_name) {
    scrape_name = F
    warning("The option scrape_name has been set to FALSE. If you require the code to scrape variable names, let the varnames option set to default value NULL.\n")
  }

  # Create data table where to store the original variable names wave by wave  
  vardf <- as.data.table(matrix(0, ncol = length(years), nrow = 0))
  colnames(vardf) <- as.character(years)
  
  # Go through each variable
  for (varindex in vars) {

    # Get html code from psidonline
    h <- htmlParse(paste0("http://simba.isr.umich.edu/cb.aspx?vList=",varindex))
    
    # Recover and edit the explicit variable name (replace special characters as they might be complicated to handle with)
    # NB: the list of special characters might need to be expanded
    if (scrape_name) {
      h.varname <- xmlValue(getNodeSet(h, "//td")[[3]])
      h.varname <- gsub("\\s+", " ", h.varname)
      h.varname <- gsub(" ", "_", h.varname)
      h.varname <- gsub("\\\\|\\-|/|\\?|\\(|\\)", "_", h.varname)
      h.varname <- gsub("\\#", "n", h.varname)
      h.varname <- tolower(h.varname)
      if (grepl("^[[:digit:]]*$",substr(h.varname,nchar(h.varname)-1,nchar(h.varname))) & substr(h.varname,nchar(h.varname)-2,nchar(h.varname)-2) == "_") {
        h.varname = substr(h.varname,1,nchar(h.varname)-3)
      }
    }

    # Locate a list of type [99]AB1234, then build a matrix with column of type c("99", "AB1234")
    for (node in rev(getNodeSet(h, "//td"))) {

      if (substr(xmlValue(node),1,1) == "[" & substr(xmlValue(node),4,4) == "]") {

        # Get the list
        h.list <- xmlValue(node)
        h.list <- strsplit(gsub("\\[", "", h.list, perl=T), " ")
        h.list <- lapply(h.list[[1]], function(x) strsplit(x,"]"))
        h.list <- as.matrix(matrix(unlist(h.list), nrow=length(unlist(h.list[1]))))

        # Change year from a two-digit character to a four-digit    
        sel <- which(nchar(h.list) == 2 & (substr(h.list,1,1) < 6))
        h.list[sel] <- paste0("20",h.list[sel])
        sel <- which(nchar(h.list) == 2 & (substr(h.list,1,1) > 5))
        h.list[sel] <- paste0("19",h.list[sel])
        h.df <- as.data.table(t(h.list[2,]), stringsAsFactors = F)
        colnames(h.df) <- h.list[1,]
        
        # Trim years that are available but not required
        for (y in h.list[1,]) {
          if (!(y %in% years)) h.df[[y]] = NULL
        }
        
        # Add "NA" if year is required but variable not available
        for (y in years) {
          if (typeof(y) != "character") y = as.character(y)
          if (is.null(h.df[[y]])) h.df[[y]] = NA
        }
        
        break

      }

      # Some variables, such as "sex", are available all years with the same code: instead of a list, the online dictionary
      # simply says they are available "All Years".
      if (xmlValue(node) == "All Years") {
        
        h.df <- data.frame(t(rep(varindex,length(years))))
        colnames(h.df) <- years
        break

      }
    
    }

    # Add variable to main data table
    if (scrape_name) rownames(h.df) <- h.varname
    vardf <- rbind(vardf, h.df)

  }
  
  # Edit data table to make it ready to use as an input in build.panel.
  # Name columns with variable names, whether supplied or scraped. Name first column as "year".
  vardf <- as.data.table(t(vardf))
  if (using_names) colnames(vardf) <- varnames 
  vardf <- vardf[, "year":= as.numeric(years)]
  setcolorder(vardf, c(ncol(vardf), 1:(ncol(vardf)-1)))
  vardf <- vardf[order(vardf$year),]

  return(vardf)
  
}
  
