
#' This function helps the user save a lot of time in building a data table containing the original PSID variable codes (e.g. ER30046) for 
#' each desired year and checks on the psidonline web site in which years a variable is available or not. Information is obtained on 
#' https://simba.isr.umich.edu/VS/s.aspx.
#'
#' The function build.panel requires the user to input a table containing the original PSID variable codes. However, this requires looking up
#' the unique variable identifier for each different wave: for instance, ER30046 stands for age of individual in 1970, while ER30070 for age of
#' individual in 1971. If the user wants to build a long panel, looking up EACH variable code for EACH wave manually takes a long time. This 
#' function simplifies the user's life by requiring to provide a list with the original variable codes from one wave only. It then 
#' automatically looks up the corresponding variable codes for the same variable in other waves, if available. For instance, if the user  
#' includes the code ER30046 (age in 1970) in the list and declares that the panel should span from 1970 to 1971, then the function will 
#' return a table containing the original PSID codes for the variable "age of individual" for all desired years, in this case ER30046 (1970)
#' and ER30070 (1971).
#'
#' Arguments:
#' @param vars character vector contains the original PSID variable codes of type ER30046.
#' @param years numeric vector contains the waves for which the user wants to retrieve the original variable codes.
#' @param scrape_name logical indicating whether to retrieve the original variable names (T) or not (F) (e.g. age_of_individual for ER30046).
#' @param varnames character vector contains customized variable names; applied if scrape_name set to F and must be of same length as vars.
#'
#' Note that, in the list of variables provided, one could include a variable code from year y and a variable code from year y'. This is 
#' allowed and irrelevant to the functioning of the function. Only remark: if the codes refer to the same variable, then the same variable will
#' show up twice in the output table.
#'
#' - If scrape_name equals TRUE, then the function recovers the variable names from the online dictionary and edits them to avoid the presence
#' of special characters such as "/", "#" or "-". The dictionary names suggested by psidonline might be long and cryptic and are NOT always
#' consistent over years (even if the underlying variable is consistent).
#' - If the vector varnames is provided and scrape_name equals FALSE, then the function names the variables after the user-provided labels;
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
#' > head(ind)
#'      year age_of_individual grade_finished individual_weight
#' 1968 1968           ER30004        ER30010           ER30019
#' 1969 1969           ER30023             NA           ER30042
#' 1970 1970           ER30046        ER30052           ER30066
#' 1971 1971           ER30070        ER30076           ER30090
#' 1972 1972           ER30094        ER30100           ER30116
#' 1973 1973           ER30120        ER30126           ER30137
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

get.variable.names <- function(vars, years, scrape_name = F, varnames = NULL){
  
  # Define certicificate file
  cafile <- system.file("CurlSSL", "cacert.pem", package = "RCurl")
  
  # Check if the function is required to scrape variable names from the dictionary or not
  using_names <- !(is.null(varnames))
  
  if (using_names & scrape_name) {
    scrape_name <- F
    warning("The option scrape_name has been set to FALSE. If you require the code to scrape variable names, let the varnames option set to default value NULL.\n")
  }

  # Create data table where to store the original variable names wave by wave  
  vardf <- as.data.table(matrix(0, ncol = length(years), nrow = 0))
  colnames(vardf) <- as.character(years)
  
  # Go through each variable
  for (varindex in vars) {

    # Get html code from psidonline
    h.source <- GET(paste0("https://simba.isr.umich.edu/cb.aspx?vList=", varindex),
      config(cainfo = cafile))
    h.page <- htmlParse(h.source)
    
    # Recover and edit the explicit variable name (replace special characters as they might be complicated to handle with)
    if (scrape_name) {
      h.varname <- xmlValue(getNodeSet(h.source, "//td")[[3]])
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
    for (node in rev(getNodeSet(h.page, "//td"))) {

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

  
