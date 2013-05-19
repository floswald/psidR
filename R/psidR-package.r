

#' psidR
#'
#' psidR is a package that helps the task of building longitudinal datasets from the Panel Study of Income Dynamics (PSID). 
#' The user must supply the PSID variable names that correspond to the variables of interest in each desired wave. The data may be in .dta, .csv format on disk. Creation
#' of .dta or .csv datasets requires access to Stata or SAS software. There is an option to bypass this requirement by directly downloading the data from the server into a
#' data.frame.
#' @import data.table
#' @import SAScii 
#' @import RCurl
#' @suggest survey
#' @doctype package
#' @name psidR
#' 
NULL
