

#' psidR
#'
#' psidR is a package that helps the task of building longitudinal datasets from the Panel Study of Income Dynamics (PSID). 
#' The user must supply the PSID variable names that correspond to the variables of interest in each desired wave. Data can be supplied via Stata, or directly downloaded from PSID servers without any need for STATA.
#' data.frame.
#' @import data.table
#' @import SAScii 
#' @import RCurl
#' @import foreign
#' @import futile.logger
#' @importFrom openxlsx read.xlsx
#' @importFrom stats rlnorm rnorm runif
#' @importFrom utils object.size tail unzip head
#' @docType package
#' @name psidR
#' 
NULL
