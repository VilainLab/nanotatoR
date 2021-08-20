#' nanotatoR: Annotation package for Bionano Data
#'
#' Annotation of Bionano data using available databases
#'
#' @import knitr
#' @import testthat
#' @import utils
#' @examples
#' path <- system.file("extdata", "Bionano_config/", package = "nanotatoR")
#' pattern <- "_hg19.txt"
#' mergedSmap <- makeInternalBNDatabase(path = path, 
#'    pattern = pattern, dbOutput = "dataframe")
#' mergedSmap[1,]
#' @docType package
#' @name nanotatoR
NULL
# > [1] '_PACKAGE' 



