#' nanotatoR: Annotation package for Bionano Data
#'
#' Annotation of Bionano data using available databases
#'
#' @import knitr
#' @examples
#' path <- system.file("extdata", "Bionano_config/", package = "nanotatoR")
#' pattern <- "_hg19.txt"
#' mergedSmap <- makeInternalBNDatabase(path = path, 
#'    pattern = pattern, dbOutput = "dataframe")
#' mergedSmap[1,]
#' @docType package
#' @name nanotatoR
#' @import testthat
#' @import graphics
#' @import grDevices
#' @import knitr
#' @import utils



