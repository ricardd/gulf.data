#' @title Snow Crab Observer Data Class:
#'
#' @description Class containers for snow crab observer data.
#'
#' @param x A \code{data.frame} object containing data fields to be mapped onto a corresponding snow crab observer data object.
#'
#' @param ... Other parameters (not used).
#'

#' @export
scobs <- function(x, ...) UseMethod("scobs")

#' @describeIn scobs Create an \code{scobs} object.
#' @export
scobs.default <- function(x, ...){
   gulf.metadata::key(x) <- c("year", "month", "day", "trip.number", "trap.number", "crab.number")

   # Define class:
   class(x) <- unique(c("scobs", class(x)))

   return(x)
}
