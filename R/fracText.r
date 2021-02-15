#' Place text on a plot, with locations given as fraction from a corner rather than with absolute coordinates
#'
#' @param fx A fraction (or a vector of fractions) of the plot interval across from the left to be. These should be numeric and range from 0.0 to 1.0. 
#' @param fy A fraction (or a vector of fractions) of the plot interval across from the bottom to be. These should be numeric and range from 0.0 to 1.0. It should have the same length as \code{fx}.
#' @param labels A vector of strings. If \code{fx} and \code{fy} have length greater than 1, then \code{labels} should have the same length.
#' @param ... Further arguments passed to \code{text()}.
#' @return Nothing is returned
fracText <- function(fx,fy,labels,...){
  text(sum(par()$usr[1:2]*c(1-fx,fx)),
       sum(par()$usr[3:4]*c(1-fy,fy)),
       labels, ...)
}
