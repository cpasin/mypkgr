#' Plot of the mvnpdf function
#'
#' @param x an object of class \code{mvnpdf} resulting from a call of
#' \code{mnvpdf()} function.
#' @param ... graphical parameters passed to \code{plot()} function.
#'
#' @return Nothing is returned, only a plot is given.
#' @export
#' @importFrom graphics plot
#'
#' @examples
#' pdfvalues <- mvnpdf(x=matrix(seq(-3, 3, by = 0.1), nrow = 1), Log=FALSE)
#' plot(pdfvalues)
plot.mvnpdf <- function(x, ...) {
  graphics::plot(x$x, x$y, type = "l", ...)
}
