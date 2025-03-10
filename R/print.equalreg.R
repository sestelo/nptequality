#' @export
print.equalreg <- function(x, ...) {

  if (inherits(x, "equalreg")) {
    cat("\n")
    cat("Testing the equality of", length(unique(x$f)), "nonparametric regression functions.", "\n")

    cat("\n")
    cat("Null hypothesis: there is no difference between the ", length(unique(x$f)), " curves", sep="", ".\n")

    cat("T1 = ", formatC(x$t[1], digits = 4), "   ", "p-value = ", formatC(x$pvalue[1], digits = 4), "\n")
    cat("T2 = ", formatC(x$t[2], digits = 4), "   ", "p-value = ", formatC(x$pvalue[2], digits = 4), "\n")


    invisible(x)
  }
}
