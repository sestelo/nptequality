#' Summarizing fits of \code{equalreg} class produced by \code{equalreg}
#'
#' @description Takes a  equalreg object
#' and produces various useful summaries from it.
#'
#' @param object a equalreg object as producted by \code{equalreg}
#' @param digits number of digits to print.
#' @param \ldots additional arguments.
#'
#' @details \code{print.equalreg} tries to be smart about \code{summary.equalreg}.
#'
#' @return \code{summary.equalreg} computes and returns a list of summary
#' information for an \code{equalreg} object.
#'
#' \item{nboot}{Number of bootstrap repeats.}
#' \item{n}{Total sample size.}
#' \item{nj}{Sample size of each group.}
#' \item{test.statistic}{Test statistics used.}
#' \item{p-values}{p-values of the test.}
#' \item{components}{A collection of components of the equalreg object.}
#'
#'
#'@author Marta Sestelo, Nora M. Villanueva and Juan Carlos Pardo-Fernández.
#'
#'@examples
#'
#' set.seed(300716)
#' n1 <- 100
#' n2 <- 200
#' x1 <- runif(n1)
#' y1 <- x1 + 0.25 * rnorm(n1)
#' x2 <- runif(n2)
#' y2 <- x2 + 0.25 * rnorm(n2)
#'
#' x <- c(x1, x2)
#' y <- c(y1, y2)
#' f <- c(rep(1, n1), rep(2, n2))
#' out <- equalreg(x = x, y = y, f = f, tstat = "KS", hmin = 0.05, hmax = 0.5,
#' nh = 46, nboot = 100, cluster = FALSE, seed = 300716)
#' summary(out)
#'
#'
#'
#' @export







summary.equalreg <- function(object, digits = 4, ...){

  if (object$tstat == "KS"){
    test <- "Kolmogorov-Smirnov"
  }
  if (object$tstat == "CM"){
    test <- "Cramer-von Mises"
  }

  cat("\n")
  cat("Testing the equality of", length(unique(object$f)), "nonparametric regression functions.", "\n")
  cat("Bootstrap algorithm is applied to obtain the null distribution.", "\n")
  cat(test,"type statistic used.", "\n")

  cat("\n")
  cat("Number of bootstrap samples: ", object$nboot, "\n")
  cat("Number of observations: ",sum(unlist(object$n)), "\n")
  cat("Number of observations for each population: ", sapply(object$f, length), "\n")
  cat("\n")
  cat("Null hypothesis: there is no difference between the ", length(unique(object$f)), " curves", sep="", ".\n")

  cat("T1 = ", formatC(object$t[1], digits = digits), "   ", "p-value = ", formatC(object$pvalue[1], digits = digits), "\n")
  cat("T2 = ", formatC(object$t[2], digits = digits), "   ", "p-value = ", formatC(object$pvalue[2], digits = digits), "\n")



  cat("\n")
  cat("\nCall: ","\n")
  print(object$call)
  cat("\n")
  cat("\nAvailable components:", sep = "\n")
  print(names(object))

  invisible(object)
}






