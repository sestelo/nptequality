#' Testing the equality of k regression functions
#'
#' This function tests the hypothesis of equality of k regression functions
#' using a bootstrap-based approach.  The test is based on the comparison of
#' two estimators of the distribution of the errors in each population.
#'
#'
#' @param y Response variable.
#' @param x Independent variable.
#' @param f Categorical variable indicating the population to which
#' the observations belong.
#' @param m_bandwidths Bandwidth smoothing parameters for the estimation of the
#' regression function for each population. Default is "cv" (cross-validation).
#' @param sigma_bandwidths Bandwidth smoothing parameters for
#' the  estimation of the conditional variance function for each population.
#' Default is "cv".
#' @param hmin Minimum value for the bandwidths. Defaults to 0.
#' @param hmax Maximum value for the bandwidths. Defaults to 1.
#' @param nh Number of points in the grid for the bandwidths. Defaults to 30.
#' @param nboot Number of bootstrap repeats. Defaults to 100.
#' @param tstat Test statistic to be used. Options include "KS"
#' (Kolmogorov-Smirnov) or "CM" (Cramer-von Mises) type. Default is "KS".
#' @param cluster A logical value. If  \code{TRUE}, the
#'  testing procedure is  parallelized. Note that there are cases
#'  (e.g., a low number of bootstrap repetitions) that R will gain in
#'  performance through serial computation. R takes time to distribute tasks
#'  across the processors also it will need time for binding them all together
#'  later on. Therefore, if the time for distributing and gathering pieces
#'  together is greater than the time need for single-thread computing, it does
#'  not worth parallelize. Defaults to FALSE.
#' @param ncores An integer value specifying the number of cores to be used
#' in the parallelized procedure. If \code{NULL} (default), the number of cores
#' to be used is equal to the number of cores of the machine - 1.
#' @param seed Seed to be used in the procedure. If \code{NULL} (default), the seed is not set.
#'
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom doRNG registerDoRNG %dorng%
#' @importFrom parallel detectCores stopCluster
#' @importFrom stats ks.test bw.nrd0 rnorm sd
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @import npregfast
#'
#' @return A list containing test results,
#' including the test statistic and p-value.
#'
#' @export
#'
#' @author Marta Sestelo, Nora M. Villanueva and Juan Carlos Pardo-Fernández.
#' @examples
#'
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
#' out
#'
#'\dontrun{
#'library(npregfast)
#'out2 <- equalreg(x = barnacle$RC, y = barnacle$DW, f = barnacle$F,
#'tstat = "KS", hmin = 1, hmax = 5, nboot = 2, nh = 15)
#'}
#'

equalreg <- function(x, y, f, m_bandwidths = "cv", sigma_bandwidths = "cv",
                     hmin = 0, hmax = 1, nh = 30, nboot = 100, tstat = "KS",
                     cluster = FALSE, ncores = NULL, seed = NULL){




  # Checking cluster ------------------------------------------------

  if (isTRUE(cluster)) {
    if (is.null(ncores)) {
      num_cores <- parallel::detectCores() - 1
    }else{
      num_cores <- ncores
    }
    doParallel::registerDoParallel(cores = num_cores)
    if (!is.null(seed)) doRNG::registerDoRNG(seed)
    on.exit(doParallel::stopImplicitCluster())
  }


  # Organizing data -------------------------------------------------

  dat <- data.frame(x,y,f)
  x <- split(dat$x, dat$f)
  y <- split(dat$y, dat$f)
  f <- split(dat$f, dat$f)
  n <- list()
  n <- lapply(x, length)
  j <- length(n)





  # Obtain h ------------------------------------------------------


    h_m <- vector("list", j)
    h_sigma <- vector("list", j)

    if (m_bandwidths[1] == "cv") {
      h_m <- lapply(1:j, function(i) h.crossvalidation.ll(n[[i]], x[[i]], y[[i]], hmin, hmax, nh))
    } else {
      h_m <- as.list(m_bandwidths[1:j])
    }

    if (sigma_bandwidths[1] == "cv") {
      h_sigma <- lapply(1:j, function(i) h.crossvalidation.nw(n[[i]], x[[i]], y[[i]], hmin, hmax, nh))
    } else {
      h_sigma <- as.list(sigma_bandwidths[1:j])
    }



  # JC!! ver si entran como argumento como las hs anteriores
  h_f <- lapply(x, bw.nrd0)
  h.fmix <- bw.nrd0(unlist(x))





  # Estimation of m1, m2, etc, m0, sigma1, sigma2, etc -----------------

  sigma_hat <- list()
  for (k in 1:j) {
    nombre_vector <- paste0("sigma", k, "x", k, "_hat")
    sigma_hat[[nombre_vector]] <- sqrt(abs(
      nadarayawatson(n[[k]], x[[k]], y[[k]]^2, n[[k]], x[[k]], h_sigma[[k]]) -
        nadarayawatson(n[[k]], x[[k]], y[[k]], n[[k]], x[[k]], h_sigma[[k]])^2
    ))
  }


  m_hat <- list()
  f_hat <- list()
  for (k in 1:j) {
    for (i in 1:j) {
      nombre_vector_m <- paste0("m", k, "x", i, "_hat")
      nombre_vector_f <- paste0("f", k, "x", i, "_hat")
      m_hat[[nombre_vector_m]] <- locallinear(n[[k]], x[[k]], y[[k]], n[[i]], x[[i]], h_m[[k]])
      f_hat[[nombre_vector_f]]  <- kerneldensity(n[[k]],x[[k]],x[[i]],h_f[[k]])
    }
  }


  fmix_hat <- list()
  for (k in 1:j) {
    nombre_fmixx <- paste0("fmixx", k, "_hat")
    fmix_hat[[nombre_fmixx]] <- kerneldensity(sum(unlist(n)),unlist(x),x[[k]],h.fmix)
  }


  p <- list()
  for (i in 1:j) {
    p[[i]] <- n[[i]]/sum(unlist(n))
  }



  # Estimation of m0 ---------------------------------------------------
  m0_hat <- list()

  for (i in 1:j) {
    numerador <- 0
    for (k in 1:j) {
      f_kxi <- paste0("f", k, "x", i, "_hat")
      m_kxi <- paste0("m", k, "x", i, "_hat")
      numerador <- numerador + p[[k]] * f_hat[[f_kxi]] * m_hat[[m_kxi]]
    }

    fmixx_i <- paste0("fmixx", i, "_hat")
    nombre_m0xi <- paste0("m0x", i, "_hat")
    m0_hat[[nombre_m0xi]] <- numerador / fmix_hat[[fmixx_i]]
  }


  # Estimation of residuals -------------------------------------------------

  eps_hat <- list()
  eps0_hat <- list()

  for (k in 1:j) {
    m_kxk_hat <- paste0("m", k, "x", k, "_hat")
    m0xk_hat <- paste0("m0x", k, "_hat")
    sigma_kxk_hat <- paste0("sigma", k, "x", k, "_hat")

    nombre_eps <- paste0("eps", k, "_hat")
    eps_hat[[nombre_eps]] <- (y[[k]] - m_hat[[m_kxk_hat]]) / sigma_hat[[sigma_kxk_hat]]

    nombre_eps0 <- paste0("eps0", k, "_hat")
    eps0_hat[[nombre_eps0]] <- (y[[k]] - m0_hat[[m0xk_hat]]) / sigma_hat[[sigma_kxk_hat]]
  }



  # Test statistics -------------------------------------------------


  if(tstat == "KS"){
    KS_1 <- 0
    KS_2 <- 0
    for (k in 1:j) {
      KS_1 <- KS_1 + sqrt(n[[k]])*ks.test(eps0_hat[[k]], eps_hat[[k]])$statistic
    }
    KS_2 <- sqrt(sum(unlist(n))) * ks.test(unlist(eps0_hat), unlist(eps_hat))$statistic
    Tsample <- rbind(as.vector(KS_1), as.vector(KS_2))
  }


  if(tstat == "CM"){
    CM_1 <- 0
    CM_2 <- 0
    for (k in 1:j) {
      CM_1 <- CM_1 + cramer.vonmises(n[[k]], eps0_hat[[k]], n[[k]], eps_hat[[k]])
    }
    CM_2 <- cramer.vonmises(sum(unlist(n)), unlist(eps0_hat), sum(unlist(n)), unlist(eps_hat))
    Tsample <- rbind(as.vector(CM_1), as.vector(CM_2))
  }



  # Bootstrap-------------------


  a_smooth <- list(j)
  for (k in 1:j) {
    #a_smooth[[k]] <- 2 * n[[k]]**(-3/10)
    a_smooth[[k]] <- 0
  }

  eps_stand <- list()
  for (k in 1:j) {
    nombre_eps_hat <- paste0("eps", k, "_hat")
    eps_k_hat <- eps_hat[[nombre_eps_hat]]

    nombre_eps_stand <- paste0("eps", k, "_stand")
    eps_stand[[nombre_eps_stand]] <- (eps_k_hat - mean(eps_k_hat)) / sd(eps_k_hat)
  }

  eps_boot_matrix <- list()
  for (k in 1:j) {
    nombre_matrix <- paste0("eps", k, "_boot_matrix")
    set.seed(1)
    eps_boot_matrix[[nombre_matrix]] <- matrix(
      sqrt(1 - a_smooth[[k]]^2) * sample(eps_stand[[k]], size = nboot * n[[k]], replace = TRUE) +
        a_smooth[[k]] * rnorm(nboot * n[[k]]),
      nrow = nboot,
      ncol = n[[k]]
    )
  }

  iboot <- NULL
  if (isTRUE(cluster)) {
    Tboot <- foreach::foreach(iboot = 1:nboot, .combine = cbind) %dorng% {


      # Generation of yboot
      y_boot <- list()
      for (k in 1:j) {
        m0x_hat_k <- paste0("m0x", k, "_hat")
        sigma_hat_k <- paste0("sigma", k, "x", k, "_hat")
        eps_boot_matrix_k <- paste0("eps", k, "_boot_matrix")

        # Calcular y_k.boot para la iteración k
        nombre_y_boot <- paste0("y", k, "_boot")
        y_boot[[nombre_y_boot]] <- m0_hat[[m0x_hat_k]] + sigma_hat[[sigma_hat_k]] * eps_boot_matrix[[eps_boot_matrix_k]][iboot, ]
      }



      sigma_hat_boot <- list()
      for (k in 1:j) {
        nombre_vector <- paste0("sigma", k, "x", k, "_hat")
        sigma_hat_boot[[nombre_vector]] <- sqrt(abs(
          nadarayawatson(n[[k]], x[[k]], y_boot[[k]]^2, n[[k]], x[[k]], h_sigma[[k]]) -
            nadarayawatson(n[[k]], x[[k]], y_boot[[k]], n[[k]], x[[k]], h_sigma[[k]])^2
        ))
      }


      m_hat_boot <- list()
      for (k in 1:j) {
        for (i in 1:j) {
          nombre_vector_m <- paste0("m", k, "x", i, "_hat")
          m_hat_boot[[nombre_vector_m]] <- locallinear(n[[k]], x[[k]], y_boot[[k]], n[[i]], x[[i]], h_m[[k]])
        }
      }

      m0_hat_boot <- list()
      for (i in 1:j) {
        numerador <- 0
        for (k in 1:j) {
          f_kxi <- paste0("f", k, "x", i, "_hat")
          m_kxi <- paste0("m", k, "x", i, "_hat")

          numerador <- numerador + p[[k]] * f_hat[[f_kxi]] * m_hat_boot[[m_kxi]]
        }
        fmixx_i <- paste0("fmixx", i, "_hat")
        nombre_m0xi <- paste0("m0x", i, "_hat")
        m0_hat_boot[[nombre_m0xi]] <- numerador / fmix_hat[[fmixx_i]]
      }



      eps_hat_boot <- list()
      eps0_hat_boot <- list()

      for (k in 1:j) {
        m_kxk_hat <- paste0("m", k, "x", k, "_hat")
        m0xk_hat <- paste0("m0x", k, "_hat")
        sigma_kxk_hat <- paste0("sigma", k, "x", k, "_hat")

        nombre_eps <- paste0("eps", k, "_hat")
        eps_hat_boot[[nombre_eps]] <- (y_boot[[k]] - m_hat_boot[[m_kxk_hat]]) / sigma_hat_boot[[sigma_kxk_hat]]

        nombre_eps0 <- paste0("eps0", k, "_hat")
        eps0_hat_boot[[nombre_eps0]] <- (y_boot[[k]] - m0_hat_boot[[m0xk_hat]]) / sigma_hat_boot[[sigma_kxk_hat]]
      }



      if(tstat == "KS"){
        KS_1_boot <- 0
        KS_2_boot <- 0
        for (k in 1:j) {
          KS_1_boot <- KS_1_boot + sqrt(n[[k]])*ks.test(eps0_hat_boot[[k]], eps_hat_boot[[k]])$statistic
        }
        KS_2_boot <- sqrt(sum(unlist(n))) * ks.test(unlist(eps0_hat_boot), unlist(eps_hat_boot))$statistic
        return(c(as.vector(KS_1_boot), as.vector(KS_2_boot)))
      }

      if(tstat == "CM"){
        CM_1_boot <- 0
        CM_2_boot <- 0
        for (k in 1:j) {
          CM_1_boot <- CM_1_boot + cramer.vonmises(n[[k]], eps0_hat_boot[[k]], n[[k]], eps_hat_boot[[k]])
        }
        CM_2_boot <- cramer.vonmises(sum(unlist(n)), unlist(eps0_hat_boot), sum(unlist(n)), unlist(eps_hat_boot))
        return(c(as.vector(CM_1_boot), as.vector(CM_2_boot)))
      }

    } # close Tboot of %dorng%


  }else{

    pb <- txtProgressBar(min = 0, max = nboot, style = 3)
    Tboot <- foreach(iboot = 1:nboot, .combine = cbind) %do% {

      setTxtProgressBar(pb, iboot)
      # Generation of yboot
      y_boot <- list()
      for (k in 1:j) {
        m0x_hat_k <- paste0("m0x", k, "_hat")
        sigma_hat_k <- paste0("sigma", k, "x", k, "_hat")
        eps_boot_matrix_k <- paste0("eps", k, "_boot_matrix")

        # Calcular y_k.boot para la iteración k
        nombre_y_boot <- paste0("y", k, "_boot")
        y_boot[[nombre_y_boot]] <- m0_hat[[m0x_hat_k]] + sigma_hat[[sigma_hat_k]] * eps_boot_matrix[[eps_boot_matrix_k]][iboot, ]
      }



      sigma_hat_boot <- list()
      for (k in 1:j) {
        nombre_vector <- paste0("sigma", k, "x", k, "_hat")
        sigma_hat_boot[[nombre_vector]] <- sqrt(abs(
          nadarayawatson(n[[k]], x[[k]], y_boot[[k]]^2, n[[k]], x[[k]], h_sigma[[k]]) -
            nadarayawatson(n[[k]], x[[k]], y_boot[[k]], n[[k]], x[[k]], h_sigma[[k]])^2
        ))
      }


      m_hat_boot <- list()
      for (k in 1:j) {
        for (i in 1:j) {
          nombre_vector_m <- paste0("m", k, "x", i, "_hat")
          m_hat_boot[[nombre_vector_m]] <- locallinear(n[[k]], x[[k]], y_boot[[k]], n[[i]], x[[i]], h_m[[k]])
        }
      }

      m0_hat_boot <- list()
      for (i in 1:j) {
        numerador <- 0
        for (k in 1:j) {
          f_kxi <- paste0("f", k, "x", i, "_hat")
          m_kxi <- paste0("m", k, "x", i, "_hat")

          numerador <- numerador + p[[k]] * f_hat[[f_kxi]] * m_hat_boot[[m_kxi]]
        }
        fmixx_i <- paste0("fmixx", i, "_hat")
        nombre_m0xi <- paste0("m0x", i, "_hat")
        m0_hat_boot[[nombre_m0xi]] <- numerador / fmix_hat[[fmixx_i]]
      }



      eps_hat_boot <- list()
      eps0_hat_boot <- list()

      for (k in 1:j) {
        m_kxk_hat <- paste0("m", k, "x", k, "_hat")
        m0xk_hat <- paste0("m0x", k, "_hat")
        sigma_kxk_hat <- paste0("sigma", k, "x", k, "_hat")

        nombre_eps <- paste0("eps", k, "_hat")
        eps_hat_boot[[nombre_eps]] <- (y_boot[[k]] - m_hat_boot[[m_kxk_hat]]) / sigma_hat_boot[[sigma_kxk_hat]]

        nombre_eps0 <- paste0("eps0", k, "_hat")
        eps0_hat_boot[[nombre_eps0]] <- (y_boot[[k]] - m0_hat_boot[[m0xk_hat]]) / sigma_hat_boot[[sigma_kxk_hat]]
      }



      if(tstat == "KS"){
        KS_1_boot <- 0
        KS_2_boot <- 0
        for (k in 1:j) {
          KS_1_boot <- KS_1_boot + sqrt(n[[k]])*ks.test(eps0_hat_boot[[k]], eps_hat_boot[[k]])$statistic
        }
        KS_2_boot <- sqrt(sum(unlist(n))) * ks.test(unlist(eps0_hat_boot), unlist(eps_hat_boot))$statistic
        return(c(as.vector(KS_1_boot), as.vector(KS_2_boot)))
      }

      if(tstat == "CM"){
        CM_1_boot <- 0
        CM_2_boot <- 0
        for (k in 1:j) {
          CM_1_boot <- CM_1_boot + cramer.vonmises(n[[k]], eps0_hat_boot[[k]], n[[k]], eps_hat_boot[[k]])
        }
        CM_2_boot <- cramer.vonmises(sum(unlist(n)), unlist(eps0_hat_boot), sum(unlist(n)), unlist(eps_hat_boot))
        return(c(as.vector(CM_1_boot), as.vector(CM_2_boot)))
      }

    } # close Tboot of %do%

    close(pb)

  }


  if(tstat == "KS"){
    pvalue1 <- mean(Tboot[1, ] > Tsample[1, ])
    pvalue2 <- mean(Tboot[2, ] > Tsample[2, ])
  }

  if(tstat == "CM"){
    pvalue1 <- mean(Tboot[1, ] > Tsample[1, ])
    pvalue2 <- mean(Tboot[2, ] > Tsample[2, ])
  }

res <- list(pvalue = c(pvalue1, pvalue2), t = as.vector(Tsample), x = x, y = y,
            f = f, h_sigma = h_sigma, h_m = h_m, f_hat = f_hat,
            fmix_hat = fmix_hat, m0_hat = m0_hat, m_hat = m_hat,
            sigma_hat = sigma_hat, eps_hat = eps_hat, eps0_hat = eps0_hat,
            n = n, p = p, tstat = tstat, nboot = nboot, call = match.call())

  class(res) <- c("equalreg")

  return(res)
} # close function













