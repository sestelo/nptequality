#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
#' Visualization of \code{equalreg} objects with ggplot2 graphics
#'
#' @description Useful for drawing the estimated functions (regression and
#' conditional variance functions) for each group or population. Additionally, the empirical
#' cumulative distribution functions for the residuals of each group are displayed.
#'
#' @param object Object of \code{equalreg} class.
#' @param interactive Logical flag indicating if an interactive plot with
#' plotly is produced.
#' @param type Character string indicating the type of plot to be produced. Options
#' include "mean" (mean regression function), "var" (conditional variances function),
#' "eps" (empirical distributions functions) or "all" (it includes all) type. Default is "all".
#' @param \ldots Other options.
#'
#'
#' @return A ggplot object, so you can use common features from
#' ggplot2 package to manipulate the plot.
#'
#' @author Marta Sestelo, Nora M. Villanueva and Juan Carlos Pardo-Fernández.
#' @examples
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
#' nh = 46, nboot = 100, parallel = FALSE, seed = 300716)
#' autoplot(out)
#' autoplot(out, interactive = TRUE)
#'
#'
#' @import ggplot2
#' @importFrom scales hue_pal
#' @importFrom patchwork wrap_plots
#' @export

autoplot.equalreg <- function(object = object, interactive = FALSE, type = "all",...){


  if(!inherits(object, "equalreg")){
    stop("The object must be of class equalreg")
  }

  # if(k < 3){
  #   colgr <- brewer.pal(n = 3, name = "Dark2")
  # }else if(k<9){
  #   colgr <- brewer.pal(n = k, name = "Dark2")
  # }else{
  #   colgr <- colorRampPalette(brewer.pal(n = 8, name = "Dark2"))(k)
  # }
  #
    j <- length(object$x)

    x <- NULL; y <- NULL; f <- NULL; eps <- NULL; eps0 <- NULL # for solving the problem in check()
    # autoplot.equalreg: no visible binding for global variable ‘x’, etc.

    m_hat <- NULL
    m0_hat <- NULL
    sigma_hat <- NULL
    for (i in 1:j) {
      nombre_vector_m <- paste0("m", i, "x", i, "_hat")
      m_hat <- c(m_hat, object$m_hat[[nombre_vector_m]])

      nombre_vector_m <- paste0("m0x", i, "_hat")
      m0_hat <- c(m0_hat, object$m0_hat[[nombre_vector_m]])

      nombre_vector_m <- paste0("sigma", i, "x", i, "_hat")
      sigma_hat <- c(sigma_hat, object$sigma_hat[[nombre_vector_m]])
    }





      data <- data.frame(x = unlist(object$x), y = unlist(object$y),
                       f = as.factor(unlist(object$f)),
                       m_hat = m_hat, m0_hat = m0_hat, sigma_hat = sigma_hat)

      pal <- scales::hue_pal()(j)

    plots <- list()
    plots[[1]] <- ggplot2::ggplot(data, aes(x = x, y = y, colour = f)) + geom_point() +
      geom_line(aes(x = x, y = m_hat, colour = f)) +
      geom_line(aes(x = x, y = m0_hat, colour = f), linetype = "dashed") +
      ggtitle("Estimated regression functions") +
      xlab("covariate") +
      ylab("response") +
      scale_colour_manual(values = pal)


    plots[[2]] <- ggplot2::ggplot(data, aes(x = x, y = sigma_hat, colour = f)) + geom_line() +
      ggtitle("Conditional variance functions") +
      xlab("covariate") +
      ylab("sigma") +
      scale_y_continuous(limits = c(min(data$sigma_hat), max(data$sigma_hat))) +
      scale_colour_manual(values = pal)



    for(i in 1:j){
      dataeps <- data.frame(eps = object$eps_hat[[i]], eps0 = object$eps0_hat[[i]])
      plots[[i+2]] <- ggplot2::ggplot(dataeps, aes(eps)) +
        stat_ecdf(geom = "step", colour = pal[i]) +
        stat_ecdf(geom = "step", aes(eps0)) +
        ggtitle(paste("F_eps", i))
    }



     if(type == "mean"){
      #combined_plot <- plots[[1]]
       combined_plot <- ggplot2::ggplot(data, aes(x = x, y = y, colour = f)) +
         geom_line(aes(x = x, y = m_hat, colour = f)) +
         ggtitle("Estimated regression functions") +
         xlab("covariate") +
         ylab("response") +
         scale_colour_manual(values = pal)
    } else if(type == "var"){
      combined_plot <- plots[[2]]
    } else if(type == "eps"){
      combined_plot <- patchwork::wrap_plots(plots[-(1:2)], ncol = 2, nrow = ceiling((j) / 2  ) )
    } else if(type == "all"){
      combined_plot <- patchwork::wrap_plots(plots, ncol = 2, nrow = ceiling((j + 2) / 2  ) )
    } else {
      stop("Invalid type argument. Choose from 'mean', 'var', 'eps', or 'all'.")
    }

    #combined_plot <- patchwork::wrap_plots(plots, ncol = 2, nrow = ceiling((j + 2) / 2  ) )



    # if(interactive == TRUE){
    #
    #   interactive_plots <- lapply(plots, function(x){x + ggtitle("")})
    #   interactive_plots <- lapply(interactive_plots, plotly::ggplotly)
    #   combined_plot <- plotly::subplot(interactive_plots, nrows = ceiling((j + 2) / 2  ), titleX = TRUE,
    #                                    titleY = TRUE,  margin = 0.05)
    #
    # }


    if (interactive == TRUE) {

      if (type == "mean") {
        selected_plots <- list(
          plota <- ggplot2::ggplot(data, aes(x = x, y = y, colour = f)) +
            geom_line(aes(x = x, y = m_hat, colour = f)) +
            ggtitle("Estimated regression functions") +
            xlab("covariate") +
            ylab("response") +
            scale_colour_manual(values = pal)
        )

      } else if (type == "var") {
        selected_plots <- list(plots[[2]])
      } else if (type == "eps") {
        selected_plots <- plots[-(1:2)]
      } else if (type == "all") {
        selected_plots <- plots
      }


      interactive_plots <- lapply(selected_plots, function(x) x + ggtitle(""))
      interactive_plots <- lapply(interactive_plots, plotly::ggplotly)

      titles <- sapply(selected_plots, function(p) p$labels$title) # no funciona

      combined_plot <- plotly::subplot(
        interactive_plots,
        nrows = ceiling(length(interactive_plots) / 2),
        titleX = TRUE,
        titleY = TRUE,
        margin = 0.05
      )

}


    combined_plot

    #return(invisible(object))
  }








