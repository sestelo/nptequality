# Kernel Epanechnikov (density)
kernel = function(u) {
  (0.75 * (1 - u ^ 2)) * (u < 1) * (u > -1)
}

# Kernel density estimator. This function calculates the kernel density
# estimator on the collection of "npoints" points "points" based on the
# "ndata" observations xdata and the bandwidth "h".
kerneldensity = function(ndata, data, points, h) {
  rowSums(kernel(outer(points, data, "-") / h)) / (ndata * h)
}

# Nadaraya-Watson estimator. This function calculates the N-W estimator on the
# collection of "npoints" points "points" based on the "ndata"
# observations (xdata,ydata) and the bandwidth "h".
nadarayawatson = function(ndata, xdata, ydata, npoints, points, h)
{
  as.vector({
    matk = kernel((points %*% t(rep(1, ndata)) - t(xdata %*% t(rep(1, npoints)))) / h)

   aux <- (matk %*% ydata) / (matk %*% rep(1, ndata))
   ifelse(is.nan(aux), NA, aux)
  })

}

# Cross-validation bandwidth selector (Nadaraya-Watson)
h.crossvalidation.nw = function(n, x, y, hmin, hmax, hngrid) {
  crossvalue = rep(0, hngrid)
  h = seq(hmin, hmax, len = hngrid)
  for (j in 1:hngrid) {
    for (i in 1:n) {
      nw <- nadarayawatson(n - 1, x[-i], y[-i], 1, x[i], h[j])
      #print(nw)
      if (is.na(nw)) {
        crossvalue[j] <- 99999999
        break
      }else{
      crossvalue[j] = crossvalue[j] + (y[i] - nw) ^ 2
      }
      #print(i)
      #print(crossvalue[j])
    }
  }

  if (all(crossvalue == 99999999)) {
    stop("\nChoose a larger bandwidth range.")
  }

  crossvalue = crossvalue / n
  #plot(h,crossvalue)
  h.crossvalidation.nw = h[which.min(crossvalue)]
}

# Local-linear regression estimator
locallinear = function(ndata, xdata, ydata, npoints, points, h) {
  # See Wand and Jones (1995). "Kernel Smoothing" (page 119).
  ii <- 1
  while (sum(ii) > 0) {
    mat.x = outer(points, xdata, "-")
    mat.k = kernel(mat.x / h)
    mat.y = matrix(rep(ydata, npoints), nrow = npoints, byrow = T)
    s0 = matrix(rep(rowSums(mat.k) / ndata, ndata), ncol = ndata)
    s1 = matrix(rep(rowSums(mat.x * mat.k) / ndata, ndata), ncol = ndata)
    s2 = matrix(rep(rowSums((mat.x ^ 2) * mat.k) / ndata, ndata), ncol =
                  ndata)
    res <- (rowSums((s2 - s1 * mat.x) * mat.k * mat.y) / ndata) / (s2[, 1] * s0[, 1] - s1[, 1] ^ 2)
    ii <- is.nan(res) | is.infinite(res)
    h <- h + 0.05
  }

  # ii <- is.nan(res)
  #  res[ii] <- points[ii]
  #ii <- is.infinite(res)
  #res[ii] <- 1-48*points[ii]+218*points[ii]**2-315*points[ii]**3+145*points[ii]**4
  return(res)
}

# Cross-validation bandwidth selector (local linear)
h.crossvalidation.ll = function(n, x, y, hmin, hmax, hngrid) {
  crossvalue = rep(0, hngrid)
  h = seq(hmin, hmax, len = hngrid)
  for (j in 1:hngrid) {
    for (i in 1:n) {
      crossvalue[j] = crossvalue[j] + (y[i] - locallinear(n - 1, x[-i], y[-i], 1, x[i], h[j])) ^ 2
    }
  }
  crossvalue = crossvalue / n
  #plot(h,crossvalue)
  h.crossvalidation.ll = h[which.min(crossvalue)]
  return(h.crossvalidation.ll)
}


################################################
# - Pardo-Fernández, J. C.; Van Keilegom, I.; González-Manteiga, W. (2007).
# Testing for the equality of k regression curves. Statistica Sinica, 17, 1115-1137.

# Cramér-von Mises statistic
emp.distr = function(ndata, data, npoints, points) {
  sapply(data, function(data, points) {
    data <= points
  }, points = points) %*% rep(1, ndata) / ndata
}
cramer.vonmises = function(nx, x, ny, y) {
  sum((emp.distr(nx, x, nx + ny, c(x, y)) - emp.distr(ny, y, nx + ny, c(x, y))) ^ 2)
}

