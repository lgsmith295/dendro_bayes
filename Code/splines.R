# from Sam Clifford: http://bragqut.github.io/2016/05/24/samclifford-splines/

N <- 100
x <- sort(runif(n=N, min = 0, max = 1))
y <- sin(2*pi*x) - 5*x^2 + 3*x + rnorm(n=N, 10, sd=0.25)

dat <- data.frame(y=y, x=x)

d <- max(4, floor(N/35))
K <- floor(N/d - 1)

i = 3
y_i <- y_orig[i, f[i]:l[i]]
x <- seq(1:length(y_i))

# Set up basis function
B <- bs(x, knots = c(min(x, na.rm = T), floor(length(x) * 0.67), max(x, na.rm = T)))
K = 6 # one internal knot on cubic spline - I think this is supposed to be deg in bs


# set up penalty matrix
makeQ = function(degree, K, epsilon=1e-3){
  x <- diag(K)
  E <- diff(x, differences=degree)
  return( t(E) %*% E + x*epsilon)
} 

Q <- makeQ(2, K=K)


sink("Code/JAGS/p_spline.txt")
cat("model{
    for (i in 1:n){
    y[i] ~ dnorm(mu[i], tau.y)
    mu[i] = X[i, ] %*% beta[1:K]
    }
    
    beta[1:K] ~ dmnorm(beta.0[1:K,1] + beta.00, lambda*Q[1:K,1:K]) # 
    
    gamma <- beta - beta.00
    
    for(k in 1:K) {
    beta.0[k, 1] ~ dnorm(0, 0.0001)
    }
    
    tau.y ~ dgamma(0.001, 0.001)
    lambda ~ dgamma(0.001, 0.001)
    
    beta.00 ~ dnorm(0, 1e-6)
    
    for (j in 1:m){
    y.rep[j] ~ dnorm(mu.rep[j], tau.y)
    mu.rep[j] = X.pred[j, ] %*% gamma[1:K] + beta.00
    }
    
    }
    ", fill=TRUE)
sink()

M_spline = jags.model('Code/JAGS/p_spline.txt', data = list(y = y_i, X.pred = B, m = length(x), n = length(x), X = B, K = K, Q = Q), n.chains = 3, n.adapt = 500) # , beta.0 = as.matrix(rep(0, K))
outM_spline = coda.samples(M_spline,c("y.rep", "mu.rep", "lambda", "beta.0"), 1000) 

# plot(outM_spline)
par(mfrow = c(1,1))

gamidx = which(substr(varnames(outM_spline),1,6)=="beta.0") # finds the indices of the gamma variables
xyplot(outM_spline[,gamidx]) # traceplot of the gamma variables
summary(outM_spline[,gamidx]) # summary of the gamma variables

xidx = which(substr(varnames(outM_spline),1,3)=="mu.")
postxm = colMeans(as.matrix(outM_spline[,xidx])) 
plot(postxm,type="l")
plot(year,postxm*x_sd + x_mean,type="l")

plot(x, y_i)
lines(x, postxm)












nyrs2 <- floor(nY2 * 0.67)

library(splines)
i = 5
foo <- ffcsaps(y = y_orig[i, f[i]:l[i]])
plot(foo, type = "l", ylim = c(100, 2000))

y_i <- y_orig[i, f[i]:l[i]]
t <- seq(1:length(y_i))
B <- bs(t, knots = floor(length(y_i) * 0.67)) # quantile(y_i, probs = c(0.1, 0.5, 0.9))) # floor(length(y_i) * 0.67))
matplot(t, B,type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)

spl1 <- lm(y_i ~ B)
summary(spl1)
# plot(seq(1:length(y_i)), predict(spl1, newdata = seq(1:length(y_i))))
plot(predict(spl1), type = "l", ylim = c(100, 2000))


B_ns <- ns(t, knots = floor(length(y_i) * 0.67)) # quantile(y_i, probs = c(0.1, 0.5, 0.9))) # floor(length(y_i) * 0.67))
matplot(t, B_ns ,type="l",xlab="Time",ylab="Basis function, Bj(X)",cex.lab=1.5,cex.axis=1.5,lwd=2)

spl1_ns <- lm(y_i ~ B_ns)
summary(spl1_ns)
# plot(seq(1:length(y_i)), predict(spl1, newdata = seq(1:length(y_i))))
plot(predict(spl1_ns), type = "l", ylim = c(100, 2000))



bar <- ffcsaps(y = log_y[1, 1:485])
plot(bar, type = "l")


function (y, x = seq_along(y), nyrs = length(y)/2, f = 0.5) 
{
  ffppual <- function(breaks, c1, c2, c3, c4, x, left) {
    if (left) {
      ix <- order(x)
      x2 <- x[ix]
    }
    else {
      x2 <- x
    }
    n.breaks <- length(breaks)
    if (left) {
      index <- pmax(ffsorted(breaks[-n.breaks], x2), 1)
    }
    else {
      index <- ffsorted2(breaks[-1], x2)
    }
    x2 <- x2 - breaks[index]
    v <- x2 * (x2 * (x2 * c1[index] + c2[index]) + c3[index]) + 
      c4[index]
    if (left) 
      v[ix] <- v
    v
  }
  ffsorted <- function(meshsites, sites) {
    index <- order(c(meshsites, sites))
    which(index > length(meshsites)) - seq_along(sites)
  }
  ffsorted2 <- function(meshsites, sites) {
    index <- order(c(sites, meshsites))
    which(index <= length(sites)) - seq(from = 0, to = length(sites) - 
                                          1)
  }
  spdiags <- function(B, d, n) {
    n.d <- length(d)
    A <- matrix(0, n.d * n, 3)
    count <- 0
    for (k in seq_len(n.d)) {
      this.diag <- d[k]
      i <- inc(max(1, 1 - this.diag), min(n, n - this.diag))
      n.i <- length(i)
      if (n.i > 0) {
        j <- i + this.diag
        row.idx <- seq(from = count + 1, by = 1, length.out = n.i)
        A[row.idx, 1] <- i
        A[row.idx, 2] <- j
        A[row.idx, 3] <- B[j, k]
        count <- count + n.i
      }
    }
    A <- A[A[, 3] != 0, , drop = FALSE]
    A[order(A[, 2], A[, 1]), , drop = FALSE]
  }
  y2 <- as.numeric(y)
  if (!is.numeric(y2)) 
    stop("'y' must be coercible to a numeric vector")
  x2 <- as.numeric(x)
  if (!is.numeric(x2)) 
    stop("'x' must be coercible to a numeric vector")
  n <- length(x2)
  if (n < 3) 
    stop("there must be at least 3 data points")
  if (!is.numeric(f) || length(f) != 1 || f < 0 || f > 1) 
    stop("'f' must be a number between 0 and 1")
  if (!is.numeric(nyrs) || length(nyrs) != 1 || nyrs <= 1) 
    stop("'nyrs' must be a number greater than 1")
  ix <- order(x2)
  zz1 <- n - 1
  xi <- x2[ix]
  zz2 <- n - 2
  diff.xi <- diff(xi)
  if (any(diff.xi == 0)) 
    stop("the data abscissae must be distinct")
  yn <- length(y2)
  if (n != yn) 
    stop("abscissa and ordinate vector must be of the same length")
  arg2 <- -1:1
  odx <- 1/diff.xi
  R <- spdiags(cbind(c(diff.xi[-c(1, zz1)], 0), 2 * (diff.xi[-1] + 
                                                       diff.xi[-zz1]), c(0, diff.xi[-c(1, 2)])), arg2, zz2)
  R2 <- spdiags(cbind(c(odx[-zz1], 0, 0), c(0, -(odx[-1] + 
                                                   odx[-zz1]), 0), c(0, 0, odx[-1])), arg2, n)
  R2[, 1] <- R2[, 1] - 1
  forR <- Matrix(0, zz2, zz2, sparse = TRUE)
  forR2 <- Matrix(0, zz2, n, sparse = TRUE)
  forR[R[, 1:2, drop = FALSE]] <- R[, 3]
  forR2[R2[, 1:2, drop = FALSE]] <- R2[, 3]
  p.inv <- (1 - f) * (cos(2 * pi/nyrs) + 2)/(12 * (cos(2 * 
                                                         pi/nyrs) - 1)^2)/f + 1
  yi <- y2[ix]
  p <- 1/p.inv
  mplier <- 6 - 6/p.inv
  u <- as.numeric(solve(mplier * tcrossprod(forR2) + forR * 
                          p, diff(diff(yi)/diff.xi)))
  yi <- yi - mplier * diff(c(0, diff(c(0, u, 0))/diff.xi, 0))
  test0 <- xi[-c(1, n)]
  c3 <- c(0, u/p.inv, 0)
  x3 <- c(test0, seq(from = xi[1], to = xi[n], length = 101))
  cc.1 <- diff(c3)/diff.xi
  cc.2 <- 3 * c3[-n]
  cc.3 <- diff(yi)/diff.xi - diff.xi * (2 * c3[-n] + c3[-1])
  cc.4 <- yi[-n]
  to.sort <- c(test0, x3)
  ix.final <- order(to.sort)
  sorted.final <- to.sort[ix.final]
  tmp <- unique(data.frame(sorted.final, c(ffppual(xi, cc.1, 
                                                   cc.2, cc.3, cc.4, test0, FALSE), ffppual(xi, cc.1, cc.2, 
                                                                                            cc.3, cc.4, x3, TRUE))[ix.final]))
  tmp2 <- tmp
  tmp2[[1]] <- round(tmp2[[1]], 5)
  res <- tmp2[[2]][tmp2[[1]] %in% x2]
  if (length(res) != n) 
    res <- approx(x = tmp[[1]], y = tmp[[2]], xout = x2)$y
  res
}