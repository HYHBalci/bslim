args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else getwd(), winslash = "/")
source(file.path(root, "research/bayes-factor-pilot-2026-10-16/bayes_factor_sidecar.R"))
set.seed(20261016)
n <- 14; W <- rnorm(n); D <- matrix(rnorm(n*3), n, 3); r <- rnorm(n)
sigma2 <- 1.3; va <- 2.1; s <- .65
log_density <- function(r, V) {
  L <- chol(V)
  -.5*length(r)*log(2*pi) - sum(log(diag(L))) - .5*sum(forwardsolve(t(L),r)^2)
}
close <- function(a,b,tol=1e-9) stopifnot(max(abs(a-b)) < tol)
V0 <- sigma2*diag(n) + va*tcrossprod(W)
V1 <- V0 + s^2*tcrossprod(D)
b <- conditional_log_bf(r,D,W,s,sigma2,va)
close(b, log_density(r,V1)-log_density(r,V0))
perm <- sample.int(n)
close(b,conditional_log_bf(r[perm],D[perm,],W[perm],s,sigma2,va))
close(b,conditional_log_bf(r,D[,3:1],W,s,sigma2,va))
close(b,conditional_log_bf(r,-D,W,s,sigma2,va))
close(conditional_log_bf(r,D,W,0,sigma2,va),0)
close(conditional_log_bf(r,D*0,W,s,sigma2,va),0)
a <- drop(crossprod(D[,1],solve(V0,D[,1])))
z <- drop(crossprod(D[,1],solve(V0,r)))
close(conditional_log_bf(r,D[,1,drop=FALSE],W,s,sigma2,va),
      -.5*log1p(s^2*a)+.5*s^2*z^2/(1+s^2*a))
p <- subset_prior(3); close(sum(exp(p$log_weight)),1)
bf <- evaluate_bf_priors(r,D,W,s,sigma2,va)
direct <- sum(vapply(seq_along(p$subsets), function(j) {
  Ds <- D[,p$subsets[[j]],drop=FALSE]
  exp(p$log_weight[j] + log_density(r,V0+s^2*tcrossprod(Ds))-log_density(r,V0))
},numeric(1)))
close(exp(bf$sparse),direct)
close(score_experiment_log_bf(6,3,0),0)
close(score_experiment_log_bf(6,3,2), -.5*3*log(3)+2)
# Independent quadrature check of the null-posterior expectation identity.
v <- rnorm(n); vf <- .4
post_var <- 1/(1/vf + drop(crossprod(v,solve(V0,v))))
post_mean <- post_var*drop(crossprod(v,solve(V0,r)))
exact <- exp(log_density(r,V0+vf*tcrossprod(v)+s^2*tcrossprod(D))-
             log_density(r,V0+vf*tcrossprod(v)))
integrand <- function(t) vapply(t,function(tt)
  exp(conditional_log_bf(r-v*tt,D,W,s,sigma2,va))*
    dnorm(tt,post_mean,sqrt(post_var)), numeric(1))
quad <- integrate(integrand,post_mean-12*sqrt(post_var),
                  post_mean+12*sqrt(post_var),rel.tol=1e-10)$value
close(exact,quad,1e-8)
# Invalid finite inputs must fail rather than silently manufacture evidence.
stopifnot(inherits(try(conditional_log_bf(r,D,W,-1),silent=TRUE),"try-error"))
cat("PASS: covariance-density equality; observation/column/sign invariance; zero-slab and zero-design; 1D formula; sparse mixture normalization and direct mixture density; score g=0; independent quadrature of forest-nuisance analogue identity; invalid input rejection.\n")
