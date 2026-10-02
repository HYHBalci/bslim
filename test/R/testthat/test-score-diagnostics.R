# R-only tests: no forest fit or C++ rebuild is needed.
make_score_diagnostic_fixture <- function(unequal=FALSE,M=160L,C=2L) {
 set.seed(93930);n<-80L;H<-matrix(rnorm(n*2),n,2,dimnames=list(NULL,c('H1','H2')));H<-scale(H,TRUE,FALSE)
 e<-if(unequal)plogis(.3+.7*H[,1])else rep(.5,n);z<-rbinom(n,1,e);y<-sin(H[,1])+.6*H[,2]+rnorm(n)
 mu<-matrix(NA_real_,n,M*C)
 for(k in seq_len(C)){state<-rnorm(n,sd=.1);for(m in seq_len(M)){state<-.7*state+rnorm(n,sd=.06);mu[,(k-1L)*M+m]<-.3*y+state+(k-1.5)*.3*H[,1]}}
 a<-matrix(rnorm(C*M,.7,.15),C,M);sig<-matrix(exp(rnorm(C*M,0,.15)),C,M)
 fit<-list(test_config=list(test_type='mean',use_rao_blackwell=TRUE,null_method='gaussian',standardize_residuals=FALSE,p_valid=ncol(H),propensity_clip=NULL),alpha_samples=a,sigma2_samples=sig,mu_hat_train=mu,y_std_train=2.3)
 list(fit=fit,y=y,z=z,e=e,H=H,M=M,C=C)
}
manual_score_chain <- function(d,k,standardize=FALSE,ridge=1e-8) {
 ix<-(k-1L)*d$M+seq_len(d$M);W<-d$z-d$e
 r<-sweep(d$fit$mu_hat_train[,ix,drop=FALSE],1,d$y,function(mm,yy)yy-mm)-outer(W,d$fit$alpha_samples[k,])
 r<-if(standardize)sweep(r,2,sqrt(d$fit$sigma2_samples[k,]),'/')else r/d$fit$y_std_train
 s<-rowMeans(r)*W;s<-s-mean(s);U<-drop(crossprod(d$H,s));V<-crossprod(d$H*s)
 list(U=U,V=V,stat=drop(crossprod(U,solve(V+diag(ridge,ncol(d$H)),U))),p=pchisq(drop(crossprod(U,solve(V+diag(ridge,ncol(d$H)),U))),ncol(d$H),lower.tail=FALSE))
}
call_score_diag <- function(d,...)bcf_score_diagnostics(d$fit,d$y,d$z,d$e,d$H,batch_size=8L,...)

testthat::test_that('chain-major ordering and original-unit ridge reconstruction are exact',{
 d<-make_score_diagnostic_fixture();manual<-lapply(seq_len(d$C),function(k)manual_score_chain(d,k));d$fit$rb_pval_quad<-vapply(manual,function(x)x$p,numeric(1));f<-call_score_diag(d)
 testthat::expect_lt(f$package_quad_reconstruction_error,1e-12)
 for(k in seq_len(d$C)){testthat::expect_equal(f$per_chain[[k]]$U,manual[[k]]$U,tolerance=1e-12);testthat::expect_equal(f$per_chain[[k]]$V,manual[[k]]$V,tolerance=1e-12)}
 testthat::expect_equal(f$chain_count,2L);testthat::expect_equal(f$draws_per_chain,d$M)
})

testthat::test_that('pooling contributions differs from averaging p-values and covariance',{
 d<-make_score_diagnostic_fixture();n<-length(d$y)
 for(k in seq_len(d$C)){ix<-(k-1L)*d$M+seq_len(d$M);d$fit$mu_hat_train[,ix]<-d$fit$mu_hat_train[,ix]+if(k==1).9*d$H[,1]else -.9*d$H[,1]}
 d$fit$rb_pval_quad<-vapply(seq_len(d$C),function(k)manual_score_chain(d,k)$p,numeric(1));f<-call_score_diag(d)
 testthat::expect_gt(abs(f$pooled$p-mean(d$fit$rb_pval_quad)),.005)
 testthat::expect_gt(max(abs(f$pooled$V-Reduce('+',lapply(f$per_chain,'[[','V'))/d$C)),.01)
 testthat::expect_equal(f$pooled$U,Reduce('+',lapply(f$per_chain,'[[','U'))/d$C,tolerance=1e-12)
})

testthat::test_that('per-draw sigma weighting is reconstructed for the mean test',{
 d<-make_score_diagnostic_fixture();d$fit$test_config$standardize_residuals<-TRUE
 manual<-lapply(seq_len(d$C),function(k)manual_score_chain(d,k,TRUE));d$fit$rb_pval_quad<-vapply(manual,function(x)x$p,numeric(1));f<-call_score_diag(d)
 testthat::expect_lt(f$package_quad_reconstruction_error,1e-12)
 for(k in seq_len(d$C))testthat::expect_equal(f$per_chain[[k]]$U,manual[[k]]$U,tolerance=1e-12)
})

testthat::test_that('unequal propensity retains alpha contribution and order',{
 d<-make_score_diagnostic_fixture(TRUE);d$fit$alpha_samples[1,]<-1+seq_len(d$M)/d$M;d$fit$alpha_samples[2,]<- -1+2*seq_len(d$M)/d$M
 manual<-lapply(seq_len(d$C),function(k)manual_score_chain(d,k));d$fit$rb_pval_quad<-vapply(manual,function(x)x$p,numeric(1));f<-call_score_diag(d)
 testthat::expect_lt(f$package_quad_reconstruction_error,1e-12)
 for(k in seq_len(d$C))testthat::expect_equal(f$per_chain[[k]]$U,manual[[k]]$U,tolerance=1e-12)
 d$fit$alpha_samples[]<-0;g<-call_score_diag(d);testthat::expect_gt(max(abs(f$pooled$U-g$pooled$U)),.01)
})

testthat::test_that('full delta derivative includes covariance change',{
 d<-make_score_diagnostic_fixture();W<-d$z-d$e;n<-length(W);H<-d$H;s<-W*(d$y-rowMeans(d$fit$mu_hat_train));s<-s-mean(s)
 set.seed(192);ds<-rnorm(n);ds<-ds-mean(ds);point<-.bcf_score_diag_quad(H,s);a<-drop(solve(point$V/n,point$U/sqrt(n)))
 analytic<-2*sum(drop(H%*%a)*ds)/sqrt(n)-2*sum(drop(H%*%a)^2*s*ds)/n;epsilon<-1e-6
 finite<-(.bcf_score_diag_quad(H,s+epsilon*ds)$statistic-.bcf_score_diag_quad(H,s-epsilon*ds)$statistic)/(2*epsilon)
 testthat::expect_equal(analytic,finite,tolerance=1e-7)
 f<-call_score_diag(d);testthat::expect_true(f$mc$includes_covariance_change);testthat::expect_true(is.finite(f$mc$p_mcse));testthat::expect_gte(f$mc$p_mcse,0)
})

testthat::test_that('singular score covariance refuses primary ridge-calibrated p-value',{
 d<-make_score_diagnostic_fixture();d$H<-cbind(d$H,d$H[,1]);d$fit$test_config$p_valid<-3L;f<-call_score_diag(d)
 testthat::expect_false(f$pooled$covariance_usable);testthat::expect_true(is.na(f$pooled$p));testthat::expect_equal(f$pooled$covariance_rank,2L);testthat::expect_true(is.finite(f$pooled$legacy_ridge_p));testthat::expect_true(is.na(f$mc$p_mcse))
})

testthat::test_that('invalid dimensions, marks and budgets fail before diagnosis',{
 d<-make_score_diagnostic_fixture();bad<-d;bad$fit$mu_hat_train<-bad$fit$mu_hat_train[,-1];testthat::expect_error(call_score_diag(bad),'draw shapes')
 bad<-d;bad$fit$sigma2_samples[1,1]<-0;testthat::expect_error(call_score_diag(bad),'sigma2')
 testthat::expect_error(bcf_score_diagnostics(d$fit,d$y,d$z,d$e,d$H,batch_size=40),'eight complete')
 testthat::expect_error(call_score_diag(d,alpha=2),'alpha')
 bad<-d;bad$e[1]<-0;testthat::expect_error(call_score_diag(bad),'propensity marks')
 bad<-d;bad$fit$test_config$test_type<-'variance';testthat::expect_error(call_score_diag(bad),'mean-test')
 bad<-d;bad$fit$test_config$null_method<-'permutation';testthat::expect_error(call_score_diag(bad),'Gaussian')
})

testthat::test_that('diagnostics preserve caller RNG state',{
 d<-make_score_diagnostic_fixture();set.seed(438);state<-.Random.seed;invisible(call_score_diag(d));testthat::expect_identical(.Random.seed,state)
})
testthat::test_that('blank and duplicate prepared-design names remain diagnosable',{
 d<-make_score_diagnostic_fixture();colnames(d$H)<-c('','');f<-call_score_diag(d,keep_draws=TRUE)
 testthat::expect_length(unique(dimnames(f$linearized_draws)[[3L]]),dim(f$linearized_draws)[3L]);testthat::expect_true(f$score_diagnostics$available)
 colnames(d$H)<-c('X','X');g<-call_score_diag(d);testthat::expect_equal(g$pooled$p,f$pooled$p)
})