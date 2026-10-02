# Tests exercise the R-only score core; BART integration has a separate smoke runner.
make_orthogonal_data <- function(n=240L) {
  set.seed(823)
  X <- matrix(rnorm(n*3),n,3); colnames(X)<-paste0('X',1:3)
  e <- plogis(-.8+1.2*X[,1]); Z <- rbinom(n,1,e)
  m <- sin(X[,1])+.5*X[,2]^2
  Y <- m+.7*(Z-e)+rnorm(n)
  list(X=X,Z=Z,Y=Y,m=m,e=e)
}
fit_oracle <- function(d,...) bcf_orthogonal_score_test(d$X,d$Z,d$Y,nuisance_method='supplied',m_hat=d$m,e_hat=d$e,nuisance_source='oracle',num_null_draws=99,seed=91,...)

testthat::test_that('projected score equals the known-nuisance analytic formula', {
 d<-make_orthogonal_data();f<-fit_oracle(d,variance_correction='HC0')
 W<-d$Z-d$e; alpha<-sum(W*(d$Y-d$m))/sum(W^2)
 H<-scale(d$X,TRUE,FALSE);H<-sweep(H,2,sqrt(colMeans(H^2)),'/')
 c<-colSums(H*W^2)/sum(W^2);Hp<-sweep(H,2,c,'-')
 psi<-Hp*(W*(d$Y-d$m-alpha*W));U<-colSums(psi);V<-crossprod(scale(psi,TRUE,FALSE))
 testthat::expect_equal(f$alpha,alpha,tolerance=1e-12)
 testthat::expect_equal(f$score_sum,U,tolerance=1e-12)
 testthat::expect_equal(f$covariance,V,tolerance=1e-12)
 testthat::expect_equal(f$statistic_quad,drop(crossprod(U,solve(V,U))),tolerance=1e-10)
 testthat::expect_gt(abs(c[1]),.05)
 testthat::expect_equal(colSums(f$projected_basis*W^2),rep(0,3),ignore_attr=TRUE,tolerance=1e-10)
})

testthat::test_that('treatment relabeling and affine outcome changes preserve evidence', {
 d<-make_orthogonal_data();f<-fit_oracle(d)
 reverse<-d;reverse$Z<-1-d$Z;reverse$e<-1-d$e;g<-fit_oracle(reverse)
 affine<-d;affine$Y<-11-2*d$Y;affine$m<-11-2*d$m;h<-fit_oracle(affine)
 testthat::expect_equal(g$alpha,-f$alpha,tolerance=1e-12)
 testthat::expect_equal(h$alpha,-2*f$alpha,tolerance=1e-12)
 for(z in list(g,h)){
  testthat::expect_equal(z$statistic_quad,f$statistic_quad,tolerance=1e-10)
  testthat::expect_equal(z$p_value_max,f$p_value_max)
 }
 A<-matrix(c(1,.2,0,0,2,.3,.1,0,1),3,3)
 b<-bcf_orthogonal_score_test(d$X,d$Z,d$Y,score_basis=d$X%*%A,nuisance_method='supplied',m_hat=d$m,e_hat=d$e,nuisance_source='oracle',num_null_draws=0)
 testthat::expect_equal(b$statistic_quad,f$statistic_quad,tolerance=1e-9)
})

testthat::test_that('OOF provenance rejects fold leakage and unlabelled predictions', {
 d<-make_orthogonal_data();fold<-rep(1:3,length.out=nrow(d$X));tr<-lapply(unique(fold),function(k)which(fold!=k))
 args<-list(X_train=d$X,Z_train=d$Z,y_train=d$Y,nuisance_method='supplied',m_hat=d$m,e_hat=d$e,num_null_draws=0)
 testthat::expect_error(do.call(bcf_orthogonal_score_test,args),'fold_id')
 args$fold_id<-fold;args$training_indices<-tr
 testthat::expect_s3_class(do.call(bcf_orthogonal_score_test,args),'bcforthogonalscore')
 args$training_indices[[1]]<-c(tr[[1]],which(fold==1)[1])
 testthat::expect_error(do.call(bcf_orthogonal_score_test,args),'leakage')
})

testthat::test_that('custom learners receive only training outcomes', {
 d<-make_orthogonal_data(90);X<-cbind(id=1:90,d$X);fold<-rep(1:3,30)
 calls<-list()
 learner<-function(X_train,Z_train,y_train,X_test){
  testthat::expect_length(intersect(X_train[,1],X_test[,1]),0)
  testthat::expect_equal(y_train,d$Y[X_train[,1]])
  calls[[length(calls)+1]]<<-list(train=X_train[,1],test=X_test[,1])
  list(m_hat=rep(mean(y_train),nrow(X_test)),e_hat=rep(.5,nrow(X_test)))
 }
 fit<-bcf_orthogonal_score_test(X,d$Z,d$Y,score_basis=d$X,nuisance_method='custom',fold_id=fold,nuisance_fitter=learner,num_null_draws=0)
 testthat::expect_length(calls,3)
 testthat::expect_equal(sort(unlist(lapply(calls,'[[','test'))),1:90,ignore_attr=TRUE)
})

testthat::test_that('rank, range and missingness safeguards fail clearly', {
 d<-make_orthogonal_data()
 testthat::expect_error(fit_oracle(d,score_basis=cbind(d$X,1)),'rank deficient')
 f<-fit_oracle(d,score_basis=cbind(d$X,d$X[,1]),rank_action='drop')
 testthat::expect_equal(f$df,3)
 testthat::expect_length(f$diagnostics$dropped_columns,1)
 bad<-d;bad$e[1]<-0;testthat::expect_error(fit_oracle(bad),'strictly inside')
 bad<-d;bad$m[1]<-NA_real_;testthat::expect_error(fit_oracle(bad),'finite numeric')
 bad<-d;bad$Z[1]<-2;testthat::expect_error(fit_oracle(bad),'binary')
 testthat::expect_warning(g<-fit_oracle(d,propensity_clip=.3),'clipped')
 testthat::expect_gt(g$diagnostics$clipped,0)
 testthat::expect_gte(min(g$nuisance_predictions$e_hat),.3)
})

testthat::test_that('multiplier probabilities are reproducible, nonzero and RNG-local', {
 d<-make_orthogonal_data();set.seed(732);state<-.Random.seed
 f<-fit_oracle(d);testthat::expect_identical(.Random.seed,state)
 g<-fit_oracle(d)
 testthat::expect_equal(f$p_value_max,g$p_value_max)
 testthat::expect_gte(f$p_value_max,.01)
 testthat::expect_gte(f$p_value_quad_multiplier,.01)
 testthat::expect_lte(f$p_value_max,1)
})

testthat::test_that('BART budgets and learner parameter lists are checked before fitting', {
 d<-make_orthogonal_data()
 testthat::expect_error(bcf_orthogonal_score_test(d$X,d$Z,d$Y,num_mcmc=0),'retained draw')
 testthat::expect_error(bcf_orthogonal_score_test(d$X,d$Z,d$Y,bart_params=3),'must be lists')
})

testthat::test_that('automatic folds are independent of treatment values', {
 d<-make_orthogonal_data()
 learner<-function(X_train,Z_train,y_train,X_test)list(m_hat=rep(0,nrow(X_test)),e_hat=rep(.5,nrow(X_test)))
 a<-bcf_orthogonal_score_test(d$X,d$Z,d$Y,nuisance_method='custom',nuisance_fitter=learner,num_folds=3,num_null_draws=0,seed=44)
 b<-bcf_orthogonal_score_test(d$X,1-d$Z,d$Y,nuisance_method='custom',nuisance_fitter=learner,num_folds=3,num_null_draws=0,seed=44)
 testthat::expect_identical(a$fold_id,b$fold_id)
 testthat::expect_true(all(vapply(a$training_indices,function(tr)length(unique(d$Z[tr]))==2,logical(1))))
})
