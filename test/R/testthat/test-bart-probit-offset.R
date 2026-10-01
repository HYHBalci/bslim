# Regression checks for the latent-probit offset in both sampler paths.
testthat::test_that('probit GFR and MCMC recover imbalanced constant probabilities', {
 for(p in c(.2,.5,.8))for(mode in c('gfr','mcmc')){
  n<-300L;X<-matrix(0,n,1);Y<-as.numeric(seq_len(n)<=round(n*p));set.seed(6610+100*p)
  f<-bart(X,Y,X_test=X,num_gfr=if(mode=='gfr')100L else 5L,num_burnin=if(mode=='gfr')0L else 100L,num_mcmc=if(mode=='gfr')0L else 250L,general_params=list(verbose=FALSE,random_seed=6610+100*p,sample_sigma2_global=FALSE),mean_forest_params=list(num_trees=1L,probit_outcome_model=TRUE))
  testthat::expect_lt(abs(mean(pnorm(f$y_hat_test))-p),.05)
  testthat::expect_equal(f$model_params$outcome_mean,qnorm(mean(Y)),tolerance=1e-12)
  testthat::expect_equal(predict(f,X)$y_hat,f$y_hat_test,tolerance=1e-10)
 }
})
