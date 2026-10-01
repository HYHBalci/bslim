# Analytic, provenance and integration checks for the experimental AIPW API.
aipw_data <- function(n = 180L) {
  set.seed(623); X <- matrix(rnorm(n * 3), n, 3); colnames(X) <- paste0('X', 1:3)
  e <- .2 + .6 * plogis(X[,1]); Z <- rbinom(n, 1, e)
  m0 <- sin(X[,1]) + .5 * X[,2]^2; m1 <- m0 + .7 + .25 * X[,1]
  Y <- (1-Z)*m0 + Z*m1 + rnorm(n)
  list(X=X,Z=Z,Y=Y,m0=m0,m1=m1,e=e)
}
aipw_oracle <- function(d, ...) bcf_aipw_score_test(d$X,d$Z,d$Y,nuisance_method='supplied',m0_hat=d$m0,m1_hat=d$m1,e_hat=d$e,nuisance_source='oracle',num_null_draws=99,seed=74,...)

testthat::test_that('weighted AIPW moments equal the independent analytic construction', {
 d<-aipw_data()
 for(w in c('overlap','uniform'))for(v in c('HC0','HC3')){
  f<-aipw_oracle(d,weighting=w,variance_correction=v)
  pseudo<-d$m1-d$m0+d$Z/d$e*(d$Y-d$m1)-(1-d$Z)/(1-d$e)*(d$Y-d$m0)
  omega<-if(w=='overlap')d$e*(1-d$e)else rep(1,length(d$Y))
  alpha<-weighted.mean(pseudo,omega);proj<-colSums(d$X*omega)/sum(omega);B<-sweep(d$X,2,proj,'-')
  psi<-B*(omega*(pseudo-alpha));lev<-hat(cbind(1,B)*sqrt(omega),intercept=FALSE)
  a<-if(v=='HC3')psi/(1-lev)else psi;V<-crossprod(scale(a,TRUE,FALSE));U<-colSums(psi)
  testthat::expect_equal(f$t,pseudo*omega,tolerance=1e-12)
  testthat::expect_equal(f$alpha,alpha,tolerance=1e-12)
  testthat::expect_equal(f$covariance,V,tolerance=1e-10)
  testthat::expect_equal(f$statistic_quad,drop(crossprod(U,solve(V,U))),tolerance=1e-10)
  testthat::expect_equal(colSums(f$projected_basis*omega),rep(0,3),ignore_attr=TRUE,tolerance=1e-10)
 }
})

testthat::test_that('conditional mean cancellation holds with correct propensity and wrong outcome models', {
 set.seed(821);base<-matrix(rnorm(120),40,3);X<-rbind(base,base);e0<-.15+.7*plogis(base[,1]);e<-rep(e0,2);Z<-rep(c(0,1),each=40)
 true0<-sin(base[,1])+base[,2]^2;true1<-true0+.7+.3*base[,1];Y<-c(true0,true1)
 fold<-rep(1:2,40);tr<-lapply(unique(fold),function(k)which(fold!=k))
 f<-bcf_aipw_score_test(X,Z,Y,nuisance_method='supplied',m0_hat=rep(0,80),m1_hat=X[,3],e_hat=e,nuisance_source=c(m0='cross_fitted',m1='cross_fitted',e='oracle'),fold_id=fold,training_indices=tr,num_null_draws=0)
 # Fixed misspecified prediction rules are independent of held-out labels.
 conditional_t<-(1-e0)*f$t[1:40]+e0*f$t[41:80]
 testthat::expect_equal(conditional_t,e0*(1-e0)*(true1-true0),tolerance=1e-12)
})

testthat::test_that('true outcome means cancel propensity misspecification under constant CATE', {
 set.seed(821);base<-matrix(rnorm(120),40,3);X<-do.call(rbind,rep(list(base),4));Z<-rep(c(0,0,1,1),each=40)
 true_e<-.2+.6*plogis(base[,1]);wrong_e<-.12+.74*plogis(-base[,1]+.8*base[,2]);m0<-sin(base[,1])+base[,2]^2;m1<-m0+.7;noise<-1+.1*base[,3]
 Y<-c(m0+noise,m0-noise,m1+noise,m1-noise);fold<-rep(1:2,80);tr<-lapply(unique(fold),function(k)which(fold!=k))
 for(w in c('overlap','uniform')){
  f<-bcf_aipw_score_test(X,Z,Y,nuisance_method='supplied',m0_hat=rep(m0,4),m1_hat=rep(m1,4),e_hat=rep(wrong_e,4),nuisance_source=c(m0='oracle',m1='oracle',e='cross_fitted'),fold_id=fold,training_indices=tr,weighting=w,num_null_draws=0)
  t0<-(f$t[1:40]+f$t[41:80])/2;t1<-(f$t[81:120]+f$t[121:160])/2
  omega<-if(w=='overlap')wrong_e*(1-wrong_e)else rep(1,40)
  testthat::expect_equal((1-true_e)*t0+true_e*t1,omega*.7,tolerance=1e-12)
  testthat::expect_equal(f$alpha,.7,tolerance=1e-12)
 }
})

testthat::test_that('relabeling, affine outcomes and nonsingular basis changes preserve quadratic evidence', {
 d<-aipw_data()
 for(w in c('overlap','uniform'))for(v in c('HC0','HC3')){
  f<-aipw_oracle(d,weighting=w,variance_correction=v)
  r<-d;r$Z<-1-d$Z;r$e<-1-d$e;r$m0<-d$m1;r$m1<-d$m0;g<-aipw_oracle(r,weighting=w,variance_correction=v)
  a<-d;a$Y<-7-2*d$Y;a$m0<-7-2*d$m0;a$m1<-7-2*d$m1;h<-aipw_oracle(a,weighting=w,variance_correction=v)
  testthat::expect_equal(g$alpha,-f$alpha,tolerance=1e-12);testthat::expect_equal(h$alpha,-2*f$alpha,tolerance=1e-12)
  for(z in list(g,h)){testthat::expect_equal(z$statistic_quad,f$statistic_quad,tolerance=1e-10);testthat::expect_equal(z$p_value_max,f$p_value_max)}
  A<-matrix(c(1,.2,0,0,2,.3,.1,0,1),3,3)
  b<-aipw_oracle(d,score_basis=sweep(d$X%*%A,2,c(2,-4,.5),'+'),weighting=w,variance_correction=v)
  testthat::expect_equal(b$statistic_quad,f$statistic_quad,tolerance=1e-9)
 }
})
testthat::test_that('mixed per-nuisance provenance audits only predictions declared cross-fitted', {
 d<-aipw_data();n<-nrow(d$X);f0<-rep(c('b','a','c'),length.out=n);f1<-rep(1:2,length.out=n)
 tr0<-lapply(unique(f0),function(k)which(f0!=k&d$Z==0));tr1<-lapply(unique(f1),function(k)which(f1!=k&d$Z==1))
 args<-list(X_train=d$X,Z_train=d$Z,y_train=d$Y,nuisance_method='supplied',m0_hat=d$m0,m1_hat=d$m1,known_propensity=d$e,num_null_draws=0,nuisance_source=c(m0='oracle',m1='cross_fitted'),fold_id=list(m1=f1),training_indices=list(m1=tr1))
 f<-do.call(bcf_aipw_score_test,args)
 testthat::expect_identical(f$diagnostics$nuisance_source,c(m0='oracle',m1='cross_fitted',e='known'))
 testthat::expect_null(f$nuisance_provenance$m0$fold_id);testthat::expect_identical(f$nuisance_provenance$m1$fold_id,f1)
 args$nuisance_source<-c(m0='cross_fitted',m1='cross_fitted');args$fold_id<-list(m0=f0,m1=f1);args$training_indices<-list(m0=tr0,m1=tr1)
 testthat::expect_s3_class(do.call(bcf_aipw_score_test,args),'bcfaipwscore')
 args$training_indices$m1[[1]]<-c(tr1[[1]],which(f1==unique(f1)[1])[1]);testthat::expect_error(do.call(bcf_aipw_score_test,args),'leakage')
 args$training_indices<-NULL;testthat::expect_error(do.call(bcf_aipw_score_test,args),'training_indices')
 args$fold_id<-NULL;testthat::expect_error(do.call(bcf_aipw_score_test,args),'fold_id')
})

testthat::test_that('custom callbacks see only training labels and have complete held-out coverage', {
 d<-aipw_data(90);X<-cbind(id=seq_len(90),d$X);calls<-list()
 learner<-function(X_train,Z_train,y_train,X_test){
  testthat::expect_length(intersect(X_train[,1],X_test[,1]),0)
  testthat::expect_equal(y_train,d$Y[X_train[,1]])
  testthat::expect_equal(Z_train,d$Z[X_train[,1]])
  calls[[length(calls)+1L]]<<-list(train=X_train[,1],test=X_test[,1])
  list(m0_hat=rep(mean(y_train[Z_train==0]),nrow(X_test)),m1_hat=rep(mean(y_train[Z_train==1]),nrow(X_test)))
 }
 f<-bcf_aipw_score_test(X,d$Z,d$Y,score_basis=d$X,nuisance_method='custom',nuisance_fitter=learner,known_propensity=d$e,num_folds=3,num_null_draws=0,seed=43)
 testthat::expect_length(calls,3);testthat::expect_equal(sort(unlist(lapply(calls,'[[','test'))),1:90,ignore_attr=TRUE)
 testthat::expect_identical(f$nuisance_provenance$e$source,'known')
 for(k in seq_along(calls))testthat::expect_identical(f$nuisance_provenance$m0$training_indices[[k]],as.integer(calls[[k]]$train))
})

testthat::test_that('automatic folds do not depend on treatment or outcome values', {
 d<-aipw_data();learner<-function(X_train,Z_train,y_train,X_test)list(m0_hat=rep(0,nrow(X_test)),m1_hat=rep(1,nrow(X_test)),e_hat=rep(.5,nrow(X_test)))
 f<-bcf_aipw_score_test(d$X,d$Z,d$Y,nuisance_method='custom',nuisance_fitter=learner,num_folds=3,num_null_draws=0,seed=48)
 g<-bcf_aipw_score_test(d$X,1-d$Z,-d$Y,nuisance_method='custom',nuisance_fitter=learner,num_folds=3,num_null_draws=0,seed=48)
 testthat::expect_identical(f$fold_id,g$fold_id)
})

testthat::test_that('GLM models are arm-specific and their saved predictions reproduce independently', {
 d<-aipw_data();f<-bcf_aipw_score_test(d$X,d$Z,d$Y,nuisance_method='glm',num_folds=3,num_null_draws=0,seed=83)
 for(k in seq_along(unique(f$fold_id))){te<-which(f$fold_id==unique(f$fold_id)[k]);Dte<-cbind(1,d$X[te,,drop=FALSE])
  for(a in 0:1){key<-paste0('m',a);tr<-f$nuisance_provenance[[key]]$training_indices[[k]];testthat::expect_true(all(d$Z[tr]==a));testthat::expect_length(intersect(tr,te),0);ref<-lm.fit(cbind(1,d$X[tr,,drop=FALSE]),d$Y[tr]);testthat::expect_equal(f$nuisance_predictions[[paste0(key,'_hat')]][te],drop(Dte%*%ref$coefficients),ignore_attr=TRUE,tolerance=1e-12)}
 }
 testthat::expect_true(all(is.finite(f$nuisance_predictions$e_hat)))
 testthat::expect_true(all(f$nuisance_predictions$e_hat>0&f$nuisance_predictions$e_hat<1))
})

testthat::test_that('BART predicts both treatment arms and probit propensity with deterministic seeds', {
 d<-aipw_data(80)
 args<-list(X_train=d$X,Z_train=d$Z,y_train=d$Y,nuisance_method='bart',num_folds=2,num_gfr=2,num_burnin=3,num_mcmc=5,mean_forest_params=list(num_trees=5),num_null_draws=0,seed=91)
 set.seed(512);state<-.Random.seed;f<-do.call(bcf_aipw_score_test,args);testthat::expect_identical(.Random.seed,state)
 g<-do.call(bcf_aipw_score_test,args);testthat::expect_identical(f$nuisance_predictions,g$nuisance_predictions)
 testthat::expect_true(is.finite(f$p_value_quad));testthat::expect_true(all(f$nuisance_predictions$e_hat>0&f$nuisance_predictions$e_hat<1))
 for(key in c('m0','m1','e'))testthat::expect_length(f$nuisance_provenance[[key]]$seeds,2)
 args$known_propensity<-d$e;h<-do.call(bcf_aipw_score_test,args);testthat::expect_equal(h$nuisance_predictions$e_hat,d$e);testthat::expect_length(h$nuisance_provenance$e$training_indices,0)
})

testthat::test_that('RNG state is restored after success, failure, and absent initial state', {
 d<-aipw_data();set.seed(732);state<-.Random.seed;f<-aipw_oracle(d);testthat::expect_identical(.Random.seed,state);g<-aipw_oracle(d);testthat::expect_equal(f$p_value_max,g$p_value_max);testthat::expect_gte(f$p_value_max,.01);testthat::expect_gte(f$p_value_quad_multiplier,.01)
 fail<-function(X_train,Z_train,y_train,X_test){runif(3);stop('deliberate callback error')}
 testthat::expect_error(bcf_aipw_score_test(d$X,d$Z,d$Y,nuisance_method='custom',nuisance_fitter=fail,num_null_draws=0),'deliberate');testthat::expect_identical(.Random.seed,state)
 rm('.Random.seed',envir=.GlobalEnv);aipw_oracle(d);testthat::expect_false(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE));assign('.Random.seed',state,envir=.GlobalEnv)
})

testthat::test_that('rank, overlap, source and numeric safeguards reject unsupported inputs', {
 d<-aipw_data();testthat::expect_error(aipw_oracle(d,score_basis=cbind(d$X,1)),'rank deficient');f<-aipw_oracle(d,score_basis=cbind(d$X,d$X[,1]),rank_action='drop');testthat::expect_equal(f$df,3)
 bad<-d;bad$e[1]<-0;testthat::expect_error(aipw_oracle(bad),'strictly inside');testthat::expect_warning(aipw_oracle(bad,propensity_clip=.1),'clipped')
 bad$e[1]<- -1;testthat::expect_error(aipw_oracle(bad,propensity_clip=.1),'probabilities')
 bad<-d;bad$m0[1]<-NA_real_;testthat::expect_error(aipw_oracle(bad),'finite numeric')
 testthat::expect_error(bcf_aipw_score_test(d$X,d$Z,d$Y,num_mcmc=0),'retained draw')
 testthat::expect_error(bcf_aipw_score_test(d$X,d$Z,d$Y,nuisance_method='supplied',m0_hat=d$m0,m1_hat=d$m1,e_hat=d$e,nuisance_source=c(m0='oracle',m1='oracle'),num_null_draws=0),'declare')
 testthat::expect_error(bcf_aipw_score_test(d$X,d$Z,d$Y,nuisance_method='supplied',m0_hat=d$m0,m1_hat=d$m1,e_hat=d$e,known_propensity=d$e,nuisance_source='oracle',num_null_draws=0),'either')
 zero<-d;zero$Y<-zero$m0<-zero$m1<-rep(0,length(d$Y));testthat::expect_error(aipw_oracle(zero),'singular')
})
testthat::test_that('explicitly clipped known propensities remain visibly modified design inputs', {
 d<-aipw_data()
 testthat::expect_warning(f<-bcf_aipw_score_test(d$X,d$Z,d$Y,nuisance_method='supplied',m0_hat=d$m0,m1_hat=d$m1,nuisance_source=c(m0='oracle',m1='oracle'),known_propensity=d$e,propensity_clip=.4,num_null_draws=0),'clipped')
 testthat::expect_true(f$diagnostics$known_propensity)
 testthat::expect_true(f$diagnostics$known_propensity_modified)
 testthat::expect_identical(f$diagnostics$used_propensity_source,'clipped_known')
 testthat::expect_identical(f$nuisance_provenance$e$learner,'design')
 testthat::expect_equal(f$nuisance_predictions$raw_e_hat,d$e)
 testthat::expect_gt(length(f$nuisance_provenance$e$prediction_transform$clipped_rows),0)
})