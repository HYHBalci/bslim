# Run from the repository root. All result artifacts go to the requested D: path.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(getwd(), winslash = "/")
out <- if (length(args)) args[1] else "D:/heterogeneity-bayes-research-2026-10-16/pilot-results"
mode <- if (length(args) > 1L) args[2] else "full"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Output directory could not be created: ",out)
source(file.path(root,"research/bayes-factor-pilot-2026-10-16/bayes_factor_sidecar.R"))
.libPaths(c(file.path(root,".research-library/fresh-2026-10-01"),.libPaths()))
suppressPackageStartupMessages(library(stochtree))
# Source present-day R wrappers without changing the loaded package namespace.
api <- new.env(parent = asNamespace("stochtree"))
sys.source(file.path(root,"R/bcf_score_test.R"),api)
sys.source(file.path(root,"R/bcf_orthogonal_score_test.R"),api)
api$bart <- get("bart",asNamespace("stochtree"))
log_file <- file.path(out,paste0("run_",mode,".log"))
sink(log_file,split=TRUE); on.exit(sink(),add=TRUE)
cat("Start",format(Sys.time(),tz="Europe/Amsterdam"),"\n")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))

smoke <- mode == "smoke"
cfg <- list(n=300L,p=6L,oracle_reps=if(smoke)2L else 300L,
            forest_reps=if(smoke)2L else 40L,
            crossfit_reps=if(smoke)1L else 25L,
            burnin=if(smoke)50L else 300L,mcmc=if(smoke)100L else 750L,
            chains=2L,num_trees=25L,slab_sd=c(.25,.5,1),alpha_variance=4,
            sigma2=1,leaf_variance=.04,base_seed=261016L)
dput(cfg,file.path(out,paste0("config_",mode,".R")))
warning_rows <- list(); results <- list(); fit_checks <- list()
add_warning <- function(stage,scenario,design,rep,message) {
  warning_rows[[length(warning_rows)+1L]] <<- data.frame(stage,scenario,design,rep,message)
}
capture_case <- function(expr,stage,scenario,design,rep) {
  tryCatch(withCallingHandlers(expr,warning=function(w) {
    add_warning(stage,scenario,design,rep,conditionMessage(w)); invokeRestart("muffleWarning")
  }),error=function(e) {
    add_warning(stage,scenario,design,rep,paste("ERROR",conditionMessage(e))); NULL
  })
}
record <- function(stage,scenario,design,rep,seed,basis,method,prior="none",scale=NA_real_,
                   p=NA_real_,log_bf=NA_real_,diag=NULL,runtime=NA_real_) {
  d <- data.frame(stage,scenario,design,rep,seed,n=cfg$n,basis,method,prior,scale,p,log_bf,runtime,
                  importance_ess=NA_real_,importance_ess_fraction=NA_real_,
                  largest_weight_fraction=NA_real_,top_one_percent_weight_fraction=NA_real_,
                  chain_log_bf_range=NA_real_,batch_relative_se=NA_real_,diagnostic_warning="none")
  if (!is.null(diag)) {
    for (nm in intersect(names(d),names(diag))) d[[nm]] <- diag[[nm]]
    d$diagnostic_warning <- diag$warning
  }
  results[[length(results)+1L]] <<- d
}
flush_results <- function() {
  if(length(results)) write.csv(do.call(rbind,results),file.path(out,paste0("results_",mode,".csv")),row.names=FALSE)
  if(length(warning_rows)) write.csv(do.call(rbind,warning_rows),file.path(out,paste0("warnings_",mode,".csv")),row.names=FALSE)
  if(length(fit_checks)) write.csv(do.call(rbind,fit_checks),file.path(out,paste0("fit_diagnostics_",mode,".csv")),row.names=FALSE)
}
make_data <- function(scenario,design,seed) {
  set.seed(seed)
  X <- matrix(rnorm(cfg$n*cfg$p),cfg$n,cfg$p)
  colnames(X) <- paste0("X",seq_len(cfg$p))
  e <- if(design=="randomized") rep(.5,cfg$n) else
    .15+.7*plogis(.6*X[,1]-.6*X[,2]+.4*X[,3]^2)
  m <- sin(X[,1])+.5*(X[,2]^2-1)+.5*X[,3]*X[,4]
  h <- switch(scenario,null=rep(0,cfg$n),sparse=.5*X[,1],
              dense=.5*rowSums(X)/sqrt(cfg$p),quadratic=.6*(X[,1]^2-1)/sqrt(2))
  Z <- as.numeric(rbinom(cfg$n,1,e)); W <- Z-e
  Y <- m+W*(.5+h)+rnorm(cfg$n)
  H <- list(linear=standardize_modifier_basis(X))
  if(scenario=="quadratic") H$augmented <- standardize_modifier_basis(cbind(X,quad=(X[,1]^2-1)/sqrt(2)))
  list(X=X,e=e,m=m,Z=Z,W=W,Y=Y,H=H)
}
split_rhat <- function(mat) {
  m <- floor(nrow(mat)/2)
  chains <- cbind(mat[seq_len(m),,drop=FALSE],mat[nrow(mat)-m+seq_len(m),,drop=FALSE])
  W <- mean(apply(chains,2,var)); B <- m*var(colMeans(chains))
  sqrt(((m-1)/m*W+B/m)/W)
}
scenarios <- c("null","sparse","dense","quadratic")
designs <- c("randomized","confounded_known")

for (dd in seq_along(designs)) for(ss in seq_along(scenarios)) {
  design<-designs[dd]; scenario<-scenarios[ss]
  cat("Oracle",design,scenario,"\n")
  for(rr in seq_len(cfg$oracle_reps)) {
    seed <- cfg$base_seed+dd*100000L+ss*10000L+rr
    dat<-make_data(scenario,design,seed)
    for(basis in names(dat$H)) {
      H<-dat$H[[basis]]; D<-H*dat$W
      t0<-proc.time()[3]
      score<-capture_case(api$bcf_orthogonal_score_test(dat$X,dat$Z,dat$Y,score_basis=H,
                nuisance_method="supplied",m_hat=dat$m,known_propensity=dat$e,
                nuisance_source="oracle",num_null_draws=0,seed=seed),"oracle",scenario,design,rr)
      if(!is.null(score)) {
        record("oracle",scenario,design,rr,seed,basis,"orthogonal_score",p=score$p_value_quad,
               runtime=proc.time()[3]-t0)
        for(g in c(1,10,cfg$n)) record("oracle",scenario,design,rr,seed,basis,
          "approximate_score_bf",prior="normal_score",scale=g,
          log_bf=score_experiment_log_bf(score$statistic_quad,score$df,g))
      }
      for(s in cfg$slab_sd) {
        t0<-proc.time()[3]
        bfs<-evaluate_bf_priors(dat$Y-dat$m,D,dat$W,s,cfg$sigma2,cfg$alpha_variance)
        for(pr in names(bfs)) record("oracle",scenario,design,rr,seed,basis,
          "exact_conditional_gaussian_bf",prior=pr,scale=s,log_bf=bfs[[pr]],runtime=proc.time()[3]-t0)
      }
    }
  }
  flush_results()
}

for(dd in seq_along(designs)) for(ss in seq_along(scenarios)) {
  design<-designs[dd];scenario<-scenarios[ss]
  cat("Forest",design,scenario,"\n")
  for(rr in seq_len(cfg$forest_reps)) {
    seed<-cfg$base_seed+dd*100000L+ss*10000L+rr
    dat<-make_data(scenario,design,seed)
    t0<-proc.time()[3]
    fit<-capture_case(api$bcf_restricted_score_test(dat$X,dat$Z,dat$Y,
      propensity_train=dat$e,num_gfr=5,num_burnin=cfg$burnin,num_mcmc=cfg$mcmc,
      num_null_draws=0,general_params=list(verbose=FALSE,standardize=FALSE,
        sample_sigma2_global=FALSE,sigma2_global_init=cfg$sigma2,alpha_prior_var=cfg$alpha_variance,
        interaction_rule="none",random_seed=seed+1000L,num_chains=cfg$chains),
      prognostic_forest_params=list(num_trees=cfg$num_trees,sample_sigma2_leaf=FALSE,
        sigma2_leaf_init=cfg$leaf_variance)),"forest",scenario,design,rr)
    if(is.null(fit)) next
    fit_seconds<-proc.time()[3]-t0
    stopifnot(all(fit$sigma2_samples==cfg$sigma2),fit$y_std_train==1,fit$y_bar_train==0)
    chain_id<-rep(seq_len(cfg$chains),each=cfg$mcmc)
    residual<-dat$Y-fit$mu_hat_train
    rhat_alpha<-split_rhat(t(fit$alpha_samples))
    rhat_mu_mean<-split_rhat(matrix(colMeans(fit$mu_hat_train),cfg$mcmc,cfg$chains))
    rhat_mu_w<-split_rhat(matrix(drop(crossprod(dat$W,fit$mu_hat_train)),cfg$mcmc,cfg$chains))
    fit_checks[[length(fit_checks)+1L]]<-data.frame(design,scenario,rep=rr,seed,
      n=cfg$n,chains=cfg$chains,burnin=cfg$burnin,retained_per_chain=cfg$mcmc,
      alpha_split_rhat=rhat_alpha,mu_mean_split_rhat=rhat_mu_mean,
      mu_treatment_projection_split_rhat=rhat_mu_w,fit_seconds,
      forest_mse=mean((rowMeans(fit$mu_hat_train)-dat$m)^2),
      alpha_posterior_mean=mean(fit$alpha_samples),
      internal_score_p_chain1=fit$rb_pval_quad[1],internal_score_p_chain2=fit$rb_pval_quad[2])
    if(rr==1) saveRDS(list(data=dat,fit=fit,config=cfg),
      file.path(out,paste0("representative_",mode,"_",design,"_",scenario,".rds")),compress="xz")
    for(basis in names(dat$H)) {
      H<-dat$H[[basis]];D<-H*dat$W
      ms<-capture_case(restricted_moment_test(dat$Y,rowMeans(fit$mu_hat_train),
                   mean(fit$alpha_samples),dat$W,H),"forest_score",scenario,design,rr)
      if(!is.null(ms)) {
        record("forest",scenario,design,rr,seed,basis,"restricted_external_score",p=ms$p,runtime=fit_seconds)
        for(g in c(1,10,cfg$n)) record("forest",scenario,design,rr,seed,basis,
          "approximate_score_bf",prior="normal_score",scale=g,
          log_bf=score_experiment_log_bf(ms$Q,ms$df,g))
      }
      for(s in cfg$slab_sd) {
        t0<-proc.time()[3]
        bfs<-evaluate_bf_priors(residual,D,dat$W,s,cfg$sigma2,cfg$alpha_variance)
        plugin<-evaluate_bf_priors(dat$Y-rowMeans(fit$mu_hat_train),D,dat$W,s,cfg$sigma2,cfg$alpha_variance)
        for(pr in names(bfs)) {
          dg<-importance_summary(bfs[[pr]],chain_id)
          record("forest",scenario,design,rr,seed,basis,"null_posterior_integrated_bf",
                 prior=pr,scale=s,log_bf=dg$log_bf,diag=dg,runtime=proc.time()[3]-t0)
          record("forest",scenario,design,rr,seed,basis,"plugin_conditional_bf",
                 prior=pr,scale=s,log_bf=plugin[[pr]])
        }
        if(rr==1 && s==.5) write.csv(data.frame(chain=chain_id,draw=seq_along(chain_id),
          dense_log_ratio=bfs$dense,sparse_log_ratio=bfs$sparse),
          file.path(out,paste0("ratios_",mode,"_",design,"_",scenario,"_",basis,".csv")),row.names=FALSE)
      }
    }
    if(rr %% 5L==0L) {cat(" completed",rr,"of",cfg$forest_reps,"\n");flush_results()}
  }
  flush_results()
}

# A separate nuisance-learning sensitivity, using the existing cross-fitting API.
# Cross-fitted residual Gaussian evidence below is a plug-in diagnostic, not a full BF.
for(scenario in c("null","sparse")) for(known in c(TRUE,FALSE)) {
  design<-if(known)"confounded_known" else "confounded_learned"
  cat("Crossfit",design,scenario,"\n")
  for(rr in seq_len(cfg$crossfit_reps)) {
    ss<-match(scenario,scenarios);seed<-cfg$base_seed+200000L+ss*10000L+rr
    dat<-make_data(scenario,"confounded_known",seed)
    t0<-proc.time()[3]
    score<-capture_case(api$bcf_orthogonal_score_test(dat$X,dat$Z,dat$Y,score_basis=dat$H$linear,
        nuisance_method="bart",known_propensity=if(known)dat$e else NULL,num_folds=3L,
        num_gfr=5,num_burnin=75,num_mcmc=150,num_null_draws=0,seed=seed+2000L,
        mean_forest_params=list(num_trees=25L)),"crossfit",scenario,design,rr)
    if(is.null(score)) next
    record("crossfit",scenario,design,rr,seed,"linear","orthogonal_score",p=score$p_value_quad,
           runtime=proc.time()[3]-t0)
    pred<-score$nuisance_predictions;W<-dat$Z-pred$e_hat;H<-dat$H$linear
    bfs<-evaluate_bf_priors(dat$Y-pred$m_hat,H*W,W,.5,cfg$sigma2,cfg$alpha_variance)
    for(pr in names(bfs)) record("crossfit",scenario,design,rr,seed,"linear",
      "crossfit_plugin_bf",prior=pr,scale=.5,log_bf=bfs[[pr]])
    for(g in c(1,10,cfg$n)) record("crossfit",scenario,design,rr,seed,"linear",
      "approximate_score_bf",prior="normal_score",scale=g,
      log_bf=score_experiment_log_bf(score$statistic_quad,score$df,g))
  }
  flush_results()
}
cat("Completed",length(results),"result rows;",length(warning_rows),"warnings/errors\n")
cat("End",format(Sys.time(),tz="Europe/Amsterdam"),"\n")
if(!length(warning_rows))write.csv(data.frame(stage=character(),scenario=character(),design=character(),
   rep=integer(),message=character()),file.path(out,paste0("warnings_",mode,".csv")),row.names=FALSE)
flush_results()
source(file.path(root,"research/bayes-factor-pilot-2026-10-16/summarize_pilot.R"),local=TRUE)
