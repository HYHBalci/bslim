# Diagnostic-selected follow-up, not an additional simulation performance sample.
args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1] else "D:/heterogeneity-bayes-research-2026-10-16/pilot-results"
source("research/bayes-factor-pilot-2026-10-16/bayes_factor_sidecar.R")
.libPaths(c(".research-library/fresh-2026-10-01",.libPaths()))
suppressPackageStartupMessages(library(stochtree))
api<-new.env(parent=asNamespace("stochtree"));sys.source("R/bcf_score_test.R",api)
cfg<-dget(file.path(out,"config_full.R"))
baseline_file<-if(length(args)>1L)args[2] else "results_full.csv"
results<-read.csv(file.path(out,baseline_file),stringsAsFactors=FALSE)
sink(file.path(out,"long_chain_diagnostics.log"),split=TRUE)
cat("Diagnostic-selected follow-up start",format(Sys.time(),tz="Europe/Amsterdam"),"\n")
chosen<-subset(results,method=="null_posterior_integrated_bf" & prior=="dense" &
                 scale==.5 & basis=="linear" & scenario %in% c("null","sparse"))
stopifnot(length(unique(chosen$design))==2L,
          all(table(interaction(chosen$design,chosen$scenario))==cfg$forest_reps))
chosen<-do.call(rbind,lapply(split(chosen,interaction(chosen$design,chosen$scenario)),function(z)
  z[which.max(z$chain_log_bf_range),]))
write.csv(chosen,file.path(out,"diagnostic_selected_cases.csv"),row.names=FALSE)
long_cfg<-list(chains=4L,burnin=1000L,mcmc=3000L,slab_sd=.5,seed_offset=5000L,
  selection="Largest dense-prior between-chain log-BF range in each of 4 design x null/sparse cells")
dput(long_cfg,file.path(out,"long_chain_config.R"))
rows<-list()
for(j in seq_len(nrow(chosen))) {
  z<-chosen[j,];set.seed(z$seed)
  # Same generator as run_pilot.R, including all draws before Y.
  X<-matrix(rnorm(cfg$n*cfg$p),cfg$n,cfg$p);colnames(X)<-paste0("X",seq_len(cfg$p))
  e<-if(z$design=="randomized")rep(.5,cfg$n) else .15+.7*plogis(.6*X[,1]-.6*X[,2]+.4*X[,3]^2)
  m<-sin(X[,1])+.5*(X[,2]^2-1)+.5*X[,3]*X[,4]
  h<-if(z$scenario=="null")rep(0,cfg$n) else .5*X[,1]
  Z<-as.numeric(rbinom(cfg$n,1,e));W<-Z-e;Y<-m+W*(.5+h)+rnorm(cfg$n)
  H<-standardize_modifier_basis(X);D<-H*W
  cat("Refitting",z$design,z$scenario,"rep",z$rep,"\n")
  tt<-proc.time()[3]
  fit<-api$bcf_restricted_score_test(X,Z,Y,propensity_train=e,num_gfr=5,
    num_burnin=long_cfg$burnin,num_mcmc=long_cfg$mcmc,num_null_draws=0,
    general_params=list(verbose=FALSE,standardize=FALSE,sample_sigma2_global=FALSE,
      sigma2_global_init=1,alpha_prior_var=4,interaction_rule="none",
      random_seed=as.integer(z$seed+long_cfg$seed_offset),num_chains=long_cfg$chains),
    prognostic_forest_params=list(num_trees=25,sample_sigma2_leaf=FALSE,sigma2_leaf_init=.04))
  seconds<-proc.time()[3]-tt
  bf<-evaluate_bf_priors(Y-fit$mu_hat_train,D,W,.5,1,4)
  chain<-rep(seq_len(long_cfg$chains),each=long_cfg$mcmc)
  for(pr in names(bf)) {
    s<-importance_summary(bf[[pr]],chain)
    base<-subset(results,method=="null_posterior_integrated_bf" & basis=="linear" & scale==.5 &
      prior==pr & seed==z$seed)
    rows[[length(rows)+1L]]<-data.frame(design=z$design,scenario=z$scenario,rep=z$rep,seed=z$seed,
       prior=pr,baseline_log_bf=base$log_bf,long_log_bf=s$log_bf,
       baseline_weight_ess=base$importance_ess,long_weight_ess=s$importance_ess,
       baseline_chain_range=base$chain_log_bf_range,long_chain_range=s$chain_log_bf_range,
       baseline_max_weight=base$largest_weight_fraction,long_max_weight=s$largest_weight_fraction,
       long_flag=s$warning,fit_seconds=seconds,log_bf_chain1=s$chain_log_bf[1],
       log_bf_chain2=s$chain_log_bf[2],log_bf_chain3=s$chain_log_bf[3],log_bf_chain4=s$chain_log_bf[4])
  }
  trace<-data.frame(chain,draw=rep(seq_len(long_cfg$mcmc),long_cfg$chains),
    alpha=as.vector(t(fit$alpha_samples)),mu_mean=colMeans(fit$mu_hat_train),
    mu_projection_W=drop(crossprod(W,fit$mu_hat_train)),dense_log_ratio=bf$dense,sparse_log_ratio=bf$sparse)
  write.csv(trace,file.path(out,paste0("long_chain_trace_",z$design,"_",z$scenario,".csv")),row.names=FALSE)
  write.csv(do.call(rbind,rows),file.path(out,"long_chain_comparison.csv"),row.names=FALSE)
}
cat("Diagnostic follow-up end",format(Sys.time(),tz="Europe/Amsterdam"),"\n")
sink()
