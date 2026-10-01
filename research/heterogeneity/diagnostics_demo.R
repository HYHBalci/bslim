# Fresh synthetic smoke for a built current source package; no external simulation files.
args <- commandArgs(TRUE)
if(length(setdiff(args,"--source-overlay"))) stop("Only --source-overlay is supported; set STOCHTREE_R_LIBRARY if needed.")
overlay <- "--source-overlay" %in% args
file_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
stopifnot(length(file_arg)==1L)
here <- dirname(normalizePath(sub("^--file=","",file_arg),winslash="/",mustWork=TRUE))
repo <- normalizePath(file.path(here,"../.."),winslash="/",mustWork=TRUE)
lib <- Sys.getenv("STOCHTREE_R_LIBRARY",file.path(repo,".research-library",paste0("R-",getRversion(),"-",R.version$platform)))
if(dir.exists(lib)) .libPaths(c(normalizePath(lib,winslash="/"),.libPaths()))
if(!requireNamespace("stochtree",quietly=TRUE)) stop("Build the source package using bootstrap.R --install first.")
if(!overlay && !"bcf_score_diagnostics" %in% getNamespaceExports("stochtree")) stop("Installed package is older than the current helper; build current source first or explicitly request --source-overlay.")
set.seed(26100102L);n<-80L
X<-matrix(runif(n*2L,-1,1),n,2L);z<-as.numeric(rbinom(n,1,.5))
y<-sin(pi*X[,1])+X[,2]^2+.5*(z-.5)+rnorm(n)
ns<-asNamespace("stochtree")
functions <- ns
if(overlay) {
  functions <- new.env(parent=ns)
  sys.source(file.path(repo,"R","bcf_score_test.R"),envir=functions)
  sys.source(file.path(repo,"R","bcf_score_diagnostics.R"),envir=functions)
}
handled<-get("standardize_X_by_index",ns)(X,process_data=FALSE,interaction_rule="continuous",cat_coding_method="difference")
H<-get("prepare_score_test_design",functions,inherits=TRUE)(handled$X_final,handled$X_final_var_info,"continuous")$X_cen
fit<-get("bcf_restricted_score_test",functions,inherits=TRUE)(X,z,y,propensity_train=rep(.5,n),
  num_gfr=5L,num_burnin=64L,num_mcmc=128L,num_null_draws=99L,
  test_type="mean",use_rao_blackwell=TRUE,standardize_residuals=FALSE,
  general_params=list(verbose=FALSE,num_chains=2L,random_seed=26100103L,standardize=TRUE,
                      standardize_cov=FALSE,interaction_rule="continuous",propensity_covariate="mu"))
report<-get("bcf_score_diagnostics",functions,inherits=TRUE)(fit,y,z,.5,H,batch_size=16L)
stopifnot(report$package_quad_reconstruction_error<=1e-8,
          report$pooled$covariance_rank==ncol(H),is.finite(report$pooled$p))
out<-file.path(here,"local-output");dir.create(out,showWarnings=FALSE,recursive=TRUE)
summary<-list(status="PASS",purpose="Fresh functionality smoke; short MCMC is not a power or calibration result.",
  package_version=as.character(packageVersion("stochtree")),
  runtime_mode=if(overlay)"Current R source overlay on explicitly selected installed binary" else "Installed package namespace",
  p=report$pooled$p,
  reconstruction_error=report$package_quad_reconstruction_error,rank=report$pooled$covariance_rank,
  p_mcse=report$mc$p_mcse,mcse_scope="Conditional on stationarity; diagnostic warnings remain relevant.")
writeLines(capture.output(print(summary)),file.path(out,"diagnostics_smoke.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"diagnostics_session.txt"))
print(summary)
