args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1] else "D:/heterogeneity-bayes-research-2026-10-16/pilot-results"
source("research/bayes-factor-pilot-2026-10-16/bayes_factor_sidecar.R")
res<-read.csv(file.path(out,"results_full.csv"),stringsAsFactors=FALSE)
sumtab<-read.csv(file.path(out,"summary_full.csv"),stringsAsFactors=FALSE)
fit<-read.csv(file.path(out,"fit_diagnostics_full.csv"),stringsAsFactors=FALSE)
warn<-read.csv(file.path(out,"warnings_full.csv"),stringsAsFactors=FALSE)
config<-dget(file.path(out,"config_full.R"))
counts<-data.frame(stage=c("oracle_datasets","oracle_basis_score_calls","forest_fits","crossfit_fits"),
  attempted=c(8*config$oracle_reps,10*config$oracle_reps,8*config$forest_reps,4*config$crossfit_reps),
  successful=c(length(unique(res$seed[res$stage=="oracle"])),
    sum(res$stage=="oracle" & res$method=="orthogonal_score"),nrow(fit),
    sum(res$stage=="crossfit" & res$method=="orthogonal_score")))
counts$unsuccessful<-counts$attempted-counts$successful
write.csv(counts,file.path(out,"attempted_and_successful.csv"),row.names=FALSE)
fmt<-function(x) sprintf("%.1f",100*x)
rate<-function(stage,scenario,design,method,prior="none",scale=NA_real_,basis="linear",event="p_lt_0.05") {
  z<-sumtab[sumtab$stage==stage & sumtab$scenario==scenario & sumtab$design==design &
      sumtab$method==method & sumtab$prior==prior & sumtab$basis==basis & sumtab$event==event &
      if(is.na(scale))is.na(sumtab$scale) else !is.na(sumtab$scale)&sumtab$scale==scale,]
  if(nrow(z)!=1L)return("unavailable")
  paste0(z$k,"/",z$n,"; ",fmt(z$rate),"% [",fmt(z$ci_low),", ",fmt(z$ci_high),"]")
}
lines<-c("# Pilot findings: actual runs, 16 October 2026","",
  "These are exploratory simulation results, not a validation study. All intervals below are 95% Wilson binomial intervals. The BF and score thresholds are unequal decision rules; their frequencies must not be interpreted as a matched-size power ranking.","",
  "The forest BF estimates integrate the constant effect analytically and average conditional evidence over the existing restricted forest's null posterior. The target is the fixed-prior Gaussian subset model in METHODS_AND_REPRODUCTION.md, not the original horseshoe-integrated model. Weight concentration and finite-chain error remain material.","",
  "## Primary forest-linked experiment (n=300, slab SD=0.5)","",
  "| Design | Truth | Restricted score p<.05 | Dense BF10>10 | Sparse BF10>10 | Sparse BF01>3 |",
  "|---|---|---|---|---|---|")
for(design in c("randomized","confounded_known"))for(sc in c("null","sparse","dense","quadratic")) {
  row<-c(design,sc,rate("forest",sc,design,"restricted_external_score"),
    rate("forest",sc,design,"null_posterior_integrated_bf","dense",.5,event="BF10_gt_10"),
    rate("forest",sc,design,"null_posterior_integrated_bf","sparse",.5,event="BF10_gt_10"),
    rate("forest",sc,design,"null_posterior_integrated_bf","sparse",.5,event="BF01_gt_3"))
  lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
}
lines<-c(lines,"","The last column is evidence for the point null **within the chosen finite modifier model**. Under the quadratic truth, the linear basis can miss heterogeneity; apparent support for its zero coefficients does not establish a globally constant CATE.","",
  "## Basis sensitivity for the quadratic alternative","",
  "| Design | Basis | Restricted score p<.05 | Sparse BF10>10 | Dense BF10>10 |",
  "|---|---|---|---|---|")
for(design in c("randomized","confounded_known"))for(basis in c("linear","augmented")) {
  row<-c(design,basis,rate("forest","quadratic",design,"restricted_external_score",basis=basis),
    rate("forest","quadratic",design,"null_posterior_integrated_bf","sparse",.5,basis,event="BF10_gt_10"),
    rate("forest","quadratic",design,"null_posterior_integrated_bf","dense",.5,basis,event="BF10_gt_10"))
  lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
}
lines<-c(lines,"","## Prior sensitivity under randomization","",
  "Equal slab standard deviations do not match total prior signal energy. With q=6 and SD=.5, E||delta||² is 1.5 for the dense alternative and approximately .37588 for the sparse mixture. Dense SD=.25 has energy .375 and is almost, but not exactly, matched to sparse SD=.5. The comparisons below therefore reflect prior-scale and sparsity choices jointly.","",
  "| Truth | Slab SD | Dense BF10>10 | Sparse BF10>10 | Dense BF01>3 | Sparse BF01>3 |",
  "|---|---|---|---|---|---|")
for(sc in c("null","sparse","dense"))for(s in c(.25,.5,1)) {
  row<-c(sc,s,rate("forest",sc,"randomized","null_posterior_integrated_bf","dense",s,event="BF10_gt_10"),
    rate("forest",sc,"randomized","null_posterior_integrated_bf","sparse",s,event="BF10_gt_10"),
    rate("forest",sc,"randomized","null_posterior_integrated_bf","dense",s,event="BF01_gt_3"),
    rate("forest",sc,"randomized","null_posterior_integrated_bf","sparse",s,event="BF01_gt_3"))
  lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
}
lines<-c(lines,"","## Approximate score-experiment evidence sensitivity","",
  "This calculation recycles the restricted external quadratic score. It is exact only for the stipulated Gaussian score experiment and is not a forest marginal likelihood. Changing g changes the alternative's prior scale; it does not merely relabel a p-value.","",
  "| Truth under randomization | g | Approximate BF10>10 | Approximate BF01>3 |",
  "|---|---|---|---|")
for(sc in c("null","sparse","dense"))for(g in c(1,10,300)) {
  row<-c(sc,g,rate("forest",sc,"randomized","approximate_score_bf","normal_score",g,event="BF10_gt_10"),
    rate("forest",sc,"randomized","approximate_score_bf","normal_score",g,event="BF01_gt_3"))
  lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
}
lines<-c(lines,"","## Oracle benchmark (300 replicates per cell)","",
  "| Design | Truth | Orthogonal score p<.05 | Dense BF10>10 | Sparse BF10>10 |",
  "|---|---|---|---|---|")
for(design in c("randomized","confounded_known"))for(sc in c("null","sparse","dense","quadratic")) {
  row<-c(design,sc,rate("oracle",sc,design,"orthogonal_score"),
    rate("oracle",sc,design,"exact_conditional_gaussian_bf","dense",.5,event="BF10_gt_10"),
    rate("oracle",sc,design,"exact_conditional_gaussian_bf","sparse",.5,event="BF10_gt_10"))
  lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
}
lines<-c(lines,"","## Cross-fitted nuisance sensitivity (25 replicates per cell)","",
  "These residual-likelihood values are plug-in diagnostics. Neither sample splitting nor a Gaussian formula makes them full-model Bayes factors.","",
  "| Propensity | Truth | Orthogonal score p<.05 | Sparse plug-in BF10>10 |",
  "|---|---|---|---|")
for(design in c("confounded_known","confounded_learned"))for(sc in c("null","sparse")) {
  row<-c(design,sc,rate("crossfit",sc,design,"orthogonal_score"),
    rate("crossfit",sc,design,"crossfit_plugin_bf","sparse",.5,event="BF10_gt_10"))
  lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
}
lines<-c(lines,"","## Evidence integration diagnostics","",
  "Weight ESS below measures concentration of importance ratios; it is **not** MCMC ESS. A seemingly acceptable value cannot reveal important tail regions never visited by either null chain. Flags mean ESS<100, a single weight>10%, or chain BFs differing by more than a factor of 3.","",
  "| Design | Truth | Prior | Median weight ESS / 1500 | Flagged estimates | Median chain log-BF range |",
  "|---|---|---|---|---|---|")
z<-subset(res,method=="null_posterior_integrated_bf" & scale==.5 & basis=="linear")
for(design in c("randomized","confounded_known"))for(sc in c("null","sparse","dense","quadratic"))for(pr in c("dense","sparse")) {
  a<-z[z$design==design & z$scenario==sc & z$prior==pr,]
  row<-c(design,sc,pr,round(median(a$importance_ess)),
      paste0(sum(a$diagnostic_warning!="none"),"/",nrow(a)),sprintf("%.2f",median(a$chain_log_bf_range)))
  lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
}
lines<-c(lines,"",sprintf("Across %d forest fits, median basic split-Rhat was %.3f for alpha, %.3f for the average forest prediction, and %.3f for the forest projection onto treatment. The largest values were %.3f, %.3f, and %.3f, respectively. These are basic split diagnostics, not rank-normalized Rhat.",
  nrow(fit),median(fit$alpha_split_rhat),median(fit$mu_mean_split_rhat),median(fit$mu_treatment_projection_split_rhat),
  max(fit$alpha_split_rhat),max(fit$mu_mean_split_rhat),max(fit$mu_treatment_projection_split_rhat)),"",
  sprintf("Median forest fitting time was %.2f seconds (range %.2f–%.2f). Post-processing runtime on BF rows covers both priors and both plug-in comparisons together and is repeated on the two rows, so those timings must not be summed.",median(fit$fit_seconds),min(fit$fit_seconds),max(fit$fit_seconds)),"",
  sprintf("The run returned %d result rows and %d captured model warnings/errors. Four R locale warnings occurred before the run's logging sink; they are documented in the methods file and test log.",nrow(res),nrow(warn)),"",
  "## Execution counts","",
  "| Stage | Attempted | Successful | Unsuccessful |","|---|---:|---:|---:|")
for(i in seq_len(nrow(counts)))lines<-c(lines,paste0("| ",paste(counts[i,],collapse=" | ")," |"))
lines<-c(lines,"",paste0("Captured model errors: ",sum(grepl("^ERROR",warn$message)),
  "; captured model warnings: ",sum(!grepl("^ERROR",warn$message)),". Oracle datasets are counted separately from basis-score evaluations because the quadratic scenario also uses the augmented basis."),"",
  "## Interpretation and next experiment","",
  "The reduced-model identity supplies a practical route from one restricted forest fit to several proper alternative priors. It also exposes the central computational risk: alternative evidence may depend on rare null-posterior draws. Plug-in evidence should not replace that integration merely because it is numerically smoother.","",
  "Prior sensitivity, uncertainty in the nuisance model, and the effect-modifier basis are separate issues. A successful coefficient test cannot rule out heterogeneity outside its basis. Before treating this as a method for publication, use matched-size operating-characteristic comparisons, substantially more replications, long independent chains, a native-target audit, unknown-noise and shrinkage-scale integration, and an alternative evidence estimator such as bridge sampling for difficult overlap cases.","",
  "The tests in mathematical_test_log.txt validate the Gaussian algebra and the integration identity in a jointly integrated Gaussian nuisance special case. They do not validate the forest sampler's invariant distribution, prove finite-chain convergence, or validate Bayes-factor use with estimated propensities.")
writeLines(lines,file.path(out,"PILOT_FINDINGS.md"))
jsonlite::write_json(list(scope="Exploratory reduced-model Gaussian-slab BF pilot; not horseshoe-integrated evidence",
  unequal_thresholds=TRUE,noise_variance_known=TRUE,full_result_rows=nrow(res),forest_fits=nrow(fit),
  captured_model_errors=sum(grepl("^ERROR",warn$message)),captured_model_warnings=sum(!grepl("^ERROR",warn$message)),
  execution_counts=counts,configuration=config,
  primary_summary=sumtab[sumtab$basis=="linear" & sumtab$event %in% c("p_lt_0.05","BF10_gt_10","BF01_gt_3") &
    (is.na(sumtab$scale) | sumtab$scale==.5),],
  basis_summary=sumtab[sumtab$scenario=="quadratic" & sumtab$event %in% c("p_lt_0.05","BF10_gt_10") &
    (is.na(sumtab$scale) | sumtab$scale==.5),]),
  file.path(out,"pilot_summary.json"),pretty=TRUE,auto_unbox=TRUE,na="null",digits=12)
cat("Wrote PILOT_FINDINGS.md\n")
