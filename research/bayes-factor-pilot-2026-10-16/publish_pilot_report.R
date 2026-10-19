args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1] else "D:/heterogeneity-bayes-research-2026-10-16/pilot-results"
report_path<-file.path(dirname(out),"reports","04_pilot_results.md")
body<-readLines(file.path(out,"PILOT_FINDINGS.md"),warn=FALSE)
tab<-read.csv(file.path(out,"summary_full.csv"),stringsAsFactors=FALSE)
null_case<-subset(tab,stage=="forest" & design=="confounded_known" & scenario=="null" &
                  basis=="linear" & method=="null_posterior_integrated_bf" & prior=="sparse" &
                  scale==.5 & event=="BF10_gt_10")
stopifnot(nrow(null_case)==1L)
intro<-c("# 4. Bayes-factor pilot in the existing semiparametric forest model","",
 "**The algebraic route is promising; the simplest null-posterior averaging estimator is not ready to use as a dependable default.** The existing restricted forest fit can supply the nuisance draws for a properly derived evidence identity. In the simulations below, however, many alternative-case estimates depend strongly on a small number of draws. The sparse prior, total prior scale, and chosen modifier basis also materially affect the answer.","",
 "All new code is isolated in `C:/Users/Yusuf Balci/Documents/stochtree/research/bayes-factor-pilot-2026-10-16`. Existing package functions, exports, C++ code, and other functionality were left unchanged. The output artifacts are in `../pilot-results/`.","",
 "## What was actually implemented","",
 "The working model is Y=f(X,e)+(Z-e){alpha+H(X)delta}+error. The forest is fitted only under delta=0; conditional Gaussian integration then gives alternative-to-null likelihood ratios that are averaged over its null posterior. Alpha is integrated analytically under the same N(0,4) prior in both models. The noise variance is known and fixed at 1, and the 25-tree forest has fixed Gaussian leaf variance .04 per tree. Outcome standardization and both variance updates are disabled. The selected modifier coefficients have either a dense Gaussian prior or a proper nonempty sparse subset mixture. This is a reduced-model experiment, not a marginal-likelihood implementation for the original horseshoe model.","",
 "The original sparse extension has tau(X)=alpha+X beta+X_products beta_int. Both beta and beta_int are heterogeneity coefficients; setting only beta_int to zero would leave linear heterogeneity. The centered treatment coding here follows the restricted score model. Changing between mu(X)+Z tau(X) and m(X)+(Z-e)tau(X) also changes the nuisance function by e tau, so independent priors cannot simply be carried over and declared equivalent. The full project must preserve compatible priors under the tested models.","",
 "The actual run reuses the installed native stochtree library with isolated current-source R wrappers. The evidence identity is exact for the stipulated compatible models, while its numerical estimate is conditional on the sampler approximating their null posterior. Current source and installed DLL are separately fingerprinted; their hashes do not establish a matching native build.","",
 "The [methods and reproduction notes](../pilot-results/METHODS_AND_REPRODUCTION.md) give the formula, priors, generator, commands, diagnostic definitions, and limitations. Thirteen algebraic checks passed, including independent quadrature against a fully integrated Gaussian nuisance benchmark. A smoke run exercised every code path with 544 result rows and zero captured model errors or warnings.","",
 "![Exploratory decision frequencies with Wilson intervals](../pilot-results/pilot_evidence_rates_print.png)","",
 "[Original full-size landscape figure](../pilot-results/pilot_evidence_rates.png).", "",
 "The plot uses unequal BF and p-value decision thresholds and should not be read as a matched-size power comparison.","")
intro<-c(intro,sprintf("One caution is visible even under the generated constant-CATE null: in the confounded-known-propensity case, sparse BF10 exceeded 10 in %d/%d datasets (%.1f%%; Wilson interval %.1f–%.1f%%). This is a fixed-data-generating-process false-evidence frequency under the chosen forest and known-noise assumptions. BF10>10 was not calibrated as a 5%% frequentist test, and this result does not contradict a prior-predictive BF identity. It does rule out describing this pilot as evidence of uniformly reliable null behavior.",
  null_case$k,null_case$n,100*null_case$rate,100*null_case$ci_low,100*null_case$ci_high),"")
lines<-c(intro,body[-1])
longfile<-file.path(out,"long_chain_comparison.csv")
if(file.exists(longfile)) {
  z<-read.csv(longfile,stringsAsFactors=FALSE)
  lines<-c(lines,"","## Diagnostic-selected longer-chain follow-up","",
    "After the four relevant null/sparse cells had completed, the dataset with the largest dense-prior between-chain log-BF disagreement was selected separately in each null/sparse-by-design cell. Each of these four datasets was refitted with four chains under a fresh call seed, 1,000 burn-in and 3,000 retained iterations per chain. The original runs used two chains with 750 retained draws each. This deliberately selects problematic cases; it is not a random extra performance sample. Changes may reflect chain length, fresh random streams, and regression to the mean. Even agreement would not establish general numerical stability.","",
    "The table reports log(BF10), so 0 means equal evidence and 2.303 corresponds to BF10=10. Weight ESS is importance concentration, not chain ESS.","",
    "| Design | Truth | Replicate | Prior | Original log BF10 | Longer log BF10 | Original weight ESS | Longer weight ESS | Longer chain log-BF range |",
    "|---|---|---:|---|---:|---:|---:|---:|---:|")
  for(i in seq_len(nrow(z))) {
    row<-c(z$design[i],z$scenario[i],z$rep[i],z$prior[i],sprintf("%.2f",z$baseline_log_bf[i]),
       sprintf("%.2f",z$long_log_bf[i]),round(z$baseline_weight_ess[i]),round(z$long_weight_ess[i]),
       sprintf("%.2f",z$long_chain_range[i]))
    lines<-c(lines,paste0("| ",paste(row,collapse=" | ")," |"))
  }
  lines<-c(lines,"",sprintf("The longer runs attempted and completed four selected fits, producing %d prior-specific evidence comparisons. %d/%d longer estimates still triggered at least one concentration/disagreement diagnostic. The largest absolute change in log evidence was %.2f; this is descriptive of the selected cases only.",
      nrow(z),sum(z$long_flag!="none"),nrow(z),max(abs(z$long_log_bf-z$baseline_log_bf))),"",
    "Selection records, all four chain-specific evidence values, and long-chain parameter/ratio traces are saved in `diagnostic_selected_cases.csv`, `long_chain_comparison.csv`, and `long_chain_trace_*.csv`. The follow-up has no separate operating-characteristic interpretation. The main-run frequencies above are deliberately left unchanged: selectively replacing only these problematic datasets and recomputing a headline rate would mix protocols and introduce selection bias.")
}
lines<-c(lines,"","## Files to retain","",
 "* `results_full.csv`: every primary per-replication endpoint, seed, prior, and numerical diagnostic.",
 "* `summary_full.csv` and `pilot_summary.json`: event frequencies, uncertainty intervals, and run metadata.",
 "* `attempted_and_successful.csv`, `warnings_full.csv`, and `fit_diagnostics_full.csv`: execution counts, captured model conditions, and posterior checks.",
 "* `config_full.R`, `long_chain_config.R`, `session_info.txt`, and run logs: exact budgets and runtime environment.",
 "* `representative_full_*.rds`, `ratios_full_*.csv`, and long-chain traces: selected saved fits and evidence diagnostics.",
 "* `mathematical_test_log.txt`, `provenance.json`, `source_and_binary_sha256.csv`, and `artifact_sha256.csv`: algebra validation and provenance.","",
 "The next substantive step is to retain this analytic integration as a benchmark while developing a more stable full-model evidence estimator, then incorporate unknown noise and the actual local/global shrinkage hierarchy. In parallel, use exact energy-matched prior sensitivity and larger matched-error-rate comparisons. The current pilot supports that research direction; it does not establish that the present estimator is calibrated, stable, or superior to the existing heterogeneity test.")
dir.create(dirname(report_path),recursive=TRUE,showWarnings=FALSE)
writeLines(lines,report_path)
cat("Wrote",report_path,"\n")
