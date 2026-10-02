# Rscript test/R/run-aipw-score-tests.R [optional isolated package-library path]
script<-sub('^--file=','',commandArgs(FALSE)[startsWith(commandArgs(FALSE),'--file=')][1])
repository<-normalizePath(file.path(dirname(script),'../..'),winslash='/',mustWork=TRUE)
args<-commandArgs(trailingOnly=TRUE)
if(length(args)) .libPaths(c(normalizePath(args[1],winslash='/',mustWork=TRUE),.libPaths()))
loadNamespace('stochtree')
cat('Compiled runtime:',find.package('stochtree'),'\n')
test_runtime<-new.env(parent=asNamespace('stochtree'))
sys.source(file.path(repository,'R/bart.R'),test_runtime)
sys.source(file.path(repository,'R/bcf_aipw_score_test.R'),test_runtime)
bcf_aipw_score_test<-test_runtime$bcf_aipw_score_test
bart<-test_runtime$bart
results<-testthat::test_file(file.path(repository,'test/R/testthat/test-aipw-score.R'),reporter='summary',stop_on_failure=TRUE)
summary<-as.data.frame(results)
cat('Assertions:',sum(summary$nb),'; failed:',sum(summary$failed),'; warnings:',sum(summary$warning),'\n')