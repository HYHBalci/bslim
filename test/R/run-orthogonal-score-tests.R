# Run from any working directory with Rscript test/R/run-orthogonal-score-tests.R.
script<-sub('^--file=','',commandArgs(FALSE)[startsWith(commandArgs(FALSE),'--file=')][1])
repository<-normalizePath(file.path(dirname(script),'../..'),winslash='/',mustWork=TRUE)
loadNamespace('spbart')
test_runtime<-new.env(parent=asNamespace('spbart'))
sys.source(file.path(repository,'R/bart.R'),test_runtime)
sys.source(file.path(repository,'R/bcf_orthogonal_score_test.R'),test_runtime)
bart<-test_runtime$bart
bcf_orthogonal_score_test<-test_runtime$bcf_orthogonal_score_test
# Registered predict method lives in the installed runtime, with unchanged format.
testthat::test_file(file.path(repository,'test/R/testthat/test-orthogonal-score.R'),reporter='summary',stop_on_failure=TRUE)
testthat::test_file(file.path(repository,'test/R/testthat/test-bart-probit-offset.R'),reporter='summary',stop_on_failure=TRUE)
