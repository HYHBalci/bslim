# Rscript test/R/run-score-diagnostics-tests.R [optional R library directory]
script<-sub('^--file=','',commandArgs(FALSE)[startsWith(commandArgs(FALSE),'--file=')][1])
repository<-normalizePath(file.path(dirname(script),'../..'),winslash='/',mustWork=TRUE)
args<-commandArgs(TRUE);if(length(args)) .libPaths(c(args[1],.libPaths()))
source(file.path(repository,'R/bcf_score_diagnostics.R'))
testthat::test_file(file.path(repository,'test/R/testthat/test-score-diagnostics.R'),reporter='summary',stop_on_failure=TRUE)