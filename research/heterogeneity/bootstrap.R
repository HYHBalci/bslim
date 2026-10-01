# Rscript research/heterogeneity/bootstrap.R [--install] [--check]
# Check by default. --install explicitly builds THIS source, never the CRAN stochtree.
args <- commandArgs(TRUE)
unknown <- setdiff(args,c("--install","--check"))
if(length(unknown)) stop("Unknown option(s): ",paste(unknown,collapse=", "))
file_arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
if(length(file_arg)!=1L) stop("Run with Rscript.")
here <- dirname(normalizePath(sub("^--file=","",file_arg),winslash="/",mustWork=TRUE))
repo <- normalizePath(file.path(here,"../.."),winslash="/",mustWork=TRUE)
desc <- read.dcf(file.path(repo,"DESCRIPTION"))
stopifnot(desc[1,"Package"]=="stochtree")
lib <- Sys.getenv("STOCHTREE_R_LIBRARY",file.path(repo,".research-library",paste0("R-",getRversion(),"-",R.version$platform)))
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
lib <- normalizePath(lib,winslash="/",mustWork=TRUE)
.libPaths(c(lib,.libPaths()))
fields <- intersect(c("Imports","LinkingTo"),colnames(desc))
deps <- unique(trimws(sub("\\s*\\(.*$","",unlist(strsplit(paste(desc[1,fields],collapse=","),",")))))
# posterior supports new diagnostics; the other packages reproduce the pilot comparisons.
deps <- setdiff(unique(c(deps,"posterior","digest","grf","sandwich","lmtest","testthat")),c("R","stats"))
header_paths <- c("deps/eigen/Eigen/Dense","deps/fmt/include/fmt/format.h","deps/fast_double_parser/include/fast_double_parser.h")
missing_headers <- header_paths[!file.exists(file.path(repo,header_paths))]
if("--install" %in% args) {
  if(length(missing_headers)) stop("Missing source submodules: ",paste(missing_headers,collapse=", "),". Initialize the checkout submodules first.")
  needed <- deps[!vapply(deps,requireNamespace,logical(1),quietly=TRUE)]
  if(length(needed)) install.packages(needed,lib=lib,repos=Sys.getenv("R_REPOSITORY","https://cloud.r-project.org"))
  install.packages(repo,lib=lib,repos=NULL,type="source",INSTALL_opts="--preclean")
}
packages <- c("stochtree",deps)
available <- vapply(packages,requireNamespace,logical(1),quietly=TRUE)
versions <- vapply(packages,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p)) else NA_character_,character(1))
report <- data.frame(package=packages,available=available,version=versions,row.names=NULL)
out <- file.path(repo,"research","heterogeneity","local-output")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.csv(report,file.path(out,"runtime_packages.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"runtime_session.txt"))
cat("Source:",repo,"\nLibrary:",lib,"\n");print(report,row.names=FALSE)
if(length(missing_headers)) cat("Missing source headers:",paste(missing_headers,collapse=", "),"\n")
helper <- available[1] && "bcf_score_diagnostics" %in% getNamespaceExports("stochtree")
cat("Installed package exports bcf_score_diagnostics:",helper,"\n")
if(any(!available)||length(missing_headers)||!helper) {
  cat("Runtime is not ready for the current source-package demo. See README.md.\n")
  quit(status=2L)
}
cat("PASS: package/dependency/source-header checks. This does not certify historical DLL identity.\n")
