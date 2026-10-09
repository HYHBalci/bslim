# Standalone regression check of the exact global-prior/slice helpers in live C++.
# Uses Rcpp to compile extracted helper text; does not install or mutate spbart.
args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=',args)])
if(length(script)!=1L)stop('Run this file using Rscript.')
repo <- normalizePath(file.path(dirname(script),'../..'),mustWork=TRUE)
root <- tempfile('hchs-global-prior-'); dir.create(root)
source_file <- file.path(repo,'src/horseshoe_samplers_old.cpp')
Sys.setenv(PATH=paste('C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/bin','C:/RBuildTools/4.4/usr/bin',Sys.getenv('PATH'),sep=.Platform$path.sep))
s<-readLines(source_file);start<-grep('^double logPosteriorTauGlob_old',s);end<-grep('^// SECTION 3:',s)-2L;stopifnot(length(start)==1L,length(end)==1L)
utility<-s[grep('^inline double safe_log',s):(grep('^double rinvgamma',s)-1L)]
wrapper<-c('#include <Rcpp.h>','#include <vector>','#include <cmath>','#include <limits>',utility,s[start:end],
'// [[Rcpp::export]]',
'double audited_log_global(double t, std::string prior="hc-hs", double A=1.0) {',
' std::vector<double>b{.3,-.2},g{.7},l{.9,1.1,.8};std::vector<std::pair<int,int>>pairs{{0,1}};',
' return logPosteriorTauGlob_old(t,b,g,l,pairs,1.0,1.3,true,prior,A,0);}',
'// [[Rcpp::export]]',
'Rcpp::NumericVector audited_slice_draws(int n, int burn, double initial, double width=.5, int steps=50, double A=1.0) {',
' std::vector<double>b{.3,-.2},g{.7},l{.9,1.1,.8};std::vector<std::pair<int,int>>pairs{{0,1}};Rcpp::NumericVector out(n);double t=initial;',
' for(int i=0;i<n+burn;i++){t=sample_tau_global_slice_old(t,b,g,l,pairs,1.0,1.3,true,width,steps,"hc-hs",A,0);if(i>=burn)out[i-burn]=t;}return out;}',
'// [[Rcpp::export]]',
'Rcpp::NumericVector audited_aux_draws(int n, int burn, double initial, double A=1.0) {Rcpp::NumericVector out(n);double t=initial;for(int i=0;i<n+burn;i++){double xi=1/R::rgamma(1,1/(1/(A*A)+1/(t*t)));double v=1/R::rgamma(1,1/(1/xi+.7*.7/(2*1.3*1.3*.8*.8)));t=std::sqrt(v);if(i>=burn)out[i-burn]=t;}return out;}')
writeLines(wrapper,file.path(root,'extracted_slice_helpers.cpp'))
Rcpp::sourceCpp(file.path(root,'extracted_slice_helpers.cpp'),cacheDir=file.path(root,'cache'),rebuild=TRUE,showOutput=FALSE)
analytic<-function(t,prior,A){q<-if(prior%in%c('hc-hs','hybrid'))1 else 3;S<-.7^2/.8^2;if(q==3)S<-S+.3^2/.9^2+.2^2/1.1^2;lp<-if(prior%in%c('hc-hs','half-cauchy'))-log(A*A+t*t)else if(prior=='half-normal')-t*t/(2*A*A)else 0;-q*log(t)-S/(2*1.3^2*t*t)+lp}
rows<-list();ii<-0L
for(prior in c('hc-hs','half-cauchy','half-normal','hybrid'))for(A in c(.5,1,2.5)){grid<-exp(seq(log(.02),log(100),length.out=101));actual<-vapply(grid,audited_log_global,numeric(1),prior=prior,A=A);ref<-analytic(grid,prior,A);err<-max(abs((actual-actual[51])-(ref-ref[51])));ii<-ii+1L;rows[[ii]]<-data.frame(prior,A,max_constant_free_error=err,pass=err<1e-10)}
analytic_rows<-do.call(rbind,rows);stopifnot(all(analytic_rows$pass));write.csv(analytic_rows,file.path(root,'analytic_log_density.csv'),row.names=FALSE)
set.seed(88147);check<-list();ii<-0L
for(A in c(1,2.5)){
 density<-function(t)ifelse(t>0,exp(analytic(t,'hc-hs',A)),0)
 norm<-integrate(density,0,Inf,rel.tol=1e-10,subdivisions=1000)$value
 cdf<-function(t)integrate(density,0,t,rel.tol=1e-10,subdivisions=1000)$value/norm
 probs<-c(.1,.25,.5,.75,.9);quant<-vapply(probs,function(p)uniroot(function(t)cdf(t)-p,c(.001,100),tol=1e-10)$root,numeric(1))
 for(method in c('slice','auxiliary')){
  x<-if(method=='slice')audited_slice_draws(200000L,10000L,1,A=A)else audited_aux_draws(200000L,10000L,1,A=A)
  for(j in seq_along(probs)){z<-as.numeric(x<=quant[j]);batches<-colMeans(matrix(z,nrow=1000L));se<-sd(batches)/sqrt(length(batches));ii<-ii+1L;check[[ii]]<-data.frame(A,method,p=probs[j],quantile=quant[j],empirical=mean(z),batch_mcse=se,z=(mean(z)-probs[j])/se,pass=abs(mean(z)-probs[j])<max(.004,5*se))}
  saveRDS(list(A=A,method=method,draws=x),file.path(root,paste0(method,'_A',A,'.rds')))
 }
}
validation<-do.call(rbind,check);write.csv(validation,file.path(root,'conditional_distribution.csv'),row.names=FALSE)
write.csv(data.frame(file=source_file,sha256=digest::digest(file=source_file,algo='sha256')),file.path(root,'source_identity.csv'),row.names=FALSE)
print(analytic_rows);print(validation);cat('ALL CHECKS PASS:',all(validation$pass),'\n')

stopifnot(all(validation$pass))
