.bcf_score_diag_controls <- function(alpha=.05,rank_tol=1e-10,ridge=0) {
 if(length(alpha)!=1L||!is.finite(alpha)||alpha<=0||alpha>=1)stop('alpha must be in (0,1).')
 if(length(rank_tol)!=1L||!is.finite(rank_tol)||rank_tol<=0||rank_tol>=1)stop('rank_tol must be in (0,1).')
 if(length(ridge)!=1L||!is.finite(ridge)||ridge<0)stop('ridge must be finite and nonnegative.')
 invisible(TRUE)
}
# Standalone same-sample restricted-forest score diagnostics, 2026-09-30.
# No calibration or stationarity guarantee is established by this helper.
.bcf_score_diag_basis <- function(H,n) {
 H<-as.matrix(H)
 if(!is.numeric(H)||any(!is.finite(H))||nrow(H)!=n||ncol(H)<1L)stop('H must be a finite n by q numeric matrix.')
 H<-sweep(H,2L,colMeans(H),'-')
 names<-colnames(H);if(is.null(names))names<-paste0('H',seq_len(ncol(H)));empty<-is.na(names)|!nzchar(names);names[empty]<-paste0('H',which(empty));colnames(H)<-make.unique(names)
 H
}
.bcf_score_diag_quad <- function(H,s,rank_tol=1e-10,ridge=0) {
 n<-nrow(H);q<-ncol(H);.bcf_score_diag_controls(rank_tol=rank_tol,ridge=ridge);if(!is.numeric(s)||length(s)!=n||any(!is.finite(s)))stop('s must be finite numeric of length n.');s<-as.numeric(s)-mean(s);U<-drop(crossprod(H,s));V<-crossprod(H*s)
 ev<-eigen(V,symmetric=TRUE,only.values=TRUE)$values;vmax<-max(ev);vmin<-min(ev)
 rank<-if(vmax>0)sum(ev>rank_tol*vmax)else 0L;condition<-if(vmin>0)vmax/vmin else Inf
 usable<-rank==q&&vmax>0&&is.finite(condition);statistic<-p<-NA_real_
 if(usable){statistic<-drop(crossprod(U,solve(V+diag(ridge,q),U)));p<-pchisq(statistic,q,lower.tail=FALSE)}
 legacy_statistic<-legacy_p<-NA_real_
 if(ridge>0){legacy_statistic<-tryCatch(drop(crossprod(U,solve(V+diag(ridge,q),U))),error=function(e)NA_real_);legacy_p<-pchisq(legacy_statistic,q,lower.tail=FALSE)}
 list(U=U,V=V,statistic=statistic,p=p,q=q,n=n,score=s,covariance_rank=rank,covariance_condition=condition,eigenvalues=ev,covariance_usable=usable,legacy_ridge_statistic=legacy_statistic,legacy_ridge_p=legacy_p)
}
.bcf_score_diag_batch_mc_cov <- function(draws,batch_size=NULL) {
 if(length(dim(draws))!=3L||any(!is.finite(draws)))stop('draws must be a finite iteration by chain by variable array.')
 M<-dim(draws)[1L];C<-dim(draws)[2L];D<-dim(draws)[3L]
 if(is.null(batch_size))batch_size<-max(2L,floor(sqrt(M)))
 if(length(batch_size)!=1L||!is.finite(batch_size)||batch_size!=floor(batch_size)||batch_size<2L)stop('batch_size must be an integer >=2.')
 B<-as.integer(batch_size);L<-floor(M/B);if(L<8L)stop('At least eight complete batches per chain are required.')
 covs<-vector('list',C);bm<-array(NA_real_,c(L,C,D))
 for(k in seq_len(C)){
  xx<-matrix(draws[,k,],M,D)
  for(b in seq_len(L))bm[b,k,]<-colMeans(xx[(b-1L)*B+seq_len(B),,drop=FALSE])
  covs[[k]]<-B*stats::cov(matrix(bm[,k,],L,D))/M
 }
 covariance<-Reduce('+',covs)/C^2;dimnames(covariance)<-rep(list(dimnames(draws)[[3L]]),2L)
 list(covariance=covariance,mcse=sqrt(pmax(0,diag(covariance))),batch_size=B,batches_per_chain=L,unused_tail_per_chain=M-L*B,chain_covariance=covs,batch_means=bm,
 interpretation='Approximate within-chain Monte Carlo error conditional on the data; assumes stationary independent chains and adequate batch length. It excludes burn-in bias and between-chain nonstationarity.')
}
.bcf_score_diag_posterior_diagnostics <- function(draws) {
 if(!requireNamespace('posterior',quietly=TRUE))return(list(available=FALSE,reason='posterior package is unavailable.'))
 tab<-as.data.frame(posterior::summarise_draws(posterior::as_draws_array(draws),'rhat','ess_bulk','ess_tail','mcse_mean'))
 list(available=TRUE,table=tab,max_rhat=if(all(!is.finite(tab$rhat)))NA_real_ else max(tab$rhat,na.rm=TRUE),min_bulk_ess=if(all(!is.finite(tab$ess_bulk)))NA_real_ else min(tab$ess_bulk,na.rm=TRUE))
}
.bcf_score_diag_delta_draws <- function(H,score_draws,M,C,ridge=0,rank_tol=1e-10) {
 n<-nrow(H);q<-ncol(H);sbar<-rowMeans(score_draws);pooled<-.bcf_score_diag_quad(H,sbar,rank_tol,ridge)
 pair<-which(upper.tri(matrix(0,q,q),diag=TRUE),arr.ind=TRUE);D<-q+nrow(pair)
 lin<-matrix(NA_real_,M*C,D);lin[,seq_len(q)]<-t(crossprod(H,score_draws))/sqrt(n)
 change<-sweep(score_draws,1L,sbar,'-')
 for(j in seq_len(nrow(pair)))lin[,q+j]<-drop(crossprod(2*sbar*H[,pair[j,1L]]*H[,pair[j,2L]]/n,change))
 colnames(lin)<-c(paste0('U_',colnames(H)),paste0('V',pair[,1L],'_',pair[,2L]))
 draws<-array(NA_real_,c(M,C,D),dimnames=list(NULL,NULL,colnames(lin)))
 for(k in seq_len(C))draws[,k,]<-lin[(k-1L)*M+seq_len(M),,drop=FALSE]
 grad<-rep(NA_real_,D)
 if(pooled$covariance_usable){a<-drop(solve((pooled$V+diag(ridge,q))/n,pooled$U/sqrt(n)));grad[seq_len(q)]<-2*a;grad[q+seq_len(nrow(pair))]<- -a[pair[,1L]]*a[pair[,2L]]*ifelse(pair[,1L]==pair[,2L],1,2)}
 list(pooled=pooled,draws=draws,gradient=grad,pairs=pair)
}
.bcf_score_diag_delta_mc <- function(delta,mc,alpha=.05) {
 q<-delta$pooled$q;n<-delta$pooled$n;mc$U_mcse<-sqrt(n)*mc$mcse[seq_len(q)];mc$V_linearized_mcse<-n*mc$mcse[q+seq_len(nrow(delta$pairs))]
 mc$statistic_mcse<-mc$p_mcse<-NA_real_
 if(delta$pooled$covariance_usable){mc$statistic_mcse<-sqrt(max(0,drop(crossprod(delta$gradient,mc$covariance%*%delta$gradient))));density<-dchisq(delta$pooled$statistic,q);if(is.finite(density))mc$p_mcse<-density*mc$statistic_mcse}
 mc$delta_gradient<-delta$gradient;mc$includes_covariance_change<-TRUE;mc$centering_correction_available<-TRUE
 mc$near_decision_threshold<-if(is.finite(mc$p_mcse))abs(delta$pooled$p-alpha)<=2*mc$p_mcse else NA
 mc
}
#' Diagnose pooled restricted-BCF mean scores and Monte Carlo precision
#'
#' Reconstructs one score test after pooling retained chain contributions, and
#' provides covariance-rank and score-specific Monte Carlo diagnostics. It does
#' not change the fitted sampler or the legacy per-chain p-values.
#' @param fit A returned posterior-averaged mean-test bcfscoretest object.
#' @param y,z Original numeric outcome and binary treatment vectors.
#' @param propensity The actual propensity marks used by the fit, scalar or length n.
#'   Internally estimated marks are not currently retained by the sampler and must
#'   be saved externally; this argument does not assert a known true propensity.
#' @param H Prespecified finite numeric score matrix with n rows. Columns are
#'   recentered. Supply the exact prepared package score matrix to reconstruct
#'   its per-chain ridge p-values; using another basis defines another test.
#' @param batch_size Contiguous draws per Monte Carlo batch. Defaults to
#'   floor(sqrt(retained draws)); at least eight complete batches are required.
#' @param alpha Decision threshold for descriptive chain/MC proximity flags.
#' @param rank_tol Relative eigenvalue tolerance for unregularized covariance rank.
#' @param ridge Nonnegative fixed ridge in sampler units, default 1e-8 matching
#'   the sampler. A ridge is never accepted as repairing a singular covariance.
#' @param keep_draws Whether to return linearized score/covariance trajectories.
#' @return A list with pooled U,V,quadratic statistic and p-value, per-chain tests,
#'   unregularized covariance rank/condition, batch Monte Carlo covariance and
#'   delta MCSE of the statistic/p-value, and optional posterior Rhat/ESS.
#' @details Returned posterior quantities are already on original outcome units.
#'   Contributions are pooled before computing covariance; averaging p-values
#'   is a different procedure. Per-draw residual standardization is reconstructed
#'   when enabled. The quadratic reference is chi-squared; permutation references
#'   and other test types are not supported. Singular or ill-conditioned declared
#'   covariance returns an NA primary p-value, even when a separate legacy ridge
#'   number is computable. A matching H is required for package agreement.
#'
#'   Delta Monte Carlo error includes changes in both score and covariance.
#'   Batch estimates assume stationary independent chains and adequate batch
#'   length; they exclude burn-in bias and between-chain nonstationarity. The
#'   optional posterior package supplies rank-normalized Rhat and ESS. None of
#'   these numerical diagnostics establishes frequentist calibration, nuisance
#'   convergence, a same-sample inferential guarantee, or observational validity.
#' @export
bcf_score_diagnostics <- function(fit,y,z,propensity,H,batch_size=NULL,alpha=.05,rank_tol=1e-10,ridge=1e-8,keep_draws=FALSE) {
 .bcf_score_diag_controls(alpha,rank_tol,ridge)
 cfg<-fit$test_config
 if(is.null(cfg)||!identical(cfg$test_type,'mean')||!isTRUE(cfg$use_rao_blackwell))stop('A posterior-averaged mean-test bcfscoretest object is required.')
 if(!cfg$null_method%in%c('gaussian','multiplier'))stop('Only the Gaussian quadratic reference is reconstructed.')
 if(!is.logical(keep_draws)||length(keep_draws)!=1L||is.na(keep_draws))stop('keep_draws must be TRUE or FALSE.');n<-length(y);H<-.bcf_score_diag_basis(H,n)
 if(!is.numeric(y)||!is.null(dim(y))||any(!is.finite(y))||!is.numeric(z)||!is.null(dim(z))||length(z)!=n||any(!z%in%c(0,1)))stop('Finite numeric y and binary z of equal length are required.')
 if(length(propensity)==1L)propensity<-rep(propensity,n)
 if(!is.numeric(propensity)||!is.null(dim(propensity))||length(propensity)!=n||any(!is.finite(propensity))||any(propensity<=0|propensity>=1))stop('Actual propensity marks must be finite, length n, and inside (0,1).')
 if(!is.null(cfg$propensity_clip))propensity<-pmin(pmax(propensity,cfg$propensity_clip),1-cfg$propensity_clip)
 a<-as.matrix(fit$alpha_samples);sig<-as.matrix(fit$sigma2_samples);mu<-as.matrix(fit$mu_hat_train);C<-nrow(a);M<-ncol(a)
 if(!is.numeric(a)||!is.numeric(sig)||!is.numeric(mu)||!identical(dim(sig),dim(a))||!identical(dim(mu),c(n,C*M))||any(!is.finite(c(a,sig,mu)))||any(sig<=0))stop('Returned draw shapes or sigma2 values are invalid.')
 W<-z-propensity;resid<-sweep(mu,1L,y,function(mm,yy)yy-mm)-outer(W,as.numeric(t(a)))
 if(isTRUE(cfg$standardize_residuals))resid<-sweep(resid,2L,sqrt(as.numeric(t(sig))),'/') else {
  ys<-fit$y_std_train;if(length(ys)!=1L||!is.finite(ys)||ys<=0)stop('Positive y_std_train is required.');resid<-resid/ys
 }
 sdraw<-resid*W;sdraw<-sweep(sdraw,2L,colMeans(sdraw),'-')
 delta<-.bcf_score_diag_delta_draws(H,sdraw,M,C,ridge,rank_tol);mc<-.bcf_score_diag_delta_mc(delta,.bcf_score_diag_batch_mc_cov(delta$draws,batch_size),alpha)
 chains<-lapply(seq_len(C),function(k).bcf_score_diag_quad(H,rowMeans(sdraw[,(k-1L)*M+seq_len(M),drop=FALSE]),rank_tol,ridge))
 pchain<-vapply(chains,function(x)x$legacy_ridge_p,numeric(1));agreement<-if(length(fit$rb_pval_quad)==C&&cfg$p_valid==ncol(H))max(abs(pchain-fit$rb_pval_quad))else NA_real_
 out<-list(pooled=delta$pooled,per_chain=chains,mc=mc,score_diagnostics=.bcf_score_diag_posterior_diagnostics(delta$draws),package_quad_reconstruction_error=agreement,package_reconstruction_requires_exact_H=TRUE,basis_rank=qr(H)$rank,chain_count=C,draws_per_chain=M,chain_decision_disagreement=length(unique(vapply(chains,function(x)x$p<=alpha,logical(1))))>1L,declaration='Numerical/MCMC diagnostics do not establish calibration, nuisance contraction, or chain stationarity.')
 if(keep_draws)out$linearized_draws<-delta$draws
 out
}
