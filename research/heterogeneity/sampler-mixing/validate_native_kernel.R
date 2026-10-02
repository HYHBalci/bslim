# Independently enumerable fixed-variance posterior versus the actual tree kernel.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
repo <- normalizePath(args[1], winslash='/')
baseline_lib <- normalizePath(args[2], winslash='/')
out <- normalizePath(args[3], winslash='/', mustWork=FALSE)
dir.create(out, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(baseline_lib, .libPaths()))
library(stochtree); library(digest)
dll <- getLoadedDLLs()[['stochtree']][['path']]
stopifnot(digest(file=dll, algo='sha256') == '33ba6750fda9193cf5d5d54d1a97f1716b8dab95f65c3407f91e637824068297')
source_cpp <- file.path(repo, 'research/heterogeneity/sampler-mixing/native_kernel_probe.cpp')
old_header <- file.path(out, 'baseline/include/stochtree/tree_sampler.h')
dir.create(dirname(old_header), recursive=TRUE, showWarnings=FALSE)
status <- system2('git', c('-c', 'core.fsmonitor=false', '-C', shQuote(repo), 'show',
  'cad3d6c9f7fd6219798fc021fcec8f5182fe9c2e:include/stochtree/tree_sampler.h'), stdout=old_header)
stopifnot(status == 0L)
original_env <- Sys.getenv(c('PKG_CPPFLAGS','PKG_LIBS','PKG_CXXFLAGS','PATH'), unset=NA_character_)
build <- function(version) {
  stage <- file.path(out, version); dir.create(stage, recursive=TRUE, showWarnings=FALSE)
  name <- paste0('native_probe_',version)
  stopifnot(file.copy(source_cpp, file.path(stage,paste0(name,'.cpp')), overwrite=TRUE))
  inc <- c(if(version=='baseline')file.path(out,'baseline/include'), system.file('include',package='Rcpp'),
           file.path(repo,c('include','deps/eigen','deps/fmt/include','deps/fast_double_parser/include')))
  Sys.setenv(PKG_CPPFLAGS=paste(c('-DSTOCHTREE_R_BUILD',paste0('-I',shQuote(inc))),collapse=' '),
    PKG_LIBS=shQuote(dll), PKG_CXXFLAGS='-std=gnu++17 -O2',
    PATH=paste('C:/rtools44/x86_64-w64-mingw32.static.posix/bin','C:/rtools44/usr/bin',original_env[['PATH']],sep=';'))
  wd <- getwd(); setwd(stage)
  status <- system2(file.path(R.home('bin'),'R.exe'), c('CMD','SHLIB',paste0(name,'.cpp')),
                    stdout='build.log', stderr='build.log')
  setwd(wd); if(status!=0L)stop('Native probe failed; inspect ',stage,'/build.log')
  path <- file.path(stage,paste0(name,'.dll')); dyn.load(path)
  list(name=name, dll=path, sha256=digest(file=path,algo='sha256'))
}
bridges <- setNames(lapply(c('baseline','repaired'),build),c('baseline','repaired'))
probe <- function(version,case) .Call('native_kernel_probe',case, PACKAGE=bridges[[version]]$name)

# Independent Gaussian integration: y ~ N(0, I + 11'), one prior N(0,1)
# per leaf. Include the ordinary branching stop factor at the depth cap.
exact <- function(case) {
  x <- case$x; y <- case$y
  enumerate <- function(lo,hi,depth) {
    idx <- lo:hi; nn <- length(idx); p <- case$alpha/(1+depth)^case$beta
    ll <- -.5*log1p(nn)-.5*(sum(y[idx]^2)-sum(y[idx])^2/(1+nn))
    states <- data.frame(state='L',log_mass=log1p(-p)+ll)
    if(depth<case$maxdepth && nn>=2*case$minleaf) for(k in lo:(hi-1L)) {
      if(k-lo+1<case$minleaf || hi-k<case$minleaf || x[k+1L]<=x[k])next
      a <- enumerate(lo,k,depth+1L); b <- enumerate(k+1L,hi,depth+1L)
      z <- expand.grid(a=seq_len(nrow(a)),b=seq_len(nrow(b)))
      states <- rbind(states,data.frame(state=paste0('(',k,':',a$state[z$a],':',b$state[z$b],')'),
        log_mass=log(p)+log((x[k+1L]-x[k])/(x[hi]-x[lo]))+a$log_mass[z$a]+b$log_mass[z$b]))
    }
    states
  }
  a <- enumerate(1L,length(x),0L); a$probability <- exp(a$log_mass-max(a$log_mass));
  a$probability <- a$probability/sum(a$probability); a
}
cases <- list(binary=list(x=rep(c(0,1),each=10),y=rep(0,20),alpha=.3,beta=2,minleaf=2L,maxdepth=1L),
              deeper=list(x=0:5,y=c(-1.2,-.8,-.3,.6,.9,1.4),alpha=.65,beta=2,minleaf=1L,maxdepth=2L))
tables <- list(); summary <- list(); traces <- list()
for (cn in names(cases)) {
  oracle <- exact(cases[[cn]])
  for(version in names(bridges)) {
    case <- c(cases[[cn]],list(draws=500000L,burn=10000L,seed=730001L+match(cn,names(cases))))
    a <- probe(version,case); stopifnot(all(a$states %in% oracle$state))
    tab <- table(factor(a$states, levels=oracle$state)); obs <- as.numeric(tab)/length(a$states)
    batches <- split(seq_along(a$states), rep(1:100,each=5000))
    mcse <- vapply(oracle$state,function(s)sd(vapply(batches,function(ix)mean(a$states[ix]==s),numeric(1)))/10,numeric(1))
    z <- (obs-oracle$probability)/pmax(mcse,1/length(a$states))
    tables[[paste(cn,version)]] <- data.frame(case=cn,version=version,oracle,observed=obs,mcse=mcse,z=z)
    summary[[paste(cn,version)]] <- data.frame(case=cn,version=version,states=nrow(oracle),draws=length(a$states),
      total_variation=.5*sum(abs(obs-oracle$probability)),max_abs_mcse_units=max(abs(z)),
      pass=max(abs(z))<6 && .5*sum(abs(obs-oracle$probability))<.015)
    traces[[paste(cn,version)]] <- a$states
  }
}
range_cases <- list(descending=c(3,2,1),negative=c(-3,-2,-1),mixed=c(2,-3,1,-1),
                    constant_negative=rep(-2,4),constant_positive=rep(2,4),zero=rep(0,4))
ranges <- do.call(rbind,lapply(names(bridges),function(version)do.call(rbind,lapply(names(range_cases),function(cn){
  x <- range_cases[[cn]]; a <- probe(version,list(x=x,y=rep(0,length(x)),alpha=.3,beta=2,minleaf=1L,
    maxdepth=1L,draws=0L,burn=0L,seed=1L))
  data.frame(version=version,case=cn,min=a$range[1],max=a$range[2],true_min=min(x),true_max=max(x),
    nonconstant=a$nonconstant,true_nonconstant=diff(range(x))>0,
    pass=identical(as.numeric(a$range),range(x)) && a$nonconstant==(diff(range(x))>0))
}))))
write.csv(do.call(rbind,tables),file.path(out,'oracle_states.csv'),row.names=FALSE)
write.csv(do.call(rbind,summary),file.path(out,'oracle_summary.csv'),row.names=FALSE)
write.csv(ranges,file.path(out,'range_checks.csv'),row.names=FALSE)
boundary <- probe('repaired',list(x=c(0,0,1,1),y=rep(0,4),alpha=.3,beta=2,minleaf=2L,
  maxdepth=1L,draws=20000L,burn=1000L,seed=731001L))
small <- probe('repaired',list(x=0,y=0,alpha=.3,beta=2,minleaf=2L,
  maxdepth=1L,draws=100L,burn=10L,seed=731002L))
stopifnot(length(unique(boundary$states))==2L,all(small$states=='L'))
write.csv(data.frame(check=c('exact_2_minleaf_boundary_can_grow_and_prune','small_root_is_valid_self_loop'),
  pass=c(length(unique(boundary$states))==2L,all(small$states=='L'))),
  file.path(out,'boundary_checks.csv'),row.names=FALSE)
saveRDS(list(cases=cases,bridges=bridges,traces=traces,baseline_dll=dll,
  baseline_dll_sha256=digest(file=dll,algo='sha256'),header_sha256=digest(file=file.path(repo,'include/stochtree/tree_sampler.h'),algo='sha256'),
  boundary=boundary,small_root=small),
  file.path(out,'native_probe_raw.rds'))
stopifnot(all(ranges$pass[ranges$version=='repaired']),all(do.call(rbind,summary)$pass[do.call(rbind,summary)$version=='repaired']))
print(do.call(rbind,summary)); print(ranges)
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'))
