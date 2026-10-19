# May be sourced by run_pilot.R or run standalone from the repository root.
if (!exists("out",inherits=FALSE)) out <- "D:/heterogeneity-bayes-research-2026-10-16/pilot-results"
if (!exists("mode",inherits=FALSE)) mode <- "full"
if (!exists("wilson_interval")) source("research/bayes-factor-pilot-2026-10-16/bayes_factor_sidecar.R")
res <- read.csv(file.path(out,paste0("results_",mode,".csv")),stringsAsFactors=FALSE)
res$scale_group <- ifelse(is.na(res$scale),"NA",as.character(res$scale))
grp <- interaction(res[,c("stage","scenario","design","basis","method","prior","scale_group")],drop=TRUE)
summaries <- lapply(split(res,grp),function(x) {
  b <- x[1,c("stage","scenario","design","basis","method","prior","scale")]
  events <- if(all(is.na(x$log_bf)))list(p_lt_0.05=x$p<.05) else list(
    BF10_gt_3=x$log_bf>log(3),BF10_gt_10=x$log_bf>log(10),
    BF01_gt_3=x$log_bf< -log(3),BF01_gt_10=x$log_bf< -log(10),
    inconclusive_1over3_to_3=abs(x$log_bf)<=log(3))
  do.call(rbind,lapply(names(events),function(event) {
    v<-events[[event]];n<-sum(!is.na(v));k<-sum(v,na.rm=TRUE);ci<-wilson_interval(k,n)
    cbind(b,event=event,n=n,k=k,rate=k/n,ci_low=ci[1],ci_high=ci[2],
          median_log_bf=if(all(is.na(x$log_bf)))NA_real_ else median(x$log_bf,na.rm=TRUE),
          median_importance_ess=if(all(is.na(x$importance_ess)))NA_real_ else median(x$importance_ess,na.rm=TRUE),
          concentration_warning_rate=mean(x$diagnostic_warning!="none"))
  }))
})
summary<-do.call(rbind,summaries)
write.csv(summary,file.path(out,paste0("summary_",mode,".csv")),row.names=FALSE)
if(mode=="full") {
  png(file.path(out,"pilot_evidence_rates.png"),width=1800,height=1100,res=160)
  par(mfrow=c(2,4),mar=c(7,4,3,1),oma=c(0,0,3,0))
  for(design in c("randomized","confounded_known"))for(scenario in c("null","sparse","dense","quadratic")) {
    # Explicit indexing avoids data-mask variable capture.
    z<-summary[summary$stage=="forest" & summary$design==design & summary$scenario==scenario &
      summary$basis=="linear" & ((summary$method=="restricted_external_score" & summary$event=="p_lt_0.05") |
      (summary$method=="null_posterior_integrated_bf" & summary$scale==.5 & summary$event=="BF10_gt_10")),]
    z<-z[order(match(z$method,c("restricted_external_score","null_posterior_integrated_bf")),z$prior),]
    labels<-ifelse(z$method=="restricted_external_score","Score p<.05",paste(z$prior,"BF>10"))
    bx<-barplot(z$rate,names.arg=labels,las=2,ylim=c(0,1),col=c("#537CA6","#E4A653","#7EAB94"),
                main=paste(design,scenario,sep="\n"),ylab="Decision frequency")
    arrows(bx,z$ci_low,bx,z$ci_high,angle=90,code=3,length=.04)
  }
  mtext("Exploratory pilot: unequal decision thresholds; 95% Wilson intervals, 40 replications/cell",outer=TRUE,cex=1)
  dev.off()
  # Portrait layout for a full page in the printed research report. Keep the
  # original wide figure above for screen use and reproducibility.
  png(file.path(out,"pilot_evidence_rates_print.png"),width=1800,height=2600,res=200)
  par(mfrow=c(4,2),mar=c(4.2,4.5,3.5,1),oma=c(0,0,4,0),cex=1.05)
  scenario_names<-c(null="Constant effect",sparse="Sparse linear heterogeneity",
                    dense="Dense linear heterogeneity",quadratic="Quadratic heterogeneity")
  design_names<-c(randomized="Randomized treatment",confounded_known="Known confounded propensity")
  for(scenario in c("null","sparse","dense","quadratic"))for(design in c("randomized","confounded_known")) {
    z<-summary[summary$stage=="forest" & summary$design==design & summary$scenario==scenario &
      summary$basis=="linear" & ((summary$method=="restricted_external_score" & summary$event=="p_lt_0.05") |
      (summary$method=="null_posterior_integrated_bf" & summary$scale==.5 & summary$event=="BF10_gt_10")),]
    z<-z[order(match(z$method,c("restricted_external_score","null_posterior_integrated_bf")),z$prior),]
    labels<-ifelse(z$method=="restricted_external_score","Score\np < .05",
                   ifelse(z$prior=="dense","Dense\nBF > 10","Sparse\nBF > 10"))
    bx<-barplot(z$rate,names.arg=labels,las=1,ylim=c(0,1),
      col=c("#537CA6","#E4A653","#7EAB94"),
      main=paste(design_names[design],scenario_names[scenario],sep="\n"),ylab="Decision frequency")
    arrows(bx,z$ci_low,bx,z$ci_high,angle=90,code=3,length=.045)
  }
  mtext("Exploratory decision frequencies",outer=TRUE,line=2.1,cex=1.25,font=2)
  mtext("Unequal thresholds; 95% Wilson intervals; 40 replications per cell",outer=TRUE,line=.7,cex=1)
  dev.off()
}
cat("Wrote",nrow(summary),"summary rows\n")
