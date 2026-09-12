# Patient-cluster bootstrap design stability; no endpoints.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
set.seed(42)
root<-'/root/projects/linezolid_platelet_recovery'
mods<-readRDS(file.path(root,'cache/phase4_overlap_models.rds'))
out<-list();B<-200L
for(nm in c('core_baseline','extended_baseline','baseline_emar_covered')) {
  obj<-mods[[nm]];d<-obj$data;f<-formula(obj$fit)
  groups<-split(seq_len(nrow(d)),d$subject_id)
  # A patient exposed to both drugs is resampled as a unit, not split across arms.
  patterns<-vapply(groups,function(ii)paste(sort(unique(d$A[ii])),collapse=','),character(1))
  strata<-split(names(groups),patterns)
  reps<-vector('list',B)
  for(b in seq_len(B)) {
    sampled<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE)
    ix<-unlist(groups[sampled],use.names=FALSE);x<-d[ix,,drop=FALSE]
    warns<-character()
    fit<-tryCatch(withCallingHandlers(glm(f,data=x,family=binomial(),control=glm.control(maxit=100,epsilon=1e-10)),
      warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
    if(inherits(fit,'error')) {reps[[b]]<-data.frame(replicate=b,error=TRUE,converged=FALSE,full_rank=FALSE,warning=TRUE,large_coefficient=NA,ESS_LZD=NA,ESS_VAN=NA);next}
    e<-fitted(fit);w<-ifelse(x$A==1,1-e,e)
    # Large coefficients are a numerical-instability flag, not a formal separation test.
    reps[[b]]<-data.frame(replicate=b,error=FALSE,converged=fit$converged,
      full_rank=fit$rank==length(coef(fit)),warning=length(warns)>0,
      large_coefficient=any(abs(coef(fit))>15,na.rm=TRUE),
      ESS_LZD=sum(w[x$A==1])^2/sum(w[x$A==1]^2),ESS_VAN=sum(w[x$A==0])^2/sum(w[x$A==0]^2))
  }
  z<-rbindlist(reps);fwrite(z,file.path(root,paste0('cache/internal_reports/bootstrap_',nm,'.csv')))
  out[[nm]]<-list(replicates=B,errors=sum(z$error),not_converged=sum(!z$converged),
    rank_deficient=sum(!z$full_rank),warning_replicates=sum(z$warning),
    abs_coefficient_gt15_replicates=sum(z$large_coefficient,na.rm=TRUE),
    ESS_LZD_quantiles_025_50_975=unname(quantile(z$ESS_LZD,c(.025,.5,.975),na.rm=TRUE)),
    ESS_VAN_quantiles_025_50_975=unname(quantile(z$ESS_VAN,c(.025,.5,.975),na.rm=TRUE)))
}
report<-list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),seed=42,
 method='200 patient-cluster resamples within patient treatment-pattern strata. Refit PS each time. Baseline median fills fixed from design snapshot; diagnostic bootstrap, not outcome inference.',
 models=out,interpretation='Convergence and exact mean balance do not prove adequate information. Large coefficient flags are not a formal separation test. Bootstrap ESS is not an effect confidence interval.')
write_json(report,file.path(root,'reports/WEIGHT_MODEL_STABILITY_v0.4.json'),pretty=TRUE,auto_unbox=TRUE,digits=8)
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE,digits=5))
