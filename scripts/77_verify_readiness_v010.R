# Reproduce baseline-only membership, weights and cluster draws; export aggregates only.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery'
for(e in parse(file.path(R,'scripts/76_baseline_model_readiness_v010.R'))){
 if(is.call(e)&&identical(e[[1]],as.name('for'))&&identical(e[[2]],as.name('pop')))break
 eval(e)
}
res<-list()
for(pop in cfg$populations){
 saved<-readRDS(file.path(R,'cache',paste0('phase10_readiness_',pop,'.rds')))
 r<-saved$raw;x<-prep(r);point<-fitone(x)
 stopifnot(identical(r$hadm_id,d$hadm_id[which(d[[pop]])]),!anyDuplicated(r$hadm_id))
 stopifnot(!any(grepl('recovery|disch|death|outcome|coded_sepsis',names(r),ignore.case=TRUE)))
 err<-max(abs(point$w-saved$weights));stopifnot(point$valid,err<1e-10)
 a<-x$A==1;es<-c(ess(point$w[a]),ess(point$w[!a]))
 stopifnot(max(abs(es-c(saved$report$ESS_LZD,saved$report$ESS_VAN)))<1e-8)
 cl<-split(seq_len(nrow(r)),r$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat)
 set.seed(cfg$bootstrap$seed);checked<-0L;refiterr<-0;refitn<-0L
 for(i in seq_along(saved$bootstrap)){
  nd<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE)
  ix<-unlist(cl[nd],use.names=FALSE);bb<-saved$bootstrap[[i]]
  if(is.null(bb)){stopifnot(!saved$diagnostics$valid[i]);next}
  stopifnot(identical(ix,bb$ix),length(bb$weights)==length(ix),all(is.finite(bb$weights)),all(bb$weights>0&bb$weights<1))
  checked<-checked+1L
  if(i%in%c(1,100,200)){re<-fitone(prep(r[ix,,drop=FALSE]));stopifnot(re$valid);refiterr<-max(refiterr,max(abs(re$w-bb$weights)));refitn<-refitn+1L}
 }
 stopifnot(refiterr<1e-10,checked==saved$report$valid_draws)
 res[[pop]]<-list(membership_exact=TRUE,no_outcome_fields_in_model=TRUE,point_weight_reproduction_max_error=err,ESS_reproduced=TRUE,valid_draws_seed_and_cluster_indices_verified=checked,selected_refits=c(1,100,200),selected_valid_refits=refitn,selected_refit_weight_max_error=refiterr,all_retained_weights_finite_positive=TRUE)
}
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),patient_rows_exported=FALSE,populations=res),file.path(R,'reports/READINESS_REFIT_CHECKS_v0.10.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
print(res)
