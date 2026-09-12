# Internal-only weight/refit checks; export aggregate errors, never patient rows or draws.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';ee<-parse(file.path(R,'scripts/61_comparability_models_v09.R'))
for(e in ee){if(is.call(e)&&identical(e[[1]],as.name('for'))&&identical(e[[2]],as.name('pop')))break;eval(e)}
res<-list()
for(pop in cfg$populations){
 saved<-readRDS(file.path(R,'cache',paste0('phase9_',pop,'.rds')));r<-saved$raw;x<-prep(r);point<-fitone(x);err<-max(abs(point$w-saved$weights));stopifnot(point$valid,err<1e-10)
 cl<-split(seq_len(nrow(r)),r$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat);set.seed(cfg$bootstrap$seed)
 checked<-0L;drawerr<-0;refiterr<-0
 for(i in seq_along(saved$bootstrap)){
  nd<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE);ix<-unlist(cl[nd],use.names=FALSE);bb<-saved$bootstrap[[i]]
  if(is.null(bb)){stopifnot(!saved$bootstrap_diagnostics$valid[i]);next}
  stopifnot(identical(ix,bb$ix),length(bb$weights)==length(ix),all(is.finite(bb$weights)),all(bb$weights>0&bb$weights<1));checked<-checked+1L
  if(i%in%c(1,500,2000)){re<-fitone(prep(r[ix,,drop=FALSE]));stopifnot(re$valid);refiterr<-max(refiterr,max(abs(re$w-bb$weights)))}
 }
 stopifnot(refiterr<1e-10,checked==saved$report$bootstrap_valid)
 res[[pop]]<-list(point_weight_reproduction_max_error=err,valid_draws_seed_and_cluster_indices_verified=checked,selected_refits=c(1,500,2000),selected_refit_weight_max_error=refiterr,all_retained_weights_finite_positive=TRUE)
}
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),patient_rows_exported=FALSE,populations=res),file.path(R,'reports/CLINICAL_REFIT_CHECKS_v0.9.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
print(res)
