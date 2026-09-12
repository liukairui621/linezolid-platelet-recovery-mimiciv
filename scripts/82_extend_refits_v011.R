# Frozen two-population patient-cluster refits. No outcomes loaded.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery'
for(e in parse(file.path(R,'scripts/76_baseline_model_readiness_v010.R'))){
 if(is.call(e)&&identical(e[[1]],as.name('for'))&&identical(e[[2]],as.name('pop')))break
 eval(e)
}
plan<-fromJSON(file.path(R,'config/analysis_plan_v0.11.json'));res<-list()
for(pop in c(plan$populations$primary,plan$populations$sensitivity)){
 m<-readRDS(file.path(R,'cache',paste0('phase10_readiness_',pop,'.rds')));raw<-m$raw;x<-prep(raw);point<-fitone(x)
 stopifnot(point$valid,max(abs(point$w-m$weights))<1e-10,length(m$bootstrap)==200)
 cl<-split(seq_len(nrow(raw)),raw$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat)
 set.seed(plan$bootstrap$seed);B<-plan$bootstrap$total_B;boot<-vector('list',B);diag<-vector('list',B)
 for(i in seq_len(B)){
  nd<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE);ix<-unlist(cl[nd],use.names=FALSE)
  if(i<=200){stopifnot(identical(ix,m$bootstrap[[i]]$ix));boot[[i]]<-m$bootstrap[[i]];diag[[i]]<-m$diagnostics[i,];diag[[i]]$reused<-TRUE
  }else{
   z<-tryCatch(fitone(prep(raw[ix,,drop=FALSE])),error=function(e)list(valid=FALSE,warnings=conditionMessage(e)))
   if(z$valid)boot[[i]]<-list(ix=ix,weights=z$w)
   diag[[i]]<-data.table(replicate=i,valid=z$valid,warnings=length(z$warnings),fallback=isTRUE(z$numerical_fallback),reused=FALSE)
  }
  if(i%%250==0){cat(pop,i,'/',B,'\n');flush.console()}
 }
 dg<-rbindlist(diag);rp<-m$report;rp$valid_draws<-sum(dg$valid);rp$requested_draws<-B;rp$numerical_fallbacks<-sum(dg$fallback);rp$warning_draws<-sum(dg$warnings>0);rp$readiness_gate<-point$valid&&rp$max_all_level_SMD<=.1&&mean(dg$valid)>=.98
 m$bootstrap<-boot;m$diagnostics<-dg;m$report<-rp
 saveRDS(m,file.path(R,'cache',paste0('phase11_',pop,'.rds')))
 res[[pop]]<-rp;cat(toJSON(rp,auto_unbox=TRUE),'\n');flush.console()
}
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),no_outcomes_loaded=TRUE,seed=plan$bootstrap$seed,first200_indices_identical=TRUE,first200_weights_reused=TRUE,point_weights_unchanged=TRUE,models=res),file.path(R,'reports/CLINICAL_MODELS_v0.11.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.11.txt'))
