suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';cfg<-fromJSON(file.path(R,'config/analysis_plan_v0.9.json'))$clinical
m<-readRDS(file.path(R,'cache/phase61_overlap.rds'));raw<-m$raw
ct<-as.data.frame(fread(file.path(R,'cache/phase9_context.csv')));ct<-ct[match(raw$hadm_id,ct$hadm_id),]
stopifnot(!anyNA(ct$hadm_id),identical(raw$drug,ct$drug));raw$rrt_documented_prior24h<-ct$combined_rrt24h;raw$invasive_vent_documented_prior24h<-ct$combined_vent24h
raw$pressors_documented_prior6h<-ct$pressor_available6h;raw$culture72h<-ct$nonscreen_culture72h;raw$respiratory72h<-ct$nonscreen_respiratory72h
extra<-c('nonscreen_blood72h','nonscreen_urine72h','known_named_nonscreen7d','known_enterococcus7d','observed_creatinine_aki','ufh_active_available7d','lmwh_available7d','marrowdrug_available7d','beta_lactam_available24h')
raw[extra]<-ct[extra];raw$icu_unit<-ct$icu_units
f<-update(formula(m$fit),. ~ .+icu_unit+nonscreen_blood72h+nonscreen_urine72h+known_named_nonscreen7d+known_enterococcus7d+observed_creatinine_aki+ufh_active_available7d+lmwh_available7d+marrowdrug_available7d+beta_lactam_available24h)
expr<-parse(file.path(R,'scripts/52_clinical_subset_models_v08.R'))
for(e in expr)if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%c('prep','fitone','ess'))eval(e)
balall<-list();reports<-list()
for(pop in cfg$populations){
 ix0<-which(ct[[pop]]);r<-raw[ix0,,drop=FALSE];r$icu_unit<-factor(r$icu_unit);x<-prep(r);fit<-fitone(x);stopifnot(fit$valid)
 mm<-fit$mm[,colnames(fit$mm)!='(Intercept)',drop=FALSE];a<-x$A==1
 bal<-rbindlist(lapply(colnames(mm),function(k){v<-mm[,k];s<-sqrt((var(v[a])+var(v[!a]))/2);data.table(population=pop,term=k,SMD_before=(mean(v[a])-mean(v[!a]))/s,SMD_after=(weighted.mean(v[a],fit$w[a])-weighted.mean(v[!a],fit$w[!a]))/s)}))
 cl<-split(seq_len(nrow(r)),r$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat)
 B<-cfg$bootstrap$B;set.seed(cfg$bootstrap$seed);boot<-vector('list',B);diag<-vector('list',B)
 for(i in 1:B){
  nd<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE);ix<-unlist(cl[nd],use.names=FALSE)
  ob<-tryCatch(fitone(prep(r[ix,,drop=FALSE])),error=function(e)list(valid=FALSE,warnings=conditionMessage(e)))
  if(ob$valid)boot[[i]]<-list(ix=ix,weights=ob$w)
  diag[[i]]<-data.table(replicate=i,valid=ob$valid,warnings=length(ob$warnings),fallback=isTRUE(ob$numerical_fallback))
  if(i%%500==0){cat(pop,i,'/2000\n');flush.console()}
 }
 dg<-rbindlist(diag);mx<-max(abs(bal$SMD_after),na.rm=TRUE)
 rp<-list(population=pop,n_LZD=sum(a),n_VAN=sum(!a),ESS_LZD=ess(fit$w[a]),ESS_VAN=ess(fit$w[!a]),nonintercept_df=ncol(mm),max_included_SMD=mx,
 omitted=fit$omitted,bootstrap_B=B,bootstrap_valid=sum(dg$valid),bootstrap_warnings=sum(dg$warnings>0),bootstrap_numerical_fallbacks=sum(dg$fallback),
 diagnostic_gate=mx<=.1&&mean(dg$valid)>=.98,formula=paste(deparse(f),collapse=' '),no_outcomes_loaded=TRUE)
 saveRDS(list(raw=r,data=x,context=ct[ix0,],fit=fit$fit,weights=fit$w,bootstrap=boot,bootstrap_diagnostics=dg,report=rp),file.path(R,'cache',paste0('phase9_',pop,'.rds')))
 reports[[pop]]<-rp;balall[[pop]]<-bal;cat(toJSON(rp,auto_unbox=TRUE),'\n')
}
fwrite(rbindlist(balall),file.path(R,'reports/COVARIATE_BALANCE_v0.9.csv'))
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),no_outcomes_loaded=TRUE,populations=reports),file.path(R,'reports/CLINICAL_MODELS_v0.9.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.9_clinical.txt'))
