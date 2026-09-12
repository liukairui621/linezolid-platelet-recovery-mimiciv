# Outcome-free point weighting and200 patient-cluster refits for clinical readiness.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';cfg<-fromJSON(file.path(R,'config/readiness_addendum_v0.10.json'));m<-readRDS(file.path(R,'cache/phase61_overlap.rds'))
d<-as.data.frame(fread(file.path(R,'cache/phase10_baseline.csv')))
d$routine_culture_no_recorded_lmwh_shared<-d$routine_bacterial_culture72h_shared & !d$lmwh_available7d
shared<-with(d[d$routine_culture_no_recorded_lmwh_shared,],table(icu_unit,drug));units<-rownames(shared)[apply(shared>0,1,all)];d$routine_culture_no_recorded_lmwh_shared<-d$routine_culture_no_recorded_lmwh_shared & d$icu_unit%in%units
# Source substitutions specified before this diagnostic run.
d$culture72h<-d$routine_bacterial_culture72h;d$respiratory72h<-d$routine_respiratory72h;d$nonscreen_blood72h<-d$routine_blood72h;d$nonscreen_urine72h<-d$routine_urine72h;d$known_named_nonscreen7d<-d$known_named_routine7d
f<-update(formula(m$fit),. ~ .+icu_unit+nonscreen_blood72h+nonscreen_urine72h+known_named_nonscreen7d+known_enterococcus7d+observed_creatinine_aki+ufh_active_available7d+lmwh_available7d+marrowdrug_available7d+beta_lactam_available24h)
for(e in parse(file.path(R,'scripts/52_clinical_subset_models_v08.R')))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])&&as.character(e[[2]])%in%c('prep','fitone','ess'))eval(e)
cols<-c('subject_id','hadm_id','drug','age','gender','known_platelets','known_creatinine','rrt_documented_prior24h','invasive_vent_documented_prior24h','pressors_documented_prior6h','platelet_change','prior_recorded_marrow','days_since_admission','culture72h','respiratory72h','baseline_inr','baseline_bilirubin','mapped_era','any_emar_available_before','icu_unit','nonscreen_blood72h','nonscreen_urine72h','known_named_nonscreen7d','known_enterococcus7d','observed_creatinine_aki','ufh_active_available7d','lmwh_available7d','marrowdrug_available7d','beta_lactam_available24h')
allbal<-list();report<-list();flow<-list()
for(pop in cfg$populations){
 ix0<-which(d[[pop]]);raw<-d[ix0,cols];raw$icu_unit<-factor(raw$icu_unit);x<-prep(raw);ob<-fitone(x);stopifnot(ob$valid)
 mm<-cbind(model.matrix(~icu_unit-1,x),ob$mm[,!grepl('^icu_unit|Intercept',colnames(ob$mm)),drop=FALSE]);a<-x$A==1
 bal<-rbindlist(lapply(colnames(mm),function(k){v<-mm[,k];sd0<-sqrt((var(v[a])+var(v[!a]))/2);data.table(population=pop,term=k,mean_LZD_before=mean(v[a]),mean_VAN_before=mean(v[!a]),mean_LZD_after=weighted.mean(v[a],ob$w[a]),mean_VAN_after=weighted.mean(v[!a],ob$w[!a]),SMD_before=(mean(v[a])-mean(v[!a]))/sd0,SMD_after=(weighted.mean(v[a],ob$w[a])-weighted.mean(v[!a],ob$w[!a]))/sd0)}))
 mx<-max(abs(bal$SMD_after),na.rm=TRUE);cat(pop,'point N',sum(a),sum(!a),'maxSMD',mx,'\n');flush.console()
 cl<-split(seq_len(nrow(raw)),raw$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat);set.seed(cfg$bootstrap$seed);B<-cfg$bootstrap$B;boot<-vector('list',B);dg<-list()
 for(i in seq_len(B)){
  nd<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE);ix<-unlist(cl[nd],use.names=FALSE);bb<-tryCatch(fitone(prep(raw[ix,,drop=FALSE])),error=function(e)list(valid=FALSE,warnings=conditionMessage(e)))
  if(bb$valid)boot[[i]]<-list(ix=ix,weights=bb$w)
  dg[[i]]<-data.table(replicate=i,valid=bb$valid,warnings=length(bb$warnings),fallback=isTRUE(bb$numerical_fallback))
 }
 diag<-rbindlist(dg);rp<-list(population=pop,n_LZD=sum(a),n_VAN=sum(!a),patients_LZD=length(unique(raw$subject_id[a])),patients_VAN=length(unique(raw$subject_id[!a])),ESS_LZD=ess(ob$w[a]),ESS_VAN=ess(ob$w[!a]),nonintercept_df=ncol(ob$mm)-1,omitted=ob$omitted,max_all_level_SMD=mx,valid_draws=sum(diag$valid),requested_draws=B,numerical_fallbacks=sum(diag$fallback),warning_draws=sum(diag$warnings>0),readiness_gate=ob$valid&&mx<=.1&&mean(diag$valid)>=.98)
 saveRDS(list(raw=raw,data=x,weights=ob$w,fit=ob$fit,bootstrap=boot,diagnostics=diag,report=rp),file.path(R,'cache',paste0('phase10_readiness_',pop,'.rds')))
 report[[pop]]<-rp;allbal[[pop]]<-bal;cat(toJSON(rp,auto_unbox=TRUE),'\n');flush.console()
 if(pop==cfg$primary_candidate){
  pr<-d[ix0,];flow[[pop]]<-rbindlist(lapply(c('linezolid','vancomycin'),function(dr){z<-pr[pr$drug==dr,];data.table(drug=dr,n=nrow(z),patients=length(unique(z$subject_id)),old_cohort=sum(z$in_original_cohort),no_discharge_sepsis_code=sum(!z$coded_sepsis))}))
 }
}
fwrite(rbindlist(allbal),file.path(R,'reports/READINESS_BALANCE_v0.10.csv'));fwrite(rbindlist(flow),file.path(R,'reports/SELECTED_BASELINE_COHORT_v0.10.csv'))
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),no_recovery_outcomes_loaded=TRUE,primary_candidate=cfg$primary_candidate,models=report,limits='Readiness only.200 resamples are not endpoint confidence intervals. Measured overlap does not establish same clinical indication or biological no-LMWH exposure.'),file.path(R,'reports/BASELINE_MODEL_READINESS_v0.10.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.10_readiness.txt'))
