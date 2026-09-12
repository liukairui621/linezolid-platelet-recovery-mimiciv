# Population and baseline-source corrections frozen before these associations.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';cfg<-fromJSON(file.path(R,'config/analysis_plan_v0.8.json'))$clinical
m<-readRDS(file.path(R,'cache/phase61_overlap.rds'));raw<-m$raw;f<-formula(m$fit)
ct<-as.data.frame(fread(file.path(R,'cache/phase8_context_with_chart.csv')))
ct<-ct[match(raw$hadm_id,ct$hadm_id),];stopifnot(!anyNA(ct$hadm_id),identical(as.character(raw$hadm_id),as.character(ct$hadm_id)),identical(raw$drug,ct$drug))
raw$rrt_documented_prior24h<-ct$combined_rrt24h
raw$invasive_vent_documented_prior24h<-ct$combined_vent24h
raw$pressors_documented_prior6h<-ct$pressor_available6h
raw$culture72h<-ct$nonscreen_culture72h;raw$respiratory72h<-ct$nonscreen_respiratory72h
prep<-function(x){
 lc<-names(x)[vapply(x,is.logical,logical(1))];x[lc]<-lapply(x[lc],as.integer)
 x$A<-as.integer(x$drug=='linezolid');x$gender<-factor(x$gender,levels=levels(m$data$gender))
 x$mapped_era<-factor(x$mapped_era,levels=levels(m$data$mapped_era))
 fill<-function(z,fallback=0){z[!is.finite(z)]<-NA_real_;v<-if(all(is.na(z)))fallback else median(z,na.rm=TRUE);z[is.na(z)]<-v;z}
 x$creatinine_missing<-as.integer(is.na(x$known_creatinine)|x$known_creatinine<=0)
 x$log_creatinine<-log(fill(ifelse(x$known_creatinine>0,x$known_creatinine,NA_real_),1))
 x$trend_missing<-as.integer(is.na(x$platelet_change));x$platelet_change_filled<-fill(x$platelet_change)
 x$inr_missing<-as.integer(is.na(x$baseline_inr)|x$baseline_inr<=0);x$log_inr<-log(fill(ifelse(x$baseline_inr>0,x$baseline_inr,NA_real_),1))
 x$bilirubin_missing<-as.integer(is.na(x$baseline_bilirubin)|x$baseline_bilirubin<0)
 x$log1p_bilirubin<-log1p(fill(ifelse(x$baseline_bilirubin>=0,x$baseline_bilirubin,NA_real_)))
 x$log1p_admission_days<-log1p(x$days_since_admission);x
}
fitone<-function(x){
 mm<-model.matrix(f,x);qrmm<-qr(mm);cols<-sort(qrmm$pivot[seq_len(qrmm$rank)]);MM<-mm[,cols,drop=FALSE]
 attempt<-function(maxit,slowit=1){warnings<-character();fit<-withCallingHandlers(brglmFit(x=MM,y=x$A,family=binomial('logit'),
 control=brglmControl(type='AS_mean',maxit=maxit,epsilon=1e-8,slowit=slowit)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')});list(fit=fit,warnings=warnings)}
 initial<-attempt(200);fit<-initial$fit;warnings<-initial$warnings;fallback<-!isTRUE(fit$converged)||any(!is.finite(coef(fit)))
 if(fallback){ob<-attempt(2000,.1);fit<-ob$fit;warnings<-ob$warnings}
 e<-fit$fitted.values;w<-ifelse(x$A==1,1-e,e)
 list(fit=fit,w=w,warnings=warnings,initial_warnings=initial$warnings,numerical_fallback=fallback,mm=MM,omitted=setdiff(colnames(mm),colnames(MM)),valid=isTRUE(fit$converged)&&all(is.finite(w))&&all(w>0&w<1)&&all(is.finite(coef(fit))))
}
# Check direct matrix fitting reproduces the previous formula implementation.
check<-fitone(m$data);err<-max(abs(check$w-m$weights));stopifnot(check$valid,err<1e-8)
ess<-function(w)sum(w)^2/sum(w^2)
allbalance<-list();reports<-list();models<-list()
for(pop in cfg$populations){
 ix0<-which(ct[[pop]]);r<-raw[ix0,,drop=FALSE];x<-prep(r);con<-ct[ix0,,drop=FALSE]
 fit<-fitone(x);stopifnot(fit$valid,all(table(x$A)>0))
 mm<-fit$mm[,colnames(fit$mm)!='(Intercept)',drop=FALSE]
 extras<-cbind(as.matrix(x[,c('ufh_documented_prior7d','lmwh_documented_prior7d','selected_other_drug_prior7d')]),
 as.matrix(con[,c('known_named_nonscreen7d','known_enterococcus7d','nonscreen_blood72h','nonscreen_urine72h','observed_creatinine_aki')]),
 model.matrix(~factor(icu_units)-1,transform(con,icu_units=ifelse(is.na(icu_units),'not_explicit_icu',icu_units))))
 audit<-cbind(mm,extras);storage.mode(audit)<-'double'
 bal<-rbindlist(lapply(seq_len(ncol(audit)),function(j){v<-audit[,j];a<-x$A==1;b<-!a;s<-sqrt((var(v[a])+var(v[b]))/2)
  data.table(population=pop,term=colnames(audit)[j],included=j<=ncol(mm),SMD_before=if(s>0)(mean(v[a])-mean(v[b]))/s else NA_real_,
   SMD_after=if(s>0)(weighted.mean(v[a],fit$w[a])-weighted.mean(v[b],fit$w[b]))/s else NA_real_)}))
 cl<-split(seq_len(nrow(r)),r$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat)
 set.seed(cfg$bootstrap$seed);B<-cfg$bootstrap$B;boot<-vector('list',B);diag<-vector('list',B)
 priorpath<-file.path(R,'cache',paste0('phase8_',pop,'.rds'));prior<-if(file.exists(priorpath))readRDS(priorpath) else NULL
 if(!is.null(prior))stopifnot(identical(prior$raw,r),max(abs(prior$weights-fit$w))<1e-10,prior$report$bootstrap_valid==B,prior$report$bootstrap_warnings==0)
 for(i in seq_len(B)){
  namesdraw<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE);ix<-unlist(cl[namesdraw],use.names=FALSE)
  if(!is.null(prior)){
   stopifnot(identical(ix,prior$bootstrap[[i]]$ix));boot[[i]]<-prior$bootstrap[[i]];diag[[i]]<-data.table(replicate=i,valid=TRUE,warnings=0L,numerical_fallback=FALSE,reused=TRUE);next
  }
  ob<-tryCatch(fitone(prep(r[ix,,drop=FALSE])),error=function(e)list(valid=FALSE,warnings=conditionMessage(e)))
  if(ob$valid)boot[[i]]<-list(ix=ix,weights=ob$w)
  diag[[i]]<-data.table(replicate=i,valid=ob$valid,warnings=length(ob$warnings),numerical_fallback=isTRUE(ob$numerical_fallback),reused=FALSE)
  if(i%%500==0){cat(pop,i,'/',B,'refits\n');flush.console()}
 }
 dg<-rbindlist(diag);maxsmd<-max(abs(bal[included==TRUE]$SMD_after),na.rm=TRUE)
 rp<-list(population=pop,n_LZD=sum(x$A==1),n_VAN=sum(x$A==0),ESS_LZD=ess(fit$w[x$A==1]),ESS_VAN=ess(fit$w[x$A==0]),
 nonintercept_df=ncol(fit$mm)-1,omitted_columns=fit$omitted,max_included_SMD=maxsmd,max_audited_not_included_SMD=max(abs(bal[included==FALSE]$SMD_after),na.rm=TRUE),
 bootstrap_B=B,bootstrap_valid=sum(dg$valid),bootstrap_warnings=sum(dg$warnings>0),bootstrap_numerical_fallbacks=sum(dg$numerical_fallback),point_numerical_fallback=fit$numerical_fallback,reused_draws=sum(dg$reused),diagnostic_gate=fit$valid&&maxsmd<=.1&&mean(dg$valid)>=.98,
 interpretation='Population-specific overlap association. Gate addresses numeric stability/included balance only; not common indication or causal exchangeability.')
 models[[pop]]<-list(raw=r,data=x,context=con,fit=fit$fit,weights=fit$w,bootstrap=boot,bootstrap_diagnostics=dg,point_initial_warnings=fit$initial_warnings,report=rp)
 reports[[pop]]<-rp;allbalance[[pop]]<-bal
 saveRDS(models[[pop]],file.path(R,'cache',paste0('phase8_',pop,'.rds')))
 write_json(rp,file.path(R,'reports',paste0('MODEL_',pop,'_v0.8.json')),pretty=TRUE,auto_unbox=TRUE,digits=10)
 cat(toJSON(rp,auto_unbox=TRUE),'\n')
}
fwrite(rbindlist(allbalance),file.path(R,'reports/COVARIATE_BALANCE_v0.8.csv'))
write_json(list(completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),no_outcomes_loaded=TRUE,original_formula_vs_matrix_weight_error=err,populations=reports),file.path(R,'reports/CLINICAL_MODELS_v0.8.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.8_clinical.txt'))
