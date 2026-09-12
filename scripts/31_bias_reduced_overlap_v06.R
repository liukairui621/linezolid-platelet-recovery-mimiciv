# Outcome-free mean-bias-reduced overlap model and patient-cluster resamples.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';start<-format(Sys.time(),tz='UTC',usetz=TRUE)
planfile<-file.path(R,'config/analysis_plan_v0.6.json');plan<-fromJSON(planfile)
stopifnot(plan$version=='0.6')
raw<-as.data.frame(fread(file.path(R,'cache/phase4_baseline_only.csv')))
cal<-fread(file.path(R,'cache/phase5_calendar.csv'))
raw$mapped_era<-cal$mapped_era[match(raw$hadm_id,cal$hadm_id)]
stopifnot(nrow(raw)==1576,!anyNA(raw$mapped_era),!anyDuplicated(raw$hadm_id))
old<-readRDS(file.path(R,'cache/phase5_weight_models.rds'))$calendar_source
f<-formula(old$fit)
prep<-function(x){
 lc<-names(x)[vapply(x,is.logical,logical(1))];x[lc]<-lapply(x[lc],as.integer)
 x$A<-as.integer(x$drug=='linezolid');x$gender<-factor(x$gender,levels=levels(old$data$gender))
 x$mapped_era<-factor(x$mapped_era,levels=c('mapped_early','boundary_uncertain','mapped_late'))
 fill<-function(z){z[!is.finite(z)]<-NA_real_;stopifnot(any(!is.na(z)));z[is.na(z)]<-median(z,na.rm=TRUE);z}
 x$creatinine_missing<-as.integer(is.na(x$known_creatinine)|x$known_creatinine<=0)
 x$log_creatinine<-log(fill(ifelse(x$known_creatinine>0,x$known_creatinine,NA_real_)))
 x$trend_missing<-as.integer(is.na(x$platelet_change));x$platelet_change_filled<-fill(x$platelet_change)
 x$inr_missing<-as.integer(is.na(x$baseline_inr)|x$baseline_inr<=0);x$log_inr<-log(fill(ifelse(x$baseline_inr>0,x$baseline_inr,NA_real_)))
 x$bilirubin_missing<-as.integer(is.na(x$baseline_bilirubin)|x$baseline_bilirubin<0)
 x$log1p_bilirubin<-log1p(fill(ifelse(x$baseline_bilirubin>=0,x$baseline_bilirubin,NA_real_)))
 x$log1p_admission_days<-log1p(x$days_since_admission);x
}
fitone<-function(x){
 warnings<-character();fit<-withCallingHandlers(glm(f,data=x,family=binomial('logit'),method=brglmFit,
   type='AS_mean',control=brglmControl(maxit=200,epsilon=1e-8)),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 e<-fitted(fit);w<-ifelse(x$A==1,1-e,e)
 list(fit=fit,w=w,warnings=warnings,valid=isTRUE(fit$converged)&&all(is.finite(w))&&all(w>0&w<1)&&all(is.finite(coef(fit)[!is.na(coef(fit))])))
}
# A separated synthetic case must give finite predictions and converge.
smoke<-glm(a~z,data=data.frame(a=c(0,0,0,1,1,1),z=1:6),family=binomial(),method=brglmFit,type='AS_mean')
stopifnot(smoke$converged,all(is.finite(coef(smoke))),all(fitted(smoke)>0&fitted(smoke)<1))
x<-prep(raw);fit<-fitone(x);stopifnot(fit$valid,fit$fit$rank-1==21)
mm<-model.matrix(fit$fit)[,-1,drop=FALSE]
audit<-cbind(mm,as.matrix(x[,c('ufh_documented_prior7d','lmwh_documented_prior7d','selected_other_drug_prior7d')]))
bal<-rbindlist(lapply(colnames(audit),function(k){v<-audit[,k];a<-x$A==1;b<-!a;s<-sqrt((var(v[a])+var(v[b]))/2)
 data.table(term=k,included=k%in%colnames(mm),SMD_before=(mean(v[a])-mean(v[b]))/s,
 SMD_after=(weighted.mean(v[a],fit$w[a])-weighted.mean(v[b],fit$w[b]))/s)}))
ess<-function(w)sum(w)^2/sum(w^2)
arm<-lapply(0:1,function(a){ii<-x$A==a;w<-fit$w[ii];list(drug=if(a)'linezolid' else 'vancomycin',n=sum(ii),
 ESS=ess(w),patient_weight_ESS=ess(rowsum(w,x$subject_id[ii])[,1]),max_normalized_weight=max(w)/sum(w),
 weighted_late=weighted.mean(x$mapped_era[ii]=='mapped_late',w),weighted_early=weighted.mean(x$mapped_era[ii]=='mapped_early',w))})
set.seed(plan$clinical$phase6_execution$bootstrap$seed);B<-plan$clinical$phase6_execution$bootstrap$B
cl<-split(seq_len(nrow(raw)),raw$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),'');strata<-split(names(cl),pat)
boot<-vector('list',B);stat<-vector('list',B)
for(i in seq_len(B)){
 namesdraw<-unlist(lapply(strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE)
 ix<-unlist(cl[namesdraw],use.names=FALSE);z<-prep(raw[ix,,drop=FALSE])
 ob<-tryCatch(fitone(z),error=function(e)list(valid=FALSE,error=conditionMessage(e),warnings=character()))
 if(ob$valid)boot[[i]]<-list(ix=ix,weights=ob$w)
 stat[[i]]<-data.table(replicate=i,valid=ob$valid,error=if(is.null(ob$error))NA_character_ else ob$error,
  warnings=length(ob$warnings),converged=if(is.null(ob$fit))FALSE else isTRUE(ob$fit$converged),
  max_abs_estimable_coefficient=if(ob$valid)max(abs(coef(ob$fit)),na.rm=TRUE) else NA_real_,
  ESS_LZD=if(ob$valid)ess(ob$w[z$A==1]) else NA_real_,ESS_VAN=if(ob$valid)ess(ob$w[z$A==0]) else NA_real_)
 if(i%%50==0){cat(i,'/',B,'completed\n');flush.console()}
}
bs<-rbindlist(stat);maxsmd<-max(abs(bal$SMD_after[bal$included]));pass<-fit$valid&&maxsmd<=.1&&mean(bs$valid)>=.98
report<-list(started_utc=start,completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),plan_version='0.6',
 formula=paste(deparse(f),collapse=' '),estimator='brglm2::brglmFit,AS_mean',package_version=as.character(packageVersion('brglm2')),
 nonintercept_df=fit$fit$rank-1,converged=fit$fit$converged,max_abs_SMD_included=maxsmd,
 eMAR_SMD=bal[term=='any_emar_available_before']$SMD_after,arms=arm,
 bootstrap=list(B=B,valid=sum(bs$valid),failed=sum(!bs$valid),nonconverged=sum(!bs$converged),draws_with_warnings=sum(bs$warnings>0),
  ESS_LZD_quantiles=unname(quantile(bs$ESS_LZD,c(.025,.5,.975),na.rm=TRUE))),
 outcome_free_reporting_gate=pass,gate_scope='Allows exploratory weighted association only; not proof of exchangeability, completeness or causal validity.',
 limits=c('Bias reduction prevents infinite estimates under supported conditions; it does not remove separation in the observed design or add treatment information.',
 'Median fills and missingness indicators do not solve missing-data bias; fills are recomputed in each bootstrap.',
 'No outcomes loaded during this model selection and diagnostic step. Same resamples later support joint endpoint uncertainty.'))
write_json(report,file.path(R,'reports/BIAS_REDUCED_OVERLAP_v0.6.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
fwrite(bal,file.path(R,'reports/COVARIATE_BALANCE_v0.6.csv'));fwrite(bs,file.path(R,'cache/internal_reports/phase6_bootstrap_diagnostics.csv'))
saveRDS(list(raw=raw,data=x,fit=fit$fit,weights=fit$w,bootstrap=boot,report=report),file.path(R,'cache/phase6_overlap.rds'))
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.6.txt'))
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE,digits=8),'\n')
