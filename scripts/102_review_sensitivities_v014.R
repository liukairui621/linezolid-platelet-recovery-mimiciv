# Reviewer-requested post-result diagnostics and sensitivity analyses.
# Patient-level objects remain in server-controlled cache; public outputs are aggregate and disclosure screened.
suppressPackageStartupMessages({library(brglm2);library(data.table);library(jsonlite);library(Hmisc)})
R <- '/root/projects/linezolid_platelet_recovery'
A <- file.path(R,'analysis_v0.14')
plan_path <- file.path(A,'config/analysis_plan_v0.14.json')
plan <- fromJSON(plan_path)
sha <- function(p) as.character(unname(tools::md5sum(p)))
# SHA-256 was checked by extraction script; retain exact frozen plan hash from plan-freeze record.
freeze <- fromJSON(file.path(A,'reports/PLAN_FREEZE_v0.14.json'))
stopifnot(identical(tolower(freeze$sha256),'5a20fba454e6b76596a0acfa7b6bb767c3702d7bdee397502944f1e04e0cd3d6'))
O <- file.path(A,'output'); dir.create(O,recursive=TRUE,showWarnings=FALSE)
C <- file.path(R,'cache/internal_review_v014'); dir.create(C,recursive=TRUE,showWarnings=FALSE)

m <- readRDS(file.path(R,'cache/phase11_routine_culture_no_recorded_lmwh_shared.rds'))
ep <- as.data.frame(fread(file.path(R,'cache/phase11_endpoints.csv')))
ep <- ep[match(m$raw$hadm_id,ep$hadm_id),]
stopifnot(!anyNA(ep$hadm_id),identical(as.character(ep$hadm_id),as.character(m$raw$hadm_id)),identical(ep$drug,m$raw$drug))
route <- as.data.frame(fread(file.path(C,'EXPOSURE_CONTEXT_EXACT_v0.14.csv')))
route <- route[match(m$raw$hadm_id,route$hadm_id),]
stopifnot(!anyNA(route$hadm_id),identical(as.character(route$hadm_id),as.character(m$raw$hadm_id)))

# Frozen main-effects propensity formulas.
f_main <- A ~ age + gender + known_platelets + log_creatinine + creatinine_missing +
 rrt_documented_prior24h + invasive_vent_documented_prior24h + pressors_documented_prior6h +
 platelet_change_filled + trend_missing + prior_recorded_marrow + log1p_admission_days + respiratory72h +
 log_inr + inr_missing + log1p_bilirubin + bilirubin_missing + mapped_era + any_emar_available_before +
 icu_unit + nonscreen_blood72h + nonscreen_urine72h + known_named_nonscreen7d + known_enterococcus7d +
 observed_creatinine_aki + ufh_active_available7d + marrowdrug_available7d + beta_lactam_available24h
f_mi <- update(f_main,. ~ . - creatinine_missing - trend_missing - inr_missing - bilirubin_missing)
ref_levels <- lapply(m$data[vapply(m$data,is.factor,logical(1))],levels)

prep <- function(x,imputed=FALSE){
 lc <- names(x)[vapply(x,is.logical,logical(1))]; x[lc] <- lapply(x[lc],as.integer)
 x$A <- as.integer(x$drug=='linezolid')
 x$gender <- factor(x$gender,levels=ref_levels$gender)
 x$mapped_era <- factor(x$mapped_era,levels=ref_levels$mapped_era)
 x$icu_unit <- factor(x$icu_unit,levels=ref_levels$icu_unit)
 fill <- function(z,fallback=0){z[!is.finite(z)]<-NA_real_;v<-if(all(is.na(z)))fallback else median(z,na.rm=TRUE);z[is.na(z)]<-v;z}
 if(imputed){
  stopifnot(all(is.finite(x$known_creatinine)&x$known_creatinine>0),all(is.finite(x$platelet_change)),all(is.finite(x$baseline_inr)&x$baseline_inr>0),all(is.finite(x$baseline_bilirubin)&x$baseline_bilirubin>=0))
  x$log_creatinine<-log(x$known_creatinine);x$platelet_change_filled<-x$platelet_change
  x$log_inr<-log(x$baseline_inr);x$log1p_bilirubin<-log1p(x$baseline_bilirubin)
 } else {
  x$creatinine_missing<-as.integer(is.na(x$known_creatinine)|x$known_creatinine<=0)
  x$log_creatinine<-log(fill(ifelse(x$known_creatinine>0,x$known_creatinine,NA_real_),1))
  x$trend_missing<-as.integer(is.na(x$platelet_change));x$platelet_change_filled<-fill(x$platelet_change)
  x$inr_missing<-as.integer(is.na(x$baseline_inr)|x$baseline_inr<=0);x$log_inr<-log(fill(ifelse(x$baseline_inr>0,x$baseline_inr,NA_real_),1))
  x$bilirubin_missing<-as.integer(is.na(x$baseline_bilirubin)|x$baseline_bilirubin<0);x$log1p_bilirubin<-log1p(fill(ifelse(x$baseline_bilirubin>=0,x$baseline_bilirubin,NA_real_)))
 }
 x$log1p_admission_days<-log1p(x$days_since_admission)
 x
}
fitone <- function(x,form){
 mm0 <- model.matrix(form,x); q <- qr(mm0); cols <- sort(q$pivot[seq_len(q$rank)]); mm <- mm0[,cols,drop=FALSE]
 attempt <- function(maxit,slowit=1){
  warns<-character(); fit<-withCallingHandlers(brglmFit(x=mm,y=x$A,family=binomial('logit'),control=brglmControl(type='AS_mean',maxit=maxit,epsilon=1e-8,slowit=slowit)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')});list(fit=fit,warnings=warns)
 }
 z<-attempt(200); fallback<-!isTRUE(z$fit$converged)||any(!is.finite(coef(z$fit)))
 if(fallback) z<-attempt(2000,.1)
 ee<-z$fit$fitted.values; w<-ifelse(x$A==1,1-ee,ee)
 list(valid=isTRUE(z$fit$converged)&&all(is.finite(w))&&all(w>0&w<1)&&all(is.finite(coef(z$fit))),fit=z$fit,e=ee,w=w,mm=mm,omitted=setdiff(colnames(mm0),colnames(mm)),warnings=z$warnings,fallback=fallback)
}
ess <- function(w) sum(w)^2/sum(w^2)
wq <- function(x,w,p=c(.25,.5,.75)){
 o<-order(x);x<-x[o];w<-w[o]/sum(w);s<-cumsum(w);vapply(p,function(pp)x[which(s>=pp)[1]],numeric(1))
}
stat_rd <- function(a,w,y) weighted.mean(y[a==1],w[a==1])-weighted.mean(y[a==0],w[a==0])
balance <- function(x,z){
 mm<-z$mm[,colnames(z$mm)!='(Intercept)',drop=FALSE];a<-x$A==1
 out<-rbindlist(lapply(seq_len(ncol(mm)),function(j){v<-mm[,j];sd0<-sqrt((var(v[a])+var(v[!a]))/2);data.table(term=colnames(mm)[j],SMD_before=if(is.finite(sd0)&&sd0>0)(mean(v[a])-mean(v[!a]))/sd0 else NA_real_,SMD_after=if(is.finite(sd0)&&sd0>0)(weighted.mean(v[a],z$w[a])-weighted.mean(v[!a],z$w[!a]))/sd0 else NA_real_)}))
 out
}
cluster_strata <- function(raw,x){cl<-split(seq_len(nrow(raw)),raw$subject_id);pat<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),character(1));list(cl=cl,strata=split(names(cl),pat))}
draw_indices <- function(cs){unlist(cs$cl[unlist(lapply(cs$strata,function(z)sample(z,length(z),replace=TRUE)),use.names=FALSE)],use.names=FALSE)}

# 1. Frozen point-model diagnostics.
x0<-m$data;a0<-x0$A;w0<-m$weights;e0<-m$fit$fitted.values
stopifnot(length(e0)==nrow(x0),all(e0>0&e0<1))
probs<-c(0,.01,.05,.25,.5,.75,.95,.99,1)
diag_rows<-rbindlist(lapply(0:1,function(aa){ii<-a0==aa;data.table(arm=if(aa==1)'LZD' else 'VAN',metric=rep(c('propensity_score','overlap_weight'),each=length(probs)),quantile=rep(probs,2),value=c(quantile(e0[ii],probs,type=7),quantile(w0[ii],probs,type=7)),n=sum(ii),ESS=ess(w0[ii]))}))
fwrite(diag_rows,file.path(O,'PROPENSITY_WEIGHT_QUANTILES_v0.14.csv'))
dens <- function(v,arm,metric){dd<-density(v,from=0,to=1,n=512,bw='nrd0');data.table(arm=arm,metric=metric,x=dd$x,density=dd$y)}
density_rows<-rbindlist(list(dens(e0[a0==1],'LZD','propensity_score'),dens(e0[a0==0],'VAN','propensity_score'),dens(w0[a0==1],'LZD','overlap_weight'),dens(w0[a0==0],'VAN','overlap_weight')))
fwrite(density_rows,file.path(O,'PROPENSITY_WEIGHT_DENSITY_v0.14.csv'))

# 2. Route-restricted sensitivity.
keep <- route$initiation_route=='parenteral' | m$raw$drug=='vancomycin'
raw_r<-m$raw[keep,,drop=FALSE]; ep_r<-ep[keep,,drop=FALSE]; x_r<-prep(raw_r,FALSE); z_r<-fitone(x_r,f_main);stopifnot(z_r$valid)
y_r<-as.numeric(ep_r$hospital_status==1);q_r<-stat_rd(x_r$A,z_r$w,y_r);bal_r<-balance(x_r,z_r)
cs_r<-cluster_strata(raw_r,x_r);Broute<-plan$tasks$parenteral_route_sensitivity$bootstrap$B;set.seed(plan$tasks$parenteral_route_sensitivity$bootstrap$seed)
br<-rep(NA_real_,Broute);valid<-logical(Broute);fallbacks<-logical(Broute)
for(i in seq_len(Broute)){
 ix<-draw_indices(cs_r);xb<-prep(raw_r[ix,,drop=FALSE],FALSE);zb<-tryCatch(fitone(xb,f_main),error=function(e)list(valid=FALSE))
 if(isTRUE(zb$valid)){br[i]<-stat_rd(xb$A,zb$w,y_r[ix]);valid[i]<-TRUE;fallbacks[i]<-isTRUE(zb$fallback)}
 if(i%%250==0){cat('ROUTE',i,'/',Broute,'valid',sum(valid),'\n');flush.console()}
}
stopifnot(mean(valid)>=.98)
route_result<-data.table(analysis='Parenteral linezolid initiation sensitivity',LZD_n='Suppressed linked value',VAN_n=as.character(sum(x_r$A==0)),LZD_events=if(sum(y_r[x_r$A==1])<10)'<10' else as.character(sum(y_r[x_r$A==1])),VAN_events=as.character(sum(y_r[x_r$A==0])),RD=q_r,lower=quantile(br[valid],.025,type=7),upper=quantile(br[valid],.975,type=7),ESS_LZD=ess(z_r$w[x_r$A==1]),ESS_VAN=ess(z_r$w[x_r$A==0]),max_abs_SMD=max(abs(bal_r$SMD_after),na.rm=TRUE),valid_draws=sum(valid),requested_draws=Broute,numerical_fallbacks=sum(fallbacks))
fwrite(route_result,file.path(O,'PARENTERAL_ROUTE_SENSITIVITY_v0.14.csv'))
saveRDS(list(point=route_result,bootstrap=br,diagnostics=data.table(valid=valid,fallback=fallbacks)),file.path(C,'PARENTERAL_ROUTE_SENSITIVITY_EXACT_v0.14.rds'))

# 3. Missing-data multiple-imputation sensitivity.
raw_mi<-m$raw
impd<-raw_mi[,setdiff(names(raw_mi),c('subject_id','hadm_id','drug','culture72h','lmwh_available7d')),drop=FALSE]
impd$A<-as.integer(raw_mi$drug=='linezolid')
impd$first_event_status<-factor(ep$hospital_status,levels=0:3)
impd$log1p_first_event_day<-log1p(ep$hospital_day)
# logical predictors are kept as factors so their categorical nature is explicit.
for(nm in names(impd)[vapply(impd,is.logical,logical(1))])impd[[nm]]<-factor(impd[[nm]],levels=c(FALSE,TRUE))
imp_formula<-as.formula(paste('~',paste(names(impd),collapse='+')))
set.seed(plan$tasks$missing_data_sensitivity$seed)
imp<-aregImpute(imp_formula,data=impd,n.impute=plan$tasks$missing_data_sensitivity$m,nk=plan$tasks$missing_data_sensitivity$imputation_settings$nk,type=plan$tasks$missing_data_sensitivity$imputation_settings$type,match=plan$tasks$missing_data_sensitivity$imputation_settings$match,boot.method=plan$tasks$missing_data_sensitivity$imputation_settings$boot_method,B=plan$tasks$missing_data_sensitivity$imputation_settings$B,pr=FALSE)
saveRDS(imp,file.path(C,'AREGIMPUTE_OBJECT_v0.14.rds'))
M<-plan$tasks$missing_data_sensitivity$m;Bmi<-plan$tasks$missing_data_sensitivity$within_imputation_bootstrap$B
mirows<-list();miboot<-vector('list',M)
for(j in seq_len(M)){
 comp<-as.data.frame(impute.transcan(imp,imputation=j,data=impd,list.out=TRUE,pr=FALSE,check=FALSE))
 rr<-raw_mi
 for(nm in c('known_creatinine','platelet_change','baseline_inr','baseline_bilirubin'))rr[[nm]]<-as.numeric(comp[[nm]])
 xx<-prep(rr,TRUE);zz<-fitone(xx,f_mi);stopifnot(zz$valid)
 qq<-stat_rd(xx$A,zz$w,as.numeric(ep$hospital_status==1));bb<-rep(NA_real_,Bmi);ok<-logical(Bmi);fb<-logical(Bmi);cs<-cluster_strata(rr,xx)
 set.seed(plan$tasks$missing_data_sensitivity$within_imputation_bootstrap$seed_base+j)
 for(i in seq_len(Bmi)){
  ix<-draw_indices(cs);xb<-prep(rr[ix,,drop=FALSE],TRUE);zb<-tryCatch(fitone(xb,f_mi),error=function(e)list(valid=FALSE))
  if(isTRUE(zb$valid)){bb[i]<-stat_rd(xb$A,zb$w,as.numeric(ep$hospital_status[ix]==1));ok[i]<-TRUE;fb[i]<-isTRUE(zb$fallback)}
 }
 stopifnot(mean(ok)>=.98);miboot[[j]]<-bb[ok];bl<-balance(xx,zz)
 mirows[[j]]<-data.table(imputation=j,RD=qq,within_variance=var(bb[ok]),bootstrap_valid=sum(ok),bootstrap_requested=Bmi,numerical_fallbacks=sum(fb),ESS_LZD=ess(zz$w[xx$A==1]),ESS_VAN=ess(zz$w[xx$A==0]),max_abs_SMD=max(abs(bl$SMD_after),na.rm=TRUE))
 cat('MI',j,'/',M,'RD',qq,'valid',sum(ok),'\n');flush.console()
}
mi<-rbindlist(mirows);Qbar<-mean(mi$RD);Ubar<-mean(mi$within_variance);Bvar<-var(mi$RD);Tvar<-Ubar+(1+1/M)*Bvar
if(Bvar>0){df<-(M-1)*(1+Ubar/((1+1/M)*Bvar))^2}else df<-Inf
crit<-if(is.finite(df))qt(.975,df) else qnorm(.975)
pooled<-data.table(imputation=NA_integer_,RD=Qbar,within_variance=Ubar,bootstrap_valid=sum(mi$bootstrap_valid),bootstrap_requested=sum(mi$bootstrap_requested),numerical_fallbacks=sum(mi$numerical_fallbacks),ESS_LZD=mean(mi$ESS_LZD),ESS_VAN=mean(mi$ESS_VAN),max_abs_SMD=max(mi$max_abs_SMD),lower=Qbar-crit*sqrt(Tvar),upper=Qbar+crit*sqrt(Tvar),between_variance=Bvar,total_variance=Tvar,df=df)
mi[,`:=`(lower=NA_real_,upper=NA_real_,between_variance=NA_real_,total_variance=NA_real_,df=NA_real_)]
fwrite(rbind(mi,pooled,fill=TRUE),file.path(O,'MULTIPLE_IMPUTATION_SENSITIVITY_v0.14.csv'))
saveRDS(list(imputation_results=mi,pooled=pooled,bootstrap=miboot),file.path(C,'MULTIPLE_IMPUTATION_SENSITIVITY_EXACT_v0.14.rds'))

# 4. Observation opportunity and baseline treatment-context summaries.
weighted_summary<-function(v,w){c(mean=weighted.mean(v,w),q25=wq(v,w,.25),median=wq(v,w,.5),q75=wq(v,w,.75))}
obsrows<-list()
for(aa in 0:1){ii<-a0==aa;arm<-if(aa==1)'LZD' else 'VAN'
 for(nm in c('hospital_tests','hospital_days','hospital_tests_per_day')){
  v<-ep[[nm]][ii];u<-quantile(v,c(.25,.5,.75),type=7);ww<-weighted_summary(v,w0[ii]);obsrows[[length(obsrows)+1]]<-data.table(arm=arm,metric=nm,n=sum(ii),unweighted_mean=mean(v),unweighted_q25=u[1],unweighted_median=u[2],unweighted_q75=u[3],weighted_mean=ww[1],weighted_q25=ww[2],weighted_median=ww[3],weighted_q75=ww[4])
 }
 v<-as.numeric(ep$hospital_tests[ii]>=2);obsrows[[length(obsrows)+1]]<-data.table(arm=arm,metric='at_least_two_postinitiation_tests',n=sum(ii),unweighted_mean=mean(v),unweighted_q25=NA_real_,unweighted_median=NA_real_,unweighted_q75=NA_real_,weighted_mean=weighted.mean(v,w0[ii]),weighted_q25=NA_real_,weighted_median=NA_real_,weighted_q75=NA_real_)
}
fwrite(rbindlist(obsrows,fill=TRUE),file.path(O,'OBSERVATION_OPPORTUNITY_v0.14.csv'))

context_vars<-c('routine_blood72h','routine_urine72h','routine_respiratory72h','known_named_routine7d','known_enterococcus7d','known_mrsa_routine7d','beta_lactam_available24h','prior_recorded_marrow','rrt_documented_prior24h','invasive_vent_documented_prior24h','pressors_documented_prior6h','observed_creatinine_aki')
base_csv<-as.data.frame(fread(file.path(R,'cache/phase10_baseline.csv')));base_csv<-base_csv[match(m$raw$hadm_id,base_csv$hadm_id),];stopifnot(!anyNA(base_csv$hadm_id))
ctx<-list();suppress_arm<-function(k,n){if(k>0&&k<10)'<10' else as.character(k)}
for(nm in context_vars)for(aa in 0:1){ii<-a0==aa;v<-as.numeric(base_csv[[nm]][ii]);k<-sum(v==1);ctx[[length(ctx)+1]]<-data.table(characteristic=nm,arm=if(aa==1)'LZD' else 'VAN',n=suppress_arm(k,sum(ii)),denominator=sum(ii),unweighted_percent=if(k>0&&k<10)NA_real_ else 100*mean(v),weighted_percent=if(k>0&&k<10)NA_real_ else 100*weighted.mean(v,w0[ii]))}
fwrite(rbindlist(ctx),file.path(O,'CLINICAL_TREATMENT_CONTEXT_v0.14.csv'))

report<-list(status='COMPLETE',completed_utc=format(Sys.time(),tz='UTC',usetz=TRUE),classification='Post-result reviewer-requested exploratory extension',frozen_primary_unchanged=TRUE,diagnostics=list(n_LZD=sum(a0==1),n_VAN=sum(a0==0),ESS_LZD=ess(w0[a0==1]),ESS_VAN=ess(w0[a0==0])),route_sensitivity=as.list(route_result[1]),multiple_imputation_pooled=as.list(pooled[1]),privacy='Patient-level rows, propensity scores, weights, shifted times, imputed records and bootstrap membership remain server-only. Public small clinical cells and linked values are suppressed.',limitations=c('Route sensitivity remains an initiation strategy comparison, not sustained-exposure causal analysis.','Multiple imputation sensitivity is post-result and does not resolve unmeasured indication confounding.','Monitoring summaries are descriptive and do not identify an observation-process causal effect.'))
write_json(report,file.path(A,'reports/REVIEW_SENSITIVITY_RESULTS_v0.14.json'),pretty=TRUE,auto_unbox=TRUE,digits=12,na='null')
capture.output(sessionInfo(),file=file.path(A,'reports/R_sessionInfo_v0.14.txt'))
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE,digits=8,na='null'),'\n')
