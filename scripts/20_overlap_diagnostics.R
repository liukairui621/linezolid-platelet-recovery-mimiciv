# Baseline-only design diagnostics under v0.4. No outcomes or treatment effects.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(ggplot2)})
set.seed(42)
root <- '/root/projects/linezolid_platelet_recovery'
plan <- fromJSON(file.path(root,'config/analysis_plan_v0.4.json'))
stopifnot(plan$version=='0.4')
d <- as.data.frame(fread(file.path(root,'cache/phase4_baseline_only.csv')))
stopifnot(!any(c('deathtime','dischtime','n_platelet_14d','observed_treatment_days') %in% names(d)))
figdir <- file.path(root,'figures/v0.4');dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(root,'cache/internal_reports'),showWarnings=FALSE)
logical_cols <- names(d)[vapply(d,is.logical,logical(1))]
d[logical_cols] <- lapply(d[logical_cols],as.integer)
d$A <- as.integer(d$drug=='linezolid')
d$gender <- factor(d$gender);d$anchor_year_group <- factor(d$anchor_year_group)
fill <- function(x) {x[!is.finite(x)]<-NA_real_;x[is.na(x)]<-median(x,na.rm=TRUE);x}
d$creatinine_missing <- as.integer(is.na(d$known_creatinine)|d$known_creatinine<=0)
d$log_creatinine <- log(fill(ifelse(d$known_creatinine>0,d$known_creatinine,NA_real_)))
d$trend_missing <- as.integer(is.na(d$platelet_change));d$platelet_change_filled<-fill(d$platelet_change)
d$inr_missing<-as.integer(is.na(d$baseline_inr)|d$baseline_inr<=0)
d$log_inr<-log(fill(ifelse(d$baseline_inr>0,d$baseline_inr,NA_real_)))
d$bilirubin_missing<-as.integer(is.na(d$baseline_bilirubin)|d$baseline_bilirubin<0)
d$log1p_bilirubin<-log1p(fill(ifelse(d$baseline_bilirubin>=0,d$baseline_bilirubin,NA_real_)))
stopifnot(all(d$days_since_admission>=0))
d$log1p_admission_days<-log1p(d$days_since_admission)
core <- c('age','gender','known_platelets','log_creatinine','creatinine_missing',
 'rrt_documented_prior24h','invasive_vent_documented_prior24h','pressors_documented_prior6h',
 'anchor_year_group','platelet_change_filled','trend_missing','prior_recorded_marrow')
extra <- c('log1p_admission_days','culture72h','respiratory72h','log_inr','inr_missing',
 'log1p_bilirubin','bilirubin_missing','decision_year_minus_anchor_year')
audit_terms <- unique(c(core,extra,'ufh_documented_prior7d','lmwh_documented_prior7d','selected_other_drug_prior7d','any_emar_available_before'))
getformula<-function(terms,x) reformulate(terms[vapply(terms,function(v) length(unique(x[[v]]))>1,logical(1))],response='A')
ess<-function(w) sum(w)^2/sum(w^2)
wmean<-function(x,w) sum(x*w)/sum(w)
wvar<-function(x,w) sum(w*(x-wmean(x,w))^2)/sum(w)
wks<-function(x,g,w) {
 ok<-is.finite(x);x<-x[ok];g<-g[ok];w<-w[ok]
 z<-data.table(x=x,g=g,w=w)[,.(w=sum(w)),by=.(x,g)]
 grid<-sort(unique(x));f<-lapply(0:1,function(a){zz<-z[g==a][order(x)];v<-c(0,cumsum(zz$w)/sum(zz$w));v[findInterval(grid,zz$x)+1]})
 max(abs(f[[1]]-f[[2]]))
}
model_report<-list();balance_all<-list();curves<-list();model_objects<-list()
run_model<-function(name,x,terms=NULL,spline=FALSE) {
 x<-droplevels(x)
 f<-if(spline) A~splines::ns(log_creatinine,df=3) else getformula(terms,x)
 warns<-character()
 fit<-withCallingHandlers(glm(f,data=x,family=binomial(),control=glm.control(maxit=100,epsilon=1e-10)),
   warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 e<-fitted(fit);w<-ifelse(x$A==1,1-e,e)
 stopifnot(all(is.finite(w)),all(w>=0),all(w<=1))
 mm<-model.matrix(fit)
 balance_matrix<-model.matrix(getformula(audit_terms,x),data=x)[,-1,drop=FALSE]
 extra_columns<-setdiff(colnames(mm),c('(Intercept)',colnames(balance_matrix)))
 if(length(extra_columns)) balance_matrix<-cbind(balance_matrix,mm[,extra_columns,drop=FALSE])
 bal<-rbindlist(lapply(colnames(balance_matrix),function(nm) {
   y<-balance_matrix[,nm];a<-x$A==1;b<-!a
   sd0<-sqrt((var(y[a])+var(y[b]))/2)
   if(!is.finite(sd0)||sd0==0) return(NULL)
   data.table(model=name,term=nm,in_fitted_matrix=nm %in% colnames(mm),
    SMD_before=(mean(y[a])-mean(y[b]))/sd0,
    SMD_after=(wmean(y[a],w[a])-wmean(y[b],w[b]))/sd0)
 }))
 max_balance<-max(abs(bal$SMD_after[bal$in_fitted_matrix]),na.rm=TRUE)
 # Actual score-equation check, not an assertion that omitted variables are balanced.
 score_error<-max(abs(drop(crossprod(mm,x$A-e))))
 if(fit$converged && all(is.finite(coef(fit)))) stopifnot(score_error<1e-4)
 perarm<-lapply(0:1,function(a) {
   ok<-x$A==a; ww<-w[ok]; cp<-rowsum(ww,x$subject_id[ok],reorder=FALSE)[,1]
   cr<-x$known_creatinine[ok];valid<-is.finite(cr)&cr>0
   list(drug=if(a==1)'linezolid' else 'vancomycin',n=sum(ok),patients=length(unique(x$subject_id[ok])),
    ESS_admission=ess(ww),ESS_patient_weight_concentration=ess(cp),
    weight_sum=sum(ww),maximum_weight=max(ww),largest_normalized_weight=max(ww)/sum(ww),
    PS_quantiles=unname(quantile(e[ok],c(0,.025,.25,.5,.75,.975,1))),
    weighted_mean_creatinine=wmean(cr[valid],ww[valid]),
    weighted_logcreatinine_mean=wmean(x$log_creatinine[ok],ww),
    weighted_fraction_creatinine_1_5_to_2_5=wmean(as.numeric(cr[valid]>=1.5&cr[valid]<=2.5),ww[valid]))
 })
 report<-list(model=name,formula=paste(deparse(f),collapse=' '),converged=fit$converged,
  matrix_columns=ncol(mm),rank=fit$rank,nonintercept_df=fit$rank-1,
  exposure_count_per_nonintercept_df=sum(x$A)/(fit$rank-1),warnings=warns,
  aliased_coefficients=sum(is.na(coef(fit))),score_equation_max_abs_error=score_error,
  max_abs_SMD_included_terms=max_balance,creatinine_weighted_KS=wks(x$known_creatinine,x$A,w),arms=perarm,
  interpretation='Design diagnostic only. Kish ESS reflects weight concentration, not outcome-specific information or proof of exchangeability.')
 model_report[[name]]<<-report;balance_all[[name]]<<-bal
 model_objects[[name]]<<-list(fit=fit,data=x,weights=w)
 for(a in 0:1) {
   ix<-x$A==a; den<-density(e[ix],weights=w[ix]/sum(w[ix]),from=0,to=max(e)*1.05,n=512)
   curves[[paste(name,a)]]<<-data.table(model=name,drug=if(a==1)'Linezolid' else 'Vancomycin',x=den$x,density=den$y)
 }
}
renal <- d[!d$creatinine_missing,]
run_model('renal_log_only',renal,'log_creatinine')
run_model('renal_spline_3df',renal,spline=TRUE)
run_model('core_baseline',d,core)
run_model('extended_baseline',d,c(core,extra))
run_model('baseline_emar_covered',d[d$any_emar_available_before==1,],
 c(core,extra,'ufh_documented_prior7d','lmwh_documented_prior7d','selected_other_drug_prior7d'))

van<-renal$known_creatinine[renal$A==0];lzd<-renal$known_creatinine[renal$A==1]
q<-unname(quantile(van,c(.025,.975)));rng<-range(van)
renal_report<-list(comparator_central95_creatinine=q,
 linezolid_observed_n=length(lzd),linezolid_inside_central95=sum(lzd>=q[1]&lzd<=q[2]),
 linezolid_outside_central95=sum(lzd<q[1]|lzd>q[2]),
 linezolid_outside_comparator_full_range=sum(lzd<rng[1]|lzd>rng[2]),
 note='Central95 is not a positivity boundary. No trimming performed. Univariate ESS is not a bound on multivariable ESS.')
report<-list(plan_version=plan$version,created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
 renal=renal_report,models=model_report,
 limitations=c('No clinical outcomes or effects computed. Shared infection indication remains unvalidated.',
 'Median filling plus missingness indicators is a provisional design diagnostic, not proof of absence of missing-data bias.',
 'Finite-sample model complexity and positivity assessed by actual matrices; no automatic6-coefficient rule.',
 'The eMAR subset is baseline-record availability, not guaranteed complete drug ascertainment.',
 'Patient-level weights/models and exact small counts remain server-only.'))
write_json(report,file.path(root,'cache/internal_reports/OVERLAP_DIAGNOSTICS_v0.4_exact.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
saveRDS(model_objects,file.path(root,'cache/phase4_overlap_models.rds'))
pub<-report
hide<-function(x) if(x>0&&x<10) '<10' else x
for(nm in c('linezolid_observed_n','linezolid_inside_central95','linezolid_outside_central95','linezolid_outside_comparator_full_range')) pub$renal[[nm]]<-hide(pub$renal[[nm]])
for(nm in names(pub$models)) for(j in seq_along(pub$models[[nm]]$arms)) {
  a<-pub$models[[nm]]$arms[[j]];a$n<-hide(a$n);a$patients<-hide(a$patients);pub$models[[nm]]$arms[[j]]<-a
}
write_json(pub,file.path(root,'reports/OVERLAP_DIAGNOSTICS_v0.4.json'),pretty=TRUE,auto_unbox=TRUE,digits=8)
bal<-rbindlist(balance_all);fwrite(bal,file.path(root,'reports/COVARIATE_BALANCE_v0.4.csv'))

pal<-c(Linezolid='#B44E34',Vancomycin='#257A91')
theme_set(theme_classic(base_size=12)+theme(legend.position='top',plot.title.position='plot'))
# Smooth density and ECDF contain no patient IDs or individual rug marks.
dens<-rbindlist(lapply(0:1,function(a){de<-density(log(renal$known_creatinine[renal$A==a]),n=512)
 data.table(drug=if(a==1)'Linezolid' else 'Vancomycin',log_creatinine=de$x,density=de$y)}))
fwrite(dens,file.path(root,'reports/RENAL_DENSITY_v0.4.csv'))
p<-ggplot(dens,aes(exp(log_creatinine),density,color=drug))+geom_line(linewidth=1)+
 geom_vline(xintercept=q,linetype=3,color='grey45')+scale_x_log10(breaks=c(.25,.5,1,2,4,8,16))+
 scale_color_manual(values=pal)+labs(x='Baseline creatinine (mg/dL; log scale)',y='Density on log-creatinine scale',color=NULL,
 title='Baseline renal overlap',subtitle='Dotted lines: vancomycin empirical 2.5th and 97.5th percentiles',
 caption='Descriptive central range, not an exclusion rule. Positive observed values only.')
ggsave(file.path(figdir,'FigD1A_renal_density.png'),p,width=7.4,height=4.7,dpi=180)
ggsave(file.path(figdir,'FigD1A_renal_density.pdf'),p,width=7.4,height=4.7)
grid<-seq(log(.1),log(30),length.out=401)
ec<-rbindlist(lapply(0:1,function(a)data.table(drug=if(a==1)'Linezolid' else 'Vancomycin',creatinine=exp(grid),ECDF=ecdf(log(renal$known_creatinine[renal$A==a]))(grid))))
fwrite(ec,file.path(root,'reports/RENAL_ECDF_v0.4.csv'))
p<-ggplot(ec,aes(creatinine,ECDF,color=drug))+geom_step(linewidth=1)+scale_x_log10(breaks=c(.25,.5,1,2,4,8,16))+scale_color_manual(values=pal)+
 labs(x='Baseline creatinine (mg/dL; log scale)',y='Cumulative fraction',color=NULL,title='Observed creatinine distributions')
ggsave(file.path(figdir,'FigD1B_renal_ecdf.png'),p,width=7.4,height=4.7,dpi=180)
ggsave(file.path(figdir,'FigD1B_renal_ecdf.pdf'),p,width=7.4,height=4.7)
b<-bal[model=='extended_baseline'];bl<-melt(b,id.vars=c('model','term','in_fitted_matrix'),measure.vars=c('SMD_before','SMD_after'),variable.name='stage',value.name='SMD')
p<-ggplot(bl,aes(abs(SMD),reorder(term,abs(SMD)),color=stage,shape=in_fitted_matrix))+geom_point(size=2.1)+geom_vline(xintercept=.1,linetype=3,color='grey50')+
 scale_color_manual(values=c(SMD_before='#B44E34',SMD_after='#257A91'),labels=c('Before weighting','After overlap weighting'))+
 labs(x='Absolute standardized mean difference',y=NULL,color=NULL,shape='Included term',title='Extended baseline model: measured balance',
 caption='Exact mean balance on included terms does not establish clinical exchangeability.')+theme(axis.text.y=element_text(size=8))
ggsave(file.path(figdir,'FigD1C_balance.png'),p,width=9.3,height=9.3,dpi=160)
ggsave(file.path(figdir,'FigD1C_balance.pdf'),p,width=9.3,height=9.3)
capture.output(sessionInfo(),file=file.path(root,'reports/R_sessionInfo_v0.4.txt'))
cat(toJSON(list(renal=pub$renal,model_summary=lapply(pub$models,function(x)x[c('model','converged','nonintercept_df','max_abs_SMD_included_terms','arms')])),pretty=TRUE,auto_unbox=TRUE,digits=5))
