# Descriptive submission tables from frozen v0.11 inputs. No model/outcome changes.
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
R <- '/root/projects/linezolid_platelet_recovery'
O <- file.path(R,'reports/submission_v0.12');dir.create(O,recursive=TRUE,showWarnings=FALSE)
pops <- c('routine_culture_no_recorded_lmwh_shared','routine_bacterial_culture72h_shared')
bal <- fread(file.path(R,'reports/READINESS_BALANCE_v0.10.csv'))
continuous <- c(age='Age, years',known_platelets='Platelet count, 10^9/L',known_creatinine='Creatinine, mg/dL',platelet_change='Platelet change, 10^9/L',baseline_inr='International normalized ratio',baseline_bilirubin='Bilirubin, mg/dL',days_since_admission='Admission to initiation, days')
binary <- c(rrt_documented_prior24h='Recorded renal replacement, 24 h',invasive_vent_documented_prior24h='Recorded invasive ventilation, 24 h',pressors_documented_prior6h='Recorded vasopressor, 6 h',prior_recorded_marrow='Prior recorded marrow disorder',respiratory72h='Routine respiratory culture, 72 h',nonscreen_blood72h='Routine blood culture, 72 h',nonscreen_urine72h='Routine urine culture, 72 h',known_named_nonscreen7d='Known named routine organism, 7 d',known_enterococcus7d='Known enterococcus, 7 d',observed_creatinine_aki='Observed creatinine AKI criteria',ufh_active_available7d='Recorded unfractionated heparin, 7 d',lmwh_available7d='Recorded low-molecular-weight heparin, 7 d',marrowdrug_available7d='Selected marrow-related medication, 7 d',beta_lactam_available24h='Recorded intravenous beta-lactam, 24 h',any_emar_available_before='Prior electronic medication record')
smdmap <- c(age='age',known_platelets='known_platelets',known_creatinine='log_creatinine',platelet_change='platelet_change_filled',baseline_inr='log_inr',baseline_bilirubin='log1p_bilirubin',days_since_admission='log1p_admission_days')
hide_n <- function(n) if(n>0&&n<10) '<10' else as.character(n)
wsd <- function(v,w) {mu<-weighted.mean(v,w);sqrt(sum(w*(v-mu)^2)/sum(w))}
rows <- list();miss <- list();meta<-list()
for(pop in pops){
 m<-readRDS(file.path(R,'cache',paste0('phase11_',pop,'.rds')));x<-m$raw;w<-m$weights;a<-x$drug=='linezolid';stopifnot(length(w)==nrow(x),!anyNA(w),all(w>0),!anyDuplicated(x$hadm_id))
 add<-function(key,label,v,is_binary=FALSE,term=key){
  z<-list(population=pop,key=key,label=label,model_term=term)
  for(g in c('LZD','VAN')){ii<-if(g=='LZD')a else !a;ok<-ii&is.finite(v);n<-sum(ok)
   if(is_binary){k<-sum(v[ok]==1);supp<-k>0&&k<10;z[[paste0(g,'_before')]]<-if(supp)'<10 (suppressed)' else sprintf('%d (%.1f)',k,100*mean(v[ok]));z[[paste0(g,'_after')]]<-if(supp)'Suppressed' else sprintf('%.1f',100*weighted.mean(v[ok],w[ok]))
   }else{z[[paste0(g,'_before')]]<-if(n>=10)sprintf('%.1f (%.1f)',mean(v[ok]),sd(v[ok])) else 'Suppressed';z[[paste0(g,'_after')]]<-if(n>=10)sprintf('%.1f (%.1f)',weighted.mean(v[ok],w[ok]),wsd(v[ok],w[ok])) else 'Suppressed'}
  }
  bb<-bal[population==pop & term==z$model_term];z$SMD_before<-if(nrow(bb)==1)abs(bb$SMD_before) else NA_real_;z$SMD_after<-if(nrow(bb)==1)abs(bb$SMD_after) else NA_real_
  rows[[length(rows)+1]]<<-as.data.table(z)
 }
 for(k in names(continuous)){
  v<-as.numeric(x[[k]]);v[!is.finite(v)]<-NA_real_
  if(k%in%c('known_creatinine','baseline_inr'))v[v<=0]<-NA_real_
  if(k=='baseline_bilirubin')v[v<0]<-NA_real_
  add(k,continuous[[k]],v,term=smdmap[[k]])
  for(g in c('LZD','VAN')){ii<-if(g=='LZD')a else !a;nm<-sum(ii&is.na(v));miss[[length(miss)+1]]<-data.table(population=pop,variable=k,label=continuous[[k]],drug=g,total=sum(ii),missing=hide_n(nm),missing_percent=if(nm>0&&nm<10)NA_real_ else 100*nm/sum(ii))}
 }
 add('female','Female sex',as.numeric(x$gender=='F'),TRUE,'genderM')
 for(k in names(binary))add(k,binary[[k]],as.numeric(x[[k]]),TRUE)
 for(k in sort(unique(as.character(x$icu_unit))))add(paste0('icu_',k),k,as.numeric(x$icu_unit==k),TRUE,paste0('icu_unit',k))
 for(k in sort(unique(as.character(x$mapped_era))))add(paste0('era_',k),paste('Calendar era:',k),as.numeric(x$mapped_era==k),TRUE,paste0('mapped_era',k))
 meta[[pop]]<-list(admissions=nrow(x),patients=length(unique(x$subject_id)),n_LZD=sum(a),n_VAN=sum(!a),data_names=names(m$data),model_terms=names(coef(m$fit)),weight_sum_LZD=sum(w[a]),weight_sum_VAN=sum(w[!a]))
}
tab<-rbindlist(rows,fill=TRUE);fwrite(tab,file.path(O,'TABLE1_BASELINE_v0.12.csv'));fwrite(rbindlist(miss),file.path(O,'TABLE_S1_MISSINGNESS_v0.12.csv'))
fwrite(bal[population%in%pops],file.path(O,'TABLE_S2_MODEL_BALANCE_v0.12.csv'))
write_json(list(utc=format(Sys.time(),tz='UTC',usetz=TRUE),source='Frozen phase11 population RDS and phase10 balance export',models=meta,continuous='Observed nonmissing original units; unweighted sample SD and weighted population SD. Model SMDs may use transformations and median filling.',binary='Unweighted n (percent), weighted percent; small raw counts and associated percentages suppressed. Weighting does not create patient counts.',no_refitting=TRUE,no_outcome_changes=TRUE),file.path(O,'TABLE_PROVENANCE_v0.12.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
cat('Aggregate baseline tables exported. No individual data exported.\n')
