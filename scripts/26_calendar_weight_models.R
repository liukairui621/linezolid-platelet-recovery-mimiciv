# Calendar-aware weighting and separation checks; no outcome analysis.
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(lpSolve);library(ggplot2)})
root<-'/root/projects/linezolid_platelet_recovery';set.seed(42)
stopifnot(fromJSON(file.path(root,'config/analysis_plan_v0.5.json'))$version=='0.5')
old<-readRDS(file.path(root,'cache/phase4_overlap_models.rds'))
cal<-as.data.frame(fread(file.path(root,'cache/phase5_calendar.csv')))
d<-old$extended_baseline$data;i<-match(d$hadm_id,cal$hadm_id);stopifnot(!anyNA(i))
d$mapped_era<-factor(cal$mapped_era[i],levels=c('mapped_early','boundary_uncertain','mapped_late'))
d$mapped_year_lower<-cal$mapped_year_lower[i];d$mapped_year_upper<-cal$mapped_year_upper[i]
d$early_anchor<-as.integer(cal$early_anchor[i]);d$year_midpoint<-(d$mapped_year_lower+d$mapped_year_upper)/2
ess<-function(w)sum(w)^2/sum(w*w)
sep_lp<-function(X,A){
 # A nonzero nonnegative signed margin implies complete or quasi separation.
 # Scaling and L1<=1 avoid numerical scale and unboundedness artifacts.
 Z<-X
 for(j in seq_len(ncol(Z)))if(sd(Z[,j])>0)Z[,j]<-(Z[,j]-mean(Z[,j]))/sd(Z[,j])
 Z<-Z*(2*A-1);K<-cbind(Z,-Z)
 ans<-lp('max',objective.in=colMeans(K),const.mat=rbind(K,rep(1,ncol(K))),
    const.dir=c(rep('>=',nrow(K)),'<='),const.rhs=c(rep(0,nrow(K)),1))
 if(ans$status!=0)return(list(solver_status=ans$status,separated=NA,positive_margin_objective=NA))
 list(solver_status=0,separated=ans$objval>1e-7,positive_margin_objective=ans$objval)
}
stopifnot(sep_lp(cbind(1,1:6),c(0,0,0,1,1,1))$separated)
stopifnot(!sep_lp(cbind(1,c(0,0,1,1)),c(0,1,0,1))$separated)
stopifnot(sep_lp(cbind(1,c(0,0,0,1,1)),c(0,1,0,0,0))$separated)
terms_ext<-attr(terms(old$extended_baseline$fit),'term.labels')
terms_new<-c(setdiff(terms_ext,c('anchor_year_group','decision_year_minus_anchor_year')),'mapped_era','any_emar_available_before')
make_f<-function(ts,x)reformulate(ts[vapply(ts,function(v)length(unique(x[[v]]))>1,logical(1))],response='A')
newfit<-function(x,ts){f<-make_f(ts,x);fit<-glm(f,data=x,family=binomial(),control=glm.control(maxit=100,epsilon=1e-10));
 list(fit=fit,data=x,weights=ifelse(x$A==1,1-fitted(fit),fitted(fit)))}
mods<-old[c('core_baseline','extended_baseline','baseline_emar_covered')]
mods$calendar_source<-newfit(d,terms_new)
mods$calendar_source_creatinine_observed<-newfit(d[d$creatinine_missing==0,],terms_new)
mods$proposed_small_PS<-newfit(d,c('log_creatinine','prior_recorded_marrow','mapped_era','known_platelets','age','gender'))
# Common-sample ablations describe ESS burden only; not a causal decomposition.
mods$renal_common_sample<-newfit(d,c('log_creatinine','creatinine_missing'))
mods$renal_plus_calendar<-newfit(d,c('log_creatinine','creatinine_missing','mapped_era'))
mods$renal_plus_anchor<-newfit(d,c('log_creatinine','creatinine_missing','anchor_year_group'))
audit_terms<-unique(c(terms_ext,'mapped_era','any_emar_available_before','ufh_documented_prior7d','lmwh_documented_prior7d','selected_other_drug_prior7d'))
bal<-list();era<-list();reports<-list()
for(nm in names(mods)){
 ob<-mods[[nm]];x<-ob$data;j<-match(x$hadm_id,d$hadm_id)
 for(col in c('mapped_era','mapped_year_lower','mapped_year_upper','year_midpoint','early_anchor'))x[[col]]<-d[[col]][j]
 ob$data<-x;mods[[nm]]<-ob;w<-ob$weights
 mm<-model.matrix(ob$fit);a<-x$A==1;b<-!a
 mat<-model.matrix(make_f(audit_terms,x),x)[,-1,drop=FALSE]
 bb<-rbindlist(lapply(colnames(mat),function(k){v<-mat[,k];s<-sqrt((var(v[a])+var(v[b]))/2)
  if(!is.finite(s)||s==0)return(NULL)
  data.table(model=nm,term=k,included=k %in% colnames(mm),SMD_before=(mean(v[a])-mean(v[b]))/s,
   SMD_after=(weighted.mean(v[a],w[a])-weighted.mean(v[b],w[b]))/s)}))
 bal[[nm]]<-bb
 er<-rbindlist(lapply(c('anchor_year_group','mapped_era'),function(v){
  dt<-data.table(drug=x$drug,category=as.character(x[[v]]),weight=w)
  z<-dt[,.(n=.N,weight_sum=sum(weight)),by=.(drug,category)]
  z[,raw_fraction:=n/sum(n),by=drug];z[,weighted_fraction:=weight_sum/sum(weight_sum),by=drug]
  z[,`:=`(model=nm,classification=v)];z}))
 era[[nm]]<-er
 reports[[nm]]<-list(formula=paste(deparse(formula(ob$fit)),collapse=' '),df=ob$fit$rank-1,
  converged=ob$fit$converged,lp_separation=sep_lp(mm,x$A),
  eMAR_SMD_after=bb[term=='any_emar_available_before']$SMD_after,
  largest_abs_SMD_included=max(abs(bb$SMD_after[bb$included])),
  largest_abs_SMD_omitted=max(abs(bb$SMD_after[!bb$included])),
  arm_statistics=lapply(0:1,function(k){ok<-x$A==k;z<-w[ok]
   list(drug=if(k==1)'linezolid' else 'vancomycin',n=sum(ok),ESS=ess(z),
    weighted_mean_calendar_midpoint=weighted.mean(x$year_midpoint[ok],z),
    weighted_mapped_early=weighted.mean(as.integer(x$mapped_era[ok]=='mapped_early'),z),
    weighted_mapped_late=weighted.mean(as.integer(x$mapped_era[ok]=='mapped_late'),z))}))
}
boots<-list()
for(nm in c('core_baseline','extended_baseline','calendar_source','calendar_source_creatinine_observed')){
 set.seed(42);ob<-mods[[nm]];x<-ob$data;f<-formula(ob$fit)
 cl<-split(seq_len(nrow(x)),x$subject_id);patt<-vapply(cl,function(ii)paste(sort(unique(x$A[ii])),collapse=','),character(1));st<-split(names(cl),patt)
 rr<-vector('list',200)
 for(b in 1:200){
  selected<-unlist(lapply(st,function(s)sample(s,length(s),replace=TRUE)),use.names=FALSE)
  z<-x[unlist(cl[selected],use.names=FALSE),,drop=FALSE];warns<-character()
  fit<-tryCatch(withCallingHandlers(glm(f,data=z,family=binomial(),control=glm.control(maxit=100,epsilon=1e-10)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)e)
  if(inherits(fit,'error')){rr[[b]]<-data.table(replicate=b,error=TRUE,warning=TRUE,converged=FALSE,separated=NA,largecoef=NA,large_nonintercept=NA,ESS_LZD=NA,ESS_VAN=NA);next}
  lp<-sep_lp(model.matrix(fit),z$A);w<-ifelse(z$A==1,1-fitted(fit),fitted(fit));cf<-coef(fit)
  rr[[b]]<-data.table(replicate=b,error=FALSE,warning=length(warns)>0,converged=fit$converged,
    separated=lp$separated,largecoef=any(abs(cf)>15,na.rm=TRUE),large_nonintercept=any(abs(cf[names(cf)!='(Intercept)'])>15,na.rm=TRUE),
    ESS_LZD=ess(w[z$A==1]),ESS_VAN=ess(w[z$A==0]))
 }
 zz<-rbindlist(rr);fwrite(zz,file.path(root,paste0('cache/internal_reports/v05_boot_',nm,'.csv')))
 boots[[nm]]<-list(B=200,errors=sum(zz$error),warnings=sum(zz$warning),nonconverged=sum(!zz$converged),
  LP_separated=sum(zz$separated,na.rm=TRUE),LP_missing=sum(is.na(zz$separated)),abscoef_gt15=sum(zz$largecoef,na.rm=TRUE),
  abs_nonintercept_gt15=sum(zz$large_nonintercept,na.rm=TRUE),ESS_LZD_quantiles=unname(quantile(zz$ESS_LZD,c(.025,.5,.975),na.rm=TRUE)))
 cat(nm,'bootstrap completed\n',file=stdout());flush.console()
}
report<-list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),models=reports,bootstrap=boots,
 limitations=c('LP checks complete or quasi separation of a specified design; it does not certify model correctness or exchangeability.',
 'Original-scale large coefficient flags depend on units and intercept; no universal8-10exposures-per-df criterion imposed.',
 'Calendar intervals are approximate; boundary_uncertain stays explicit.',
 'No effect estimates or clinical endpoints computed; complete-case branch changes the target population.'))
write_json(report,file.path(root,'reports/CALENDAR_WEIGHT_MODELS_v0.5.json'),pretty=TRUE,auto_unbox=TRUE,digits=8)
saveRDS(mods,file.path(root,'cache/phase5_weight_models.rds'))
fwrite(rbindlist(bal),file.path(root,'reports/COVARIATE_BALANCE_v0.5.csv'))
er<-rbindlist(era);fwrite(er,file.path(root,'cache/internal_reports/ERA_DISTRIBUTIONS_v0.5_exact.csv'))
er[,n:=ifelse(n>0&n<10,'<10',as.character(n))]
fwrite(er,file.path(root,'reports/ERA_DISTRIBUTIONS_v0.5.csv'))
figdir<-file.path(root,'figures/v0.5');dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
ee<-er[model %in% c('extended_baseline','calendar_source') & classification=='mapped_era']
ee[,label:=factor(category,levels=c('mapped_early','boundary_uncertain','mapped_late'),labels=c('Upper bound <=2013','Crosses 2013/2014','Lower bound >=2014'))]
g<-ggplot(ee,aes(drug,weighted_fraction,fill=label))+geom_col(width=.65)+facet_wrap(~model)+
 scale_fill_manual(values=c('#7895A5','#D4B86B','#438676'))+scale_y_continuous(labels=function(z)paste0(round(z*100),'%'))+
 labs(x=NULL,y='ATO-weighted fraction',fill='Mapped decision-year interval',title='Who does the overlap population represent?',subtitle='Calendar intervals use anchor-year displacement; anchor labels alone are insufficient.')+
 theme_classic(base_size=12)+theme(legend.position='bottom',legend.direction='vertical')
ggsave(file.path(figdir,'FigD2A_weighted_calendar.png'),g,width=8.3,height=5.8,dpi=180)
ggsave(file.path(figdir,'FigD2A_weighted_calendar.pdf'),g,width=8.3,height=5.8)
capture.output(sessionInfo(),file=file.path(root,'reports/R_sessionInfo_v0.5_clinical.txt'))
cat(toJSON(list(models=reports,bootstrap=boots),auto_unbox=TRUE,pretty=TRUE,digits=4))
