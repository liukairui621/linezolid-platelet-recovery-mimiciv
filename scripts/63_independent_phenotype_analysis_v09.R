# v0.9 independent phenotype evidence; all contrasts and exact pathway sets frozen before these results.
suppressPackageStartupMessages({library(edgeR);library(limma);library(data.table);library(jsonlite)})
R<-'/root/projects/linezolid_platelet_recovery';P<-file.path(R,'public_data/phenotype_evidence_v0.9');O<-file.path(P,'analysis');dir.create(O,recursive=TRUE,showWarnings=FALSE)
cfg<-fromJSON(file.path(R,'config/pathway_sets_v0.8.json'),simplifyVector=FALSE)$sets
setids<-vapply(cfg,function(z)z$id,'');sets<-setNames(lapply(cfg,function(z)unlist(z$genes)),setids);setnames<-setNames(vapply(cfg,function(z)z$name,''),setids)
checks<-list();allres<-list();coverage<-list();gene_sets<-list()
# Reuse the immutable v0.8 pathway function only; do not rerun its old analyses.
expr<-parse(file.path(R,'scripts/51_pathways_and_perturbation_v08.R'))
fun<-Filter(function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('testsets')),as.list(expr));stopifnot(length(fun)==1);eval(fun[[1]])
analyses<-list()
runfit<-function(v,X,coef,symbol,dataset,branch){
 stopifnot(all(is.finite(v$E)),qr(X)$rank==ncol(X),nrow(X)==ncol(v$E),length(symbol)==nrow(v$E))
 fit<-eBayes(lmFit(v,X),robust=TRUE);tt<-topTable(fit,coef=coef,number=Inf,sort.by='none')
 fwrite(data.table(gene_id=rownames(tt),symbol=symbol,tt),file.path(O,paste0(dataset,'_',branch,'_DE.tsv.gz')),sep='\t',compress='gzip')
 analyses[[paste(dataset,branch,sep='__')]]<<-list(samples=nrow(X),genes=nrow(v$E),model_columns=colnames(X),rank=qr(X)$rank,residual_df=nrow(X)-ncol(X),DE_FDR05=sum(tt$adj.P.Val<.05),min_DE_FDR=min(tt$adj.P.Val),positive_contrast=colnames(X)[coef])
 testsets(v,X,coef,symbol,dataset,branch);cat(dataset,branch,'complete\n');flush.console()
}
# MARS: supplied log-scale, processed expression; fixed representative probes selected across all99 without labels.
s<-as.data.frame(fread(file.path(P,'GSE65682_phenotype_samples.tsv'),na.strings=c('NA','')))
d<-fread(file.path(P,'GSE65682_selected_expression.tsv.gz'));E<-as.matrix(d[,s$GSM,with=FALSE]);rownames(E)<-as.character(d[[1]])
a<-fread(file.path(P,'GPL13667_annotation.tsv.gz'),na.strings=c('NA','---',''));sy<-trimws(a[['symbol']][match(rownames(E),a[['probe_id']])])
stopifnot(length(sy)==nrow(E),sum(!is.na(sy))>10000)
ok<-!is.na(sy)&nzchar(sy)&!grepl('///',sy,fixed=TRUE);oo<-order(-rowMeans(E),rownames(E));oo<-oo[ok[oo]];oo<-oo[!duplicated(sy[oo])]
E<-E[oo,,drop=FALSE];sy<-sy[oo];s$age<-as.numeric(s$age);s$gender<-factor(s$gender);s$pneumonia_source<-factor(ifelse(is.na(s[['pneumonia diagnoses']]),'unknown',s[['pneumonia diagnoses']]))
s$grade<-unname(c(D_normal=0,C_low=1,B_medium_low=2,A_very_low=3)[s$thrombocytopenia]);stopifnot(nrow(s)==99,!anyDuplicated(s$GSM),!anyDuplicated(s$sample_suffix),all(complete.cases(s[,c('age','gender','pneumonia_source','grade')])))
# No paired post-admission samples selected. Metadata and source-study admission design support statistical units.
stopifnot(all(is.na(s$icu_acquired_infection_paired)|s$icu_acquired_infection_paired=='Admis_ICUAI'))
fwrite(s,file.path(O,'GSE65682_sample_design.tsv'),sep='\t');fwrite(data.table(probe=rownames(E),symbol=sy),file.path(O,'GSE65682_representative_probes.tsv'),sep='\t')
corE<-cor(E);diag(corE)<-NA
checks$GSE65682<-list(selected_samples=nrow(s),primary_cases=sum(s$grade>=2),primary_controls=sum(s$grade==0),all_low_cases=sum(s$grade>=1),supplied_probes=nrow(d),unambiguous_unique_symbols=nrow(E),expression_range=range(E),max_inter_sample_correlation=max(corE,na.rm=TRUE),pairs_correlation_above_099=sum(corE[upper.tri(corE)]>.99,na.rm=TRUE),normalization='Deposited preprocessed log-scale values as supplied. No second log or normalization.',unknown_pneumonia_source=sum(s$pneumonia_source=='unknown'))
for(branch in c('lt100_vs_normal_adjusted','lt150_vs_normal_adjusted','ordinal_grade_adjusted')){
 ix<-if(branch=='lt100_vs_normal_adjusted')which(s$grade!=1)else seq_len(nrow(s));sm<-s[ix,,drop=FALSE]
 sm$phenotype<-if(branch=='ordinal_grade_adjusted')sm$grade else as.integer(sm$grade>0)
 X<-model.matrix(~age+gender+pneumonia_source+phenotype,sm);v<-new('EList',list(E=E[,ix,drop=FALSE],design=X))
 runfit(v,X,match('phenotype',colnames(X)),sy,'GSE65682',branch)
}
# Day1 only: 52 unique patients. Blood-derived support; unresolved PBMC/whole-blood metadata retained.
s<-as.data.frame(fread(file.path(P,'GSE273700_samples.tsv')));s<-s[s$time=='Day1',];stopifnot(nrow(s)==52,!anyDuplicated(s$patient),!anyDuplicated(s$sample))
d<-fread(file.path(P,'GSE273700_counts_mapped.tsv.gz'));cc<-as.matrix(d[,s$sample,with=FALSE]);rownames(cc)<-d[[1]];stopifnot(all(cc>=0),all(cc==floor(cc)),!anyDuplicated(rownames(cc)))
h<-fread(file.path(R,'public_data/reference_sources_v0.8/hgnc_complete_set.txt'),select=c('symbol','status','ensembl_gene_id'));h<-unique(h[status=='Approved' & !is.na(ensembl_gene_id) & nzchar(ensembl_gene_id),.(symbol,ensembl_gene_id)])
amb<-h[,.(n=uniqueN(symbol)),by=ensembl_gene_id][n>1]$ensembl_gene_id;h<-h[!ensembl_gene_id%in%amb];sy<-h$symbol[match(sub('\\.[0-9]+$','',rownames(cc)),h$ensembl_gene_id)]
s$phenotype<-as.integer(s$condition=='Thrombocytopenia(+)');stopifnot(sum(s$phenotype)==22)
globins<-c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ');s$globin_fraction<-colSums(cc[sy%in%globins,,drop=FALSE])/colSums(cc)
s$globin_centered<-qlogis(pmax(1e-6,pmin(1-1e-6,s$globin_fraction)));s$globin_centered<-s$globin_centered-mean(s$globin_centered)
keep<-rowSums(t(t(cc)/colSums(cc))*1e6>=1)>=22;y<-calcNormFactors(DGEList(cc),method='TMM');y<-y[keep,,keep.lib.sizes=TRUE]
fwrite(s,file.path(O,'GSE273700_sample_design.tsv'),sep='\t')
checks$GSE273700<-list(samples=nrow(s),unique_patients=length(unique(s$patient)),cases=sum(s$phenotype),controls=sum(s$phenotype==0),supplied_genes=nrow(cc),tested_genes=sum(keep),mapped_tested_genes=sum(!is.na(sy[keep])),globin_range=range(s$globin_fraction),globin_above20pct=sum(s$globin_fraction>.2),tissue_conflict_retained=TRUE,threshold='Admission platelet <=150 x10^9/L',serial_recovery_validation=FALSE)
for(branch in c('Day1_group','Day1_group_globin')){
 X<-model.matrix(if(branch=='Day1_group')~phenotype else ~phenotype+globin_centered,s);v<-voom(y,X,plot=FALSE,normalize.method='none');runfit(v,X,match('phenotype',colnames(X)),sy[keep],'GSE273700',branch)
}
res<-rbindlist(allres);fwrite(res,file.path(R,'reports/INDEPENDENT_PHENOTYPE_PATHWAYS_v0.9.tsv'),sep='\t');fwrite(rbindlist(coverage),file.path(R,'reports/INDEPENDENT_PHENOTYPE_COVERAGE_v0.9.tsv'),sep='\t');fwrite(rbindlist(gene_sets),file.path(O,'pathway_member_statistics.tsv.gz'),sep='\t',compress='gzip')
checks$analyses<-analyses;checks$completed_utc<-format(Sys.time(),tz='UTC',usetz=TRUE);checks$limitations<-c('Outcome-informed external follow-up, not preregistration.', 'Phenotype association, not linezolid treatment or recovery validation.', 'GSE65682 pneumonia NA retained as unknown. Severity and processing batch cannot be adjusted without sample-level values.', 'GSE273700 tissue descriptions conflict; no platelet-intrinsic claim.', 'Gene-set tests have distinct null hypotheses; fixed-correlation or fry results cannot replace primary camera results.')
write_json(checks,file.path(R,'reports/INDEPENDENT_PHENOTYPE_ANALYSIS_v0.9.json'),pretty=TRUE,auto_unbox=TRUE,digits=12)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.9_omics.txt'))
print(res[,.(dataset,branch,set_id,expressed_N,camera_direction,camera_FDR,fixed01_FDR,fry_FDR,mean_log2FC)])
