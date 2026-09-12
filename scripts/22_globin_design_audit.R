# Quantify design consequences of globin adjustment/removal, without DE or enrichment.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(jsonlite)})
set.seed(42)
root<-'/root/projects/linezolid_platelet_recovery';p<-file.path(root,'public_data/GSE252275_processed')
tab<-fread(cmd=paste('gzip -dc',shQuote(file.path(p,'total_exon_counts.tsv.gz'))))
ids<-tab[[1]];counts<-as.matrix(tab[,-1]);rownames(counts)<-ids
ann<-fread(file.path(p,'gene_annotations.tsv'));ann<-ann[match(ids,ensembl_id)]
sam<-fread(file.path(p,'samples_and_qc.tsv'));sam<-sam[match(colnames(counts),sample)]
globins<-c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
isglobin<-ann$symbol %in% globins
sam$globin_fraction<-colSums(counts[isglobin,,drop=FALSE])/colSums(counts)
sam$flagged<-sam$sample %in% c('sepsis_decrease_4','sepsis_decrease_5','sepsis_decrease_6')
sam$group<-factor(sam$group,levels=c('sepsis_preserve','sepsis_decrease','control'))
res<-list();leverage<-list()
for(scope in c('all20','exclude3')) {
 keep<-sam$group!='control' & (scope=='all20'|!sam$flagged)
 x<-droplevels(as.data.frame(sam[keep,]));y<-counts[,keep,drop=FALSE]
 x$globin_logit<-qlogis(pmin(pmax(x$globin_fraction,1e-6),1-1e-6))
 x$globin_logit<-x$globin_logit-mean(x$globin_logit)
 design0<-model.matrix(~group,data=x);design1<-model.matrix(~group+globin_logit,data=x)
 rsq<-summary(lm(as.numeric(group=='sepsis_decrease')~globin_logit,data=x))$r.squared
 lev0<-hat(design0);lev1<-hat(design1)
 base<-calcNormFactors(DGEList(y),method='TMM')
 onlyrows<-base[!isglobin,,keep.lib.sizes=TRUE]
 recompute<-calcNormFactors(DGEList(y[!isglobin,,drop=FALSE]),method='TMM')
 off0<-base$samples$lib.size*base$samples$norm.factors
 off1<-recompute$samples$lib.size*recompute$samples$norm.factors
 delta<-log2(off1/off0);delta<-delta-mean(delta)
 a<-cpm(base,log=TRUE,prior.count=.5)[!isglobin,,drop=FALSE]
 b<-cpm(onlyrows,log=TRUE,prior.count=.5)
 stopifnot(identical(rownames(a),rownames(b)),max(abs(a-b))<1e-12)
 res[[scope]]<-list(n=nrow(x),group_counts=as.list(table(x$group)),
  group_only_rank=qr(design0)$rank,group_only_residual_df=nrow(x)-qr(design0)$rank,
  globin_adjusted_rank=qr(design1)$rank,globin_adjusted_residual_df=nrow(x)-qr(design1)$rank,
  group_globin_R_squared=rsq,group_coefficient_VIF=1/(1-rsq),
  max_leverage_unadjusted=max(lev0),max_leverage_adjusted=max(lev1),
  max_abs_logCPM_change_removing_only_globin_rows_fixed_offsets=max(abs(a-b)),
  max_abs_centered_log2_offset_change_if_recomputed=max(abs(delta)))
 leverage[[scope]]<-data.table(scope=scope,sample=x$sample,GSM=x$GSM,group=x$group,
  globin_fraction=x$globin_fraction,flagged=x$flagged,leverage_unadjusted=lev0,leverage_adjusted=lev1)
}
report<-list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),globin_symbols=globins,
 exact_globin_rows=sum(isglobin),results=res,
 interpretation=c('Globin adjustment is estimable but an endogenous sensitivity, not automatically the unique correct model.',
 'VIF/hat leverage are unweighted design diagnostics, not voom observation-level leverage or proof of adequate power.',
 'Removing globin target rows while keeping library offsets leaves other logCPM values unchanged; it does not decontaminate other genes.',
 'No differential expression, pathway testing or method selection by significance in this audit.'))
write_json(report,file.path(root,'reports/GLOBIN_DESIGN_AUDIT_v0.4.json'),pretty=TRUE,auto_unbox=TRUE,digits=10)
fwrite(rbindlist(leverage),file.path(p,'globin_design_leverage_v0.4.tsv'),sep='\t')
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE,digits=5))
