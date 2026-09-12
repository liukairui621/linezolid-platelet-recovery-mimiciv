# First exploratory disease contrast. v0.5 frozen before these results.
suppressPackageStartupMessages({library(edgeR);library(limma);library(data.table);library(jsonlite)})
R <- '/root/projects/linezolid_platelet_recovery'
P <- file.path(R,'public_data/GSE252275_processed')
O <- file.path(R,'public_data/GSE252275_DE_v0.5');dir.create(O,showWarnings=FALSE)
F <- file.path(R,'figures/v0.5');dir.create(F,recursive=TRUE,showWarnings=FALSE)
planfile <- file.path(R,'config/analysis_plan_v0.5.json')
stopifnot(strsplit(system2('sha256sum',planfile,stdout=TRUE),' +')[[1]][1] ==
 'c2f348bb0cec1af380a887d33be8392d03a93b79daf130b92db400ffa152ed70')
plan <- fromJSON(planfile)
ss <- as.data.frame(fread(file.path(P,'samples_and_qc.tsv')))
ss <- ss[ss$group %in% c('sepsis_preserve','sepsis_decrease'),]
ss$group <- factor(ss$group,levels=c('sepsis_preserve','sepsis_decrease'))
stopifnot(nrow(ss)==20,all(as.integer(table(ss$group))==c(11,9)),!anyDuplicated(ss$sample))
dt <- fread(file.path(P,'total_exon_counts.tsv.gz'));ids <- dt[[1]]
cts <- as.matrix(dt[,ss$sample,with=FALSE]);rownames(cts)<-ids;storage.mode(cts)<-'double'
ann <- as.data.frame(fread(file.path(P,'gene_annotations.tsv')))
ann <- ann[match(ids,ann$gene_id_as_supplied),]
stopifnot(!anyNA(ann$gene_id_as_supplied),!anyDuplicated(ids),all(is.finite(cts)),all(cts>=0),all(cts==floor(cts)))
globin <- ann$symbol %in% plan$omics$composition_flags$globin_genes
ss$globin_fraction <- colSums(cts[globin,,drop=FALSE])/colSums(cts)
ss$globin_logit <- qlogis(pmax(1e-6,pmin(1-1e-6,ss$globin_fraction)))
flag <- ss$globin_fraction>.2
stopifnot(setequal(ss$sample[flag],paste0('sepsis_decrease_',4:6)))
keep <- rowSums(t(t(cts)/colSums(cts))*1e6>=1)>=9
stopifnot(sum(keep)>1000)
gene <- data.frame(gene_id=ids[keep],symbol=ann$symbol[keep],is_globin=globin[keep])
fwrite(cbind(ss,flagged_high_globin=flag),file.path(O,'samples_used.tsv'),sep='\t')
fwrite(data.frame(gene_id=ids,retained=keep,is_globin=globin),file.path(O,'fixed_filter.tsv.gz'),sep='\t',compress='gzip')
fit_branch <- function(ix,adjust=FALSE) {
 d <- ss[ix,,drop=FALSE];d$globin_centered<-d$globin_logit-mean(d$globin_logit)
 design <- model.matrix(if(adjust)~group+globin_centered else ~group,data=d)
 stopifnot(qr(design)$rank==ncol(design))
 y <- DGEList(cts[,ix,drop=FALSE]);y<-calcNormFactors(y,method='TMM');y<-y[keep,,keep.lib.sizes=TRUE]
 v <- voom(y,design,plot=FALSE,normalize.method='none')
 fit <- eBayes(lmFit(v,design),robust=TRUE)
 tt <- topTable(fit,coef='groupsepsis_decrease',number=Inf,sort.by='none')
 stopifnot(identical(rownames(tt),gene$gene_id),all(is.finite(tt$P.Value)),all(is.finite(tt$logFC)))
 out <- cbind(gene,tt)
 out$q_non_globin_reporting_universe <- NA_real_
 out$q_non_globin_reporting_universe[!out$is_globin]<-p.adjust(out$P.Value[!out$is_globin],'BH')
 list(table=out,logCPM=v$E,design=design,normalization=y$samples,
  n_preserve=sum(d$group=='sepsis_preserve'),n_decrease=sum(d$group=='sepsis_decrease'),
  residual_df=nrow(design)-ncol(design),n_FDR05=sum(tt$adj.P.Val<.05),
  n_non_globin_FDR05=sum(out$q_non_globin_reporting_universe<.05,na.rm=TRUE))
}
specs<-list(all20_group=list(ix=seq_len(20),adjust=FALSE),
 all20_group_globin=list(ix=seq_len(20),adjust=TRUE),
 exclude3_group=list(ix=which(!flag),adjust=FALSE),
 exclude3_group_globin=list(ix=which(!flag),adjust=TRUE))
fits<-list();summ<-list()
for(nm in names(specs)){
 cat('Fitting',nm,'\n');flush.console()
 f<-do.call(fit_branch,specs[[nm]]);fits[[nm]]<-f
 fwrite(f$table,file.path(O,paste0(nm,'.tsv.gz')),sep='\t',compress='gzip')
 fwrite(cbind(sample=rownames(f$normalization),f$normalization),file.path(O,paste0(nm,'_normalization.tsv')),sep='\t')
 summ[[nm]]<-f[c('n_preserve','n_decrease','residual_df','n_FDR05','n_non_globin_FDR05')]
}
ref<-fits$all20_group$table
looFC<-looQ<-matrix(NA_real_,sum(keep),20,dimnames=list(gene$gene_id,ss$sample))
for(i in seq_len(20)){
 cat('Leave one out',ss$sample[i],'\n');flush.console()
 z<-fit_branch(setdiff(seq_len(20),i),FALSE)$table
 looFC[,i]<-z$logFC;looQ[,i]<-z$adj.P.Val
}
fwrite(data.frame(gene_id=gene$gene_id,looFC,check.names=FALSE),file.path(O,'loo_log2FC.tsv.gz'),sep='\t',compress='gzip')
fwrite(data.frame(gene_id=gene$gene_id,looQ,check.names=FALSE),file.path(O,'loo_BH_q.tsv.gz'),sep='\t',compress='gzip')
E<-fits$all20_group$logCPM
rho<-as.numeric(cor(t(E),ss$globin_fraction,method='spearman'))
ranks<-t(apply(E,1,rank,ties.method='average'))
X<-fits$all20_group$design
rgenes<-qr.resid(qr(X),t(ranks));rglobin<-qr.resid(qr(X),rank(ss$globin_fraction,ties.method='average'))
partial_rho<-as.numeric(cor(rgenes,rglobin))
out<-gene
for(nm in names(fits)){
 out[[paste0(nm,'_log2FC')]]<-fits[[nm]]$table$logFC
 out[[paste0(nm,'_q')]]<-fits[[nm]]$table$adj.P.Val
}
out$globin_rho<-rho;out$globin_partial_rank_rho_given_group<-partial_rho
out$abs_globin_rho_ge_05_annotation<-abs(rho)>=.5
out$LOO_min_log2FC<-apply(looFC,1,min);out$LOO_max_log2FC<-apply(looFC,1,max)
out$LOO_same_direction_count<-rowSums(sign(looFC)==sign(ref$logFC))
out$LOO_FDR05_count<-rowSums(looQ<.05)
fcs<-sapply(fits,function(x)x$table$logFC);qs<-sapply(fits,function(x)x$table$adj.P.Val)
out$all4_same_direction<-apply(sign(fcs),1,function(z)length(unique(z))==1)
out$all4_FDR05<-rowSums(qs<.05)==4
out$reference_candidate_FDR05<-ref$adj.P.Val<.05
out<-out[order(out$all20_group_q),]
fwrite(out,file.path(O,'gene_sensitivity_and_globin_annotations.tsv.gz'),sep='\t',compress='gzip')
fwrite(out[out$reference_candidate_FDR05,],file.path(O,'reference_candidates.tsv'),sep='\t')
fwrite(head(out,30),file.path(O,'top30_reference_ranked.tsv'),sep='\t')
cor_pairs<-sapply(fits[-1],function(x)cor(ref$logFC,x$table$logFC,method='spearman'))
candidate<-out[out$reference_candidate_FDR05,]
report<-list(created_utc=format(Sys.time(),tz='UTC',usetz=TRUE),stage='Exploratory disease-state discovery; no drug exposure or recovery outcome in this cohort',
 tested_genes=sum(keep),supplied_genes=nrow(cts),globin_genes_retained=sum(globin[keep]),
 contrast='sepsis_decrease minus sepsis_preserve',branches=summ,
 full_gene_log2FC_spearman_vs_reference=as.list(cor_pairs),
 reference_candidates=list(n=nrow(candidate),non_globin=sum(!candidate$is_globin),
  all4_same_direction=sum(candidate$all4_same_direction),all4_FDR05=sum(candidate$all4_FDR05),
  LOO_all20_same_direction=sum(candidate$LOO_same_direction_count==20),LOO_all20_FDR05=sum(candidate$LOO_FDR05_count==20),
  absolute_globin_rho_ge05=sum(candidate$abs_globin_rho_ge_05_annotation),
  median_abs_globin_rho=if(nrow(candidate))median(abs(candidate$globin_rho)) else NA_real_),
 top_reference_genes=head(out[,c('gene_id','symbol','all20_group_log2FC','all20_group_q','all20_group_globin_q',
  'exclude3_group_q','exclude3_group_globin_q','globin_rho','globin_partial_rank_rho_given_group','LOO_same_direction_count','LOO_FDR05_count')],10),
 limits=c('Sample is the deposited statistical unit; individual donor independence not explicitly documented.',
  'All comparisons and correlations reuse the same cohort; branches are sensitivity analyses, not independent validation.',
  'Globin is an endogenous composition proxy, not a verified technical contaminant; adjustment may remove biology or introduce conditioning bias.',
  'Deleting globin rows only changes the reporting multiplicity universe, not other gene fits or contamination.',
  'BH correction is within each prespecified branch, not across a selected union of branch discoveries.',
  'No pathway test, mechanistic target claim, clinical effect or drug-disease signature correlation computed.'))
write_json(report,file.path(R,'reports/DISEASE_DE_v0.5.json'),pretty=TRUE,auto_unbox=TRUE,na='null',digits=10)
capture.output(sessionInfo(),file=file.path(R,'reports/R_sessionInfo_v0.5_omics.txt'))
draw<-function(){
 par(mfrow=c(2,2),mar=c(4.5,4.5,3,1));pal<-ifelse(ref$adj.P.Val<.05,'#D55E00','#77777755')
 plot(ref$logFC,-log10(pmax(ref$adj.P.Val,1e-300)),pch=16,cex=.45,col=pal,
  xlab='log2 fold change (decrease - preserve)',ylab='-log10(BH q)',main=paste0('All 20: ',sum(ref$adj.P.Val<.05),' genes, q < 0.05'))
 abline(h=-log10(.05),lty=2,col='#333333')
 for(nm in names(fits)[-1]){
  plot(ref$logFC,fits[[nm]]$table$logFC,pch=16,cex=.4,col=pal,xlab='All 20 group-only log2FC',
   ylab='Sensitivity log2FC',main=paste(nm,'\nSpearman',round(cor_pairs[[nm]],3)))
  abline(a=0,b=1,lty=2,col='#0072B2');abline(h=0,v=0,col='#BBBBBB')
 }
}
png(file.path(F,'FigD3_disease_DE_sensitivity.png'),width=2400,height=2000,res=220);draw();dev.off()
pdf(file.path(F,'FigD3_disease_DE_sensitivity.pdf'),width=11,height=9);draw();dev.off()
inputpaths<-c(planfile,file.path(R,'scripts/27_disease_DE_v05.R'),file.path(P,c('total_exon_counts.tsv.gz','gene_annotations.tsv','samples_and_qc.tsv')))
hashes<-setNames(vapply(inputpaths,function(p)strsplit(system2('sha256sum',p,stdout=TRUE),' +')[[1]][1],''),inputpaths)
write_json(list(created_utc=report$created_utc,input_sha256=as.list(hashes),plan_verified=TRUE),file.path(R,'manifests/disease_DE_v0.5.json'),pretty=TRUE,auto_unbox=TRUE)
cat(toJSON(report,pretty=TRUE,auto_unbox=TRUE,na='null',digits=8),'\n')
