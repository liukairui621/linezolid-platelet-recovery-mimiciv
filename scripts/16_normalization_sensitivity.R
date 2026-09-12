# Descriptive normalization audit only: no differential expression or pathway test.
suppressPackageStartupMessages({library(data.table);library(edgeR);library(jsonlite)})
set.seed(42)
root <- '/root/projects/linezolid_platelet_recovery'
p <- file.path(root,'public_data/GSE252275_processed')
tab <- fread(cmd=paste('gzip -dc',shQuote(file.path(p,'total_exon_counts.tsv.gz'))))
ids <- tab[[1]];counts <- as.matrix(tab[,-1]);rownames(counts) <- ids
ann <- fread(file.path(p,'gene_annotations.tsv'));ann <- ann[match(ids,ensembl_id)]
sam <- fread(file.path(p,'samples_and_qc.tsv'));sam <- sam[match(colnames(counts),sample)]
globin <- c('HBA1','HBA2','HBB','HBD','HBE1','HBG1','HBG2','HBM','HBQ1','HBZ')
reference <- ann$chromosome!='MT' & !ann$symbol %in% globin
report <- list()
for (scope in c('all_40','sepsis_20')) {
    use <- if(scope=='all_40') rep(TRUE,ncol(counts)) else sam$group!='control'
    y <- counts[,use,drop=FALSE]
    full <- calcNormFactors(DGEList(y),method='TMM')
    ref <- calcNormFactors(DGEList(y[reference,,drop=FALSE]),method='TMM')
    offset_full <- full$samples$lib.size*full$samples$norm.factors
    offset_ref <- ref$samples$lib.size*ref$samples$norm.factors
    # Global offset multiplier has no between-sample meaning.
    offset_full <- offset_full/exp(mean(log(offset_full)))
    offset_ref <- offset_ref/exp(mean(log(offset_ref)))
    normtab <- data.frame(sample=colnames(y),group=sam$group[use],
        full_TMM_relative_offset=offset_full,nuclear_non_globin_TMM_relative_offset=offset_ref,
        log2_offset_difference=log2(offset_ref/offset_full))
    fwrite(normtab,file.path(p,paste0('TMM_sensitivity_',scope,'.tsv')),sep='\t')
    corr <- cor(log2(offset_full),log2(offset_ref),method='pearson')
    keep <- reference & rowSums(cpm(y)>=1)>=min(9,ncol(y))
    z <- cpm(y,lib.size=offset_ref*exp(mean(log(colSums(y)))),log=TRUE,prior.count=.5)[keep,,drop=FALSE]
    top <- order(apply(z,1,var),decreasing=TRUE)[seq_len(min(2000,nrow(z)))]
    pc <- prcomp(t(z[top,,drop=FALSE]),center=TRUE,scale.=FALSE)
    report[[scope]] <- list(n_samples=ncol(y),offset_log_correlation=corr,
        maximum_abs_log2_offset_difference=max(abs(normtab$log2_offset_difference)),
        nuclear_PCA_PC1_variance=pc$sdev[1]^2/sum(pc$sdev^2),
        abs_PC1_detected_features_spearman=abs(cor(pc$x[,1],sam$detected_features[use],method='spearman')))
}
report$interpretation <- 'TMM uses trimming; exclusion of mitochondrial/globin reference features is a sensitivity analysis, not an automatic universal requirement. Similar offsets do not prove that biological contamination or composition has been removed.'
write_json(report,file.path(root,'reports/NORMALIZATION_SENSITIVITY_v0.3.json'),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(root,'reports/R_sessionInfo_v0.3.txt'))
print(toJSON(report,pretty=TRUE,auto_unbox=TRUE))
