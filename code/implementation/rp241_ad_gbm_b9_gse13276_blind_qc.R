args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else '/mnt/h/SCI2/MRQI'
.libPaths(c(file.path(root,'0_admin/environment/rp241_b9_microarray/R_library'), .libPaths()))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(hgu133acdf))
suppressPackageStartupMessages(library(jsonlite))
date <- '2026-10-14'
src <- file.path(root,'3_results/rp241_ad_gbm_b9_source_audit')
out <- file.path(root,'3_results/rp241_ad_gbm_b9_blind_qc')
stopifnot(!dir.exists(out));dir.create(out)
fig <- file.path(root,'4_figures/rp241_ad_gbm_b9_blind_qc')
stopifnot(!dir.exists(fig));dir.create(fig)
meta <- read.delim(file.path(src,paste0('GSE13276_Primary12_Metadata_',date,'.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(meta)==12,length(unique(meta$donor_key))==5,all(meta$role %in% c('TUMOR','HOST_BAT')))
files <- file.path(root,'1_data/transcriptomics/glioma/rp24_a1_geo_candidate_matrices/GSE13276/raw/B9_exact_CEL',paste0(meta$gsm,'.CEL.gz'))
stopifnot(all(file.exists(files)))
headers <- lapply(files,affyio::read.celfile.header)
stopifnot(all(vapply(headers,function(h) identical(h$cdfName,'HG-U133A') && all(h[['CEL dimensions']]==712L),logical(1))))
raw <- ReadAffy(filenames=files, sampleNames=meta$gsm, cdfname='hgu133acdf')
stopifnot(identical(sampleNames(raw),meta$gsm),identical(dim(raw@assayData$exprs),c(712L*712L,12L)))
qc <- data.frame(gsm=meta$gsm,donor_key=meta$donor_key,role=meta$role,
                 raw_nonfinite=colSums(!is.finite(exprs(raw))),raw_min=apply(exprs(raw),2,min),
                 raw_q25=apply(exprs(raw),2,quantile,.25),raw_median=apply(exprs(raw),2,median),
                 raw_q75=apply(exprs(raw),2,quantile,.75),raw_max=apply(exprs(raw),2,max))
stopifnot(all(qc$raw_nonfinite==0))
norm <- rma(raw,background=TRUE,normalize=TRUE)
x <- exprs(norm)
stopifnot(nrow(x)==22283,ncol(x)==12,all(is.finite(x)),identical(colnames(x),meta$gsm),!anyDuplicated(rownames(x)))
saveRDS(x,file.path(out,paste0('GSE13276_RMA12_All22283_Probe_Expression_',date,'.rds')))
qc$rma_min <- apply(x,2,min);qc$rma_median <- apply(x,2,median);qc$rma_max <- apply(x,2,max)
qc$rma_IQR <- apply(x,2,IQR)
write.table(qc,file.path(out,paste0('All12_Raw_RMA_Array_QC_',date,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
map <- read.delim(file.path(src,paste0('GPL96_AllProbe_Exact_Entrez_HGNC_Identity_',date,'.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
map$HGNC_ID[is.na(map$HGNC_ID)] <- ''
stopifnot(setequal(map$probe_id,rownames(x)))
map <- map[match(rownames(x),map$probe_id),]
ok <- map$HGNC_ID!=''
sumx <- rowsum(x[ok,,drop=FALSE],group=map$HGNC_ID[ok],reorder=TRUE)
nprobe <- table(map$HGNC_ID[ok]);gx <- sumx/as.numeric(nprobe[rownames(sumx)])
stopifnot(nrow(gx)==12405,all(is.finite(gx)))
saveRDS(gx,file.path(out,paste0('GSE13276_RMA12_All12405_UniqueHGNC_Expression_',date,'.rds')))
agg <- list();md <- list()
for (donor in sort(unique(meta$donor_key))) for (role in c('HOST_BAT','TUMOR')) {
    index <- which(meta$donor_key==donor & meta$role==role)
    key <- paste(donor,role,sep='__')
    agg[[key]] <- rowMeans(gx[,index,drop=FALSE])
    md[[key]] <- data.frame(analysis_sample_id=key,donor_key=donor,role=role,source_gsms=paste(meta$gsm[index],collapse=';'),source_arrays=length(index))
}
y <- do.call(cbind,agg);am <- do.call(rbind,md);rownames(am)<-NULL
saveRDS(y,file.path(out,paste0('GSE13276_RMA_Aggregated10_All12405_Expression_',date,'.rds')))
write.table(am,file.path(out,paste0('GSE13276_5Pairs_Aggregated10_Analysis_Metadata_',date,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
am$patient <- factor(am$donor_key);am$compartment <- factor(am$role,levels=c('HOST_BAT','TUMOR'))
design <- model.matrix(~patient+compartment,am)
stopifnot(nrow(design)==10,ncol(design)==6,qr(design)$rank==6)
write.table(design,file.path(out,paste0('GSE13276_Primary_Paired_Design_Rank6_df4_',date,'.tsv')),sep='\t',quote=FALSE,col.names=NA)
pc <- prcomp(t(x),scale.=FALSE,center=TRUE)
pcdf <- data.frame(meta[,c('gsm','title','donor_key','role')],PC1=pc$x[,1],PC2=pc$x[,2],PC3=pc$x[,3])
write.table(pcdf,file.path(out,paste0('All12_AllProbe_Blind_PCA_',date,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
v <- 100*pc$sdev^2/sum(pc$sdev^2)
png(file.path(fig,paste0('GSE13276_All12_Array_Blind_QC_',date,'.png')),width=1950,height=850,res=150)
par(mfrow=c(1,2),mar=c(9,4,3,1))
color <- ifelse(meta$role=='TUMOR','#BB5544','#36708D')
boxplot(x,col=color,outline=FALSE,las=2,names=meta$gsm,ylab='RMA log2 probe expression',main='12 arrays / 5 patients; no candidate tests')
par(mar=c(5,4,3,2))
plot(pc$x[,1],pc$x[,2],col=color,pch=19,cex=1.2,xlab=paste0('PC1 (',round(v[1],1),'%)'),ylab=paste0('PC2 (',round(v[2],1),'%)'),main='All 22,283 probes; patient nesting retained')
for(donor in unique(meta$donor_key)) lines(pc$x[meta$donor_key==donor,1:2],col='#999999',lty=2)
text(pc$x[,1],pc$x[,2],labels=meta$source_sample_suffix,pos=3,cex=.75)
legend('topright',legend=c('Tumor core','Host white matter'),col=c('#BB5544','#36708D'),pch=19,bty='n')
dev.off()
corx <- cor(x)
write.table(corx,file.path(out,paste0('All12_AllProbe_Array_Correlation_',date,'.tsv')),sep='\t',quote=FALSE,col.names=NA)
distances <- as.matrix(dist(t(x)))
write.table(distances,file.path(out,paste0('All12_AllProbe_Array_Distance_',date,'.tsv')),sep='\t',quote=FALSE,col.names=NA)
stopifnot(all(distances[lower.tri(distances)]>0))
receipt <- list(phase='B9_BLIND_QC_ONLY',date=date,raw_arrays=12,independent_patients=5,probe_family=22283,gene_identity_family=12405,
    raw_shape=dim(raw@assayData$exprs),raw_finite=all(qc$raw_nonfinite==0),RMA_finite=all(is.finite(x)),identical_expression_arrays=0,
    design_rows=nrow(design),design_columns=ncol(design),design_rank=qr(design)$rank,residual_df=4,
    source_replication_estimand='HISTORICAL_ADULT_TUMOR_CORE_MINUS_SAME_PATIENT_HOST_WHITE_MATTER',
    candidate_effects_evaluated=FALSE,automatic_sample_exclusion=FALSE,human_S6_status='PENDING',
    normalized_not_corrected_for_composition=TRUE,PCA_variance_percent=v[1:3],
    package_versions=sapply(c('affy','affyio','preprocessCore','hgu133acdf','limma','jsonlite'),function(p)as.character(packageVersion(p))))
write_json(receipt,file.path(out,paste0('Blind_QC_Execution_Receipt_',date,'.json')),pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file=file.path(out,paste0('R_sessionInfo_',date,'.txt')))
cat('B9 blind QC completed; no new candidate contrast tests.\n')
