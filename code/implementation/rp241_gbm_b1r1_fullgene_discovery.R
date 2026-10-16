# New full-gene exploratory module. Original A1D quantification stays immutable.
argv<-commandArgs(trailingOnly=TRUE)
value<-function(k){i<-match(k,argv);if(is.na(i)||i==length(argv))stop(k);argv[i+1]}
root<-normalizePath(value('--root'),mustWork=TRUE);stamp<-'2026-10-14'
qa<-'--qa-only'%in%argv;execute<-'--execute'%in%argv
stopifnot(!(qa&&execute))
Sys.setenv(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
suppressPackageStartupMessages(library(DESeq2));suppressPackageStartupMessages(library(jsonlite))
stopifnot(as.character(getRversion())=='4.1.2',as.character(packageVersion('DESeq2'))=='1.34.0')
putjson<-function(x,p){if(file.exists(p))stop('Refuse overwrite: ',p);writeLines(toJSON(x,auto_unbox=TRUE,pretty=TRUE,na='null',digits=16),p,useBytes=TRUE)}
puttable<-function(x,p){if(file.exists(p))stop('Refuse overwrite: ',p);write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')}
hash<-function(p){stopifnot(file.exists(p));v<-system2('sha256sum',shQuote(p),stdout=TRUE);stopifnot(length(v)==1);substr(v,1,64)}
full_adjust<-function(p,method){answer<-rep(NA_real_,length(p));ok<-is.finite(p);answer[ok]<-p.adjust(p[ok],method=method,n=length(p));answer}
fitmodel<-function(k,m){
 stopifnot(!anyNA(k),all(k>=0),all(k==floor(k)),!anyDuplicated(rownames(k)),identical(colnames(k),rownames(m)))
 m$patient_id<-droplevels(factor(m$patient_id));m$tissue<-factor(m$tissue,levels=c('Normal','Tumor'))
 d<-DESeqDataSetFromMatrix(k,m,~patient_id+tissue)
 d<-DESeq(d,test='Wald',fitType='parametric',betaPrior=FALSE,minReplicatesForReplace=Inf,parallel=FALSE,quiet=TRUE)
 actual<-attr(dispersionFunction(d),'fitType')
 if(!identical(actual,'parametric'))stop('HOLD_DISPERSION_FIT_DRIFT: ',actual)
 stopifnot(identical(counts(d),k))
 z<-results(d,contrast=c('tissue','Tumor','Normal'),alpha=.05,independentFiltering=FALSE)
 converged<-mcols(d)$betaConv;if(is.null(converged))stop('beta convergence field unavailable')
 maxcook<-mcols(d)$maxCooks;if(is.null(maxcook))maxcook<-rep(NA_real_,nrow(k))
 valid<-is.finite(z$pvalue)&!is.na(converged)&converged
 rawp<-z$pvalue;p<-rawp;p[!valid]<-NA_real_
 status<-ifelse(is.na(converged)|!converged,'BETA_NOT_CONVERGED',ifelse(!is.finite(rawp),'P_MISSING_COOKS_OR_MODEL_DIAGNOSTIC','TESTED_FINITE_P'))
 rr<-data.frame(gene_id=rownames(k),baseMean=z$baseMean,log2FoldChange_Tumor_minus_Normal=z$log2FoldChange,lfcSE=z$lfcSE,Wald_statistic=z$stat,
  model_raw_pvalue=rawp,pvalue=p,BH_full_primary_family=full_adjust(p,'BH'),BY_full_primary_family=full_adjust(p,'BY'),
  CI95_lower_log2=z$log2FoldChange-qnorm(.975)*z$lfcSE,CI95_upper_log2=z$log2FoldChange+qnorm(.975)*z$lfcSE,
  beta_converged=converged,maxCooks=maxcook,p_value_status=status)
 list(dds=d,result=rr,fit=actual)
}
models<-function(m){
 answer<-list(M0_PAIRED15_FULL30=seq_len(nrow(m)),S1_GRADE4_PAIRED11=which(m$who_grade==4))
 for(pid in sort(unique(m$patient_id)))answer[[paste0('LOPO_',pid)]]<-which(m$patient_id!=pid)
 stopifnot(length(answer)==17);answer
}

if(qa){
 dest<-file.path(root,'0_admin/audit/RP241_AD_B6P_Synthetic_QA_2026-10-14');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 marker<-file.path(root,'0_admin/validation/RP241_AD_B6P_Synthetic_Statistical_QA_2026-10-14.json');stopifnot(!file.exists(marker))
 set.seed(20261014);g<-1200L;n<-30L
 m<-data.frame(patient_id=rep(sprintf('SYNTHETIC_P%02d',1:15),each=2),tissue=rep(c('Normal','Tumor'),15),who_grade=rep(c(rep(4,11),rep(3,4)),each=2))
 rownames(m)<-sprintf('SYNTHETIC_LIBRARY%02d',1:n)
 base<-exp(seq(log(50),log(1500),length.out=g));patient<-rep(exp(rnorm(15,0,.2)),each=2)
 mu<-outer(base,patient);mu[1:40,m$tissue=='Tumor']<-mu[1:40,m$tissue=='Tumor']*4;mu[41:80,m$tissue=='Tumor']<-mu[41:80,m$tissue=='Tumor']/4
 dispersion<-.03+10/base
 k<-matrix(rnbinom(g*n,mu=as.vector(mu),size=rep(1/dispersion,n)),g,n)
 low<-matrix(sample(0:3,20*n,replace=TRUE),20,n);boundary15<-matrix(c(rep(10L,15),rep(0L,15)),1,n);boundary14<-matrix(c(rep(10L,14),rep(0L,16)),1,n)
 k<-rbind(k,low,boundary15,boundary14);storage.mode(k)<-'integer';rownames(k)<-sprintf('SYNTHETIC_GENE%04d',seq_len(nrow(k)));colnames(k)<-rownames(m)
 keep<-rowSums(k>=10)>=15;stopifnot(all(keep[1:g]),!any(keep[(g+1):(g+20)]),keep[g+21],!keep[g+22]);universe<-rownames(k)[keep];k<-k[keep,,drop=FALSE]
 stopifnot(length(universe)==1201)
 fixture<-c(.001,.01,NA_real_);stopifnot(isTRUE(all.equal(full_adjust(fixture,'BH'),c(.003,.015,NA_real_))),isTRUE(all.equal(full_adjust(fixture,'BY'),c(.0055,.0275,NA_real_))))
 ms<-models(m);checks<-list();orientation<-NULL
 for(name in names(ms)){
  ii<-ms[[name]];mm<-m[ii,,drop=FALSE];x<-model.matrix(~patient_id+tissue,transform(mm,patient_id=factor(patient_id),tissue=factor(tissue,levels=c('Normal','Tumor'))))
  rank<-qr(x)$rank;expected<-if(name=='M0_PAIRED15_FULL30')16 else if(name=='S1_GRADE4_PAIRED11')12 else 15;stopifnot(rank==expected)
  cat(format(Sys.time()),'SYNTHETIC_START',name,'\n');flush.console()
  fit<-fitmodel(k[,ii,drop=FALSE],mm);z<-fit$result;stopifnot(identical(z$gene_id,universe),nrow(z)==1201)
  up<-mean(z$log2FoldChange_Tumor_minus_Normal[1:40]>0);down<-mean(z$log2FoldChange_Tumor_minus_Normal[41:80]<0);stopifnot(up>=.9,down>=.9)
  if(name=='M0_PAIRED15_FULL30'){
   rev<-results(fit$dds,contrast=c('tissue','Normal','Tumor'),alpha=.05,independentFiltering=FALSE);ok<-is.finite(rev$log2FoldChange)&is.finite(z$log2FoldChange_Tumor_minus_Normal)
   stopifnot(max(abs(rev$log2FoldChange[ok]+z$log2FoldChange_Tumor_minus_Normal[ok]))<1e-10)
   pok<-is.finite(rev$pvalue)&is.finite(z$model_raw_pvalue);stopifnot(max(abs(rev$pvalue[pok]-z$model_raw_pvalue[pok]))<1e-10)
   orientation<-list(reverse_contrast_beta_sign_pass=TRUE,reverse_contrast_P_identical=TRUE)
   puttable(z,file.path(dest,'SYNTHETIC_Primary_All1201_Result.tsv'))
  }
  checks[[name]]<-data.frame(model=name,libraries=length(ii),patients=length(unique(mm$patient_id)),rank=rank,residual_df=length(ii)-rank,full_gene_family=1201,positive_fixture_correct_fraction=up,negative_fixture_correct_fraction=down,fit_type=fit$fit,count_replacement=FALSE)
  cat(format(Sys.time()),'SYNTHETIC_COMPLETE',name,'\n');flush.console();rm(fit);gc()
 }
 puttable(do.call(rbind,checks),file.path(dest,'SYNTHETIC_All17_Models_QA.tsv'))
 putjson(list(status='PASS',time=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),R_version=as.character(getRversion()),DESeq2_version=as.character(packageVersion('DESeq2')),seed=20261014,raw_synthetic_genes=1222,primary_filtered_family=1201,models=17,filter_boundary15_pass=TRUE,filter_boundary14_excluded=TRUE,full_BH_BY_missing_denominator_manual_fixture_pass=TRUE,orientation=orientation,actual_human_gene_effect_tests=0,count_replacement=FALSE),marker)
 writeLines(capture.output(sessionInfo()),file.path(dest,'R_sessionInfo.txt'));cat('SYNTHETIC_ALL17_MODELS_QA_PASS\n');quit(status=0)
}

repair_freeze<-fromJSON(file.path(root,'0_admin/protocol/RP241/RP241_GBM_B1R1_InputGuard_Repair_Freeze_2026-10-14.json'),simplifyVector=FALSE)
for(p in names(repair_freeze$sha256))stopifnot(identical(hash(file.path(root,p)),repair_freeze$sha256[[p]]))
completion_files<-list.files(file.path(root,'3_results/rp241_a1d/quantification'),pattern='^completion.json$',full.names=TRUE,recursive=TRUE)
required<-c('3_results/rp241_a1d/GSE165595_GBM15_GENCODE27_Fragment_Counts.tsv','3_results/rp241_a1d/Gene_Filter.tsv','3_results/rp241_a1d/Full_Library_QC.tsv','3_results/rp241_a1d/Blind_VST_PCA_Coordinates.tsv','3_results/rp241_a1d/VST_Sample_Correlation.tsv','3_results/rp241_a1d/Primary_Design_Matrix.tsv','3_results/rp241_a1d/Grade4_Design_Matrix.tsv','3_results/rp241_a1d/multiqc/multiqc_report.html','0_admin/decision/RP241_A1D_PostQuantification_Decision_2026-10-12.json')
missing<-required[!file.exists(file.path(root,required))]
admission<-file.path(root,'0_admin/decision/RP241_GBM_B1_FullQC_Admission.json')
if(length(completion_files)!=30||length(missing)||!file.exists(admission)){
 attempt<-format(Sys.time(),'%Y%m%dT%H%M%S%z')
 p<-file.path(root,paste0('0_admin/validation/RP241_AD_B6P_Real_Input_Preflight_2026-10-14_',attempt,'.json'))
 putjson(list(status='HOLD_FULL30_COUNTS_AND_QC_ADMISSION_REQUIRED__NO_REAL_GBM_EFFECT_OUTPUT',completed_libraries=length(completion_files),expected=30,missing_paths=missing,full_QC_admission_present=file.exists(admission),real_effect_tests=0),p)
 cat('HOLD_FULL30_COUNTS_AND_QC_ADMISSION_REQUIRED\n');quit(status=if(execute)2 else 0)
}
q<-fromJSON(admission);stopifnot(q$status=='GO_GBM_PAIRED_EXPLORATORY_FULL_GENE_DISCOVERY',isTRUE(q$exploratory_execution_allowed),identical(hash(file.path(root,required[1])),q$count_matrix_sha256),identical(hash(file.path(root,required[2])),q$gene_filter_sha256))
post<-fromJSON(file.path(root,required[9]));stopifnot(isTRUE(post$counts_complete),post$samples==30,post$patients==15,identical(hash(file.path(root,required[1])),post$count_matrix_sha256))
fr<-fromJSON(file.path(root,'0_admin/protocol/RP241/RP241_AD_B6P_Execution_Freeze_2026-10-14.json'),simplifyVector=FALSE)
for(p in names(fr$sha256))stopifnot(identical(hash(file.path(root,p)),fr$sha256[[p]]))
old<-fromJSON(file.path(root,'0_admin/manifest/RP241_A1D_Execution_Contract_2026-10-12.json'),simplifyVector=FALSE)
for(p in names(old$sha256))stopifnot(identical(hash(file.path(root,p)),old$sha256[[p]]))
m<-read.delim(file.path(root,'0_admin/manifest/RP241_A1D_Samples_2026-10-12.tsv'),check.names=FALSE);rownames(m)<-m$run
stopifnot(nrow(m)==30,length(unique(m$patient_id))==15,setequal(basename(dirname(completion_files)),m$run))
expected_config<-hash(file.path(root,'0_admin/protocol/RP241/RP241_A1D_Execution_Config_2026-10-12.json'))
for(p in completion_files){
 d<-fromJSON(p,simplifyVector=FALSE);run<-basename(dirname(p));expected<-m[run,,drop=FALSE]
 stopifnot(d$status=='PASS',identical(d$config_sha256,expected_config),identical(d$sample$run,run),identical(d$sample$patient_id,expected$patient_id),identical(d$sample$tissue,expected$tissue))
 for(k in names(d$outputs_sha256))stopifnot(identical(hash(file.path(dirname(p),k)),d$outputs_sha256[[k]]))
 for(k in names(d$fastqc_sha256))stopifnot(identical(hash(file.path(dirname(p),'fastqc',k)),d$fastqc_sha256[[k]]))
}
all_counts<-read.delim(file.path(root,required[1]),row.names=1,check.names=FALSE)
identity<-read.delim(file.path(root,'3_results/rp241_ad_b6p_discovery_preflight/GBM_All58288_Reference_Identity_Ledger_2026-10-14.tsv'),check.names=FALSE)
stopifnot(nrow(all_counts)==58288,ncol(all_counts)==30,identical(colnames(all_counts),m$run),identical(rownames(all_counts),identity$GENCODE27_gene_ID),all(as.matrix(all_counts)>=0),all(as.matrix(all_counts)==floor(as.matrix(all_counts))))
gf<-read.delim(file.path(root,required[2]),check.names=FALSE);keep<-rowSums(all_counts>=10)>=15;stopifnot(identical(gf$gene_id,rownames(all_counts)),identical(gf$kept,unname(keep)),sum(keep)>1000)
k<-as.matrix(all_counts[keep,,drop=FALSE]);storage.mode(k)<-'integer';universe<-rownames(k)
if(!execute){cat('FULL_INPUT_QC_AND_FROZEN_MODEL_PREFLIGHT_PASS; no effect execution requested\n');quit(status=0)}
stopifnot(fromJSON(file.path(root,'0_admin/validation/RP241_AD_B6P_Synthetic_Statistical_QA_2026-10-14.json'))$status=='PASS')
dest<-file.path(root,'3_results/rp241_gbm_b1_discovery');if(dir.exists(dest)&&length(list.files(dest,all.files=TRUE,no..=TRUE)))stop('Preserve existing discovery attempt')
if(!dir.create(dest,recursive=FALSE,showWarnings=FALSE))stop('Discovery directory creation failed or another execution owns it')
puttable(data.frame(gene_id=universe,source_primary_measurable=TRUE,filter='count10_in15of30'),file.path(dest,'Primary_Fixed_Gene_Universe.tsv'))
ms<-models(m);summaries<-list()
for(name in names(ms)){
 ii<-ms[[name]];mm<-m[ii,,drop=FALSE];rank<-qr(model.matrix(~patient_id+tissue,transform(mm,patient_id=factor(patient_id),tissue=factor(tissue,levels=c('Normal','Tumor')))))$rank
 expected<-if(name=='M0_PAIRED15_FULL30')16 else if(name=='S1_GRADE4_PAIRED11')12 else 15;stopifnot(rank==expected)
 cat(format(Sys.time()),'REAL_EXPLORATORY_START',name,'\n');flush.console()
 fit<-tryCatch(fitmodel(k[,ii,drop=FALSE],mm),error=function(e){putjson(list(status='HOLD_MODEL_TECHNICAL_FAILURE',model=name,message=conditionMessage(e)),file.path(dest,paste0(name,'_Technical_Hold.json')));stop(e)})
 z<-fit$result;stopifnot(identical(z$gene_id,universe));z$evidence_scope<-'ONE_HISTORICAL_GRADE3_4_PAIRED_TUMOR_HOST_SOURCE_EXPLORATORY'
 puttable(z,file.path(dest,paste0(name,'_AllPrimaryGene_Results.tsv')))
 puttable(data.frame(run=colnames(k)[ii],size_factor=sizeFactors(fit$dds)),file.path(dest,paste0(name,'_Size_Factors.tsv')))
 saveRDS(fit$dds,file.path(dest,paste0(name,'_DESeq2_Model.rds')),compress=TRUE)
 summaries[[name]]<-data.frame(model=name,libraries=length(ii),patients=length(unique(mm$patient_id)),rank=rank,df=length(ii)-rank,family=length(universe),finite_P=sum(is.finite(z$pvalue)),missing_P=sum(!is.finite(z$pvalue)),actual_fit=fit$fit)
 cat(format(Sys.time()),'REAL_EXPLORATORY_COMPLETE',name,'\n');flush.console();rm(fit);gc()
}
puttable(do.call(rbind,summaries),file.path(dest,'All17_Models_Execution_Summary.tsv'))
putjson(list(status='FULL17_GBM_SOURCE_MODELS_COMPLETE__EXPLORATORY_ONLY',family=length(universe),models=17,joint_AD_GBM_G_frozen=FALSE,independent_AD_effect_unlock=FALSE,healthy_population_claim=FALSE),file.path(dest,'Execution_Receipt.json'))
