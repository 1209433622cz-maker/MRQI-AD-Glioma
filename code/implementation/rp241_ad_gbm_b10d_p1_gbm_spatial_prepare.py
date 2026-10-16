"""CPTAC clinical/identity/model strata and all13 spatial technical joins only."""
from pathlib import Path
import csv,json,tarfile,io,openpyxl,collections,re,gzip,h5py,numpy as np
from scipy.sparse import csc_matrix
R=Path(r'H:\SCI2\MRQI');O=R/'3_results/rp241_ad_gbm_b10d_p1_preparation';O.mkdir(parents=True,exist_ok=True);P=R/'3_results/rp241_ad_gbm_b10d_p0_source_audit'
def read(p):return list(csv.DictReader(p.open(encoding='utf-8'),delimiter='\t'))
def write(n,rows):
 with (O/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
m=read(P/'All81_Fixed_HGNC_UniProt_Identity_2026-10-15.tsv');parent=read(P/'CPTAC_All_Aliquot_Clinical_Join_2026-10-15.tsv')
with tarfile.open(R/'1_data/proteomics/GBM/CPTAC_PDC000204/GBM_Metadata_CancerCell2021.tar.gz') as t:
 z=next(x for x in t.getmembers() if x.name.endswith('Clinical_Data_Feb2021_r2.xlsx'));w=openpyxl.load_workbook(io.BytesIO(t.extractfile(z).read()),read_only=True,data_only=True)
 rr=list(w['Patient_Clinical_Attributes'].iter_rows(values_only=True));clinical={x[1]:dict(zip(rr[0],x)) for x in rr[1:] if x[1]};rr=list(w['Biospecimen_Attributes'].iter_rows(values_only=True));spec={x[2]:dict(zip(rr[0],x)) for x in rr[1:] if x[2]};w.close()
with tarfile.open(R/'1_data/proteomics/GBM/CPTAC_PDC000204/Clinical_Data_Subtypes.tar.gz') as t:
 z=next(x for x in t.getmembers() if x.name.endswith('all_subtypes.v5.1.tsv'));rows=list(csv.DictReader(io.StringIO(t.extractfile(z).read().decode('utf-8')),delimiter='\t'));sub={x['case']:x for x in rows}
cohort=[]
for x in parent:
 c=clinical.get(x['case_id'],{});mtd=c.get('ancillary_studies_IDH1/2_mutational_testing_method') or '';flag=sub.get(x['case_id'],{}).get('IDH1_alter_status','UNRESOLVED');normal=x['sample_type']=='normal';withdrawn=x['effective_withdrawn']=='True'
 sequence=bool(re.search('sequenc',mtd,re.I));isihc=bool(re.search('IHC|immunohisto',mtd,re.I));negative=x['clinical_IDH_negative_candidate']=='True';adult=bool(c.get('age') is not None and float(c['age'])>=18)
 strict=bool(not normal and negative and adult and sequence and flag=='NA' and c.get('histologic_type') in ['Glioblastoma','Small cell  glioblastoma'])
 sampleids=[s.strip() for s in x['parent_sample_ids'].split(',')];ss=[spec.get(s,{}) for s in sampleids];methodlevel='SEQUENCE_BASED_REPORTED_NEGATIVE' if negative and sequence else 'IHC_REPORTED_NEGATIVE' if negative and isihc else 'OTHER_REPORTED_NEGATIVE' if negative and mtd not in ['','Unknown'] else 'METHOD_UNKNOWN_REPORTED_NEGATIVE' if negative else 'NOT_NEGATIVE_GROUP'
 cohort.append({**x,'histologic_type':c.get('histologic_type',''),'IDH_testing_method':mtd,'IDH_method_level':methodlevel,'author_IDH1_alter_status':flag,'flag_NA_proves_WT':False,'adult':adult,'specimens_total':len(sampleids),'specimens_joined':sum(bool(s) for s in ss),'source_specimen_meet_criteria':' | '.join(str(s.get('meet_criteria')) for s in ss),'source_specimen_histology_eligible':' | '.join(str(s.get('histological_type_eligible')) for s in ss),'source_specimen_tissue_types':' | '.join(str(s.get('tissue_type')) for s in ss),'main_limited_sequence_negative_set':strict,'healthy_GTEx_control':normal and not withdrawn,'expanded_clinical_negative_set':negative,'historical_nonwithdrawn_tumor_set':not normal and not withdrawn,'treatment_naive_evidence':'publication_cohort_level_only' if x['case_id'] in sub and not normal else 'NOT_RELEVANT_OR_UNRESOLVED','fully_WHO2021_reclassified':False})
for x in cohort:
 x['sequence_negative_specimen_QC_core']=bool(x['main_limited_sequence_negative_set'] and x['specimens_joined']==x['specimens_total'] and all(z=='Yes' for z in x['source_specimen_meet_criteria'].split(' | ')) and all(z=='Yes' for z in x['source_specimen_histology_eligible'].split(' | ')))
coreplex={x['plex'] for x in cohort if x['sequence_negative_specimen_QC_core']}
for x in cohort:x['primary_shared_plex_set']=x['sequence_negative_specimen_QC_core'] or (x['healthy_GTEx_control'] and x['plex'] in coreplex)
write('CPTAC_All110_Diagnostic_Specimen_Adjudication_NoPrune_2026-10-15.tsv',cohort)
designs=[]
for key in ['primary_shared_plex_set','main_limited_sequence_negative_set','expanded_clinical_negative_set','historical_nonwithdrawn_tumor_set']:
 rr=[x for x in cohort if x[key] or (key!='primary_shared_plex_set' and x['healthy_GTEx_control'])];plex=sorted({int(x['plex']) for x in rr});X=np.asarray([[1,int(x['sample_type']=='tumor'),float(x['age']),int(x['sex']=='Female'),*[int(int(x['plex'])==p) for p in plex[1:]]] for x in rr],float)
 rank=int(np.linalg.matrix_rank(X));source=np.asarray([int(x['sample_type']=='tumor') for x in rr])
 designs.append({'model_stratum':key,'tumor_N':sum(x['sample_type']=='tumor' for x in rr),'control_N':sum(x['sample_type']=='normal' for x in rr),'N':len(rr),'columns':X.shape[1],'rank':rank,'residual_df':len(rr)-rank,'source_acquisition_equals_disease_column':True,'rank_if_add_collection_source':int(np.linalg.matrix_rank(np.column_stack([X,source]))),'collection_adjusted_disease_effect_identifiable':False,'regression_contrast_is_limited_public_framework_only':True})
write('CPTAC_ThreeStrata_Design_Rank_Acquisition_Confounding_2026-10-15.tsv',designs)
for result in designs:
 key=result['model_stratum'];rr=[x for x in cohort if x[key] or (key!='primary_shared_plex_set' and x['healthy_GTEx_control'])];plex=sorted({int(x['plex']) for x in rr});X=np.asarray([[1,int(x['sample_type']=='tumor'),float(x['age']),int(x['sex']=='Female'),*[int(int(x['plex'])==p) for p in plex[1:]]] for x in rr],float)
 leverage=np.diag(X@np.linalg.pinv(X));result['max_leverage']=float(leverage.max());result['HC3_defined']=bool(leverage.max()<1-1e-8)
write('CPTAC_ThreeStrata_Design_Rank_Acquisition_Confounding_2026-10-15.tsv',designs)
hgnc={x['hgnc_id']:x for x in read(R/'1_data/reference/gene_identity/hgnc_complete_set_2026-09-04.tsv')};genes={x['HGNC_symbol'] for x in m};found={}
with (R/'1_data/proteomics/GBM/CPTAC_PDC000204/CPTAC3_Glioblastoma_Multiforme_Proteome.tmt11.tsv').open(encoding='utf-8') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['Gene'] in genes:found[x['Gene']]=x
coverage=[];mat=[];bygene={}
for x in m:
 y=found.get(x['HGNC_symbol']);entrez=hgnc[x['HGNC_ID']]['entrez_id'];exact=bool(y is not None and str(y['NCBIGeneID'])==entrez and entrez!='');values=[];uv=[]
 for c in cohort:
  def number(col):
   v=y.get(col,'') if y else ''
   try:return float(v)
   except (ValueError,TypeError):return np.nan
  values.append(number(c['aliquot_id']+' Log Ratio'));uv.append(number(c['aliquot_id']+' Unshared Log Ratio'))
 vv=np.asarray(values);uu=np.asarray(uv);main=[i for i,c in enumerate(cohort) if c['primary_shared_plex_set']];finite=int(np.isfinite(vv[main]).sum());status='MEASURABLE_IDENTITY_EXACT' if exact and finite==len(main) else 'PARTIAL_MISSING_OR_IDENTITY_UNRESOLVED' if y else 'ABSENT_FROM_SOURCE_FEATURES'
 coverage.append({**x,'source_NCBIGeneID':y['NCBIGeneID'] if y else '', 'HGNC_EntrezID':entrez,'symbol_and_Entrez_exact':exact,'main_set_finite_ordinary_N':finite,'main_set_finite_unshared_N':int(np.isfinite(uu[main]).sum()),'status':status})
 for kind,arr in [('LogRatio',vv),('Unshared',uu)]:mat.append({'HGNC_ID':x['HGNC_ID'],'HGNC_symbol':x['HGNC_symbol'],'quantity':kind,**{c['aliquot_id']:float(arr[i]) if exact and np.isfinite(arr[i]) else '' for i,c in enumerate(cohort)}})
write('CPTAC_All81_Identity_Missingness_NoPrune_2026-10-15.tsv',coverage);write('CPTAC_All81_BothQuantities_All110_Prepared_2026-10-15.tsv',mat)

# Exact spatial joins across original samples; no candidate expression contrasts.
sp=read(P/'Spatial_All_Sample_Donor_Phenotype_Role_Ledger_2026-10-15.tsv');sp=[x for x in sp if x['source']=='GSE237183' and x['group']=='IDHwt_GBM'];qcs=[];samplecov=[]
for x in sp:
 gsm=x['GSM'];folder=R/'1_data/spatial/GSE237183/processed';hp=next(folder.glob(gsm+'_*matrix.h5'))
 with h5py.File(hp,'r') as z:
  zz=z['matrix'];bar=[v.decode() for v in zz['barcodes'][:]];fid=[v.decode().split('.')[0] for v in zz['features/id'][:]];counts=csc_matrix((zz['data'][:],zz['indices'][:],zz['indptr'][:]),shape=tuple(zz['shape'][:]));tot=np.asarray(counts.sum(axis=0)).ravel();det=np.diff(counts.indptr)
 with gzip.open(next((R/'1_data/spatial/GSE237183/features').glob(gsm+'_*features.tsv.gz')),'rt',encoding='utf-8') as f:external_features=[row[0].split('.')[0] for row in csv.reader(f,delimiter='\t')]
 assert fid==external_features,(gsm,'H5/TSV feature order')
 with gzip.open(next(folder.glob(gsm+'_*barcodes.tsv.gz')),'rt',encoding='utf-8') as f:external=[line.strip() for line in f]
 assert bar==external,(gsm,'H5/TSV barcode order')
 with gzip.open(next(folder.glob(gsm+'_*tissue_positions*.gz')),'rt',encoding='utf-8') as f:pos=list(csv.reader(f))
 positions={row[0]:row[1:] for row in pos};assert len(positions)==len(pos)
 assert len(set(bar))==len(bar) and set(bar)<=set(positions)
 assert np.isfinite(counts.data).all() and np.all(counts.data>=0) and np.all(counts.data==np.floor(counts.data))
 in_tissue=[positions[b][0]=='1' for b in bar];assert all(in_tissue),(gsm,'filtered off-tissue barcode')
 title=x['title'];region='infiltrating' if 'infiltrating' in title else 'T1_enhancing' if 'T1 contrast' in title else 'necrotic' if 'necrotic' in title else 'bulk_unspecified' if 'bulk' in title else 'unannotated'
 qcs.append({'GSM':gsm,'patient':x['donor_key'],'author_region':region,'H5_name':hp.name,'genes':counts.shape[0],'filtered_spots':counts.shape[1],'H5_TSV_barcodes_equal':True,'all_spots_have_position':True,'position_barcode_duplicates':0,'filtered_off_tissue_spots':0,'median_UMI':float(np.median(tot)),'median_detected_features':float(np.median(det)),'nonnegative_integer_counts':True,'age':'NOT_RESOLVED','therapy_status':'NOT_RESOLVED','standalone_GEO_IDHwt_tissue_label':True,'old_GBM_source_patient_bridge':'NOT_RESOLVED','spatial_zone_annotations':'specimen_region_only__no_author_spot_zone_labels_frozen'})
 for g in m:
  ix=[i for i,k in enumerate(fid) if k==g['HGNC_Ensembl']];assert len(ix)<=1
  n=int(np.count_nonzero(counts[ix[0]].data)) if ix else 0
  samplecov.append({'GSM':gsm,'patient':x['donor_key'],'HGNC_ID':g['HGNC_ID'],'HGNC_symbol':g['HGNC_symbol'],'primary43':g['primary43'],'feature_present':bool(ix),'spots_with_nonzero_counts':n,'n_spots':counts.shape[1],'gene_expression_contrast_executed':False})
write('GSE237183_All13_Counts_Barcode_Position_QC_2026-10-15.tsv',qcs);write('GSE237183_All81_All13_Technical_Detection_NoPrune_2026-10-15.tsv',samplecov)
patients=[]
for patient in sorted({x['patient'] for x in qcs}):
 pp=[x for x in qcs if x['patient']==patient];regions={x['author_region'] for x in pp};patients.append({'patient':patient,'specimens':len(pp),'filtered_spots':sum(x['filtered_spots'] for x in pp),'author_regions':';'.join(sorted(regions)),'paired_infiltrating_vs_T1':{'infiltrating','T1_enhancing'}<=regions,'paired_infiltrating_vs_necrotic':{'infiltrating','necrotic'}<=regions,'independent_clinical_N':1,'age_treatment_or_independence_bridge_closed':False})
write('GSE237183_FivePatient_NestedRegion_Contrast_Readiness_2026-10-15.tsv',patients)
summary={'CPTAC':{'clinical52_method_levels':dict(collections.Counter(x['IDH_method_level'] for x in cohort if x['expanded_clinical_negative_set'])),'sequence_negative_limited_tumor_N':sum(x['main_limited_sequence_negative_set'] for x in cohort),'control_N':10,'all81_identity_exact':sum(x['symbol_and_Entrez_exact'] for x in coverage),'primary43_identity_exact':sum(x['symbol_and_Entrez_exact'] and x['primary43']=='True' for x in coverage),'designs':designs,'new_GBM_protein_effects':False},'Spatial':{'specimens':len(qcs),'patients':len(patients),'filtered_spots':sum(x['filtered_spots'] for x in qcs),'paired_infiltrating_vs_T1_patients':sum(x['paired_infiltrating_vs_T1'] for x in patients),'age_treatment_bridge_closed':False,'spot_level_author_zone_labels_closed':False,'new_candidate_spatial_effects':False}}
(O/'CPTAC_Spatial_Clinical_Model_Technical_Closure_2026-10-15.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(summary,ensure_ascii=False,indent=2))
