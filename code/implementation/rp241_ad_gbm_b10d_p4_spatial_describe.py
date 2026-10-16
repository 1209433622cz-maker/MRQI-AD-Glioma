"""Accepted-scope descriptive patient-zone aggregation. No P or FDR."""
from pathlib import Path
from datetime import datetime,timezone,timedelta
import argparse,hashlib,json,gc
import numpy as np,pandas as pd,h5py
from scipy import sparse
R=Path('H:/SCI2/MRQI');D='2026-10-15'
O=R/'3_results/rp241_ad_gbm_b10d_p4_patient_spatial_description'
ZONES=['st_mal','dis_mal','st_norm','dis_norm','Intermediate','UNANNOTATED']
PAIRS=[('st_mal','st_norm'),('dis_mal','dis_norm'),('st_mal','dis_mal')]

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def js(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 receipt={'date':D,'phase':'RP241 B10D-P4','execute_requested':args.execute,'candidate_spatial_effects_computed':False}
 review=R/f'0_admin/audit/RP241_AD_GBM_B10D_P4_Spatial_S6_User_Reported_Review_{D}.json'
 if not review.exists():
  receipt['status']='HOLD_EXACT_SPATIAL_S6_REQUIRED';js(R/f'0_admin/validation/RP241_AD_GBM_B10D_P4_Execution_Preflight_{D}.json',receipt);return 2
 att=json.loads(review.read_text(encoding='utf-8'))
 assert att['status']=='USER_REPORTED_AUTHOR_REVIEW_ACCEPTED' and att['accepted_sources']==['GSE237183']
 assert att['reviewer_name'].strip() and att['review_date']=='2026-10-15'
 assert sha(R/att['review_memo_relative_path'])==att['review_memo_sha256']
 paths=[R/'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_2026-10-14.tsv',
 R/f'3_results/rp241_ad_gbm_b10d_p3_spatial_source/All25488_Local_Spots_Author_Join_NoDrop_{D}.tsv',
 R/f'3_results/rp241_ad_gbm_b10d_p3_spatial_source/All81_Author_MP_Gene_Overlap_Circularity_Audit_{D}.tsv',
 R/f'3_results/rp241_ad_gbm_b10d_p1_preparation/GSE237183_All13_Counts_Barcode_Position_QC_{D}.tsv',
 R/f'3_results/rp241_ad_gbm_b10d_p3_spatial_source/All13_Exact_Author_Annotation_Join_Readiness_{D}.tsv',
 R/f'0_admin/protocol/RP241/RP241_AD_GBM_B10D_P4_PatientAware_Descriptive_Spatial_Protocol_{D}.md',
 R/att['review_memo_relative_path'],review,Path(__file__).resolve()]
 old=pd.read_csv(R/f'0_admin/manifest/RP241_AD_GBM_B10D_P1_New_Source_Bytes_{D}.tsv',sep='\t')
 hh=old.loc[old.relative_path.str.endswith('.h5')];assert len(hh)==13
 for _,r in hh.iterrows():
  p=R/r.relative_path;assert p.stat().st_size==r.bytes and sha(p)==r.sha256;paths.append(p)
 ident=pd.read_csv(paths[0],sep='\t');assert len(ident)==81 and ident.primary43.sum()==43
 spots=pd.read_csv(paths[1],sep='\t',keep_default_na=False);assert len(spots)==25488
 assert spots.patient.nunique()==5 and spots.GSM.nunique()==13 and not spots.duplicated(['GSM','spot_id']).any()
 assert spots.author_metadata_present.eq(True).sum()==23358
 freeze=pd.DataFrame([dict(relative_path=p.relative_to(R).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in paths])
 man=R/f'0_admin/manifest/RP241_AD_GBM_B10D_P4_Analysis_Input_Freeze_{D}.tsv'
 if man.exists():assert pd.read_csv(man,sep='\t').to_dict('records')==freeze.to_dict('records')
 else:freeze.to_csv(man,sep='\t',index=False)
 receipt.update(status='READY_ACCEPTED_SCOPE_DESCRIPTIVE_ONLY',input_files=len(paths),input_manifest_sha256=sha(man),source_S6_review=att['user_attestation'])
 js(R/f'0_admin/validation/RP241_AD_GBM_B10D_P4_Execution_Preflight_{D}.json',receipt)
 if not args.execute:print(json.dumps(receipt,ensure_ascii=False,indent=2));return 0
 if O.exists():raise RuntimeError('Formal P4 output directory already exists; no overwrite')
 O.mkdir(parents=True)
 qc=pd.read_csv(paths[3],sep='\t');patients=sorted(spots.patient.unique());agg={}
 identrows=[];libs=[]
 for _,q in qc.iterrows():
  with h5py.File(R/'1_data/spatial/GSE237183/processed'/q.H5_name) as h:
   m=h['matrix'];names=[x.decode() for x in m['features/name'][:]];bc=[x.decode() for x in m['barcodes'][:]]
   counts=m['data'][:];assert np.issubdtype(counts.dtype,np.integer) and (counts>=0).all()
   mat=sparse.csc_matrix((counts.astype(np.int64),m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:]))
  ss=spots.loc[spots.GSM.eq(q.GSM)].set_index('spot_id').loc[bc]
  assert len(ss)==len(bc) and ss.patient.nunique()==1
  zone=ss.final_zone.replace('','UNANNOTATED').to_numpy();assert set(zone)<=set(ZONES)
  name_indices={}
  for i,name in enumerate(names):name_indices.setdefault(name,[]).append(i)
  mapping=[]
  for _,g in ident.iterrows():
   ix=name_indices.get(g.HGNC_symbol,[]);mapping.append(ix[0] if len(ix)==1 else None)
   identrows.append(dict(GSM=q.GSM,HGNC_ID=g.HGNC_ID,HGNC_symbol=g.HGNC_symbol,exact_feature_matches=len(ix),source_identity_status='EXACT_UNIQUE' if len(ix)==1 else 'NOT_TESTABLE_ABSENT' if len(ix)==0 else 'NOT_TESTABLE_AMBIGUOUS'))
  finite=[i for i,x in enumerate(mapping) if x is not None];sub=mat[[mapping[i] for i in finite],:]
  library=np.asarray(mat.sum(axis=0)).ravel().astype(np.int64)
  for z in ZONES:
   mask=zone==z;n=int(mask.sum());lib=int(library[mask].sum());c=np.zeros(81,dtype=np.int64);det=np.zeros(81,dtype=np.int64)
   if n:
    mm=sub[:,mask];c[finite]=np.asarray(mm.sum(axis=1)).ravel();det[finite]=np.asarray((mm>0).sum(axis=1)).ravel()
   key=(q.patient,z)
   if key not in agg:agg[key]=dict(n=0,lib=0,count=np.zeros(81,dtype=np.int64),det=np.zeros(81,dtype=np.int64),measured=np.ones(81,dtype=bool),specimens=[])
   v=agg[key];v['n']+=n;v['lib']+=lib;v['count']+=c;v['det']+=det
   if n:v['measured']&=np.asarray([x is not None for x in mapping]);v['specimens'].append(q.GSM)
   libs.append(dict(GSM=q.GSM,patient=q.patient,zone=z,spots=n,all_gene_UMI=lib,library_source='ALL_SOURCE_GENES_NOT_CANDIDATE_PANEL'))
  print(q.GSM,'full_count_scan_complete',len(bc),flush=True)
  del mat,sub;gc.collect()
 overlap=pd.read_csv(paths[2],sep='\t').set_index('HGNC_ID')
 raw=[]
 for patient in patients:
  for z in ZONES:
   v=agg[(patient,z)]
   for i,g in ident.iterrows():
    measured=bool(v['measured'][i]);n=v['n'];lib=v['lib'];cnt=int(v['count'][i]) if measured and n else None
    cpm=cnt/lib*1e6 if cnt is not None and lib>0 else np.nan
    raw.append(dict(patient=patient,zone=z,HGNC_ID=g.HGNC_ID,HGNC_symbol=g.HGNC_symbol,primary43=g.primary43,direction=g.direction,
      spots=n,source_specimens_with_spots=';'.join(v['specimens']),all_gene_UMI=lib,raw_gene_UMI=cnt,
      gene_detected_spots=int(v['det'][i]) if measured and n else None,
      CPM=cpm,log2_CPM1=float(np.log2(cpm+1)) if np.isfinite(cpm) else np.nan,
      source_identity_unique=measured,eligible_annotated_zone=n>=50 and z!='UNANNOTATED' and measured,
      status='ZONE_NOT_OBSERVED' if n==0 else 'NOT_TESTABLE_SOURCE_IDENTITY' if not measured else 'UNANNOTATED_SOURCE_CONTEXT_ONLY' if z=='UNANNOTATED' else 'SMALL_SPOT_DEPTH_CONTEXT_ONLY' if n<50 else 'PATIENT_ZONE_DESCRIPTIVE',
      author_MP_overlap_n=int(overlap.loc[g.HGNC_ID,'membership_n']),label_derived_from_same_RNA=True,clinical_age_treatment_unknown=True))
 raw=pd.DataFrame(raw);assert len(raw)==5*6*81
 raw.to_csv(O/f'All5Patient_All6Zone_All81_Counts_Depth_Description_{D}.tsv',sep='\t',index=False,na_rep='')
 pd.DataFrame(libs).to_csv(O/f'All13Specimen_All6Zone_AllGene_Library_Receipt_{D}.tsv',sep='\t',index=False)
 pd.DataFrame(identrows).to_csv(O/f'All13Specimen_All81_Exact_Source_Identity_{D}.tsv',sep='\t',index=False)
 pairs=[];summ=[];lopo=[]
 for za,zb in PAIRS:
  for _,g in ident.iterrows():
   vals=[];who=[]
   for patient in patients:
    ra=raw.loc[(raw.patient==patient)&(raw.zone==za)&(raw.HGNC_ID==g.HGNC_ID)].iloc[0]
    rb=raw.loc[(raw.patient==patient)&(raw.zone==zb)&(raw.HGNC_ID==g.HGNC_ID)].iloc[0]
    ok=bool(ra.eligible_annotated_zone and rb.eligible_annotated_zone)
    delta=ra.log2_CPM1-rb.log2_CPM1 if ok else np.nan
    pairs.append(dict(HGNC_ID=g.HGNC_ID,HGNC_symbol=g.HGNC_symbol,primary43=g.primary43,patient=patient,zone_a=za,zone_b=zb,
      spots_a=int(ra.spots),spots_b=int(rb.spots),eligible_pair=ok,delta_log2_CPM1=delta,status='DESCRIPTIVE_PAIR' if ok else 'NOT_ELIGIBLE_ZONE_DEPTH_OR_IDENTITY'))
    if ok:vals.append(delta);who.append(patient)
   aa=np.asarray(vals)
   summ.append(dict(HGNC_ID=g.HGNC_ID,HGNC_symbol=g.HGNC_symbol,primary43=g.primary43,direction=g.direction,zone_a=za,zone_b=zb,
     eligible_patients=len(vals),patient_ids=';'.join(who),equal_patient_mean=float(aa.mean()) if len(aa) else np.nan,
     median=float(np.median(aa)) if len(aa) else np.nan,range_min=float(aa.min()) if len(aa) else np.nan,range_max=float(aa.max()) if len(aa) else np.nan,
     positive_point_n=int((aa>0).sum()),negative_point_n=int((aa<0).sum()),zero_point_n=int((aa==0).sum()),
     status='SMALL_PATIENT_N_DESCRIPTIVE_ONLY__NO_P_OR_FDR',author_MP_overlap_n=int(overlap.loc[g.HGNC_ID,'membership_n']),independent_spatial_validation=False))
   for j,patient in enumerate(who):
    bb=np.delete(aa,j);lopo.append(dict(HGNC_ID=g.HGNC_ID,HGNC_symbol=g.HGNC_symbol,primary43=g.primary43,zone_a=za,zone_b=zb,removed_patient=patient,remaining_patient_N=len(bb),remaining_equal_patient_mean=float(bb.mean()) if len(bb) else np.nan))
 pd.DataFrame(pairs).to_csv(O/f'All81_All3Pair_All5Patient_Complete_Descriptions_{D}.tsv',sep='\t',index=False,na_rep='')
 pd.DataFrame(summ).to_csv(O/f'All81_All3Pair_EqualPatient_Descriptions_NoP_{D}.tsv',sep='\t',index=False,na_rep='')
 pd.DataFrame(lopo).to_csv(O/f'All81_All3Pair_LeaveOnePatient_Point_Influence_{D}.tsv',sep='\t',index=False,na_rep='')
 for x in freeze.to_dict('records'):assert sha(R/x['relative_path'])==x['sha256']
 receipt.update(status='COMPLETE_ACCEPTED_SCOPE_PATIENT_SPATIAL_DESCRIPTIONS__NO_NEW_P_FDR_OR_CAUSAL_UPGRADE',candidate_spatial_effects_computed=True,
  independent_patients=5,patient_zone_gene_rows=len(raw),pair_patient_gene_slots=len(pairs),pair_gene_summaries=len(summ),lopo_rows=len(lopo),
  new_P_or_FDR=False,clinical_gaps_imputed=False,label_is_pathology_truth=False,Tier1_causal_targets=0,
  all_source_UMI=sum(x['all_gene_UMI'] for x in libs),all_local_spots=sum(x['spots'] for x in libs),recorded_at_local=datetime.now(timezone(timedelta(hours=8))).isoformat())
 js(O/'Execution_Receipt.json',receipt);print(json.dumps(receipt,ensure_ascii=False,indent=2));return 0
if __name__=='__main__':raise SystemExit(main())
