"""Re-read all44 matrices once; full RNA libraries; two fixed CNV masks."""
from pathlib import Path
import csv,gzip,json,tarfile,argparse,hashlib
import numpy as np,pandas as pd
from scipy.io import mmread
from threadpoolctl import threadpool_limits
from rp241_ad_gbm_b10b_count_localization import Aggregate
from rp241_ad_gbm_b10c_common import ROOT,DATE,MAP,MASKS,check_freeze,mask_index
def read_metadata(root):
 p=root/'1_data/single_cell/GBM/GSE182109/author_metadata/Meta_Data_GBMatlas.txt'
 with p.open(encoding='utf-8',newline='') as f:
  rr=csv.reader(f,delimiter='\t');hdr=next(rr);md=pd.DataFrame(rr,columns=['cell_id']+hdr)
 assert len(md)==201986 and md.cell_id.is_unique;md['author_state']=md.Assignment;md['broad_compartment']=md.Assignment.map(MAP)
 return md
def run(root):
 check_freeze(root);out=root/'3_results/rp241_ad_gbm_b10c_cnv_count'
 if out.exists() and any(out.iterdir()):raise RuntimeError('Preserve existing outputs')
 out.mkdir(parents=True);c=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_{DATE}.tsv',sep='\t')
 assert len(c)==81 and c.primary43.sum()==43;md=read_metadata(root)
 dm=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/GSE182109_18_Patient_Phenotype_Adjudication_{DATE}.tsv',sep='\t')
 roles=dm.set_index('donor_id').B10_role.map({'GBM_TUMOR_ONLY_LOCALIZATION':'MAIN_GBM','LGG_CONTEXT_ONLY':'LGG_CONTEXT','PHENOTYPE_CONFLICT_HOLD_GBM_MAIN':'PHENOTYPE_CONFLICT'}).to_dict()
 contexts=dm.set_index('donor_id').GEO_context.to_dict();aggs={m:Aggregate(c,'GSE182109',dm.donor_id,MAP.keys(),MAP.values(),roles,contexts) for m in MASKS}
 lm=pd.read_csv(root/'3_results/phase0w_c1/Phase0W_C1_GSE182109_DonorAware_Map_2026-09-28.tsv',sep='\t')
 ftfreeze=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/Tar_Triplet_Feature_Barcode_MatrixHeader_Inventory_{DATE}.tsv',sep='\t')
 old=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10b_state_localization/GSE182109_Exact_Matrix_Integer_Join_Receipt_{DATE}.tsv',sep='\t').set_index('GSM')
 rec=[]
 with tarfile.open(root/'1_data/single_cell/GBM/GSE182109/GSE182109_RAW.tar') as tf:
  look={}
  for mm in tf:
   if mm.isfile():
    gsm=Path(mm.name).name.split('_',1)[0];kind=next(k for k in ['matrix','features','barcodes'] if k in mm.name);assert (gsm,kind) not in look;look[gsm,kind]=mm
  assert len(look)==132
  for lr in lm.sort_values('gsm').itertuples():
   gsm=lr.gsm;d=lr.donor_id
   with gzip.GzipFile(fileobj=tf.extractfile(look[gsm,'features'])) as f:payload=f.read()
   ft=[x.split('\t') for x in payload.decode().splitlines()];assert all(x[2]=='Gene Expression' for x in ft)
   fr=ftfreeze[(ftfreeze.GSM==gsm)&(ftfreeze.source=='GSE182109')].iloc[0];assert hashlib.sha256(payload).hexdigest()==fr.feature_payload_sha256
   idmap={};sy={}
   for i,x in enumerate(ft):
    eid=x[0].split('.')[0];idmap.setdefault(eid,[]).append(i)
    if eid in sy:assert sy[eid]==x[1]
    sy[eid]=x[1]
   with gzip.GzipFile(fileobj=tf.extractfile(look[gsm,'barcodes'])) as f:bar=[x.decode().strip().split('\t')[0] for x in f if x.strip()]
   assert len(bar)==len(set(bar));bi={b:i for i,b in enumerate(bar)}
   with gzip.GzipFile(fileobj=tf.extractfile(look[gsm,'matrix'])) as f:
    with threadpool_limits(limits=1):mat=mmread(f).tocsr()
   assert mat.shape==(len(ft),len(bar)) and mat.nnz==int(fr.declared_nnz)
   assert np.all(np.isfinite(mat.data)) and np.all(mat.data>=0) and np.all(mat.data==np.floor(mat.data));mat=mat.astype(np.int64)
   sub=md[md.GSMID==gsm].copy();assert set(sub.Patient)=={d} and sub.barcode.is_unique;sub=sub[sub.barcode.isin(bi)].reset_index(drop=True)
   cols=np.array([bi[b] for b in sub.barcode]);joined=mat[:,cols];lib=np.asarray(joined.sum(axis=0)).ravel();counts=np.zeros((81,len(sub)),np.int64);present=np.zeros(81,bool)
   assert len(sub)==old.loc[gsm,'joined_cells'] and lib.sum()==old.loc[gsm,'source_library_UMI']
   for i,x in enumerate(c.itertuples()):
    inds=idmap.get(x.AD_external_Ensembl_stable_ID,[])
    if inds:assert sy[x.AD_external_Ensembl_stable_ID]==x.HGNC_symbol;counts[i]=np.asarray(joined[inds].sum(axis=0)).ravel();present[i]=True
   assert present.sum()==80 and c.iloc[np.flatnonzero(~present)].HGNC_symbol.tolist()==['FAAHP1']
   for mask,agg in aggs.items():
    ix=np.flatnonzero(mask_index(sub,mask));agg.add(d,sub.iloc[ix].reset_index(drop=True),lib[ix],counts[:,ix],present,gsm)
    rec.append({'mask':mask,'GSM':gsm,'donor_id':d,'role':roles[d],'exact_joined_before_mask':len(sub),'masked_cells':len(ix),
     'masked_full_source_UMI':int(lib[ix].sum()),'all_gene_library_denominator':True,'raw_shape_integer_feature_hash_verified':True})
   print(f'B10C two-mask full-source scan {len(rec)//2}/44: {gsm}',flush=True)
   del mat,joined,counts
 for mask,agg in aggs.items():
  dd=out/mask;dd.mkdir();agg.export(dd)
 pd.DataFrame(rec).to_csv(out/f'All44_Libraries_TwoMask_Exact_Count_Receipt_{DATE}.tsv',sep='\t',index=False)
 assert len(rec)==88
 receipt={'date':DATE,'status':'RAW_CNV_MASK_COUNTS_COMPLETE__NO_DISEASE_EFFECT','masks':MASKS,'all44libraries_verified':True,
  'all18roles_retained':True,'candidate81_and43_unchanged':True,'GSE211376_corresponding_CNV_prediction_available':False,
  'scope_only_author_Glioma_Myeloid_Oligo':True,'new_gene_P_or_GBM_scores':0,
  'mask_counts':{m:sum(x['masked_cells'] for x in rec if x['mask']==m) for m in MASKS}}
 (out/f'Execution_Receipt_{DATE}.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt,indent=2))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--root',default=str(ROOT));run(Path(a.parse_args().root))
