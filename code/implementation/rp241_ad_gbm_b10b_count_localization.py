"""B10B full-source RNA denominators and donor-aware complete81 localisation.
No disease contrast, score fitting, cell-level P or candidate selection.
"""
from __future__ import annotations
import argparse,csv,gzip,hashlib,json,tarfile
from pathlib import Path
from collections import Counter
import numpy as np,pandas as pd
from scipy import sparse
from scipy.io import mmread
from threadpoolctl import threadpool_limits
DATE='2026-10-14'
MAP211={'AC-like':'malignant','MES-like':'malignant','NPC-like':'malignant','OPC-like':'malignant',
 'TAM-MG':'myeloid','TAM-BDM':'myeloid','DC':'myeloid','CD4/CD8':'lymphoid',
 'Endothelial':'vascular_mixed','Pericyte':'vascular_mixed','Astrocyte':'astrocyte_nonmalignant',
 'Neuron':'neuron_nonmalignant','Oligodendrocyte':'oligodendrocyte','OPC':'opc_nonmalignant'}
MAP182={'Glioma':'malignant','Myeloid':'myeloid','Oligo':'oligodendrocyte','Endo':'vascular_mixed',
 'Pericytes':'vascular_mixed','TCells':'lymphoid','BCells':'lymphoid','Other':'other_unresolved'}
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(16<<20),b''):h.update(b)
 return h.hexdigest()
def check_freeze(root):
 p=root/f'0_admin/protocol/RP241/RP241_AD_GBM_B10B_PreExpression_Input_Code_Freeze_{DATE}.json'
 v2=root/f'0_admin/protocol/RP241/RP241_AD_GBM_B10B_Implementation_Repair_Input_Code_Freeze_V2_{DATE}.json'
 if v2.exists():p=v2
 j=json.loads(p.read_text(encoding='utf-8'))
 for r in j['inputs']:
  p=root/r['relative_path'];assert p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],r['relative_path']
 return j

class Aggregate:
 def __init__(self,cand,source,donors,states,broad,roles,contexts):
  self.cand=cand;self.source=source;self.donors=sorted(donors);self.roles=roles;self.contexts=contexts
  self.groups=[(d,level,s) for d in self.donors for level,ss in [('author_state',sorted(set(states))),('broad_compartment',sorted(set(broad)))] for s in ss]
  assert len(self.groups)==len(set(self.groups))
  self.gi={k:i for i,k in enumerate(self.groups)};n=len(self.groups)
  self.cells=np.zeros(n,np.int64);self.lib=np.zeros(n,np.int64)
  self.count=np.zeros((81,n),np.int64);self.detect=np.zeros((81,n),np.int64)
  self.present=None;self.region_rows=[];self.lib_receipts=[]
 def add(self,donor,meta,cell_lib,cand_count,present,unit):
  if self.present is None:self.present=np.array(present,bool)
  assert np.array_equal(present,self.present),'Nonuniform identity coverage: no silent zero filling'
  assert cand_count.shape==(81,len(meta)) and len(cell_lib)==len(meta)
  assert np.all(cand_count>=0) and np.all(cand_count.sum(axis=0)<=cell_lib)
  for level,col in [('author_state','author_state'),('broad_compartment','broad_compartment')]:
   for state,inds in meta.groupby(col,sort=True).indices.items():
    inds=np.asarray(inds);j=self.gi[donor,level,state];n=len(inds);library=int(cell_lib[inds].sum())
    sums=cand_count[:,inds].sum(axis=1,dtype=np.int64);det=(cand_count[:,inds]>0).sum(axis=1)
    self.cells[j]+=n;self.lib[j]+=library;self.count[:,j]+=sums;self.detect[:,j]+=det
    for i,r in enumerate(self.cand.itertuples()):
     ok=bool(present[i]);cp=float(sums[i])/library*1e6 if library and ok else np.nan
     self.region_rows.append({'source':self.source,'unit':unit,'donor_id':donor,'role':self.roles[donor],'clinical_context':self.contexts[donor],
      'level':level,'state':state,'HGNC_ID':r.HGNC_ID,'HGNC_symbol':r.HGNC_symbol,'primary43':bool(r.primary43),
      'cells':n,'source_library_UMI':library,'identity_present':ok,'raw_UMI':int(sums[i]) if ok else np.nan,
      'CPM':cp,'detection_cells':int(det[i]) if ok else np.nan,'detection_fraction':float(det[i])/n if n and ok else np.nan})
  self.lib_receipts.append({'source':self.source,'unit':unit,'donor_id':donor,'joined_cells':len(meta),
   'source_library_UMI':int(cell_lib.sum()),'candidate81_UMI_sum':int(cand_count.sum()),
   'candidate_count_is_denominator':False})
 def export(self,out):
  rows=[];q=[]
  for j,(d,level,state) in enumerate(self.groups):
   n=int(self.cells[j]);lib=int(self.lib[j]);role=self.roles[d]
   q.append({'source':self.source,'donor_id':d,'role':role,'clinical_context':self.contexts[d],
    'level':level,'state':state,'cells':n,'source_library_UMI':lib,'eligible_50':n>=50,'exploratory_20_49':20<=n<50})
   for i,r in enumerate(self.cand.itertuples()):
    ok=bool(self.present[i]);obs=n>0 and ok and lib>0;cp=float(self.count[i,j])/lib*1e6 if obs else np.nan
    rows.append({'source':self.source,'donor_id':d,'role':role,'clinical_context':self.contexts[d],
     'level':level,'state':state,'HGNC_ID':r.HGNC_ID,'HGNC_symbol':r.HGNC_symbol,'Ensembl_stable_ID':r.AD_external_Ensembl_stable_ID,
     'direction':r.direction,'primary43':bool(r.primary43),'cells':n,'source_library_UMI':lib,
     'eligible_50':n>=50,'identity_present':ok,'raw_UMI':int(self.count[i,j]) if obs else np.nan,
     'CPM':cp,'log2_CPM_plus1':np.log2(cp+1) if obs else np.nan,
     'detection_cells':int(self.detect[i,j]) if obs else np.nan,'detection_fraction':float(self.detect[i,j])/n if obs else np.nan,
     'status':'IDENTITY_NOT_TESTABLE' if not ok else 'STATE_NOT_OBSERVED' if not n else 'OBSERVED_ZERO' if not self.count[i,j] else 'OBSERVED_POSITIVE',
     'GBM_vs_normal_effect':False})
  assert sum(self.lib[j] for j,g in enumerate(self.groups) if g[1]=='author_state')==sum(r['source_library_UMI'] for r in self.lib_receipts)
  assert sum(self.lib[j] for j,g in enumerate(self.groups) if g[1]=='broad_compartment')==sum(r['source_library_UMI'] for r in self.lib_receipts)
  def save(x,name):pd.DataFrame(x).to_csv(out/f'{self.source}_{name}_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
  save(rows,'All81_DonorState_Complete_Localization');save(q,'DonorState_Depth_Library_QC')
  save(self.region_rows,'All81_RegionState_Complete_Localization');save(self.lib_receipts,'Region_Library_Conservation')
  return pd.DataFrame(rows)

def source211(root,cand,out):
 src=root/'1_data/single_cell/GBM/GSE211376';rs=[]
 with gzip.open(src/'GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz','rt',newline='') as f:
  assert next(f).strip()=='predicted.high_hierarchy,patient'
  for r in csv.reader(f):assert len(r)==3;rs.append(r)
 md=pd.DataFrame(rs,columns=['barcode','author_state','donor_id']);assert md.barcode.is_unique and len(md)==39355
 md['broad_compartment']=md.author_state.map(MAP211);assert md.broad_compartment.notna().all()
 dm=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/GSE211376_11_Donor_Public_Clinical_Map_{DATE}.tsv',sep='\t')
 assert set(dm.donor_id)==set(md.donor_id) and len(dm)==11
 agg=Aggregate(cand,'GSE211376',dm.donor_id,MAP211.keys(),MAP211.values(),{d:'MAIN_GBM' for d in dm.donor_id},{d:'primary_IDHwt_GBM_author_GEO' for d in dm.donor_id})
 wanted={s:i for i,s in enumerate(cand.HGNC_symbol)};assert len(wanted)==81
 p=src/'GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz'
 with gzip.open(p,'rt') as f:
  bar=next(f).rstrip('\r\n').split('\t');assert len(bar)==len(set(bar)) and set(bar)==set(md.barcode)
  md=md.set_index('barcode').loc[bar].reset_index()
  lib=np.zeros(len(bar),np.int64);counts=np.zeros((81,len(bar)),np.int64);found=np.zeros(81,bool);features=set();ng=0
  for line in f:
   pos=line.find('\t');assert pos>0;gene=line[:pos];assert gene not in features;features.add(gene)
   v=np.fromstring(line[pos+1:],dtype=np.int64,sep='\t');assert len(v)==len(bar) and np.all(v>=0)
   lib+=v;ng+=1
   if gene in wanted:i=wanted[gene];counts[i]=v;found[i]=True
   if ng%5000==0:print(f'GSE211376 full RNA denominator: {ng} gene rows',flush=True)
 assert found.sum()==78 and {cand.HGNC_symbol[i] for i in np.where(~found)[0]}=={'FLCN','OR7A5','FAAHP1'}
 for d,inds in md.groupby('donor_id',sort=True).indices.items():
  inds=np.asarray(inds);agg.add(d,md.iloc[inds].reset_index(drop=True),lib[inds],counts[:,inds],found,d)
 frame=agg.export(out)
 rec={'source':'GSE211376','feature_rows':ng,'source_cells':len(bar),'joined_cells':len(bar),'candidate_identity_present':int(found.sum()),
  'total_source_UMI':int(lib.sum()),'full_source_integer_nonnegative_width_eof':True,'single_cell_or_disease_P_generated':False}
 (out/f'GSE211376_Source_Count_Receipt_{DATE}.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
 return frame,rec

def source182(root,cand,out):
 src=root/'1_data/single_cell/GBM/GSE182109'
 with open(src/'author_metadata/Meta_Data_GBMatlas.txt',encoding='utf-8',newline='') as f:
  rr=csv.reader(f,delimiter='\t');hdr=next(rr);md=pd.DataFrame(rr,columns=['cell_id']+hdr)
 assert md.cell_id.is_unique and len(md)==201986
 md['author_state']=md.Assignment;md['broad_compartment']=md.Assignment.map(MAP182);assert md.broad_compartment.notna().all()
 ds=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/GSE182109_18_Patient_Phenotype_Adjudication_{DATE}.tsv',sep='\t')
 rolemap={'GBM_TUMOR_ONLY_LOCALIZATION':'MAIN_GBM','LGG_CONTEXT_ONLY':'LGG_CONTEXT','PHENOTYPE_CONFLICT_HOLD_GBM_MAIN':'PHENOTYPE_CONFLICT'}
 roles=ds.set_index('donor_id').B10_role.map(rolemap).to_dict();contexts=ds.set_index('donor_id').GEO_context.to_dict()
 assert list(roles.values()).count('MAIN_GBM')==14
 agg=Aggregate(cand,'GSE182109',ds.donor_id,MAP182.keys(),MAP182.values(),roles,contexts)
 lm=pd.read_csv(root/'3_results/phase0w_c1/Phase0W_C1_GSE182109_DonorAware_Map_2026-09-28.tsv',sep='\t')
 frozenft=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/Tar_Triplet_Feature_Barcode_MatrixHeader_Inventory_{DATE}.tsv',sep='\t')
 receipts=[];cnv=[];substates=[]
 with tarfile.open(src/'GSE182109_RAW.tar','r') as tf:
  members=[m for m in tf if m.isfile()];assert len(members)==132;look={}
  for m in members:
   gsm=Path(m.name).name.split('_',1)[0];kind=next((k for k in ['features','matrix','barcodes'] if k in m.name),None)
   if kind:assert (gsm,kind) not in look;look[gsm,kind]=m
  for lr in lm.sort_values('gsm').itertuples():
   gsm=lr.gsm;d=lr.donor_id
   with gzip.GzipFile(fileobj=tf.extractfile(look[gsm,'features'])) as f:payload=f.read()
   ft=[x.split('\t') for x in payload.decode().splitlines()]
   assert all(len(x)>=3 and x[2]=='Gene Expression' for x in ft)
   ids=[x[0].split('.')[0] for x in ft];names=[x[1] for x in ft]
   fr=frozenft[(frozenft.source=='GSE182109')&(frozenft.GSM==gsm)].iloc[0]
   assert hashlib.sha256(payload).hexdigest()==fr.feature_payload_sha256
   idmap={};symmap={}
   for i,(eid,sym) in enumerate(zip(ids,names)):
    idmap.setdefault(eid,[]).append(i)
    if eid in symmap:assert symmap[eid]==sym,'duplicate stable-ID symbols conflict'
    symmap[eid]=sym
   with gzip.GzipFile(fileobj=tf.extractfile(look[gsm,'barcodes'])) as f:bar=[x.decode().strip().split('\t')[0] for x in f if x.strip()]
   assert len(bar)==len(set(bar))
   with gzip.GzipFile(fileobj=tf.extractfile(look[gsm,'matrix'])) as f:
    with threadpool_limits(limits=1):mat=mmread(f).tocsr()
   assert mat.shape==(len(ft),len(bar)) and mat.nnz==int(fr.declared_nnz)
   assert np.all(np.isfinite(mat.data)) and np.all(mat.data>=0) and np.all(mat.data==np.floor(mat.data))
   mat=mat.astype(np.int64)
   sub=md[md.GSMID==gsm].copy();assert sub.barcode.is_unique and set(sub.Patient)=={d}
   bi={b:i for i,b in enumerate(bar)};sub=sub[sub.barcode.isin(bi)].reset_index(drop=True)
   cols=np.array([bi[x] for x in sub.barcode]);joined=mat[:,cols]
   lib=np.asarray(joined.sum(axis=0)).ravel().astype(np.int64)
   counts=np.zeros((81,len(cols)),np.int64);present=np.zeros(81,bool)
   for i,r in enumerate(cand.itertuples()):
    ri=idmap.get(r.AD_external_Ensembl_stable_ID,[])
    if ri:
     assert symmap[r.AD_external_Ensembl_stable_ID]==r.HGNC_symbol
     counts[i]=np.asarray(joined[ri,:].sum(axis=0)).ravel();present[i]=True
   assert present.sum()==80 and cand.iloc[np.where(~present)[0]].HGNC_symbol.tolist()==['FAAHP1']
   agg.add(d,sub,lib,counts,present,gsm)
   receipts.append({'GSM':gsm,'donor_id':d,'role':roles[d],'source_features':len(ft),'matrix_cells':len(bar),'metadata_cells':int((md.GSMID==gsm).sum()),
    'joined_cells':len(sub),'raw_only_excluded':len(bar)-len(sub),'metadata_only_excluded':int((md.GSMID==gsm).sum())-len(sub),
    'source_nnz':int(mat.nnz),'source_library_UMI':int(lib.sum()),'candidate_identity_present':80})
   for (state,pred),z in sub.groupby(['Assignment','CopyKatPrediction'],dropna=False):cnv.append({'GSM':gsm,'donor_id':d,'role':roles[d],'Assignment':state,'CopyKatPrediction':pred,'cells':len(z)})
   for (state,ss),z in sub.groupby(['Assignment','SubAssignment'],dropna=False):substates.append({'GSM':gsm,'donor_id':d,'role':roles[d],'Assignment':state,'SubAssignment':ss,'cells':len(z)})
   del mat,joined,counts
   print(f'GSE182109 exact raw matrix {len(receipts)}/44: {gsm}',flush=True)
 assert sum(r['joined_cells'] for r in receipts)==201893
 assert sum(r['raw_only_excluded'] for r in receipts)==63058 and sum(r['metadata_only_excluded'] for r in receipts)==93
 frame=agg.export(out)
 for rows,name in [(receipts,'Exact_Matrix_Integer_Join_Receipt'),(cnv,'Author_Assignment_CopyKat_Audit'),(substates,'Author_SubAssignment_Nesting_Audit')]:
  pd.DataFrame(rows).to_csv(out/f'GSE182109_{name}_{DATE}.tsv',sep='\t',index=False)
 rec={'source':'GSE182109','patients_all':18,'main_GBM':14,'LGG_context':2,'phenotype_conflict':2,'libraries':44,
  'source_cells':sum(r['matrix_cells'] for r in receipts),'joined_cells':201893,'candidate_identity_present':80,
  'total_source_UMI':sum(r['source_library_UMI'] for r in receipts),'all44_matrix_numeric_integer_shape_nnzs':True,
  'single_cell_or_disease_P_generated':False}
 (out/f'GSE182109_Source_Count_Receipt_{DATE}.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
 return frame,rec

def run(root):
 freeze=check_freeze(root)
 out=root/'3_results/rp241_ad_gbm_b10b_state_localization'
 if out.exists() and any(out.iterdir()):raise RuntimeError('Preserve existing outputs; do not overwrite')
 out.mkdir(parents=True,exist_ok=True)
 c=pd.read_csv(root/f'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_{DATE}.tsv',sep='\t')
 assert len(c)==81 and c.primary43.sum()==43 and c.HGNC_ID.is_unique
 a,ra=source211(root,c,out);b,rb=source182(root,c,out)
 pd.concat([a,b],ignore_index=True).to_csv(out/f'All81_TwoSource_DonorState_Complete_Localization_{DATE}.tsv',sep='\t',index=False,na_rep='NA')
 rec={'status':'RAW_COUNTS_DONOR_STATE_LOCALIZATION_COMPLETE','date':DATE,'sources':[ra,rb],'raw_donor_state_rows':len(a)+len(b),
  'candidate81_unchanged':True,'primary43_unchanged':True,'normal_controls_present':False,'new_P_values':0,'new_programme_scores':0,
  'author_and_broad_UMI_conservation':True,'cell_state_labels_retrained':False}
 (out/f'Count_Execution_Receipt_{DATE}.json').write_text(json.dumps(rec,indent=2),encoding='utf-8')
 print(json.dumps(rec,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--root',default=r'H:\SCI2\MRQI');a=ap.parse_args();run(Path(a.root))
