"""GSE203206 metadata/scale audit and gated, bounded BA17 expression validation.

Preflight never computes disease coefficients. Execution requires a real S6 record
whose packet/input/code hashes agree with the pre-effect freeze.
"""
from pathlib import Path
import argparse, hashlib, json, zipfile, xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
from scipy import stats

ROOT=Path('H:/SCI2/MRQI'); DATE='2026-10-16'
OUT=ROOT/f'3_results/rp241_eda4_gse203206_{DATE}'
QA=ROOT/f'0_admin/validation/rp241_eda4_{DATE}'
AUDIT=ROOT/f'0_admin/source_audit/rp241_eda4_{DATE}'
PROTOCOL=ROOT/f'0_admin/protocol/RP241/RP241_EDA4_GSE203206_BA17_Model_Freeze_{DATE}.md'
MEMO=ROOT/f'0_admin/protocol/RP241/RP241_EDA4_GSE203206_S6_Review_Memo_{DATE}.md'
INPUTS={
 'matrix':ROOT/'1_data/transcriptomics/alzheimer/GSE203206/processed/Caldwell2022_AdditionalFile2_RNAcounts.tsv',
 'geo_matrix':ROOT/'1_data/transcriptomics/alzheimer/GSE203206/processed/GSE203206_Subramaniam.ADRC_brain.counts.tsv.gz',
 'metadata':ROOT/'1_data/transcriptomics/alzheimer/GSE203206/metadata/Caldwell2022_AdditionalFile2_Brain_Metadata.tsv',
 'closure':ROOT/'3_results/rp241_ad_b5_source_admission/GSE203206_Paper_GEO_Donor_Control_QC_Closure_2026-10-14.tsv',
 'identities':ROOT/'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_2026-10-14.tsv',
 'Caldwell_article':ROOT/'8.Reference/RP241_AD_B5_Source_Admission_2026-10-14/Caldwell2022_GSE203206_fulltext.xml',
 'Srinivasan_article':ROOT/'8.Reference/RP241_AD_Source_Methods_2026-10-13/Srinivasan2020_GSE125583_fulltext.xml',
 'author_zip':ROOT/'8.Reference/RP241_AD_B5_Source_Admission_2026-10-14/Caldwell2022_GSE203206_AdditionalFile2.zip',
 'donor_discovery':ROOT/'3_results/rp241_ad_source_intake/GSE125583_Sample_Metadata_2026-10-13.tsv',
 'Penn_role':ROOT/'3_results/rp241_ad_gbm_b8_preflight/Primary22_Metadata_2026-10-14.tsv',
}
MODELS={
 'M0_primary':('all',['age','RIN'],'provided_log2'),
 'M1_mapping':('all',['age','RIN','mapping_rate'],'provided_log2'),
 'M2_PMI':('all',['age','RIN','PMI'],'provided_log2'),
 'M3_age_common_support':('common_age',['age','RIN'],'provided_log2'),
 'M4_EOAD':('EOAD',['age','RIN'],'provided_log2'),
 'M5_LOAD':('LOAD',['age','RIN'],'provided_log2'),
 'M6_relative_scale':('all',['age','RIN'],'relative_log2'),
}

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def record(p):return dict(path=Path(p).relative_to(ROOT).as_posix(),bytes=Path(p).stat().st_size,sha256=sha(p))
def write(p,v):
 Path(p).parent.mkdir(parents=True,exist_ok=True)
 Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n',encoding='utf-8')
def tab(frame,name,folder=OUT):
 folder.mkdir(parents=True,exist_ok=True)
 (frame if isinstance(frame,pd.DataFrame) else pd.DataFrame(frame)).to_csv(folder/f'{name}_{DATE}.tsv',sep='\t',index=False,na_rep='NA',float_format='%.12g')

def metadata():
 m=pd.read_csv(INPUTS['metadata'],sep='\t',dtype=str).rename(columns={'SampleID':'sample_id','PATHID':'donor_id','SEX':'sex','AGE':'age','MappingRate':'mapping_rate','PMHOURS':'PMI'})
 c=pd.read_csv(INPUTS['closure'],sep='\t',dtype=str).set_index('sample_id')
 m['GEO_group']=m.sample_id.map(c.GEO_group)
 for col in ['age','RIN','mapping_rate','PMI']:m[col]=pd.to_numeric(m[col],errors='coerce')
 m['diagnosis_AD']=(m.GEO_group!='Non-demented Control').astype(float)
 assert len(m)==47 and m.donor_id.is_unique and m.sample_id.is_unique
 assert m.GEO_group.value_counts().to_dict()=={'Late-onset sporadic AD':20,'Early-onset sporadic AD':19,'Non-demented Control':8}
 assert set(m.sex)=={'M','F'}
 return m

def design(m,subset,covs):
 s=m.copy()
 if subset in ['EOAD','LOAD']:
  lab={'EOAD':'Early-onset sporadic AD','LOAD':'Late-onset sporadic AD'}[subset]
  s=s[s.GEO_group.isin([lab,'Non-demented Control'])]
 if subset=='common_age':
  a=m.loc[m.diagnosis_AD==1,'age'].dropna();c=m.loc[m.diagnosis_AD==0,'age'].dropna()
  lo=max(a.min(),c.min());hi=min(a.max(),c.max());s=s[s.age.between(lo,hi)]
 s=s.dropna(subset=covs+['sex']).copy().reset_index(drop=True)
 x=pd.DataFrame({'Intercept':np.ones(len(s)),'sex_M':(s.sex=='M').astype(float)})
 for col in covs:
  sd=s[col].std(ddof=0);assert sd>0
  x[col+'_z']=(s[col]-s[col].mean())/sd
 x['diagnosis_AD']=s.diagnosis_AD
 assert np.linalg.matrix_rank(x)==x.shape[1] and len(s)>x.shape[1]+2
 return s,x

def fit(x,y):
 """Unweighted OLS; HC3 covariance with residual-df t reference (approximate)."""
 x=np.asarray(x,float); y=np.asarray(y,float)
 if y.ndim==1:y=y[:,None]
 inv=np.linalg.inv(x.T@x);a=inv@x.T;b=a@y;e=y-x@b
 h=np.einsum('ij,ji->i',x,a); assert max(h)<1-1e-9
 var3=(a[-1,:,None]**2*(e/(1-h[:,None]))**2).sum(axis=0)
 se=np.sqrt(var3);df=x.shape[0]-x.shape[1];t=b[-1]/se
 p=2*stats.t.sf(abs(t),df);se0=np.sqrt((e**2).sum(axis=0)/df*inv[-1,-1])
 crit=stats.t.ppf(.975,df)
 return dict(beta=b[-1],SE_HC3=se,P_HC3=p,CI95_low=b[-1]-crit*se,CI95_high=b[-1]+crit*se,
  SE_OLS=se0,P_OLS=2*stats.t.sf(abs(b[-1]/se0),df),residual_df=df)

def bh_fixed(p):
 p=np.asarray(p,float);work=np.where(np.isfinite(p),p,1.);idx=np.argsort(work,kind='stable')
 q=np.empty(len(p));q[idx]=np.minimum.accumulate((work[idx]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
 q[~np.isfinite(p)]=np.nan;return q

def preflight():
 if (OUT/f'Bounded_Validation_Summary_{DATE}.json').exists():
  raise SystemExit('REFUSE_POST_EFFECT_REFREEZE__PRESERVE_ORIGINAL_PROTOCOL')
 for p in [OUT,QA,AUDIT]:p.mkdir(parents=True,exist_ok=True)
 assert sha(INPUTS['identities'])=='93dff5698e60a1b03d2a232ab2963ffe5cd8fd60699bda479d11ed7f1eec24e0'
 m=metadata();tab(m,'All47_Source_Clinical_Review')
 designs=[];members=[];vifs=[]
 for name,(sub,covs,scale) in MODELS.items():
  s,x=design(m,sub,covs);a=np.linalg.inv(x.T@x)@x.T;h=np.einsum('ij,ji->i',x,a)
  f=x.copy();f.insert(0,'sample_id',s.sample_id);f.insert(1,'donor_id',s.donor_id);f['leverage']=h
  tab(f,name+'_Frozen_Design')
  members.extend(dict(model=name,sample_id=r.sample_id,donor_id=r.donor_id,included=r.sample_id in set(s.sample_id),reason='included' if r.sample_id in set(s.sample_id) else 'subset_or_missing_required_covariate') for r in m.itertuples())
  for j,col in enumerate(x.columns[1:],1):
   target=x.iloc[:,j].to_numpy();other=np.delete(x.to_numpy(),j,axis=1);res=target-other@np.linalg.lstsq(other,target,rcond=None)[0]
   vif=np.sum((target-target.mean())**2)/np.sum(res**2)
   vifs.append(dict(model=name,term=col,VIF=vif))
  designs.append(dict(model=name,subset=sub,scale=scale,covariates='sex+'+'+'.join(covs),n=len(s),AD=int(s.diagnosis_AD.sum()),Control=int((s.diagnosis_AD==0).sum()),rank=np.linalg.matrix_rank(x),parameters=x.shape[1],residual_df=len(s)-x.shape[1],condition_number=np.linalg.cond(x),max_leverage=max(h),max_VIF=max(t['VIF'] for t in vifs if t['model']==name),disease_effect_computed=False))
 tab(designs,'Model_Design_Closure_NoDiseaseEffects');tab(members,'All47_Model_Membership');tab(vifs,'Model_VIF_NoDiseaseEffects')
 df=pd.read_csv(INPUTS['matrix'],sep='\t',index_col=0);g=pd.read_csv(INPUTS['geo_matrix'],sep='\t',index_col=0)
 assert df.shape==(24448,47) and df.index.is_unique and df.columns.is_unique
 assert set(df.columns)==set(m.sample_id)
 assert np.array_equal(df.to_numpy(),g.reindex(index=df.index,columns=df.columns).to_numpy())
 v=df.to_numpy(float);assert np.isfinite(v).all() and (v>=0).all()
 scale=dict(genes=24448,samples=47,numeric_cells=v.size,fractional_cells=int((v!=np.floor(v)).sum()),minimum=float(v.min()),maximum=float(v.max()),column_sum_min=float(v.sum(0).min()),column_sum_max=float(v.sum(0).max()),all_columns_sum_one_million=bool(np.allclose(v.sum(0),1e6)),provided_stage='UNRESOLVED_EXPORT_OBJECT',raw_count_model_permitted=False,allowed_scope='log2(provided_value+1) continuous association; relative column-composition sensitivity',fractional_counts_alone_disqualify_count_methods=False)
 write(AUDIT/f'Numeric_Scale_Fingerprint_{DATE}.json',scale)
 inventory=[]
 with zipfile.ZipFile(INPUTS['author_zip']) as z:
  assert z.testzip() is None
  for i in z.infolist():inventory.append(dict(member=i.filename,bytes=i.file_size,CRC=i.CRC))
 tab(inventory,'Author_AdditionalFile2_Member_Inventory',AUDIT)
 excerpts=[]
 for key,terms in [('Caldwell_article',['TMM','voom','tximport','cluster 3','degradation','Shiley','total of 50','quality filtering']),('Srinivasan_article',['Banner Sun Health','Whole tissue studies cohorts','289 total'])]:
  root=ET.parse(INPUTS[key]).getroot()
  for i,e in enumerate(root.iter('p')):
   t=' '.join(''.join(e.itertext()).split())
   hits=[q for q in terms if q.lower() in t.lower()]
   if hits:excerpts.append(dict(source=key,paragraph_index=i,matched_terms=';'.join(hits),text=t))
 tab(excerpts,'Primary_Paper_Provenance_Excerpts',AUDIT)
 fam=[dict(source='GSE203206',bank='UCSD Shiley-Marcos ADRC',region='BA17/occipital',family='UCSD_Shiley_Marcos_ADRC',role='new bounded validation; S6 pending',evidence='Caldwell 2022 Methods and author table'),
      dict(source='GSE125583 / GSE95587',bank='Banner Sun Health Brain and Body Donation Program',region='fusiform gyrus',family='Banner_Sun_Health_GSE125583_GSE95587',role='discovery family; 89 shared donors across the two accessions; not independent of each other',evidence='Srinivasan 2020 Experimental model/Whole tissue RNA-Seq'),
      dict(source='GSE159699 / related Nativio series',bank='Penn CNDR',region='cortex',family='Penn_CNDR_Nativio',role='existing validation family',evidence='B5/B8 source admission receipts')]
 tab(fam,'Source_Family_Independence_Review',AUDIT)
 # Reviewable input/code freeze is made only after protocol, memo and validator exist.
 deps=[Path(__file__),ROOT/'2_code/rp241/rp241_eda4_validate.py',PROTOCOL,MEMO]
 freeze=dict(date=DATE,status='PRE_EFFECT_SOURCE_MODEL_FREEZE__S6_REQUIRED',inputs={k:record(p) for k,p in INPUTS.items()},governance=[record(p) for p in deps],models=designs,source_family_independence='documented different banks; not a global identity linkage',disease_effect_computed=False)
 write(QA/f'PreEffect_Input_Code_Model_Freeze_{DATE}.json',freeze)
 print(json.dumps(dict(status=freeze['status'],models=designs,scale=scale),indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)))

def execute(review_path):
 if (OUT/f'Bounded_Validation_Summary_{DATE}.json').exists():
  raise SystemExit('EXISTING_EFFECT_OUTPUTS_PRESENT__VALIDATE_OR_USE_NEW_VERSIONED_RUN')
 f=json.loads((QA/f'PreEffect_Input_Code_Model_Freeze_{DATE}.json').read_text(encoding='utf-8'))
 for r in list(f['inputs'].values())+f['governance']:
  assert record(ROOT/r['path'])==r, f"Frozen input/code drift: {r['path']}"
 review=json.loads(Path(review_path).read_text(encoding='utf-8'))
 assert review.get('source')=='GSE203206' and review.get('accepted') is True
 assert review.get('reviewer_name') and review.get('review_date') and review.get('human_statement_verbatim')
 assert review.get('review_memo_sha256')==sha(MEMO) and review.get('model_protocol_sha256')==sha(PROTOCOL)
 assert review.get('scope')=='BA17_BOUNDED_PROVIDED_SCALE_VALIDATION_ALL81_PRIMARY43'
 b=pd.read_csv(INPUTS['identities'],sep='\t');b['primary43']=b.primary43.astype(str).str.lower().isin(['true','1'])
 assert len(b)==81 and b.primary43.sum()==43 and b.HGNC_ID.is_unique
 mat=pd.read_csv(INPUTS['matrix'],sep='\t',index_col=0)
 rel=mat.div(mat.sum(axis=0),axis=1)*1e6
 m=metadata();allrows=[];loos=[]
 assert set(mat.columns)==set(m.sample_id)
 for name,(sub,covs,scale) in MODELS.items():
  s,x=design(m,sub,covs);trans=np.log2((mat if scale=='provided_log2' else rel)+1)
  keys=b.AD_external_Ensembl_stable_ID;present=keys.isin(mat.index)
  y=trans.loc[keys[present],s.sample_id].to_numpy().T;result=fit(x,y)
  rows=b.copy();rows['model']=name;rows['scale']=scale;rows['measured']=present
  rows['n']=len(s);rows['AD']=int(s.diagnosis_AD.sum());rows['Control']=len(s)-int(s.diagnosis_AD.sum())
  for col,val in result.items():
   rows[col]=np.nan;rows.loc[present,col]=val
  rows['expected_AD_sign']=np.where(rows.direction.str.startswith('R1'),1,-1)
  rows['direction_concordant']=np.where(present,np.sign(rows.beta)==rows.expected_AD_sign,None)
  rows['q_BH_81']=bh_fixed(rows.P_HC3);rows['q_BH_43']=np.nan
  rows.loc[rows.primary43,'q_BH_43']=bh_fixed(rows.loc[rows.primary43,'P_HC3'])
  rows['primary_support']=rows.primary43 & rows.measured & rows.q_BH_43.lt(.05) & rows.direction_concordant.eq(True)
  allrows.append(rows)
  if name=='M0_primary':
   for i in range(len(s)):
    xx=np.delete(x.to_numpy(),i,axis=0);yy=np.delete(y,i,axis=0);rr=fit(xx,yy)
    zz=b.loc[present,['HGNC_ID','HGNC_symbol','primary43']].copy();zz['omitted_sample']=s.iloc[i].sample_id;zz['omitted_donor']=s.iloc[i].donor_id
    for col in ['beta','SE_HC3','P_HC3']:zz[col]=rr[col]
    loos.append(zz)
 results=pd.concat(allrows,ignore_index=True);tab(results,'All81_All7Models_Expression_Validation')
 loo=pd.concat(loos,ignore_index=True);tab(loo,'M0_AllMeasurable79_LeaveOneDonorOut')
 primary=results[(results.model=='M0_primary') & results.primary43]
 summary=dict(date=DATE,status='GSE203206_BA17_BOUNDED_EXPLORATORY_VALIDATION_EXECUTED',source='GSE203206',S6=review,review_record=record(review_path),model_freeze_sha256=sha(QA/f'PreEffect_Input_Code_Model_Freeze_{DATE}.json'),primary_denominator=43,primary_measured=int(primary.measured.sum()),primary_concordant=int(primary.direction_concordant.eq(True).sum()),primary_q43_lt_005=int(primary.q_BH_43.lt(.05).sum()),primary_directional_support=int(primary.primary_support.sum()),raw_count_reproduction=False,causal_tier_promotion=False)
 write(OUT/f'Bounded_Validation_Summary_{DATE}.json',summary);print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--preflight',action='store_true');ap.add_argument('--execute',action='store_true');ap.add_argument('--review-record',type=Path);a=ap.parse_args()
 assert a.preflight != a.execute,'Choose one mode'
 if a.preflight:preflight()
 else:
  if a.review_record is None or not a.review_record.is_file():raise SystemExit('HOLD_GSE203206_S6_NAMED_SOURCE_REVIEW_REQUIRED__NO_DISEASE_EFFECTS_COMPUTED')
  execute(a.review_record)
