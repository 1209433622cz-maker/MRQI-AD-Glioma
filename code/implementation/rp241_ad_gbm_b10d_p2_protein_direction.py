"""S6-gated, SHA-locked all81/43 protein contrast execution; default preflight."""
from pathlib import Path
import csv,json,argparse,hashlib,datetime,numpy as np
from scipy.stats import t as student
R=Path(r'H:\SCI2\MRQI');P=R/'3_results/rp241_ad_gbm_b10d_p1_preparation'
def read(p):return list(csv.DictReader(p.open(encoding='utf-8'),delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,rows):
 with p.open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def fit(y,group,age,sex,batch,pmi=None):
 y=np.asarray(y,float);mask=np.isfinite(y);y=y[mask];g=np.asarray(group,float)[mask];a=np.asarray(age,float)[mask];s=np.asarray(sex,float)[mask];b=np.asarray(batch)[mask]
 result={'n':len(y),'cases':int(sum(g==1)),'controls':int(sum(g==0)),'rank':0,'df':0,'OLS_beta':'','OLS_SE':'','OLS_CI_low':'','OLS_CI_high':'','OLS_P':'','HC3_SE':'','HC3_CI_low':'','HC3_CI_high':'','HC3_P':'','max_leverage':'','status':'NOT_ESTIMABLE'}
 if min(result['cases'],result['controls'])<5:result['status']='NOT_ESTIMABLE_GROUP_N_LT5';return result
 cols=[np.ones(len(y)),g,a-a.mean(),s];levels=sorted(set(b));cols.extend((b==level).astype(float) for level in levels[1:])
 if pmi is not None:
  v=np.asarray(pmi,float)[mask];cols.append(v-v.mean())
 X=np.column_stack(cols);rank=int(np.linalg.matrix_rank(X));df=len(y)-rank;result.update(rank=rank,df=df)
 if rank<X.shape[1] or df<5:result['status']='NOT_ESTIMABLE_RANK_OR_DF';return result
 inv=np.linalg.inv(X.T@X);beta=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@beta;se=np.sqrt(max(0,float((res@res)/df*inv[1,1])));q=student.ppf(.975,df);leverage=np.sum((X@inv)*X,axis=1)
 if se<=1e-15:result['status']='NOT_ESTIMABLE_ZERO_VARIANCE';return result
 result.update(OLS_beta=float(beta[1]),OLS_SE=float(se),OLS_CI_low=float(beta[1]-q*se),OLS_CI_high=float(beta[1]+q*se),OLS_P=float(2*student.sf(abs(beta[1]/se),df)),max_leverage=float(leverage.max()),status='ESTIMABLE')
 if leverage.max()<1-1e-8:
  middle=X.T@((res/(1-leverage))[:,None]**2*X);cov=inv@middle@inv;hs=np.sqrt(max(0,cov[1,1]));result.update(HC3_SE=float(hs),HC3_CI_low=float(beta[1]-q*hs),HC3_CI_high=float(beta[1]+q*hs),HC3_P=float(2*student.sf(abs(beta[1]/hs),df)) if hs>1e-15 else '')
 return result
def bh(rows,pkey,qkey,subset):
 ii=[i for i,r in enumerate(rows) if subset(r) and isinstance(r[pkey],float)];ii.sort(key=lambda i:rows[i][pkey]);n=len(ii);vals=np.minimum.accumulate(np.asarray([rows[i][pkey]*n/(j+1) for j,i in enumerate(ii)])[::-1])[::-1] if n else []
 for i,v in zip(ii,vals):rows[i][qkey]=float(min(v,1))
 return n
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--execute',action='store_true');args=parser.parse_args()
 mf=R/'0_admin/manifest/RP241_AD_GBM_B10D_P1_Analysis_Input_Freeze_2026-10-15.tsv';assert mf.exists(),'P1 input SHA freeze required'
 for x in read(mf):assert sha(R/x['relative_path'])==x['sha256'],'INPUT_SHA_CHANGED:'+x['relative_path']
 ap=R/'0_admin/audit/RP241_AD_GBM_B10D_P1_Protein_S6_User_Reported_Review_2026-10-15.json';accepted=[]
 if ap.exists():
  review=json.loads(ap.read_text(encoding='utf-8'));assert review['status']=='USER_REPORTED_AUTHOR_REVIEW_ACCEPTED' and review['reviewer_name'] and review['review_date'] and review['user_attestation'];assert review['origin']=='user_reply_in_this_chat';assert review['review_memo_sha256']==sha(R/'0_admin/protocol/RP241/RP241_AD_GBM_B10D_P1_Protein_S6_Source_Review_Memo_2026-10-15.md');accepted=review['accepted_sources'];assert set(accepted)<={'AD_PXD020296','GBM_CPTAC_PDC000204'}
 receipt={'date':'2026-10-15','input_manifest_sha256':sha(mf),'execute_requested':args.execute,'accepted_sources':accepted,'status':'HOLD_NEW_PROTEIN_SOURCE_S6_REQUIRED' if not accepted else 'READY_FOR_ACCEPTED_SOURCE_EXECUTION','new_gene_effects_computed':False}
 rp=R/'0_admin/validation/RP241_AD_GBM_B10D_P2_Execution_Preflight_2026-10-15.json';rp.write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 if not args.execute or not accepted:
  print(json.dumps(receipt,indent=2))
  if args.execute and not accepted:raise SystemExit(2)
  return
 O=R/'3_results/rp241_ad_gbm_b10d_p2_protein_direction';assert not O.exists(),'Preserve prior scientific execution; no overwrite';O.mkdir()
 identity=read(R/'3_results/rp241_ad_gbm_b10_source_model_freeze/All81_Immutable_Identity_Direction_Input_2026-10-14.tsv');meta={x['HGNC_ID']:x for x in identity};assert len(identity)==81
 ad=read(P/'AD_All81_Prepared_Protein_Log2_GIS_Matrix_2026-10-15.tsv');ch=read(R/'3_results/rp241_ad_gbm_b10d_p0_source_audit/AD_PXD020296_All33_Channel_Donor_Join_2026-10-15.tsv');clin=read(P/'CPTAC_All110_Diagnostic_Specimen_Adjudication_NoPrune_2026-10-15.tsv');coverage={x['HGNC_ID']:x for x in read(P/'CPTAC_All81_Identity_Missingness_NoPrune_2026-10-15.tsv')};cm=read(P/'CPTAC_All81_BothQuantities_All110_Prepared_2026-10-15.tsv')
 tasks=[]
 if 'AD_PXD020296' in accepted:
  cc=[x for x in ch if x['Diagnosis'] in ['AD','control']]
  for name,pmi in [('AD_PRIMARY',False),('AD_PMI',True)]:tasks.append((name,ad,cc,'AD',pmi))
 if 'GBM_CPTAC_PDC000204' in accepted:
  for name,key,quant in [('GBM_PRIMARY','primary_shared_plex_set','LogRatio'),('GBM_UNSHARED','primary_shared_plex_set','Unshared'),('GBM_SEQUENCE12','main_limited_sequence_negative_set','LogRatio'),('GBM_CLINICAL52','expanded_clinical_negative_set','LogRatio'),('GBM_HISTORICAL98','historical_nonwithdrawn_tumor_set','LogRatio'),('GBM_CLINICAL52_NO_NONHOTSPOT','expanded_clinical_negative_set','LogRatio')]:
   cc=[x for x in clin if x[key]=='True' or (key!='primary_shared_plex_set' and x['healthy_GTEx_control']=='True')]
   if name.endswith('NO_NONHOTSPOT'):cc=[x for x in cc if x['case_id']!='C3L-03390']
   tasks.append((name,[x for x in cm if x['quantity']==quant],cc,'GBM',False))
 for name,matrix,cc,disease,pmi in tasks:
  rows=[]
  for x in matrix:
   md=meta[x['HGNC_ID']];eligible=x.get('status')=='MEASURABLE_MULTI_PEPTIDE' if disease=='AD' else coverage[x['HGNC_ID']]['symbol_and_Entrez_exact']=='True'
   if not eligible:
    rr={k:'' for k in ['n','cases','controls','rank','df','OLS_beta','OLS_SE','OLS_CI_low','OLS_CI_high','OLS_P','HC3_SE','HC3_CI_low','HC3_CI_high','HC3_P','max_leverage']};rr['status']='NOT_MEASURABLE_IDENTITY_OR_PEPTIDE'
   else:
    labels=[f'B{c["Batch"]}_{c["Channel"]}' if disease=='AD' else c['aliquot_id'] for c in cc];y=[float(x[label]) if x[label]!='' else np.nan for label in labels];groups=[int(c['Diagnosis']=='AD') if disease=='AD' else int(c['sample_type']=='tumor') for c in cc];ages=[float(c['Age'] if disease=='AD' else c['age']) for c in cc];sex=[int(c['Sex']=='female') if disease=='AD' else int(c['sex']=='Female') for c in cc];batch=[c['Batch'] if disease=='AD' else c['plex'] for c in cc]
    rr=fit(y,groups,ages,sex,batch,[float(c['PMI']) for c in cc] if pmi else None)
   expected=1 if (md['direction'].startswith('R1') and disease=='AD') or (md['direction'].startswith('R2') and disease=='GBM') else -1
   rows.append({**md,'analysis':name,'expected_beta_sign':expected,**rr,'OLS_Q_primary43':'','OLS_Q_all81_exploratory':'','HC3_Q_primary43':'','RNA_expected_direction_matches':(rr['OLS_beta']*expected>0) if isinstance(rr['OLS_beta'],float) else ''})
  bh(rows,'OLS_P','OLS_Q_primary43',lambda x:x['primary43']=='True');bh(rows,'OLS_P','OLS_Q_all81_exploratory',lambda x:True);bh(rows,'HC3_P','HC3_Q_primary43',lambda x:x['primary43']=='True');write(O/(name+'_All81_NoPrune.tsv'),rows)
 if 'AD_PXD020296' in accepted and 'GBM_CPTAC_PDC000204' in accepted:
  aa={x['HGNC_ID']:x for x in read(O/'AD_PRIMARY_All81_NoPrune.tsv')};bb={x['HGNC_ID']:x for x in read(O/'GBM_PRIMARY_All81_NoPrune.tsv')};joint=[]
  for x in identity:
   a=aa[x['HGNC_ID']];b=bb[x['HGNC_ID']];measured=all(v['status']=='ESTIMABLE' for v in [a,b]);direction=measured and all(v['RNA_expected_direction_matches']=='True' for v in [a,b]);supported=bool(direction and x['primary43']=='True' and all(v['OLS_Q_primary43']!='' and float(v['OLS_Q_primary43'])<.05 for v in [a,b]));joint.append({**x,'AD_status':a['status'],'GBM_status':b['status'],'AD_beta':a['OLS_beta'],'GBM_beta':b['OLS_beta'],'AD_main_q':a['OLS_Q_primary43'],'GBM_main_q':b['OLS_Q_primary43'],'both_point_directions_expected':direction,'both_main_BH_limited_public_P_support':supported,'CPTAC_collection_confounding_removed':False,'G_or_F_causal_support_created':False})
  write(O/'All81_CrossDisease_Protein_NoPrune_Ledger.tsv',joint)
 receipt.update(status='COMPLETE_ACCEPTED_SOURCE_PROTEIN_EXECUTION__NO_G_OR_F_UPGRADE',new_gene_effects_computed=True,analysis_sets=[x[0] for x in tasks],recorded_at_local=datetime.datetime.now().astimezone().isoformat());(O/'Execution_Receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
