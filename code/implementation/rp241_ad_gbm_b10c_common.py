"""B10C fixed input hashes and numerical primitives; no output-driven selection."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import t
DATE='2026-10-14';ROOT=Path(r'H:\SCI2\MRQI')
COMP=['malignant','myeloid','oligodendrocyte']
MAP={'Glioma':'malignant','Myeloid':'myeloid','Oligo':'oligodendrocyte'}
MASKS=['S_CNV_complete','S_CNV_concordant']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(16<<20),b''):h.update(b)
 return h.hexdigest()
def check_freeze(root=ROOT):
 j=json.loads((root/f'0_admin/protocol/RP241/RP241_AD_GBM_B10C_PreExecution_Input_Code_Freeze_{DATE}.json').read_text(encoding='utf-8'))
 for x in j['inputs']:
  p=root/x['relative_path'];assert p.stat().st_size==x['bytes'] and sha(p)==x['sha256'],x['relative_path']
 return j
def mask_index(md,mask):
 assert md.CopyKatPrediction.isin(['aneuploid','diploid','NA']).all()
 scope=md.Assignment.isin(MAP)
 if mask=='S_CNV_complete':return (scope & (md.CopyKatPrediction!='NA')).to_numpy()
 assert mask=='S_CNV_concordant'
 return (scope & (((md.Assignment=='Glioma')&(md.CopyKatPrediction=='aneuploid')) |
                 (md.Assignment.isin(['Myeloid','Oligo'])&(md.CopyKatPrediction=='diploid')))).to_numpy()
def fixed_bh(p):
 p=np.asarray(p,float);p=np.where(np.isfinite(p),p,1);n=len(p);o=np.argsort(p)
 q=np.minimum(1,np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]);v=np.empty(n);v[o]=q;return v
def ilr(counts):
 counts=np.asarray(counts,float);assert np.all(counts>=0)
 prop=(counts+.5)/(counts+.5).sum(axis=1,keepdims=True);lp=np.log(prop);K=lp.shape[1]
 return np.column_stack([np.sqrt(j/(j+1))*(lp[:,:j].mean(axis=1)-lp[:,j]) for j in range(1,K)])
def design(md,balances,rin=False,composition=False):
 xs=[np.ones(len(md))];names=['Intercept']
 for b in sorted(md.Batch.astype(str).unique())[1:]:xs.append((md.Batch.astype(str)==b).to_numpy(float));names.append('Batch_'+b)
 xs += [md.Age.to_numpy(float), (md.Sex=='M').to_numpy(float), (md.Diagnosis=='AD').to_numpy(float)]
 names += ['Age_centered_per10','Sex_M','Diagnosis_AD']
 if rin:xs.append(md.RIN.to_numpy(float));names.append('RIN_centered')
 if composition:
  for j in range(balances.shape[1]):xs.append(balances[:,j]);names.append('ILR_'+str(j+1))
 return np.column_stack(xs),names
def hc3(X,y,j):
 n,p=X.shape;rank=int(np.linalg.matrix_rank(X));df=n-rank
 base={'n':n,'columns':p,'rank':rank,'df':df,'condition_number':float(np.linalg.cond(X))}
 if rank!=p or df<3:return {**base,'status':'NOT_TESTABLE_RANK_OR_DF'}
 inv=np.linalg.inv(X.T@X);b=inv@X.T@y;res=y-X@b;h=np.einsum('ij,ij->i',X@inv,X)
 base.update(max_leverage=float(h.max()),min_one_minus_leverage=float((1-h).min()))
 if np.any(1-h<=1e-10):return {**base,'status':'NOT_TESTABLE_UNIT_LEVERAGE_HC3'}
 V=inv@(X.T@(((res/(1-h))**2)[:,None]*X))@inv;se=float(np.sqrt(max(0,V[j,j])));beta=float(b[j])
 if not np.isfinite(se) or se<=0:return {**base,'status':'NOT_TESTABLE_HC3_VARIANCE'}
 tc=float(t.ppf(.975,df));pv=float(2*t.sf(abs(beta/se),df))
 return {**base,'status':'EXPLORATORY_FIXED_SCORE_CONDITIONAL_ASSOCIATION','beta':beta,'SE':se,
         'CI_lower':beta-tc*se,'CI_upper':beta+tc*se,'P_two_sided':pv}
