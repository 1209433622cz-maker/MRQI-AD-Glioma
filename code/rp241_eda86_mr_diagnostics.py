from rp241_eda86_common import *
import math, subprocess, gzip, struct
from itertools import product
cfg=load(P/f'RP241_EDA86_CoreGap_Diagnostic_And_Public_Disclosure_Execution_Freeze_{D}.json')
fwd=[x for x in rows(R/cfg['forward_input']) if x['analysis_status']=='INCLUDE']
rev=[x for x in rows(R/cfg['reverse_input']) if x['analysis_status']=='INCLUDE']
fs={x['rsid'] for x in fwd}; rs={x['rsid'] for x in rev}; union=fs|rs
assert len(fs)==10 and len(rs)==57
panel=R/'1_data/reference/1000G_phase3/ieu_1kg_v3/EUR'
bim={};total=0
with panel.with_suffix('.bim').open() as z:
 for i,line in enumerate(z):
  a=line.split();total+=1
  if a[1] in union:
   assert a[1] not in bim,('Duplicate panel rsID',a[1])
   bim[a[1]]=dict(chr=a[0],pos=int(a[3]),a1=a[4],a2=a[5],index=i)
comp=str.maketrans('ACGT','TGCA')
bindings=[]
for direction,records,ea,oa in [('GBM_to_strictAD',fwd,'effect_allele_exposure','other_allele_exposure'),('ADproxy_to_GBM',rev,'effect_allele_x','other_allele_x')]:
 for x in records:
  p=bim.get(x['rsid']); alleles={x[ea],x[oa]}
  status='PANEL_ABSENT' if not p else ('EXACT_ALLELE_SET' if alleles=={p['a1'],p['a2']} else 'COMPLEMENT_ALLELE_SET' if {s.translate(comp) for s in alleles}=={p['a1'],p['a2']} else 'ALLELE_MISMATCH')
  bindings.append(dict(direction=direction,rsid=x['rsid'],status=status,panel_chr=p['chr'] if p else '',panel_GRCh37_pos=p['pos'] if p else '',source_build='GRCh38 outcome/AD source; published GBM identified by rsID and alleles; no coordinate mixing'))
tab(O/'All67_Direction_IV_To_EUR_Panel_Binding.tsv',bindings)
eligible={x['rsid'] for x in bindings if x['status'] in ('EXACT_ALLELE_SET','COMPLEMENT_ALLELE_SET')}
put(O/'Allele_Qualified_Panel_Extract.txt','\n'.join(sorted(eligible))+'\n')
prefix=O/'EUR_AllUnion_AllPairs'
cmd=[str(R/'2_code/tools/plink_1.9_20260913/plink.exe'),'--bfile',str(panel),'--extract',str(O/'Allele_Qualified_Panel_Extract.txt'),'--r2','inter-chr','--ld-window-r2','0','--threads','1','--memory','1024','--out',str(prefix)]
z=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
put(V/'PLINK_Command_And_Output.txt',json.dumps(cmd,ensure_ascii=False)+'\n'+z.stdout+'\n'+z.stderr)
assert z.returncode==0,z.stderr
ld={}
with prefix.with_suffix('.ld').open() as z:
 h=z.readline().split()
 for line in z:
  a=dict(zip(h,line.split()));ld[frozenset((a['SNP_A'],a['SNP_B']))]=float(a['R2'])
assert len(ld)==len(eligible)*(len(eligible)-1)//2,(len(ld),len(eligible))
cross=[]
for a,b in product(sorted(fs),sorted(rs)):
 p,q=bim.get(a),bim.get(b)
 key=frozenset((a,b));r2=1.0 if a==b else ld.get(key)
 cross.append(dict(GBM_IV=a,AD_IV=b,exact_same=a==b,panel_pair_status='ESTIMATED' if r2 is not None else 'NOT_ESTIMABLE_PANEL_OR_ALLELE',r2=r2 if r2 is not None else '',strong_LD=r2 is not None and r2>=.8,same_panel_chr=p['chr']==q['chr'] if p and q else '',panel_distance_bp=abs(p['pos']-q['pos']) if p and q and p['chr']==q['chr'] else '',scope='EUR reference genotypic LD; no shared participant inference'))
tab(O/'All570_CrossDirection_Exact_And_LD_Overlap.tsv',cross)
exact=fs&rs
strong={x['GBM_IV'] for x in cross if x['strong_LD']}|{x['AD_IV'] for x in cross if x['strong_LD']}
overlap=dict(forward=10,reverse=57,union=len(union),exact_intersection=sorted(exact),jaccard=len(exact)/len(union),forward_exact_fraction=len(exact)/len(fs),reverse_exact_fraction=len(exact)/len(rs),cross_pairs=len(cross),panel_bound_variants=len(eligible),panel_absent_or_allele_unresolved=sorted(union-eligible),cross_pairs_estimable=sum(x['r2']!='' for x in cross),strong_LD_pairs=sum(x['strong_LD'] for x in cross),strong_LD_variant_union=sorted(strong),caution='Absence of exact/high-LD overlap does not prove exclusion restriction or absence of pleiotropy; unbound variants remain unresolved')
js(O/'CrossDirection_Overlap_Summary.json',overlap)

# Independent BED decoding validates all pairwise cross-direction r2 against PLINK.
fam=[x.split() for x in panel.with_suffix('.fam').read_text().splitlines() if x.strip()]
n=len(fam);width=(n+3)//4;dos={};code={0:2,1:None,2:1,3:0}
with panel.with_suffix('.bed').open('rb') as f:
 assert f.read(3)==bytes([0x6c,0x1b,0x01])
 for snp in eligible:
  f.seek(3+bim[snp]['index']*width);b=f.read(width)
  dos[snp]=[code[(b[i//4]>>(2*(i%4)))&3] for i in range(n)]
checks=[]
for x in cross:
 if x['r2']=='':continue
 a,b=x['GBM_IV'],x['AD_IV']; pairs=[(v,w) for v,w in zip(dos[a],dos[b]) if v is not None and w is not None]
 nv=len(pairs);sx=sum(v for v,w in pairs);sy=sum(w for v,w in pairs)
 cov=sum(v*w for v,w in pairs)-sx*sy/nv
 vx=sum(v*v for v,w in pairs)-sx*sx/nv;vy=sum(w*w for v,w in pairs)-sy*sy/nv
 r2=cov*cov/(vx*vy)
 # PLINK text output has six significant figures; do not call this historical MR tolerance.
 checks.append(dict(GBM_IV=a,AD_IV=b,complete_panel_n=nv,plink_r2=x['r2'],independent_bed_r2=r2,absolute_difference=abs(r2-x['r2']),pass_text_precision=abs(r2-x['r2'])<=5.1e-7))
assert all(x['pass_text_precision'] for x in checks)
tab(V/'All570_PLINK_Vs_Independent_BED_R2_Check.tsv',checks)

def fit(bx,by,sy,mre):
 w=[1/(s*s) for s in sy];den=sum(v*v*t for v,t in zip(bx,w));b=sum(v*u*t for v,u,t in zip(bx,by,w))/den
 q=sum(t*(u-b*v)**2 for v,u,t in zip(bx,by,w));phi=max(1,q/(len(bx)-1)) if mre and len(bx)>1 else 1
 se=math.sqrt(phi/den);p=math.erfc(abs(b/se)/math.sqrt(2))
 return dict(beta=b,se=se,ci_low=b-1.959963984540054*se,ci_high=b+1.959963984540054*se,p=p,Q=q,phi=phi)
assert abs(fit([1,2,3],[-.5,-1,-1.5],[1,1,1],True)['beta']+.5)<1e-12
sens=[]
for name,data,bx,by,sy,sx in [('GBM_to_strictAD',fwd,'beta_exposure','beta_outcome','se_outcome','se_exposure'),('ADproxy_to_GBM',rev,'beta_x','beta_y','se_y','se_x')]:
 for label,excluded in [('original_fixed_set',set()),('exact_overlap_excluded',exact),('strong_LD_both_members_excluded',strong)]:
  keep=[x for x in data if x['rsid'] not in excluded]
  assert len(keep)>=3
  for method in ['IVW_FE','IVW_MRE']:
   rr=fit([float(x[bx]) for x in keep],[float(x[by]) for x in keep],[float(x[sy]) for x in keep],method.endswith('MRE'))
   sens.append(dict(direction=name,sensitivity=label,method=method,k=len(keep),removed=';'.join(sorted({x['rsid'] for x in data}&excluded)),F_min=min((float(x[bx])/float(x[sx]))**2 for x in keep),F_max=max((float(x[bx])/float(x[sx]))**2 for x in keep),identical_to_original_set=len(keep)==len(data),**rr))
tab(O/'MR_Exact_And_StrongLD_Exclusion_Sensitivities.tsv',sens)
baseline=[]
oldf=rows(R/'3_results/phase0t_a2/Phase0T_A2_FinnGen_R13_G6_AD_STRICT_MR_Methods_2026-09-27.tsv')
oldr=rows(R/'3_results/phase0w_b1/Phase0W_B1_AD_to_C3GBM_MR_Methods_2026-09-27.tsv')
for direction,old in [('GBM_to_strictAD',oldf),('ADproxy_to_GBM',oldr)]:
 for method in ['IVW_FE','IVW_MRE']:
  a=next(x for x in sens if x['direction']==direction and x['sensitivity']=='original_fixed_set' and x['method']==method)
  b=next(x for x in old if x['method'].replace('-','_')==method)
  for key in ['beta','se','p']:
   err=abs(a[key]-float(b[key]));baseline.append(dict(direction=direction,method=method,field=key,difference=err,pass_abs_1e_12=err<=1e-12))
assert all(x['pass_abs_1e_12'] for x in baseline)
tab(V/'Original_Primary_MR_Reimplementation_Check.tsv',baseline)

# Read only the 57 frozen target rows from original raw GWAS; no exposure reselection.
raws=[(R/'1_data/gwas/alzheimer/GCST90027158/GCST90027158_buildGRCh38.tsv.gz','013134d22a2a0c2dd0cb683192d0c109da1aec83ba54d94a837a8b71bcbc4a1f','ADproxy'),(R/'1_data/gwas/glioma/finngen_R13_C3_GBM/finngen_R13_C3_GBM.gz','5b126534879804f842cad522ba473337ed962b2e5f6c3470b7786f674fbfd189','GBM_outcome')]
targetrows=[];raw_receipts=[]
for path,expected,label in raws:
 actual=sha(path);assert actual==expected,('Raw SHA mismatch',path)
 found={};scanned=0
 with gzip.open(path,'rt',encoding='utf-8') as z:
  h=z.readline().rstrip('\n\r').split('\t');ix={x:i for i,x in enumerate(h)}
  for line in z:
   scanned+=1;a=line.rstrip('\n\r').split('\t')
   if label=='ADproxy':hits=[a[ix['variant_id']]] if a[ix['variant_id']] in rs else []
   else:hits=[s for s in re.split(r'[,;|\s]+',a[ix['rsids']]) if s in rs]
   for snp in hits:
    assert snp not in found,('Duplicate exact raw',snp,label)
    record=dict(zip(h,a));found[snp]=record
   if scanned%5000000==0:print(label,scanned,'rows scanned; fixed hits',len(found),flush=True)
 for snp in sorted(rs):
  a=found.get(snp)
  targetrows.append(dict(source=label,rsid=snp,found=a is not None,n_cases=a.get('n_cases','') if a else '',n_controls=a.get('n_controls','') if a else '',effect_allele=a.get('effect_allele',a.get('alt','')) if a else '',other_allele=a.get('other_allele',a.get('ref','')) if a else '',EAF=a.get('effect_allele_frequency',a.get('af_alt','')) if a else '',raw_beta=a.get('beta','') if a else '',raw_se=a.get('standard_error',a.get('sebeta','')) if a else '',interpretation='Exact source row; proxy/clinical phenotype and liability assumptions still distinct'))
 raw_receipts.append(dict(**rec(path),scanned_rows=scanned,fixed_targets=57,exact_hits=len(found)))
tab(O/'All114_Reverse_IV_Exact_Raw_EAF_CaseCount_Receipt.tsv',targetrows)
js(V/'Two_Raw_GWAS_Full_SHA_FixedTarget_Scan_Receipt.json',raw_receipts)
metadata=[]
for name,key in [('GBM_forward_exposure','GBM'),('strict_AD_forward_outcome','AD'),('proxy_AD_reverse_exposure','ADproxy'),('GBM_reverse_outcome','GBMout')]:
 reasons=[]
 if key=='GBM':reasons=['Published instrument table has N=24352, no exact exposure EAF or source-matched case/control split; Melin overall 6191+18190=24381 cannot silently replace per-effect provenance']
 elif key=='AD':reasons=['Frozen FinnGen12545/483210 and aligned EAF available; standalone row does not close missing GBM exposure fields']
 elif key=='ADproxy':reasons=['Exact variant-level n_cases/n_controls/EAF recovered from own raw file; proxy-inclusive meta-analysis is not strict clinical AD prevalence and must not borrow accession summary counts without matching']
 else:reasons=['FinnGen467/372159 and exact outcome EAF available after fixed57 scan; reverse proxy liability interpretation remains limited']
 metadata.append(dict(component=name,status='FIELD_PARTIAL_LIABILITY_HOLD' if key=='GBM' else 'FIELD_AVAILABLE_WITH_PHENOTYPE_CAVEAT',reason='; '.join(reasons),formal_Steiger_p_produced=False,liability_grid_executed=False,scientific_no_go=False))
tab(O/'Steiger_Exact_Source_Field_And_Phenotype_Gate.tsv',metadata)
js(V/'MR_Diagnostic_Validation.json',dict(overlap=overlap,independent_LD_pairs_checked=len(checks),reference_individuals=n,MR_baseline_numeric_checks=len(baseline),liability_Steiger='TECHNICAL_HOLD_SOURCE_MATCHED_GBM_EAF_COUNTS_AND_PROXY_PREVALENCE_INTERPRETATION',new_transcriptomic_DE=False))
print(json.dumps({'overlap':overlap,'raw_receipts':raw_receipts,'MR_rows':len(sens),'Steiger':'FIELD_AND_PHENOTYPE_HOLD_NO_P_CREATED'},ensure_ascii=False,indent=2))
