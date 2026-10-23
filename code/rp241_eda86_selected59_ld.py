from rp241_eda86_common import *
import subprocess,math
from itertools import product
# Exposure-side completeness supplement; the previously frozen analysed570 stays intact.
fwd=[x for x in rows(R/'3_results/phase0t_a2/Phase0T_A2_FinnGen_R13_G6_AD_STRICT_Harmonised_2026-09-27.tsv') if x['analysis_status']=='INCLUDE']
ad=rows(R/'3_results/phase0t_a/GCST90027158_AD__strict_clumped_instruments.tsv')
assert len(ad)==59
js(P/f'RP241_EDA86_ExposureSide59_LD_Completeness_Addendum_{D}.json',dict(frozen_utc=utc(),reason='RP241 overlap concerns exposure-side selection as well as analysed rows. Add two previously frozen AD IVs lacking outcome; retain original570 and no rescue.',GBM_nonpal_exposure_IVs=10,AD_exposure_IVs=59,AD_analysed_IVs=57,threshold_r2=.8,primary_estimate_change=False,no_new_exposure_selection=True))
fs={x['rsid'] for x in fwd};rs={x['SNP'] for x in ad};union=fs|rs
panel=R/'1_data/reference/1000G_phase3/ieu_1kg_v3/EUR';bim={}
with panel.with_suffix('.bim').open() as z:
 for i,line in enumerate(z):
  a=line.split()
  if a[1] in union:
   assert a[1] not in bim
   bim[a[1]]={'chr':a[0],'pos':int(a[3]),'a1':a[4],'a2':a[5],'index':i}
comp=str.maketrans('ACGT','TGCA');binding=[]
for d,items,key,ea,oa in [('GBM',fwd,'rsid','effect_allele_exposure','other_allele_exposure'),('AD_selected',ad,'SNP','EFFECT_ALLELE','OTHER_ALLELE')]:
 for x in items:
  snp=x[key];p=bim.get(snp);aset={x[ea],x[oa]}
  ok=p is not None and (aset=={p['a1'],p['a2']} or {v.translate(comp) for v in aset}=={p['a1'],p['a2']})
  binding.append(dict(direction=d,rsid=snp,panel_allele_set_pass=ok,analysed_in_reverse=snp in {v['rsid'] for v in rows(R/'3_results/phase0w_b1/Phase0W_B1_AD_to_C3GBM_Harmonised_2026-09-27.tsv') if v['analysis_status']=='INCLUDE'},scope='Exposure side, not outcome reselection'))
assert all(x['panel_allele_set_pass'] for x in binding)
tab(O/'All69_ExposureSide_IV_Panel_Binding.tsv',binding)
put(O/'Frozen_ExposureSide_Union69.txt','\n'.join(sorted(union))+'\n')
prefix=O/'EUR_ExposureSide69_AllPairs'
cmd=[str(R/'2_code/tools/plink_1.9_20260913/plink.exe'),'--bfile',str(panel),'--extract',str(O/'Frozen_ExposureSide_Union69.txt'),'--r2','inter-chr','--ld-window-r2','0','--threads','1','--memory','1024','--out',str(prefix)]
z=subprocess.run(cmd,capture_output=True,text=True,timeout=120);put(V/'PLINK_ExposureSide69_Command_Output.txt',json.dumps(cmd)+'\n'+z.stdout+'\n'+z.stderr);assert z.returncode==0
ld={}
with prefix.with_suffix('.ld').open() as z:
 h=z.readline().split()
 for line in z:
  x=dict(zip(h,line.split()));ld[frozenset((x['SNP_A'],x['SNP_B']))]=float(x['R2'])
assert len(ld)==len(union)*(len(union)-1)//2
cross=[]
for a,b in product(sorted(fs),sorted(rs)):
 v=1.0 if a==b else ld[frozenset((a,b))]
 cross.append(dict(GBM_IV=a,AD_selected_IV=b,exact_same=a==b,r2=v,strong_LD_r2_ge_0_8=v>=.8,AD_outcome_missing=b in {'rs12151021','rs12972720'},panel_same_chr=bim[a]['chr']==bim[b]['chr'],panel_distance_bp=abs(bim[a]['pos']-bim[b]['pos']) if bim[a]['chr']==bim[b]['chr'] else ''))
tab(O/'All590_ExposureSide_CrossDirection_LD.tsv',cross)
fam=panel.with_suffix('.fam').read_text().splitlines();n=len(fam);w=(n+3)//4;dose={};decode={0:2,1:None,2:1,3:0}
with panel.with_suffix('.bed').open('rb') as z:
 assert z.read(3)==b'\x6c\x1b\x01'
 for snp in union:
  z.seek(3+bim[snp]['index']*w);b=z.read(w);dose[snp]=[decode[(b[i//4]>>(2*(i%4)))&3] for i in range(n)]
checks=[]
for x in cross:
 p=[(a,b) for a,b in zip(dose[x['GBM_IV']],dose[x['AD_selected_IV']]) if a is not None and b is not None];m=len(p);sx=sum(a for a,b in p);sy=sum(b for a,b in p)
 cov=sum(a*b for a,b in p)-sx*sy/m;vx=sum(a*a for a,b in p)-sx*sx/m;vy=sum(b*b for a,b in p)-sy*sy/m;v=cov**2/(vx*vy)
 checks.append(dict(GBM_IV=x['GBM_IV'],AD_IV=x['AD_selected_IV'],r2_abs_difference=abs(v-x['r2']),pass_text_precision=abs(v-x['r2'])<=5.1e-7))
assert all(x['pass_text_precision'] for x in checks)
tab(V/'All590_Independent_BED_Check.tsv',checks)
summary=dict(exposure_GBM=10,exposure_AD=59,analysed_AD=57,union=len(union),cross_pairs=len(cross),panel_bound=len(union),exact_intersection=sorted(fs&rs),strong_pairs=sum(x['strong_LD_r2_ge_0_8'] for x in cross),max_cross_pair=max(cross,key=lambda x:x['r2']),outcome_missing_preserved=['rs12151021','rs12972720'],independent_checks=len(checks),MR_tools_or_estimates_changed=False)
js(O/'ExposureSide59_Overlap_Summary.json',summary)
print(json.dumps(summary,ensure_ascii=False,indent=2))
