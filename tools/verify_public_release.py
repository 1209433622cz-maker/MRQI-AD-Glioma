"""Verify this public mirror without network, secrets or statistical fitting."""
import csv,hashlib,json,sys,collections,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
    errors=[];manifest=rows(ROOT/'PUBLIC_MANIFEST.tsv')
    for r in manifest:
        p=(ROOT/r['path']).resolve()
        if not p.is_relative_to(ROOT):errors.append('unsafe manifest path');continue
        if not p.is_file():errors.append('missing '+r['path']);continue
        b=p.read_bytes()
        if len(b)!=int(r['bytes']) or hashlib.sha256(b).hexdigest()!=r['sha256']:errors.append('bytes/hash '+r['path'])
    a=rows(ROOT/'results/S14a_All81_Aligned_Complete_Evidence_2026-10-16.tsv')
    b=rows(ROOT/'results/S14b_Original43_Aligned_Complete_Evidence_2026-10-16.tsv')
    if len(a)!=81 or len(b)!=43 or len(a[0])!=163 or len(b[0])!=163:errors.append('81/43 shape')
    identity=lambda r:(r['HGNC_ID'],r['gene'],r['direction'])
    if len({identity(r) for r in a})!=81 or not {identity(r) for r in b}<={identity(r) for r in a}:errors.append('immutable identity/subset')
    states=collections.Counter(r['BA17_state'] for r in b)
    if states!={'DIRECTION_ONLY_NO_FDR':29,'POINT_OPPOSITE_NO_FDR':12,'MISSING_OR_NOT_ESTIMABLE':2}:errors.append('BA17 negative/conflict/missing retention')
    if {r['gene'] for r in b if r['BA17_state']=='MISSING_OR_NOT_ESTIMABLE'}!={'OR7A5','PRPF4'}:errors.append('BA17 missing identities')
    def q(x):
        try:return float(x)
        except (ValueError,TypeError):return math.nan
    if any(q(r['BA17_q43'])<.05 for r in b):errors.append('BA17 primary FDR changed')
    def joint_expected(r):
        ad_expected=1 if r['direction'].startswith('R1') else -1
        return q(r['AD_public_protein_HC3_Q_primary43'])<.05 and q(r['GBMprotein_main_HC3_primary_q'])<.05 and q(r['AD_public_protein_OLS_beta'])*ad_expected>0 and q(r['GBMprotein_main_beta'])*ad_expected<0
    if any(joint_expected(r) for r in b):errors.append('joint expected-direction HC3 negative changed')
    for name,n in [('S02a_All16068_Discovery_Effects_And_Sensitivity.tsv',16068),('S02b_All32136_Directional_IUT_Tests.tsv',32136),('S03a_All508_AD_Validation_And_Discovery_Context.tsv',508),('S13_Source_Model_Family_Alignment_Registry_2026-10-16.tsv',124)]:
        if len(rows(ROOT/'results'/name))!=n:errors.append('row count '+name)
    print(json.dumps({'status':'FAIL' if errors else 'PASS','manifest_files':len(manifest),'complete81':len(a),'original43':len(b),'new_models':0,'errors':errors},indent=2))
    return 1 if errors else 0
if __name__=='__main__':sys.exit(main())
