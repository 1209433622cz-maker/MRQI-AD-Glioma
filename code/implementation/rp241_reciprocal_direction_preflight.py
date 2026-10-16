"""Synthetic-only QA for the RP2.4.1 reciprocal-direction intersection-union test.
No human effect file is read. This does not execute AD--GBM discovery.
"""
import json, math, pathlib
import numpy as np
import pandas as pd

ROOT = pathlib.Path('H:/SCI2/MRQI')

def one_sided(beta, p_two):
    if not (np.isfinite(beta) and np.isfinite(p_two)): return np.nan, np.nan
    if not 0 <= p_two <= 1: raise ValueError('Invalid two-sided P')
    if beta == 0:
        if p_two != 1: raise ValueError('Exact zero beta requires symmetric-test P=1')
        return .5, .5
    # Construct the small tail directly; 1-(1-p/2) would erase extreme valid P.
    return (p_two/2, 1-p_two/2) if beta > 0 else (1-p_two/2, p_two/2)

def adjust_full_family(p, method='BH'):
    """Keep missing hypotheses untestable, while retaining the full family denominator."""
    p = np.asarray(p, dtype=float); n = len(p)
    valid = np.flatnonzero(np.isfinite(p)); result = np.full(n, np.nan)
    order = valid[np.argsort(p[valid], kind='stable')]
    factor = sum(1/i for i in range(1,n+1)) if method == 'BY' else 1
    if method not in ['BH','BY']: raise ValueError('Unsupported adjustment')
    q = p[order] * n / np.arange(1, len(order)+1) * factor
    result[order] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
    return result

def reciprocal(ad_beta, ad_p, gbm_beta, gbm_p):
    ap, am = one_sided(ad_beta, ad_p); gp, gm = one_sided(gbm_beta, gbm_p)
    if not np.isfinite([ap,am,gp,gm]).all(): return np.nan, np.nan
    return max(ap, gm), max(am, gp)

def main():
    out = ROOT / '3_results/rp241_reciprocal_direction_preflight'
    receipt = ROOT / '0_admin/validation/RP241_Reciprocal_Direction_Synthetic_QA_2026-10-14.json'
    if receipt.exists(): raise RuntimeError('Frozen QA exists')
    out.mkdir(parents=True, exist_ok=True)
    cases=[('R1_known',1,.0002,-1,.0004),('R2_known',-1,.0004,1,.0002),
           ('Same_sign',1,.0002,1,.0004),('Missing_source',1,.0002,-1,np.nan),
           ('Zero_effect',0,1,-1,.001),('Extreme_R1',1,1e-200,-1,2e-200),
           ('Extreme_R2',-1,2e-200,1,1e-200)]
    rows=[];ps=[]
    for name,ab,ap,gb,gp in cases:
        p1,p2=reciprocal(ab,ap,gb,gp);ps.extend([p1,p2])
        rows.append(dict(case=name,AD_beta=ab,AD_p2=ap,GBM_beta=gb,GBM_p2=gp,R1_p=p1,R2_p=p2))
    bh=adjust_full_family(ps);by=adjust_full_family(ps,'BY')
    for i,r in enumerate(rows):
        r.update(R1_BH=bh[2*i],R2_BH=bh[2*i+1],R1_BY=by[2*i],R2_BY=by[2*i+1])
    assert rows[0]['R1_BH'] < .05 and rows[1]['R2_BH'] < .05
    assert rows[2]['R1_p'] > .5 and rows[2]['R2_p'] > .5
    assert np.isnan(rows[3]['R1_BH']) and rows[4]['R1_p']==.5
    assert rows[5]['R1_p']==1e-200 and rows[6]['R2_p']==1e-200
    for p in [1e-300,1e-200,1e-30,1e-10,.5,1]:
        assert one_sided(1,p)[0]==p/2 and one_sided(-1,p)[1]==p/2
    assert np.all(by[np.isfinite(by)]>=bh[np.isfinite(bh)])
    for r in rows:
        if r['case'] != 'Missing_source':
            swapped=reciprocal(-r['AD_beta'],r['AD_p2'],-r['GBM_beta'],r['GBM_p2'])
            assert np.allclose(swapped,[r['R2_p'],r['R1_p']])
    pd.DataFrame(rows).to_csv(out/'Known_Direction_Fixtures_2026-10-14.tsv',sep='\t',index=False,na_rep='NA')
    # Check composite-null IUT calibration in named synthetic normal scenarios.
    rng=np.random.default_rng(20261014);n=100000;sim=[]
    for name,ma,mg,rho,null in [('both_null_independent',0,0,0,True),('both_null_positive_cor',0,0,.6,True),
                              ('both_null_negative_cor',0,0,-.6,True),('AD_null_GBM_negative',0,-6,0,True),
                              ('AD_positive_GBM_null',6,0,0,True),('both_negative_wrong_R1',-6,-6,0,True),
                              ('true_R1',6,-6,0,False)]:
        a=rng.normal(size=n);g=rho*a+math.sqrt(1-rho*rho)*rng.normal(size=n)
        a+=ma;g+=mg
        plus=np.fromiter((.5*math.erfc(z/math.sqrt(2)) for z in a),float,n)
        minus=np.fromiter((.5*math.erfc(-z/math.sqrt(2)) for z in g),float,n)
        rate=float(np.mean(np.maximum(plus,minus)<.05));tol=3*math.sqrt(.05*.95/n)
        if null:assert rate<=.05+tol
        sim.append(dict(scenario=name,n=n,R1_rejection_fraction=rate,composite_null=null,three_SE_tolerance=tol))
    pd.DataFrame(sim).to_csv(out/'Composite_Null_Simulation_2026-10-14.tsv',sep='\t',index=False)
    result=dict(status='PASS_SYNTHETIC_ONLY',two_direction_family=len(ps),missing_hypotheses_retained=2,
                extreme_small_P_tail_precision='PASS_DIRECT_SMALL_TAIL_NOT_SUBTRACTION_FROM_ONE',
                input_human_effect_files=0,GBM_effects_read=False,independent_AD_effects_read=False,
                analytic_rule='MAX_OF_VALID_SIGNED_ONE_SIDED_P__INTERSECTION_UNION',
                simulations=sim,limitation='Monte Carlo fixtures are not a proof of BH under arbitrary inter-gene dependence; BY sensitivity retained.')
    receipt.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
