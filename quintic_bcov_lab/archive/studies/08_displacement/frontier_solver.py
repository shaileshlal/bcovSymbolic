"""Safe physics-gated integer charge prediction and exact first-response projection.

The default is ABSTENTION for unprotected/non-extremal sites, not hallucinated zeros.
The extremal values are derived from a known theorem, while a learned surrogate
of that theorem is archived under results/model_benchmarks.json.
"""
from pathlib import Path
import csv, json
from build_data import residual_dimension, margin10, top_extremal
ROOT=Path(__file__).resolve().parent


def classify_charge(g:int,d:int)->dict:
    if g<4 or d<1 or d>residual_dimension(g):
        return {'g':g,'d':d,'status':'outside_residual_space','delta_n':None}
    m=margin10(g,d)
    if m<0:
        return {'g':g,'d':d,'status':'certified_zero','delta_n':0}
    if m==0 and d>=20 and d%5==0:
        k=d//5
        n=top_extremal(k)[2]
        return {'g':g,'d':d,'status':'extremal_from_theorem','delta_n':n,
                'sign_phase':k%4}
    return {'g':g,'d':d,'status':'ABSTAIN_unresolved','delta_n':None,
            'untrusted_naive_zero_baseline':0}


def first_nonzero_section_taylor(g:int,d:int,delta_n:int)->dict:
    """Leading delta f_g coefficient if charges at degrees <d are known zero.

    If d=D_g it is the complete reduced-section displacement delta f_g = a*x**d.
    If d<D_g higher x powers depend on higher charges, left uncomputed.
    """
    D=residual_dimension(g)
    if not 1<=d<=D:raise ValueError('not a residual GV degree')
    exponent=4*g-4-5*d
    if exponent<0:raise ValueError('unexpected negative power in quintic residual sector')
    return {'g':g,'d':d,'D':D,'delta_n':delta_n,'taylor_degree':d,
        'coefficient_factor':f'{delta_n}*5^{exponent}',
        'coefficient_exact':str(delta_n*5**exponent),
        'complete_displacement':d==D,
        'statement':('delta_f_g(x) = coefficient*x^d' if d==D
                     else '[x^d] delta_f_g(x) = coefficient; higher powers unknown')}


def export_frontier():
    places=[(51,20),(52,20),(53,20),(54,21),(55,21),(56,22),(76,25),(76,26)]
    rows=[]
    for g,d in places:
        row=classify_charge(g,d)
        if row['delta_n'] not in (None,0):row['first_response']=first_nonzero_section_taylor(g,d,row['delta_n'])
        rows.append(row)
    report={'frontier_results':rows,'g51':first_nonzero_section_taylor(51,20,175),
        'g54_unknown_template':{'g':54,'d':21,'D':21,
            'delta_f_54':'n_(54,21) * 5^107 * x^21',
            'falsifiable_naive_minheight_baseline':0,'prediction_status':'ABSTAIN_NO_EMPIRICAL_LABELS'},
        'general_rule':'If n(g,j)=0 for j<d, the first nonzero reduced-section coefficient is [x^d]delta_f = 5^(4g-4-5d)*n(g,d), exact.'}
    (ROOT/'results'/'frontier_predictions.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(rows,indent=2)[:4000])

if __name__=='__main__':export_frontier()
