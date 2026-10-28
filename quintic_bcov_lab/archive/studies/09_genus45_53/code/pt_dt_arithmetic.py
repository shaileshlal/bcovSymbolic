"""Exact PT/GV and rank-1 DT/PT leading-coefficient neighborhood checks.

All calculations are formal Laurent expansions (the exponents are stored by
shifting the origin), with coefficients entirely in Python integers.
The rank-one DT0 factor M(-q)^{-200} is used; rank-zero D4-D2-D0 DT is a
*distinct* object and is only referenced via published modular tables.
"""
from math import comb
from pathlib import Path
import json

def conv(a,b,N):
    c=[0]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<=N:c[i+j]+=x*y
    return c

def dt_zero_coeffs(N):
    # M(-q)^{-200} = \prod_{n >= 1} (1-(-q)^n)^{200n}
    ans=[1]+[0]*N
    for n in range(1,N+1):
        b=[0]*(N+1)
        for j in range(N//n+1):
            b[j*n]=comb(200*n,j)*((-1)**((n+1)*j))
        ans=conv(ans,b,N)
    return ans

def pt_gv_top(n,g,N):
    return [n*comb(2*g-2,i) for i in range(N+1)]

def main():
    root=Path(__file__).resolve().parents[1]
    N=6
    z=dt_zero_coeffs(N)
    assert z[:3]==[1,200,19500]
    for g,n in [(51,175),(54,319125)]:
        pt=pt_gv_top(n,g,N)
        dt=conv(pt,z,N)
        if g==51:
            assert pt[:3]==[175,17500,866250]
            assert dt[:3]==[175,52500,7778750]
        print(g,'PT-linear-response',pt[:4], 'rank-one-DT-response',dt[:4])
    # PT top-degree identities using the leading two GV spin contributions.
    pt_q20_next=pt_gv_top(175,51,1)[1]+9700
    pt_q21_next=pt_gv_top(319125,54,1)[1]+7158750
    assert pt_q20_next==27200 and pt_q21_next==40986000
    # Modularity-conditional degree21 leading PT cross-check:
    dt_q21_next=pt_q21_next+z[1]*319125
    assert dt_q21_next==104811000
    out={'DT0_rank1_MacMahon_coefficients':[str(x) for x in z],
         'g51_Q20_n175':{'PT_c0_c1_c2':[175,17500,866250],
                        'DT_rank1_response_c0_c1_c2':[175,52500,7778750],
                        'PT_at_m_minus49_with_g50_input':27200},
         'g54_Q21_modularity_n319125':{'PT_m_minus53':319125,
                        'PT_m_minus52_including_n53_21':40986000,
                        'rank1_DT_m_minus52_including_DT0':dt_q21_next},
         'qualification':'PT/rank1 DT leading terms are algebraically determined by the GV input; they are not independent lattice-selection equations. Rank0 DT modularity is different.'}
    (root/'data'/'pt_dt_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('All exact PT/DT identities passed')

if __name__=='__main__':main()
