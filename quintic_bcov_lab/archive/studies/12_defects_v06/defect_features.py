"""Time-causal integer-valued finite-difference features for quintic PT boundary shells.

Training uses only Q<=20.  A first update uses Q21..25; the next block Q26..30
is scored by a separate validator.  No physics formula or modular generator is used.
"""
from __future__ import annotations
import math
from collections import defaultdict
from fractions import Fraction

DELTA0={0:0,1:1,2:2,3:2,4:1}
PEER={1:4,4:1,2:3,3:2}

def integer_rows(rows):
    return {(int(r['degree']),int(r['shell'])):int(r['pt']) for r in rows if r['pt'] not in ('','None')}

def grouped_signed_rows(rows,signs,period=4,first_degree=6):
    data=integer_rows(rows)
    groups=defaultdict(list)
    for (d,t),val in data.items():
        if d<first_degree:continue
        s=d%5;k=d//5
        if (s,k%period) not in signs:continue
        groups[(s,t)].append((k,signs[(s,k%period)]*val))
    return {key:sorted(v) for key,v in groups.items()}

def val5(n):
    n=abs(int(n))
    if n==0:return None
    out=0
    while n%5==0:n//=5;out+=1
    return out

def val2(n):
    n=abs(int(n))
    if n==0:return None
    return (n&-n).bit_length()-1

def sg(x):return 1 if x>0 else -1 if x<0 else 0

def slog(x):return sg(x)*math.log10(1+abs(x))

def flratio(n,d):
    if not d:return 0.0
    return float(Fraction(n,d))

def current_state(groups,s,t):
    seq=groups[(s,t)]
    if len(seq)<3:return None
    # In production, only the most recent three labels are available.
    x=[v for k,v in seq]
    C=x[-1]-2*x[-2]+x[-3]
    prevC=x[-2]-2*x[-3]+x[-4] if len(x)>=4 else None
    return {'last_k':seq[-1][0], 'last':x[-1], 'penult':x[-2],
            'curvature':C, 'prev_curvature':prevC,'prev_defect':C-prevC if prevC is not None else None}

def feat_from_state(groups,s,t):
    own=current_state(groups,s,t)
    if own is None:return None
    p=PEER.get(s)
    peer=current_state(groups,p,t) if p is not None and (p,t) in groups else None
    C=own['curvature']; old=own['prev_curvature']
    pc=peer['curvature'] if peer else None
    gap=pc-C if pc is not None else None
    bracket=bool(old is not None and pc is not None and min(old,C)<=pc<=max(old,C))
    # These are features of already-observed integers, NEVER of a future target.
    features=dict(residue=s,shell=t,k_last=own['last_k'],t_minus_k=t-own['last_k'],
      absolute_C_log=math.log10(1+abs(C)), signed_C_log=slog(C),
      absolute_PT_log=math.log10(1+abs(own['last'])),
      local_ratio_log=math.log10((1+abs(own['last']))/(1+abs(own['penult']))),
      v2_C=val2(C) if C else 0, v5_C=val5(C) if C else 0,
      v5_PT=val5(own['last']) if own['last'] else 0,
      relative_curvature_log=math.log10((1+abs(C))/(1+abs(own['last']))),
      peer_available=int(pc is not None),
      peer_gap_log=slog(gap) if gap is not None else 0.0,
      peer_gap_ratio=flratio(gap,1+abs(C)) if gap is not None else 0.0,
      prior_defect_available=int(old is not None),
      prior_defect_log=slog(own['prev_defect']) if old is not None else 0.0,
      prior_defect_ratio=flratio(own['prev_defect'],1+abs(C)) if old is not None else 0.0,
      curvature_stable=int(old==C) if old is not None else 0,
      mirror_bracket=int(bracket))
    return dict(state=own,peer_curvature=pc,peer_bracket=bracket,features=features)

def event_records(history, future, signs, period=4):
    groups=grouped_signed_rows(history,signs,period)
    outcomes=integer_rows(future)
    records=[]
    for (s,t),qseq in sorted(groups.items()):
        o=feat_from_state(groups,s,t)
        if o is None:continue
        own=o['state'];next_d=5*(own['last_k']+1)+s
        if (next_d,t) not in outcomes:continue
        y=signs[(s,(own['last_k']+1)%period)]*outcomes[next_d,t]
        nextC=y-2*own['last']+own['penult']
        residual=nextC-own['curvature']
        records.append(dict(degree=next_d,residue=s,shell=t,k=own['last_k']+1,
          actual_PT=outcomes[next_d,t],signed_PT=y,
          actual_defect=residual,active=int(residual!=0),
          own_curvature=own['curvature'],peer_curvature=o['peer_curvature'],
          own_prevcurvature=own['prev_curvature'],own_last=own['last'],
          own_penult=own['penult'],peer_bracket=o['peer_bracket'],features=o['features']))
    return records
