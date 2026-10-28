#!/usr/bin/env python3
"""Empirical PT-scale pulse learner (v1.0).

Pure Python, exact Fractions. NEVER reads the conditional reference/validator.
Input training is limited to the supplied Q<=40 PT observations; forecasts
are for unlabelled Q41..53 sites. Shape hypotheses are disclosed below.

The normalized variables are
 A_{s,t}(k) = PT(5k+s,t)/(sigma_s(k)*W_{s,t}(k)),
 W = 5k^2+(2s-5)k+c_s-2t,
 K(s,t) = t-1 (s=0), t+3-s (s>0), r=K-k.

The learned profile hypothesis is
 A_{s,t}(k) = B_{group(s),t} + (-1)^t [P_{s,r}+Q_{s,r}/W].
B is an observed/inferred stabilized scale, P and Q are exact fractions.

The first and second residue pairings are selected from the training table;
they are empirical hypotheses, NOT imported modular/wall-crossing identities.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import csv
import json
import sys

BASE=Path(__file__).resolve().parent
ROOT=BASE.parent
sys.path.insert(0,str(BASE))
from recursive_defects import read,write,phase_fit
from deep_envelope import W,build_amplitudes


def K(s:int,t:int)->int:
    return t-1 if s==0 else t+3-s

def group(s:int)->int:
    return 0 if s==0 else (1 if s in (1,4) else 2)

@dataclass(frozen=True)
class Profile:
    P:F
    Q:F
    min_t:int
    support:int
    total:int
    origin:str

    def predict(self,t:int,w:int)->F:
        assert t>=self.min_t and w
        return (-1 if t%2 else 1)*(self.P+self.Q/F(w))

class Learner:
    def __init__(self,root:Path=ROOT):
        self.root=Path(root)
        cal=read(self.root/'data/shallow_calibration_Q1_Q20.csv')
        self.period,self.phase=phase_fit(cal)
        assert self.period==4
        self.train=read(self.root/'data/shallow_observed_Q6_Q30.csv')+read(self.root/'data/calibration_deep_Q20_Q40.csv')
        assert all(int(z['degree'])<=40 for z in self.train)
        self.A={(s,t,k):v for (s,t),seq in build_amplitudes(self.train,self.phase,self.period).items() for k,v in seq}
        self.B={}
        self.Bsource={}
        self.profiles={}
        self.audit=[]

    def collect(self,s:int,r:int):
        pts=[]
        for (ss,t,k),v in self.A.items():
            if ss!=s or K(s,t)-k!=r or (group(s),t) not in self.B:continue
            w=W(s,t,k)
            if not w:continue
            y=(-1 if t%2 else 1)*(v-self.B[group(s),t])
            pts.append((t,k,w,y))
        return sorted(pts)

    def seed_B(self):
        candidates=defaultdict(list)
        for (s,t,k),a in self.A.items():
            if k>=K(s,t):candidates[group(s),t].append((s,k,a))
        for key,rows in candidates.items():
            cnt=Counter(v for _,_,v in rows)
            best=max(cnt,key=lambda v:(cnt[v],max(k for s,k,a in rows if a==v)))
            if key[1]<6 and cnt[best]<2:continue
            # Early frontier exceptions at t<=2 make mode more robust than
            # blindly using k=K. At t>=6 the mature scale is unambiguous.
            if key[1]>=6:
                assert len(cnt)==1,(key,cnt)
            self.B[key]=best
            self.Bsource[key]='observed_mature'

    def constants_r1(self):
        for s in range(5):
            pts=self.collect(s,1)
            if len(pts)<2:continue
            counts=Counter(y for t,k,w,y in pts)
            value,n=max(counts.items(),key=lambda q:(q[1],max(t for t,k,w,y in pts if y==q[0])))
            if n<2:continue
            matched=[t for t,k,w,y in pts if y==value]
            self.profiles[s,1]=Profile(value,F(0),min(matched),n,len(pts),'constant_r1')

    @staticmethod
    def affine_fit(a,b):
        t0,k0,w0,y0=a;t1,k1,w1,y1=b
        if w0==w1:raise ValueError('same W')
        P=(y0*w0-y1*w1)/F(w0-w1)
        Q=(y0-P)*w0
        return P,Q

    def robust_affine(self,s,r,min_support=3):
        pts=self.collect(s,r)
        choices=[]
        for i,a in enumerate(pts):
            for b in pts[i+1:]:
                if a[2]==b[2]:continue
                P,Q=self.affine_fit(a,b)
                matches=[t for t,k,w,y in pts if y==P+Q/w]
                if len(matches)<min_support:continue
                if max(matches)!=max(t for t,k,w,y in pts):continue
                # choose largest exact consensus; prefer earlier support and
                # simpler denominators only on complete ties.
                choices.append((len(matches),-min(matches),-len(str(P))-len(str(Q)),P,Q,matches))
        if not choices:return None
        count,_,_,P,Q,matched=max(choices)
        pr=Profile(P,Q,min(matched),count,len(pts),'affine_consensus')
        self.profiles[s,r]=pr
        return pr

    def infer_B(self):
        proposed=defaultdict(list)
        for (s,t,k),a in sorted(self.A.items()):
            key=(group(s),t)
            if t<6 or key in self.B or k>=K(s,t):continue
            r=K(s,t)-k
            pr=self.profiles.get((s,r))
            if pr is None or t<pr.min_t:continue
            newB=a-pr.predict(t,W(s,t,k))
            proposed[key].append(((s,t,k),newB))
        added=0
        for key,options in sorted(proposed.items()):
            # Never select between conflicting B values using future targets.
            distinct={v for site,v in options}
            if len(distinct)>1:
                self.audit.append({'type':'inconsistent_B','group':key,'proposals':len(options)})
                continue
            val=distinct.pop()
            self.B[key]=val
            self.Bsource[key]='inferred_from_'+str(options[0][0])
            self.audit.append({'type':'inferred_B','group':key,'source':options[0][0]})
            added+=1
        return added

    def pair_r2(self):
        # First test the discovered sign antisymmetry on an independently
        # overdetermined pair (s0 vs s2, deficit2).
        assert self.profiles[0,1].P == -self.profiles[2,1].P
        assert self.profiles[3,1].P == -self.profiles[4,1].P
        a0=[z for z in self.collect(0,2) if z[0]>=6][-2:]
        a2=[z for z in self.collect(2,2) if z[0]>=6][-2:]
        assert len(a0)==len(a2)==2
        p0,q0=self.affine_fit(*a0);p2,q2=self.affine_fit(*a2)
        assert (p0,q0)==(-p2,-q2),(p0,q0,p2,q2)
        for s,pp,qq,pts in ((0,p0,q0,a0),(2,p2,q2,a2)):
            self.profiles[s,2]=Profile(pp,qq,min(z[0] for z in pts),len(pts),len(self.collect(s,2)), 'r2_anti_pair_verified')
        self.audit.append({'type':'r2_anti_pair_verified','pair':'0-2','two_independent_fits_agree':True})

        # The next pair has only one usable anchor per residue, but we transfer
        # the sign-antisymmetry learned at r1 and tested at r2 above.
        a3=[z for z in self.collect(3,2) if z[0]>=6]
        a4=[z for z in self.collect(4,2) if z[0]>=6]
        assert len(a3)==len(a4)==1,(a3,a4)
        t3,k3,w3,y3=a3[0];t4,k4,w4,y4=a4[0]
        assert w3!=w4
        # y3=P+Q/w3, y4=-(P+Q/w4)
        P=(y3*w3+y4*w4)/F(w3-w4)
        Q=(y3-P)*w3
        assert y4==-P-Q/F(w4)
        self.profiles[3,2]=Profile(P,Q,t3,1,1,'antisym_pair_3_4_two_joint_points')
        self.profiles[4,2]=Profile(-P,-Q,t4,1,1,'antisym_pair_3_4_two_joint_points')
        self.audit.append({'type':'r2_anti_pair_extrapolation','pair':'3-4','total_calibration_points':2})

    def deeper_r4(self):
        # One additional two-point, forward-in-t interpolation at s1,r4.
        pts=[z for z in self.collect(1,4) if z[0]>=6]
        assert len(pts)==2,(pts,'r4 requires two independent anchors')
        P,Q=self.affine_fit(*pts)
        self.profiles[1,4]=Profile(P,Q,min(z[0] for z in pts),2,len(self.collect(1,4)), 'r4_two_point_extrapolation')
        self.audit.append({'type':'two_point_extrapolation','s':1,'r':4,'calibration_t':[x[0] for x in pts]})

    def fit(self):
        self.seed_B()
        self.constants_r1()
        # r2 and r3 at residue1 are overdetermined by earlier degrees.
        assert self.robust_affine(1,2,min_support=3)
        assert self.robust_affine(1,3,min_support=3)
        assert self.infer_B()==3
        self.pair_r2()
        assert self.infer_B()==2
        self.deeper_r4()
        return self

    def predict(self,sites):
        out=[]
        for item in sorted(sites,key=lambda x:(int(x['degree']),int(x['shell']))):
            d=int(item['degree']);s=d%5;k=d//5;t=int(item['shell']);r=K(s,t)-k
            bg=group(s),t;source='';val=None;status='abstain'
            if bg in self.B and r<=0:
                val=self.B[bg];source='mature_scale'
            elif bg in self.B and (s,r) in self.profiles:
                prof=self.profiles[s,r]
                if t>=prof.min_t:
                    val=self.B[bg]+prof.predict(t,W(s,t,k))
                    source=prof.origin
            if val is not None:
                raw_pt=self.phase[s,k%self.period]*W(s,t,k)*val
                if raw_pt.denominator!=1:
                    status='nonintegral_candidate';source='';val=None
                else:
                    status='predicted';prediction=str(raw_pt.numerator)
            if val is None:prediction=''
            out.append({'degree':d,'shell':t,'residue':s,'k':k,'front_distance':r,
                'prediction':prediction,'status':status,'source':source,
                'B_origin':self.Bsource.get(bg,''),
                'profile_support':self.profiles[s,r].support if r>0 and (s,r) in self.profiles else '',
                'profile_origin':self.profiles[s,r].origin if r>0 and (s,r) in self.profiles else ''})
        return out

    def serialize(self):
        return {
            'training_max_degree':40,'reference_in_training':False,
            'period':self.period,'profile_equation':'A=B+(-1)^t*(P+Q/W)',
            'front':'K=t-1 (s=0), K=t+3-s (s>0)',
            'B': [{'class':g,'shell':t,'rational':str(v),'source':self.Bsource[g,t]} for (g,t),v in sorted(self.B.items())],
            'profiles':[{'residue':s,'front_distance':r,'P':str(pr.P),'Q':str(pr.Q),
                'min_t':pr.min_t,'support':pr.support,'total':pr.total,'origin':pr.origin}
              for (s,r),pr in sorted(self.profiles.items())],
            'audit':self.audit,
        }

def run(root=ROOT):
    m=Learner(root).fit()
    sites=read(Path(root)/'data/forecast_sites_Q41_Q53.csv')
    out=m.predict(sites)
    assert len(out)==65
    results=Path(root)/'results';results.mkdir(parents=True,exist_ok=True)
    write(results/'pulse_predictions_Q41_Q53.csv',out)
    (results/'fitted_pulses.json').write_text(json.dumps(m.serialize(),indent=2)+'\n')
    print('PT predictions',sum(bool(z['prediction']) for z in out),
          'abstentions',sum(not z['prediction'] for z in out),'out of 65')
    return out

if __name__=='__main__':
    run()
