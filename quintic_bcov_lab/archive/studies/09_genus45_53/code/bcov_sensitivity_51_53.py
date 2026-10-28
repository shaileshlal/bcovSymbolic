import sys,json,time
from math import comb
from fractions import Fraction as F
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).resolve().parent/'engine'))
from bcov_yy import YYRing, Poly,X,v1,u, integrate_gradient,Amplitude
from exact_series import evaluate_poly
from period_frames import frame
from physical_boundaries import q_coefficients, gv_kernel
ring=YYRing(); t0=time.perf_counter()
# Use scale c=1: physical n51,20=c/5**100, avoids massive integer until final
p=Poly({(0,0,0,0,j): F((-1)**(j-20)*comb(30,j-20)) for j in range(20,51)})
P1=F(125,12)-v1/2-F(28,3)*u+F(5,12)*X
p1=ring.covariant(p,51,0);p2=ring.covariant(p1,51,1)
R52=-(p2+2*P1*p1)/2
s0,s1,s2=(R52.coefficient('u',i) for i in range(3))
an52=integrate_gradient(s0,-s1+5*X*s2,-s2)
# input n51,20= 5^-100 due to normalized p coefficient 1
# nth GV response from pureX^20 at g52 = -5^-104
n51=F(1,5**100)
target_gw52_d20=gv_kernel(52,51)*n51
A=target_gw52_d20 / -F(1,5**104)
assert A==F(5**6,6),A
# orbifold vanish t^0 .. t^30, t=psi^5; p52 has line-bundle shift -153
vals,k,_=frame('orbifold',38)
base=evaluate_poly(an52,vals)
assert all(base.coeff(j)==0 for j in range(30))
B=-base.coeff(30)-A
ap52=Poly({(0,0,0,0,j): F(((-1)**(j-20)*comb(30,j-20)))*A for j in range(20,51)}) + Poly({(0,0,0,0,j):F(((-1)**(j-21)*comb(30,j-21)))*B for j in range(21,52)})
P52=an52+ap52
orb=evaluate_poly(P52,vals)
print('A and B',A,B,'orb coeff val',orb.valuation,'nterms',len(P52),'time',time.perf_counter()-t0,flush=True)
assert all(orb.coeff(j)==0 for j in range(31))
vals,k,q=frame('mum',26)
F52=evaluate_poly(P52,vals)*k**51
coeff=q_coefficients(F52,20)
print('MUM lower',[(d,str(coeff[d])) for d in range(21) if coeff[d]],'target',target_gw52_d20,flush=True)
assert all(coeff[j]==0 for j in range(20))
assert coeff[20]==target_gw52_d20
# conifold gap difference poles must vanish
vals,k,flat=frame('conifold',112)
cf=evaluate_poly(P52,vals)*k**51
print('conifold delta valuation',cf.valuation, 'time',time.perf_counter()-t0,flush=True)
assert cf.valuation>=0
assert (P52.partial('v1')-u*P52.partial('v2')-u*(u+5*X)*P52.partial('v3'))==R52
# genus53 linearized, now with *physical* genus52 perturbation
ref=json.load(open(ROOT/'data'/'reference_P2.json'))
P2=Amplitude.from_data(ref).polynomial()
P2_1=ring.covariant(P2,2,0)
P52_1=ring.covariant(P52,52,0)
P52_2=ring.covariant(P52_1,52,1)
R53=-(P52_2+2*P1*P52_1+2*P2_1*p1)/2
s0,s1,s2=(R53.coefficient('u',i) for i in range(3))
an53=integrate_gradient(s0,-s1+5*X*s2,-s2)
print('R53 terms',len(R53),'an53 terms',len(an53),'time',time.perf_counter()-t0,flush=True)
# MUM ΔN(g=53,d=20) from original genus51 and genus52 charges (n52,20=0)
# h=51, m2 derivative of sine kernel
target_gw53_d20=gv_kernel(53,51)*n51
print('gw53 target',target_gw53_d20,flush=True)
# same approach: possible δa53 polynomial X20 (1-X)^m etc;
# probe valuations and leading coeffs at orbifold and conifold
vals,k,_=frame('orbifold',42)
base=evaluate_poly(an53,vals)
print('g53 orbifold valuation',base.valuation,'min coeffs',[(j,str(base.coeff(j))) for j in range(31,33)],flush=True)
vals,k,q=frame('mum',26)
baseM=evaluate_poly(an53,vals)*k**52
cM=q_coefficients(baseM,20)
print('g53 MUM valuation',min([i for i,v in enumerate(cM) if v] or [999]),flush=True)
vals,k,flat=frame('conifold',116)
baseC=evaluate_poly(an53,vals)*k**52
print('g53 conifold valuation',baseC.valuation,'time',time.perf_counter()-t0,flush=True)
json.dump({'scale':'deltaP51=X^20(1-X)^30; physical n51,20=5^-100','alpha_52':str(A),'beta_52':str(B),'n52_anomaly_terms':len(an52),'n52_full_perturbation_terms':len(P52),
'g52_orbifold_val':orb.valuation,'g52_conifold_val':cf.valuation,'g52_mum_first_nonzero':20,'g52_mum_coefficient':str(coeff[20]), 'g52_expected':str(target_gw52_d20),'g53_anomaly_terms':len(an53),'g53_orbifold_val':base.valuation,'g53_conifold_val':baseC.valuation,'g53_mum_coeff_20':str(cM[20]),'g53_mum_target20':str(target_gw53_d20)},open(ROOT/'data'/'physical_sensitivity.json','w'),indent=2)
# Solve 5-dimensional exact physical completion of genus 53 sensitivity,
# with support X^20(1-X)^30 sum_{r=0}^4 c_r X^r.
from exact_series import Jet
basis=[]
for r in range(5):
 basis.append(Poly({(0,0,0,0,j):F((-1)**(j-20-r)*comb(30,j-20-r)) for j in range(20+r,51+r)}))

def rows_for(kind,n,orders):
 vs,kk,qq=frame(kind,n)
 K=kk**52
 vbase=evaluate_poly(an53,vs)*K
 cols=[evaluate_poly(b,vs)*K for b in basis]
 return [(label,[col.coeff(i) for col in cols],vbase.coeff(i)) for label,i in orders]

orows=rows_for('orbifold',39,[('orbifold30',30),('orbifold31',31)])
crows=rows_for('conifold',116,[('conifold-2',-2),('conifold-1',-1)])
mrows=rows_for('mum',24,[('mum20',20)])
mat=[]
for name,row,base0 in mrows+orows+crows:
 target=(target_gw53_d20 / 3125**20) if name=='mum20' else F(0)
 mat.append([*row, target-base0])
 print('constraint',name,'rhs',str(target-base0)[:100],flush=True)
# Gaussian elimination over Q
for k in range(5):
 j=next((r for r in range(k,5) if mat[r][k]),None)
 if j is None:raise RuntimeError(f'singular constraint matrix pivot {k}')
 mat[k],mat[j]=mat[j],mat[k]
 v=mat[k][k];mat[k]=[x/v for x in mat[k]]
 for i in range(5):
  if i!=k and mat[i][k]:
   c0=mat[i][k];mat[i]=[p-c0*q for p,q in zip(mat[i],mat[k])]
sol=[row[5] for row in mat]
P53=an53+sum((s*b for s,b in zip(sol,basis)), Poly())
print('g53 completed c',list(map(str,sol)),flush=True)
assert (P53.partial('v1')-u*P53.partial('v2')-u*(u+5*X)*P53.partial('v3'))==R53
for kind,n,indices in [('mum',24,list(range(21))),('orbifold',40,list(range(32))),('conifold',116,list(range(-2,0)))]:
 vs,k,q=frame(kind,n)
 F_g=evaluate_poly(P53,vs)*k**52
 if kind=='mum':
  coeff=q_coefficients(F_g,20)
  assert all(coeff[d]==0 for d in range(20))
  assert coeff[20]==target_gw53_d20
 elif kind=='orbifold':
  assert all(F_g.coeff(i)==0 for i in indices)
 else:
  assert all(F_g.coeff(i)==0 for i in indices)
 print('verified genus53',kind,'valuation',F_g.valuation,flush=True)
report=json.load(open(ROOT/'data'/'physical_sensitivity.json'))
report.update({'genus53_ambiguity_coeffs_basis_X20(1-X)^30 * poly_4':[str(s) for s in sol],
 'genus53_exact_BCOV_source_terms':len(R53), 'genus53_exact_anomaly_terms':len(an53),
 'genus53_complete_sensitivity_terms':len(P53),'genus53_pure_X_terms':len(P53.substitute({'v1':0,'v2':0,'v3':0})),
 'genus53_conifold_gap_difference_pass':True, 'genus53_orbifold_difference_pass':True,'genus53_MUM_GW_degree0to20_difference_pass':True,
 'convention':'All displayed sensitivity polys scaled by c51=1; multiply by 175*5**100 for actual n51,20=175'})
json.dump(report,open(ROOT/'data'/'physical_sensitivity.json','w'),indent=2)
for _name,_poly in [('delta_P51',p),('delta_R52',R52),('delta_P52',P52),('delta_R53',R53),('delta_P53',P53)]:
 (ROOT/'data'/f'{_name}_exact.json').write_text(json.dumps({
  'object':_name, 'poly_terms':_poly.to_data(),
  'scale':'multiply all coefficients by 175 * 5**100 to obtain the physical n_51_20=175 perturbation',
  'YY_coordinates':['u','v1','v2','v3','X'],
  'scope':'Linear sensitivity/difference, not complete physical genus-g free energy'
 },separators=(',',':'))+'\n')
print('FINISHED physical boundary g51 -> g52 -> g53',round(time.perf_counter()-t0,3),'sec')
