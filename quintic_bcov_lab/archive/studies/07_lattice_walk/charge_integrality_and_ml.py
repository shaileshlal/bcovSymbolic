"""Two tests: integral transport to degree 26; lattice-quantized held-out ML.

No physical test-genus coefficient used in the lattice decoding. For g11-12
only, the available particular is reference-derived, flagged as hindsight.
"""
from higher_integrality import *
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
import csv

def exact_transport(N=26):
 out=[]
 for g in range(4,13):
  K,m,D=dimensions(g);H=full_response(g,N)
  inversecols=[solve(H[:D],[Q(i==j) for i in range(D)]) for j in range(D)]
  A=[[sum((H[d][r]*inversecols[j][r] for r in range(D)),Q(0)) for j in range(D)] for d in range(N)]
  bad=[(d+1,j) for d in range(N) for j in range(D) if A[d][j].denominator!=1]
  assert not bad,(g,bad)
  out.append({'g':g,'D':D,'tested_degrees':N,'all_integer':True,
     'n_integrality_entries':N*D,
     'unit_charge_response_selected':{str(j+1):[int(A[d][j]) for d in range(N)] for j in range(D)}})
 (ROOT/'results'/'integral_transport_degree26.json').write_text(json.dumps(out,indent=2)+'\n')
 return out

def project_ml():
 r=json.loads((ROOT/'results'/'exact_lattice_audit.json').read_text())
 rec={x['g']:x for x in r['records']}
 inputfile=ROOT/'Quintic_Intermediate_ML'/'results'/'predictions.csv'
 rows=list(csv.DictReader(inputfile.open()))
 scores=[]
 for split in ('train_to_8','expanding_9_to_12'):
  models=sorted(set(x['model'] for x in rows if x['split']==split))
  for model in models:
   for g in range(9,13):
    z=sorted((x for x in rows if x['split']==split and x['model']==model and int(x['genus'])==g),key=lambda x:int(x['r']))
    if not z: continue
    x=rec[g]; D=x['D']; H=[[Q(q) for q in rr] for rr in x['H']];baseline=[Q(v) for v in x['n_anchor']]
    assert len(z)==D
    chosen=[];actual=[Q(v) for v in x['physical_c']];charge=[]
    with localcontext() as ctx:
     ctx.prec=95
     for k,item in enumerate(z):
      signed=(1 if int(item['pred_sign'])>0 else -1)
      c_hat=Decimal(signed)*(Decimal(10)**Decimal(item['pred_log_mag']))
      b=baseline[k]+sum((H[k][i]*chosen[i] for i in range(k)),Q(0))
      R=H[k][k]
      candidate=(Decimal(b.numerator)/Decimal(b.denominator)+
          Decimal(R.numerator)/Decimal(R.denominator)*c_hat)
      pick=int(candidate.to_integral_value(rounding=ROUND_HALF_EVEN))
      chosen.append((Q(pick)-b)/R)
      charge.append(pick)
    scores.append({'split':split,'model':model,'g':g,
         'independent_particular':x['independent_particular_matches_reference_split'],
         'n_exact':sum(a==b for a,b in zip(chosen,actual)),
         'n_total':D,'exact_all':chosen==actual,
         'predicted_GV':charge,
         'label_GV':[0]*D})
 (ROOT/'results'/'ml_lattice_projection.json').write_text(json.dumps(scores,indent=2)+'\n')
 return scores

if __name__=='__main__':
 test=exact_transport(26)
 print('Integral transport:',len(test),'genera,',sum(x['n_integrality_entries'] for x in test),'integer response entries through degree 26')
 scores=project_ml()
 for split in ('train_to_8','expanding_9_to_12'):
  for model in ('mixed_plus_genus__majority','v01_trees__combined_v01','mixed_residual__majority'):
   z=[x for x in scores if x['split']==split and x['model']==model]
   zz=[x for x in z if x['independent_particular']]
   print(split,model,'all',sum(x['n_exact'] for x in z),'/',sum(x['n_total'] for x in z),'independent',sum(x['n_exact'] for x in zz),'/',sum(x['n_total'] for x in zz),'charges',[x['predicted_GV'] for x in z])
