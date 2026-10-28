"""Create exactly labelled/censored integer-displacement dataset for quintic residual BCOV ambiguities.

The residual lattice uses unrestricted integer GV coordinates, with min-height origin n_min=0.
Never encode unresolved n_(g,d) as zero. Mathematical sources:
  https://arxiv.org/abs/2210.13411 (vanishing & extremal theorem)
  https://scipost.org/SciPostPhys.16.3.079/pdf (equation 135, omega=kappa=5)
Prior exact residual audit bundled under source/.
"""
from __future__ import annotations
from pathlib import Path
import csv, gzip, json

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'
SOURCE=ROOT/'source'
MAX_G=1000


def residual_dimension(g:int)->int:
    return (2*g-2)//5


def margin10(g:int,d:int)->int:
    """Ten times the Castelnuovo upper-genus allowance G(d)-g."""
    return d*d+5*d+10-10*g


def top_extremal(k:int)->tuple[int,int,int]:
    """Proven extremal n[g,5k] for k>=2, Liu--Ruan; only k>=4 in this residual sector."""
    if k<2:
        raise ValueError('extremal formula branch requires k>=2')
    d=5*k
    g=1+5*k*(k+1)//2
    n=(-1 if (k*(k-1)//2)%2 else 1)*25*(1+k*(k-1)//2)
    return g,d,n


def generate(max_g:int=MAX_G):
    fields=['g','r','d','D','g_mod5','r_mod5','d_mod5','k_extremal',
            'margin10','label_status','displacement','origin_minimum','data_provenance']
    summary={'max_g':max_g, 'total':0,'certified_zero':0,
             'extremal_nonzero':0,'unresolved':0}
    known=[]
    with gzip.open(DATA/'all_lattice_sites_g2_g1000.csv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields)
        w.writeheader()
        for g in range(4,max_g+1):
            for d in range(1,residual_dimension(g)+1):
                m=margin10(g,d)
                if m<0:
                    status='certified_zero';n=0;prov='Liu_Ruan_Castelnuovo_vanishing';k=''
                elif m==0 and d%5==0 and d>=20:
                    k=d//5
                    assert top_extremal(k)[0:2]==(g,d)
                    status='extremal_nonzero';n=top_extremal(k)[2];prov='Liu_Ruan_extremal_curves';
                else:
                    status='unresolved';n='';k='';prov='unknown_do_not_impute'
                row={'g':g,'r':d-1,'d':d,'D':residual_dimension(g),
                     'g_mod5':g%5,'r_mod5':(d-1)%5,'d_mod5':d%5,
                     'k_extremal':k,'margin10':m,'label_status':status,
                     'displacement':n,'origin_minimum':0,'data_provenance':prov}
                w.writerow(row)
                summary['total']+=1;summary[status]+=1
                if status!='unresolved':known.append(row)
    with (DATA/'certified_labels.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(known)
    # Compact extrapolation slice: no invented values, include unknown charges.
    with (DATA/'frontier_g49_g80.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for g in range(49,81):
            for d in range(18,residual_dimension(g)+1):
                m=margin10(g,d)
                if m<0:status='certified_zero';n=0;prov='Castelnuovo';k=''
                elif m==0 and d%5==0 and d>=20:
                    k=d//5;status='extremal_nonzero';n=top_extremal(k)[2];prov='extremal theorem'
                else:status='unresolved';n='';prov='unknown_do_not_impute';k=''
                w.writerow({'g':g,'r':d-1,'d':d,'D':residual_dimension(g),
                    'g_mod5':g%5,'r_mod5':(d-1)%5,'d_mod5':d%5,
                    'k_extremal':k,'margin10':m,'label_status':status,
                    'displacement':n,'origin_minimum':0,'data_provenance':prov})
    top=[dict(k=k, g=top_extremal(k)[0], d=top_extremal(k)[1],r=5*k-1,
             n=top_extremal(k)[2],phase=k%4,abs_n=abs(top_extremal(k)[2]),
             split='train' if k<=11 else 'test',source='proven_extremal_theorem')
             for k in range(4,80) if top_extremal(k)[0]<=max_g]
    with (DATA/'extremal_displacements.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(top[0]));w.writeheader();w.writerows(top)

    # Reference validation: all 22 genuine residual g4..12 charges are zero.
    audit=json.loads((SOURCE/'exact_lattice_audit.json').read_text())
    cross=[]
    for item in audit['records']:
        g=int(item['g'])
        for r, (ng, nmin) in enumerate(zip(item['physical_GV'],item['minimum_GV'])):
            if margin10(g,r+1)>=0:
                raise AssertionError(('expected vanishing for g<=12',g,r))
            assert int(ng)==int(nmin)==0
            cross.append({'g':g,'r':r,'d':r+1,'physical_charge':ng,
                          'minheight_charge':nmin,'displacement':int(ng)-int(nmin),
                          'particular_independent':item['independent_particular_matches_reference_split']})
    assert len(cross)==22
    with (DATA/'exact_g4_g12_calibration.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(cross[0]));w.writeheader();w.writerows(cross)
    summary.update(extremal_train=sum(t['split']=='train' for t in top),
                    extremal_test=sum(t['split']=='test' for t in top),
                    certified_labels=len(known),reference_calibration=len(cross),
                    period4_signs=[1,1,-1,-1])
    (DATA/'dataset_manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return summary

if __name__=='__main__':
    generate()
