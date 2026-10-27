"""Extract Wolfram GridBox integer tables from quinticGVDTPT.m safely.

The .m file is a Mathematica *typeset-box* export, not ordinary plaintext
Mathematica list syntax. This parser treats it as data and never evaluates it.
Unknown boxes are represented by None; they are NOT assumed zero.
"""
from pathlib import Path
import argparse,csv,gzip,json,re

NAMES=('GVTab','PTTab','DTTab','PTTabGv','GVTabExt','GVTabExt2')

def split_top(s):
    parts=[]; quoted=False; esc=False; brace=brack=0;start=0
    for i,c in enumerate(s):
        if quoted:
            if esc:esc=False
            elif c=='\\':esc=True
            elif c=='"':quoted=False
            continue
        if c=='"':quoted=True
        elif c=='{':brace+=1
        elif c=='}':brace-=1
        elif c=='[':brack+=1
        elif c==']':brack-=1
        elif c==',' and brace==brack==0:
            parts.append(s[start:i].strip());start=i+1
    assert not quoted and brace==brack==0,'Unbalanced expression'
    parts.append(s[start:].strip())
    return parts

def parse_table(source,name):
    head=re.search(r'(?m)^'+re.escape(name)+r'\s*=\s*\\!\\\(\\\*',source)
    if not head:raise ValueError('missing '+name)
    beg=source.index('GridBox[{',head.start())+len('GridBox[{')
    end=source.index('\n},\nGridBoxAlignment',beg)
    table=source[beg:end].strip().rstrip(',')
    rows=[]
    for row in split_top(table):
        if not(row.startswith('{') and row.endswith('}')):raise ValueError('invalid row in '+name)
        data=[]
        for cell in split_top(row[1:-1]):
            if re.fullmatch(r'"-?\d+"',cell):data.append(int(cell[1:-1]))
            elif m:=re.fullmatch(r'RowBox\[\{"-",\s*"(\d+)"\}\]',cell):data.append(-int(m[1]))
            else:data.append(None)
        rows.append(data)
    width=len(rows[0]);assert all(len(row)==width for row in rows)
    return rows

def export(src,out):
    source=Path(src).read_text(encoding='utf8')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    tables={name:parse_table(source,name) for name in NAMES}
    meta={}
    for name,rows in tables.items():
        total=sum(x is not None for r in rows for x in r)
        meta[name]={'rows':len(rows),'columns':len(rows[0]),
                    'numeric_cells':total,'unknown_cells':len(rows)*len(rows[0])-total,
                    'nonzero_cells':sum(x not in (None,0) for r in rows for x in r)}
        if name in ('GVTab','GVTabExt'):
            with gzip.open(out/(name+'.csv.gz'),'wt',encoding='utf8',newline='') as fh:
                w=csv.writer(fh);w.writerow(['genus']+[f'degree_{d}' for d in range(1,len(rows[0])+1)])
                for genus,row in enumerate(rows):
                    w.writerow([genus]+[str(x) if x is not None else '' for x in row])
        else:
            with gzip.open(out/(name+'.jsonl.gz'),'wt',encoding='utf8') as fh:
                for i,row in enumerate(rows):
                    fh.write(json.dumps({'row_index_0':i,'values':[str(x) if x is not None else None for x in row]},separators=(',',':'))+'\n')
    diffs=[];equal=0
    a,b=tables['GVTab'],tables['GVTabExt']
    for g in range(min(len(a),len(b))):
        for d in range(1,min(len(a[g]),len(b[g]))+1):
            x,y=a[g][d-1],b[g][d-1]
            if y is None:continue
            if x==y:equal+=1
            else:diffs.append((g,d,x,y,y-x))
    with (out/'GVTab_vs_GVTabExt_disagreements.csv').open('w',newline='',encoding='utf8') as fh:
        w=csv.writer(fh);w.writerow(['genus','degree','GVTab','GVTabExt','GVTabExt_minus_GVTab'])
        w.writerows(diffs)
    meta['GV_overlap']={'identical':equal,'different':len(diffs)}
    meta['quintic']={'chi':-200,'kappa':5,'c2':50}
    (out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
    # Lightweight assertions on unambiguous quintic cases.
    assert a[0][:4]==[2875,609250,317206375,242467530000]
    assert a[1][2]==609250 and a[2][3]==534750 and a[3][3]==8625
    assert b[64][23]==-6881505911867623750
    assert b[64][26] is None
    return meta

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('source',nargs='?',default='quinticGVDTPT.m');ap.add_argument('--out',default='quintic_integer_tables');args=ap.parse_args()
    print(json.dumps(export(args.source,args.out),indent=2))
