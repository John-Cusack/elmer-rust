import json, sys
W='/home/john/firstmate/data/elmer-swarm-l15/refs/'
NN={808:8,706:6,605:5,504:4,404:4,303:3}
def total(code,p,sd):
    import math
    edges={808:12,706:9,605:8,504:6,404:4,303:3}[code]
    if code==808: faces=[6]; 
    nodes=NN[code]
    e=edges*(p-1)
    tri=(p-1)*(p-2)//2
    quad=(p-2)*(p-3)//2 if sd else (p-1)**2
    if code in (404,303): f=b=0
    if code==808: f=6*quad; b=((p-3)*(p-4)*(p-5)//6 if sd else (p-1)**3)
    if code==706: f=2*tri+3*quad; b=((p-2)*(p-3)*(p-4)//6 if sd else (p-1)**2*(p-2)//2)
    if code==605: f=(p-1)**2+4*tri; b=(p-2)*(p-1)*p//6
    if code==504: f=4*tri; b=(p-3)*(p-2)*(p-1)//6
    if code==404: return nodes+e+max(0,(p-2)*(p-3)//2 if sd else (p-1)**2)
    if code==303: return nodes+e+max(0,(p-2)*(p-1)//2)
    return nodes+e+f+max(0,b)
def mk(name, out, maxpts_per_elem=3):
    seen=set(); lines=[]; per={}
    for l in open(W+name+'.jsonl'):
        d=json.loads(l)
        if 'error' in d or d['code'] not in NN: continue
        if name.startswith('h_') and d['code'] in (404,303): continue
        key=(d['code'],d['u'],d['v'],d['w'])
        if key in seen: continue
        seen.add(key)
        import re
        pi=d['pinfo']
        p=int(re.search(r'_p(\d)',name).group(1)); sd='_sd' in name
        pi={'tetratype':pi.get('tetratype',0)}
        n=total(d['code'],p,sd); nn=NN[d['code']]
        per[d['code']]=per.get(d['code'],0)+1
        if per[d['code']]%max(1,STEP)!=1 and STEP>1: continue
        b=d['basis'][0][:n]; g=d['dbasisdx']
        toks=[d['code'],nn,n,pi['tetratype'],p,int(sd)]
        for i in range(nn):
            toks+= [d['nodeindexes'][i], d['xyz'][0][i], d['xyz'][1][i], d['xyz'][2][i]]
        toks+=[d['u'],d['v'],d['w'],d['detj']]
        toks+=b
        for i in range(n): toks+=[g[0][i],g[1][i],g[2][i]]
        lines.append(' '.join(repr(t) if isinstance(t,float) else str(t) for t in toks))
    open(out,'w').write('\n'.join(lines)+'\n')
    print(name,len(lines))
STEP=int(sys.argv[1]) if len(sys.argv)>1 else 1
import os
os.makedirs('/tmp/fx',exist_ok=True)
for n in ['h_p2_id','h_p2_sh','h_p3_id','h_p3_sh','h_p4_sh','h_p3_sd_sh','h2d_p2_sh','h2d_p3_sh','h2d_p4_sh','h2d_p5_sh','h2d_p4_sd_sh','h2d_p3_id']:
    mk(n,'/tmp/fx/%s.txt'%n)
