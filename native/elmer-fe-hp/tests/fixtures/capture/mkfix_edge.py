import json, sys, os
W='/home/john/firstmate/data/elmer-swarm-l15/refs/'
NN={808:8,706:6,605:5,504:4,404:4,303:3}
NF={0:{303:3,404:4,504:6,605:8,706:9,808:12},1:{303:3,404:6,504:6,605:10,706:15,808:27},2:{303:8,404:12,504:20,605:31,706:36,808:54}}
MODE={'edge':0,'edge_piola':1,'edge_quadratic':2}
STEP=int(sys.argv[1]) if len(sys.argv)>1 else 1
os.makedirs('/tmp/fx',exist_ok=True)
for src in MODE:
    for perm in ('id','sh'):
        name='e_%s_%s'%(src,perm); mode=MODE[src]
        seen=set(); lines=[]; per={}
        for l in open(W+name+'.jsonl'):
            d=json.loads(l)
            if 'error' in d: continue
            key=(d['code'],d['u'],d['v'],d['w'])
            if key in seen: continue
            seen.add(key)
            per[d['code']]=per.get(d['code'],0)+1
            if STEP>1 and per[d['code']]%STEP!=1: continue
            nn=NN[d['code']]; nf=NF[mode][d['code']]
            toks=[d['code'],nn,nf,mode]
            for i in range(nn): toks+=[d['nodeindexes'][i],d['xyz'][0][i],d['xyz'][1][i],d['xyz'][2][i]]
            toks+=[d['u'],d['v'],d['w'],d['detj']]
            e=d['edgebasis']; c=d['rotbasis']
            for i in range(nf): toks+=[e[0][i],e[1][i],e[2][i]]
            for i in range(nf): toks+=[c[0][i],c[1][i],c[2][i]]
            lines.append(' '.join(repr(t) if isinstance(t,float) else str(t) for t in toks))
        open('/tmp/fx/%s.txt'%name,'w').write('\n'.join(lines)+'\n')
        print(name,len(lines),per)
