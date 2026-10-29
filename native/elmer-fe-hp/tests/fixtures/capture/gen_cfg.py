import os, re, random, shutil, sys
W='/home/john/firstmate/data/elmer-swarm-l15/refs/'
def make(name, src, p, perm, sd=False, extra=""):
    d=W+name
    if os.path.exists(d): shutil.rmtree(d)
    shutil.copytree(W+src, d)
    if p is not None:
        sif=open(d+'/case.sif').read()
        sif=re.sub(r'(?m)^  Element = "[^"]*"\n','',sif)
        sif=sif.replace('  Serendipity P Elements = False\n','').replace('  Serendipity P Elements = True\n','')
        sif=sif.replace('  Steady State Max Iterations = 1\n','  Steady State Max Iterations = 1\n  Serendipity P Elements = %s\n'%('True' if sd else 'False'),1)
        sif=sif.replace('  Procedure = "SingleEdgeElements" "SingleEdgeElements"\n','  Procedure = "SingleEdgeElements" "SingleEdgeElements"\n  Element = "p:%d"\n'%p,1)
        open(d+'/case.sif','w').write(sif)
    # renumber nodes per element block
    el=[l.split() for l in open(d+'/unitelems/mesh.elements') if l.strip()]
    mp={}
    rng=random.Random(12345)
    for e in el:
        ids=[int(x) for x in e[3:]]
        if perm=='identity': new=ids[:]
        elif perm=='reversed': new=ids[::-1]
        else:
            new=ids[:]; rng.shuffle(new)
        for a,b in zip(ids,new): mp[a]=b
    nodes={int(l.split()[0]):l.split()[1:] for l in open(d+'/unitelems/mesh.nodes') if l.strip()}
    with open(d+'/unitelems/mesh.nodes','w') as f:
        newn={mp[k]:v for k,v in nodes.items()}
        for k in sorted(newn): f.write('%d %s\n'%(k,' '.join(newn[k])))
    with open(d+'/unitelems/mesh.elements','w') as f:
        for e in el: f.write(' '.join(e[:3]+[str(mp[int(x)]) for x in e[3:]])+' \n')
    bl=[l.split() for l in open(d+'/unitelems/mesh.boundary') if l.strip()]
    with open(d+'/unitelems/mesh.boundary','w') as f:
        for e in bl: f.write(' '.join(e[:5]+[str(mp[int(x)]) for x in e[5:]])+'\n')
if __name__=='__main__':
    for name,p,perm,sd in [('h_p2_id',2,'identity',False),('h_p2_sh',2,'shuffle',False),('h_p3_id',3,'identity',False),('h_p3_sh',3,'shuffle',False),('h_p4_sh',4,'shuffle',False),('h_p3_sd_sh',3,'shuffle',True)]:
        make(name,'ph3',p,perm,sd)
