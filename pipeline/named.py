"""Validate alignment by street IDENTITY, not just geometry.
For each named street we measure (a) the across-street offset and (b) the width of
the bright corridor found there. If a mosaic were locked onto the wrong parallel
street, an arterial's centreline would sit over a narrow residential street and the
measured width would collapse."""
import numpy as np, math, json, dtmap
import os as _os
_DATA = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'data')
def _d(p):
    return p if _os.path.isabs(p) or _os.path.exists(p) else _os.path.join(_DATA, _os.path.basename(p))


def load_named(path='named_centerlines.json'):
    raw=json.load(open(_d(path)))
    by={}
    for nm,la1,lo1,la2,lo2 in raw:
        by.setdefault(nm,[]).append(((lo1-dtmap.LON0)*dtmap.MLON,(la1-dtmap.LAT0)*dtmap.MLAT,
                                     (lo2-dtmap.LON0)*dtmap.MLON,(la2-dtmap.LAT0)*dtmap.MLAT))
    return {k:np.array(v,dtype=np.float32) for k,v in by.items()}

def profile_stats(blk,S,rng=40.0,step=1.0,min_len=50,max_n=4000):
    """returns offset list and corridor-width list for one named street"""
    d=np.arange(-rng,rng+step,step)
    offs=[];wids=[];nx=[];ny=[]
    I=blk.I; H,W=I.shape
    for k in range(len(S)):
        x0,y0,x1,y1=S[k]
        L=math.hypot(x1-x0,y1-y0)
        if L<min_len: continue
        for t in (np.arange(max(1,int(L//40)))+0.5)/max(1,int(L//40)):
            mx=x0+t*(x1-x0); my=y0+t*(y1-y0)
            px,py=blk.px(mx,my)
            if not (60<px<W-60 and 60<py<H-60): continue
            ux,uy=(x1-x0)/L,(y1-y0)/L; nxu,nyu=-uy,ux
            sx=(px+nxu*d/blk.mppx).astype(np.int32); sy=(py-nyu*d/blk.mppy).astype(np.int32)
            if sx.min()<0 or sy.min()<0 or sx.max()>=W or sy.max()>=H: continue
            v=I[sy,sx].astype(np.float32)
            if (v==0).any(): continue
            base=np.percentile(v,15); pk=v.max()
            if pk-base<18: continue
            half=base+0.5*(pk-base)
            above=v>=half
            # contiguous run containing the profile centre-most peak
            i=int(np.argmax(v)); a=i
            while a>0 and above[a-1]: a-=1
            b=i
            while b<len(v)-1 and above[b+1]: b+=1
            width=(b-a+1)*step
            if width>2*rng*0.8: continue
            centre=(d[a]+d[b])/2
            offs.append(centre); wids.append(width); nx.append(nxu); ny.append(nyu)
            if len(offs)>=max_n: break
        if len(offs)>=max_n: break
    return np.array(offs),np.array(wids),np.array(nx),np.array(ny)
