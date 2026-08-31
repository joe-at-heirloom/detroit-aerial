"""Anchor using NAMED ARTERIALS only.
Detroit's arterials lie on the mile grid (~1.6 km apart), so matching to the wrong
one would require a ~1.6 km error - unlike the residential grid, whose ~100 m
period made every correlation method ambiguous. Each arterial yields one 1-D
constraint (offset along its normal); arterials of differing bearing in a chunk
jointly pin a 2-D shift."""
import numpy as np, math

ART = ['Grand River','Michigan','Livernois','Wyoming','Warren','Joy','McNichols','Davison',
       'Outer','Fort','Vernor','Springwells','Tireman','Plymouth','Chicago','Puritan',
       'Fenkell','Schaefer','Greenfield','Evergreen','7 Mile','8 Mile','Schoolcraft',
       'Jefferson','Gratiot','Woodward','Grand','Dexter','Linwood','Hamilton','Oakman']

def profile(blk, S, y0, y1, rng=320.0, step=2.0, min_n=10):
    """Averaged across-street brightness for segments falling in rows [y0,y1)."""
    d=np.arange(-rng,rng+step,step)
    P=[];NX=[];NY=[]
    I=blk.I; H,W=I.shape
    for x0,yy0,x1,yy1 in S:
        L=math.hypot(x1-x0,yy1-yy0)
        if L<40: continue
        n=max(1,int(L//35))
        for t in (np.arange(n)+0.5)/n:
            mx=x0+t*(x1-x0); my=yy0+t*(yy1-yy0)
            px,py=blk.px(mx,my)
            if not (y0<=py<y1): continue
            ux,uy=(x1-x0)/L,(yy1-yy0)/L
            sx=px-uy*d/blk.mppx; sy=py-ux*d/blk.mppy
            if sx.min()<0 or sy.min()<0 or sx.max()>=W-1 or sy.max()>=H-1: continue
            v=I[sy.astype(np.int32),sx.astype(np.int32)].astype(np.float32)
            if (v==0).any(): continue
            P.append(v); NX.append(-uy); NY.append(ux)
    if len(P)<min_n: return None
    m=np.array(P).mean(0)
    i=int(np.argmax(m))
    med=np.median(m); mad=np.median(np.abs(m-med))+1e-6
    prom=(m[i]-med)/mad
    # reject if the peak sits at the search edge (means the corridor is outside the window)
    if i<3 or i>len(m)-4: return None
    return dict(off=float(d[i]), prom=float(prom), n=len(P),
                nx=float(np.mean(NX)), ny=float(np.mean(NY)), peak=float(m[i]), med=float(med))

def solve_chunk(cons, min_cond=0.25):
    """cons: list of dicts with nx,ny,off,prom -> robust 2-D shift"""
    if len(cons)<2: return None
    A=np.array([[c['nx'],c['ny']] for c in cons])
    b=np.array([c['off'] for c in cons])
    w=np.array([min(c['prom'],8.0) for c in cons])
    u,sv,_=np.linalg.svd(A,full_matrices=False)
    if sv[-1]/sv[0] < min_cond: return None      # all arterials parallel -> underdetermined
    d=np.linalg.lstsq(A*np.sqrt(w)[:,None], b*np.sqrt(w), rcond=None)[0]
    for _ in range(3):
        r=b-A@d; s=max(np.median(np.abs(r))*1.4826, 8.0)
        w2=w/(1+(r/(2.0*s))**2)
        d=np.linalg.lstsq(A*np.sqrt(w2)[:,None], b*np.sqrt(w2), rcond=None)[0]
    resid=b-A@d
    return dict(dE=float(d[0]), dN=float(d[1]), n=len(cons),
                resid=float(np.median(np.abs(resid))), cond=float(sv[-1]/sv[0]))
