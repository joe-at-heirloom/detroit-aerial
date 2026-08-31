"""Rubber-sheet a historical mosaic onto the modern road network.

Each sample gives ONE constraint: the offset along the road's normal (the
across-road direction is all a single road can tell you).  Roads of differing
bearing inside a neighbourhood jointly constrain both components, so we solve a
local 2-D translation per grid cell, then interpolate a smooth displacement field."""
import json, glob, math, numpy as np, rasterio, dtmap
from PIL import Image
import os as _os
_DATA = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'data')
def _d(p):
    return p if _os.path.isabs(p) or _os.path.exists(p) else _os.path.join(_DATA, _os.path.basename(p))

Image.MAX_IMAGE_PIXELS=None

MAIN={'motorway','trunk','primary','secondary','tertiary'}
def load_roads(classes=MAIN):
    segs=[]
    for p in sorted(glob.glob(_os.path.join(_DATA, 'osmtiles/r_*.json'))):
        for w in json.load(open(p)).get('elements',[]):
            if w.get('tags',{}).get('highway') not in classes: continue
            g=w.get('geometry')
            if not g or len(g)<2: continue
            pts=[dtmap.to_m(q['lat'],q['lon']) for q in g]
            for a,b in zip(pts,pts[1:]): segs.append((a[0],a[1],b[0],b[1]))
    return np.array(segs,dtype=np.float32)

class Block:
    def __init__(self,tag,ds=2):
        g=json.load(open(f'{tag}_full_geo.json')); self.g=g
        s,w,n,e=g['bbox']
        self.W=g['W']//ds; self.H=g['H']//ds
        self.mppx=(e-w)*dtmap.MLON/self.W
        self.mppy=(n-s)*dtmap.MLAT/self.H
        self.minE=(w-dtmap.LON0)*dtmap.MLON
        self.maxN=(n-dtmap.LAT0)*dtmap.MLAT
        src=rasterio.open(f'detroit_{tag}_0.63m.tif')
        self.I=src.read(1,out_shape=(self.H,self.W))
        self.tag=tag
    def px(self,ex,ey): return (ex-self.minE)/self.mppx, (self.maxN-ey)/self.mppy

def sample_offsets(blk,S,rng=26.0,step=1.0,min_len=60,stride=1,max_n=200000):
    """Returns arrays: x,y (block px), nx,ny (unit normal in world), off (metres)."""
    d=np.arange(-rng,rng+step,step)
    X=[];Y=[];NX=[];NY=[];OFF=[]
    I=blk.I; H,W=I.shape
    idx=range(0,len(S),stride)
    for k in idx:
        x0,y0,x1,y1=S[k]
        L=math.hypot(x1-x0,y1-y0)
        if L<min_len: continue
        nseg=max(1,int(L//45))
        for t in (np.arange(nseg)+0.5)/nseg:
            mx=x0+t*(x1-x0); my=y0+t*(y1-y0)
            px,py=blk.px(mx,my)
            if not (40<px<W-40 and 40<py<H-40): continue
            ux,uy=(x1-x0)/L,(y1-y0)/L
            nx,ny=-uy,ux
            sx=px+nx*d/blk.mppx; sy=py-ny*d/blk.mppy
            xi=sx.astype(np.int32); yi=sy.astype(np.int32)
            if xi.min()<0 or yi.min()<0 or xi.max()>=W or yi.max()>=H: continue
            v=I[yi,xi].astype(np.float32)
            if (v==0).any(): continue
            base=np.percentile(v,20); wgt=np.clip(v-base,0,None)
            if wgt.sum()<1e-6 or (v.max()-base)<14: continue
            c=(wgt*d).sum()/wgt.sum()
            sp=math.sqrt((wgt*(d-c)**2).sum()/wgt.sum())
            if sp>15: continue
            X.append(px);Y.append(py);NX.append(nx);NY.append(ny);OFF.append(c)
            if len(OFF)>=max_n: break
        if len(OFF)>=max_n: break
    return (np.array(X),np.array(Y),np.array(NX),np.array(NY),np.array(OFF))

def solve_cells(X,Y,NX,NY,OFF,W,H,ncx,ncy,min_pts=25,min_cond=0.12):
    """Local 2-D translation per cell. Returns (cx,cy,dE,dN,n,ok) arrays."""
    cw=W/ncx; ch=H/ncy
    out=[]
    for j in range(ncy):
        for i in range(ncx):
            m=(X>=i*cw)&(X<(i+1)*cw)&(Y>=j*ch)&(Y<(j+1)*ch)
            n=int(m.sum())
            if n<min_pts: out.append((( i+.5)*cw,(j+.5)*ch,0,0,n,False)); continue
            A=np.column_stack([NX[m],NY[m]]); b=OFF[m]
            # need bearing diversity, else the normal direction is degenerate
            u,sv,vt=np.linalg.svd(A,full_matrices=False)
            if sv[-1]/sv[0] < min_cond: out.append((((i+.5)*cw,(j+.5)*ch,0,0,n,False))); continue
            d=np.linalg.lstsq(A,b,rcond=None)[0]
            for _ in range(3):
                r=b-A@d; s=max(np.median(np.abs(r))*1.4826,1.0)
                w=1/(1+(r/(2.5*s))**2)
                d=np.linalg.lstsq(A*np.sqrt(w)[:,None],b*np.sqrt(w),rcond=None)[0]
            out.append(((i+.5)*cw,(j+.5)*ch,float(d[0]),float(d[1]),n,True))
    a=np.array([(o[0],o[1],o[2],o[3],o[4],1.0 if o[5] else 0.0) for o in out])
    return a
