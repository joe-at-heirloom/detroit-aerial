"""Second RBF pass.

The first warp leaves ~15 m median but with cells past 100 m, because the control
stations were solved against an unwarped mosaic where the search had to be wide
and was therefore ambiguous. Re-solving on the warped result lets the search be
tight (+/-90 m), which is unambiguous, so the stations are far cleaner and a
second, shorter-length RBF can take out what the first pass missed."""
import sys, json, math, os, time, numpy as np
sys.path.insert(0,'pipeline')
import absorient, dtmap, stations, rbfwarp, rbfapply, warpapply
from PIL import Image
import rasterio
Image.MAX_IMAGE_PIXELS=None

MODB=json.load(open('data/modern_west_geo.json'))['bbox']
MOD=np.asarray(Image.open('mosaics/modern_west.png').convert('L'))
def modern_for(bbox,W,H):
    ms,mw,mn,me=MODB; mh,mwd=MOD.shape
    s,w,n,e=bbox
    y0=max(0,int((mn-n)/(mn-ms)*mh)); y1=min(mh,int((mn-s)/(mn-ms)*mh))
    x0=max(0,int((w-mw)/(me-mw)*mwd)); x1=min(mwd,int((e-mw)/(me-mw)*mwd))
    if y1-y0<10 or x1-x0<10: return None
    return np.asarray(Image.fromarray(MOD[y0:y1,x0:x1]).resize((W,H),Image.LANCZOS))

def trend_at(t,cE,cN):
    a=warpapply.design(np.array([(cE-t['mx'])/t['sx']]),np.array([(cN-t['my'])/t['sy']]),t['deg'])
    return np.array([float((a@np.array(c))[0]) for c in t['C']])
def poly_trend(st,deg=1):
    P=np.array([[s['cE'],s['cN']] for s in st],float)
    V=np.array([[s['dE'],s['dN'],0.0] for s in st],float)
    W=np.array([min(s.get('ratio',1),4.0) for s in st],float)
    mx,my=P[:,0].mean(),P[:,1].mean(); sx,sy=max(P[:,0].std(),1),max(P[:,1].std(),1)
    A=warpapply.design((P[:,0]-mx)/sx,(P[:,1]-my)/sy,deg); sw=np.sqrt(W)[:,None]
    C=[np.linalg.lstsq(A*sw,V[:,k]*np.sqrt(W),rcond=None)[0] for k in range(3)]
    return dict(deg=deg,mx=mx,my=my,sx=sx,sy=sy,C=[c.tolist() for c in C])

def grid_resid(arr,mod,MPP,NY=16,NX=3,search=200.0):
    Hp,Wp=arr.shape; ch=Hp//NY; cw=Wp//NX; vals=[]
    for j in range(NY):
        for i in range(NX):
            a=arr[j*ch:(j+1)*ch, i*cw:(i+1)*cw]; b=mod[j*ch:(j+1)*ch, i*cw:(i+1)*cw]
            if (a>0).mean()<0.45: continue
            ra,va=absorient.ridge(a.astype(np.float32),MPP); rb,vb=absorient.ridge(b.astype(np.float32),MPP)
            A=ra*va; B=rb*vb
            if A.std()<1e-6 or B.std()<1e-6: continue
            A=A-A[va].mean(); A[~va]=0; B=B-B[vb].mean(); B[~vb]=0
            C=np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(A)*np.conj(np.fft.rfft2(B)),s=A.shape))
            lim=int(search/MPP); cy,cx=A.shape[0]//2,A.shape[1]//2
            sb=C[max(0,cy-lim):cy+lim,max(0,cx-lim):cx+lim]
            pk=np.unravel_index(np.argmax(sb),sb.shape)
            vals.append(math.hypot((pk[1]-lim)*MPP,(pk[0]-lim)*MPP))
    return np.array(vals)

for tag in sys.argv[1:]:
    t0=time.time()
    src=f'mosaics/detroit_{tag}_rbf.tif'
    geo=json.load(open(f'data/{tag}_rbf_geo.json'))
    bbox=geo['bbox']
    ds=rasterio.open(src)
    MPP=2.5
    Wp=int((bbox[3]-bbox[1])*dtmap.MLON/MPP); Hp=int((bbox[2]-bbox[0])*dtmap.MLAT/MPP)
    arr=ds.read(1,out_shape=(Hp,Wp))
    mod=modern_for(bbox,Wp,Hp)
    before=grid_resid(arr,mod,MPP)
    print(f"=== {tag} pass 2 ===",flush=True)
    print(f"  before: median {np.median(before):.1f}  p90 {np.percentile(before,90):.1f}  "
          f"max {before.max():.1f}  >25m {(before>25).sum()}/{len(before)}",flush=True)
    st=stations.solve_grid(arr,mod,bbox,MPP,dtmap.MLAT,dtmap.MLON,
                           win_m=1800.0, overlap=0.4, search_m=90.0, min_valid=0.45)
    good=[s for s in st if s['ratio']>=1.08]
    print(f"  stations {len(st)} -> {len(good)} confident (tight +/-90 m search)",flush=True)
    if len(good)<25:
        print("  too few stations, skipping"); continue
    tr=poly_trend(good,1)
    res=[dict(cE=s['cE'],cN=s['cN'],
              dE=s['dE']-trend_at(tr,s['cE'],s['cN'])[0],
              dN=s['dN']-trend_at(tr,s['cE'],s['cN'])[1],
              rot=0.0,ratio=s['ratio']) for s in good]
    R=rbfwarp.RBFWarp(res,length=800.0,trim=True,min_ratio=0.0)
    def warp_at(cE,cN,tr=tr,R=R):
        a=trend_at(tr,cE,cN); d=R.at(cE,cN)
        return (a[0]+d[0], a[1]+d[1], 0.0)
    g2=dict(minE=geo['minE'],maxN=geo['maxN'],W=geo['W'],H=geo['H'],mpp=geo['mpp'])
    out=rbfapply.apply(src,g2,warp_at,f'mosaics/detroit_{tag}_p2.tif',f'/tmp/p2_{tag}.raw',log=lambda s:None)
    json.dump(out,open(f'data/{tag}_p2_geo.json','w'))
    json.dump(dict(kind='poly1_rbf',trend=tr,length=800.0),open(f'data/rbffit2_{tag}.json','w'))
    R.dump(f'data/rbfres2_{tag}.json')
    ds2=rasterio.open(f'mosaics/detroit_{tag}_p2.tif')
    arr2=ds2.read(1,out_shape=(Hp,Wp))
    after=grid_resid(arr2,mod,MPP)
    print(f"  after : median {np.median(after):.1f}  p90 {np.percentile(after,90):.1f}  "
          f"max {after.max():.1f}  >25m {(after>25).sum()}/{len(after)}   [{time.time()-t0:.0f}s]",flush=True)
    del arr, arr2, mod
print("REFINE DONE",flush=True)
