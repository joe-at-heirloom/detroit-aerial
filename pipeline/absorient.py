"""ABSOLUTE ORIENTATION via the mile grid.

Detroit's arterials lie on a ~1609 m lattice. Correlating against them therefore
cannot alias inside a +/-400 m search, unlike the ~100 m residential grid that
defeated every earlier attempt. Signal is a wide-road ridge response, not raw
brightness, so tree cover and film tone matter far less."""
import numpy as np, math, json
from PIL import Image, ImageDraw
from skimage.filters import sato
from scipy.ndimage import uniform_filter
import named, dtmap
import os as _os
_DATA = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'data')
def _d(p):
    return p if _os.path.isabs(p) or _os.path.exists(p) else _os.path.join(_DATA, _os.path.basename(p))


ART = ['Grand River','Michigan','Livernois','Wyoming','Warren','Joy','McNichols','Davison',
       'Outer','Fort','Vernor','Springwells','Tireman','Plymouth','Chicago','Puritan',
       'Fenkell','Schaefer','Greenfield','Evergreen','7 Mile','8 Mile','Schoolcraft',
       'Jefferson','Gratiot','Woodward','Grand','Dexter','Linwood','Hamilton','Oakman',
       'Southfield','Telegraph','Ford','Michigan Ave','West Chicago','Lyndon','Curtis',
       'Six Mile','Seven Mile','Eight Mile','Meyers','Hubbell','Trinity','Burt','Lahser']

def arterial_segments(use_osm=True):
    """City named arterials (accurate, Detroit only) + OSM arterial classes
    (covers the suburbs where the city file stops). Both are on the mile grid."""
    BY = named.load_named()
    out=[BY[nm] for nm in ART if nm in BY]
    if use_osm:
        import glob
        KEEP={'motorway','trunk','primary','secondary'}
        segs=[]
        for p in sorted(glob.glob(_os.path.join(_DATA, 'osmtiles/r_*.json'))):
            for w in json.load(open(p)).get('elements',[]):
                if w.get('tags',{}).get('highway') not in KEEP: continue
                g=w.get('geometry')
                if not g or len(g)<2: continue
                pts=[dtmap.to_m(q['lat'],q['lon']) for q in g]
                for a,b in zip(pts,pts[1:]): segs.append((a[0],a[1],b[0],b[1]))
        if segs: out.append(np.array(segs,dtype=np.float32))
    return np.vstack(out) if out else np.zeros((0,4),np.float32)

def bearing_diversity(S, cE, cN, rad):
    """smallest singular value ratio of the segment-normal set - low means every
    arterial in view is parallel, so the across-street direction is unconstrained"""
    m=(((S[:,0]+S[:,2])/2-cE)**2+(((S[:,1]+S[:,3])/2)-cN)**2) < rad*rad
    if m.sum()<12: return 0.0, int(m.sum())
    v=S[m]; dx=v[:,2]-v[:,0]; dy=v[:,3]-v[:,1]
    L=np.hypot(dx,dy); k=L>1
    if k.sum()<12: return 0.0, int(k.sum())
    n=np.column_stack([-dy[k]/L[k], dx[k]/L[k]])
    n=n*np.sign(n[:,:1]+1e-9)
    _,sv,_=np.linalg.svd(n, full_matrices=False)
    return float(sv[-1]/sv[0]), int(k.sum())

def ridge(I, mppx, road_w=(18.0,46.0), flat_m=320.0):
    valid = I > 0
    J = I.astype(np.float32)
    J = J - uniform_filter(J, size=int(max(5, flat_m/mppx)))
    J[~valid] = 0
    sig = np.linspace(max(1.0, road_w[0]/(2*mppx)), road_w[1]/(2*mppx), 5)
    r = sato(J, sigmas=sig, black_ridges=False)
    r[~valid] = 0
    if r.max() > 0: r = r/r.max()
    return r, valid

def mask_from(S, W, H, minE, maxN, mppx, mppy, dE=0.0, dN=0.0, width=3):
    m = Image.new('L', (W, H), 0); d = ImageDraw.Draw(m)
    X0=(S[:,0]+dE-minE)/mppx; Y0=(maxN-(S[:,1]+dN))/mppy
    X1=(S[:,2]+dE-minE)/mppx; Y1=(maxN-(S[:,3]+dN))/mppy
    k=~(((X0<0)&(X1<0))|((X0>W)&(X1>W))|((Y0<0)&(Y1<0))|((Y0>H)&(Y1>H)))
    for a,b,c,e in zip(X0[k],Y0[k],X1[k],Y1[k]): d.line([(a,b),(c,e)], fill=255, width=width)
    return np.asarray(m, dtype=np.float32)/255.0

def solve(r, valid, S, minE, maxN, mppx, mppy, search_m=400.0):
    H,W = r.shape
    A = r.copy()
    A[~valid]=0
    A = A - A[valid].mean(); A[~valid]=0
    M = mask_from(S, W, H, minE, maxN, mppx, mppy)
    cov = M.mean()
    if M.sum() < 200: return None
    Mz = M - M.mean()
    C = np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(A)*np.conj(np.fft.rfft2(Mz)), s=A.shape))
    lx = int(search_m/mppx); ly = int(search_m/mppy)
    cy,cx = H//2, W//2
    sub = C[max(0,cy-ly):cy+ly, max(0,cx-lx):cx+lx]
    pk = np.unravel_index(np.argmax(sub), sub.shape)
    peak = sub[pk]
    mm = sub.copy()
    # exclusion must stay well inside the search window or "second peak" vanishes
    exm = min(180.0, search_m*0.40)
    ey = max(1,int(exm/mppy)); ex = max(1,int(exm/mppx))
    mm[max(0,pk[0]-ey):pk[0]+ey, max(0,pk[1]-ex):pk[1]+ex] = -9e9
    second = mm.max()
    dE = (pk[1]-lx)*mppx; dN = -(pk[0]-ly)*mppy
    return dict(dE=float(dE), dN=float(dN), ratio=float(peak/max(second,1e-9)),
                sigma=float((peak-sub.mean())/(sub.std()+1e-9)), coverage=float(cov))

def solve_block(tag, S, mpp=2.5, half=2400.0, nlat=11, search=400.0):
    """Stage A: wide search at many stations -> robust block-level consensus.
       Stage B: re-solve each station constrained to +/-130 m of that consensus."""
    import rasterio, json as _j
    g=_j.load(open(f'{tag}_full_geo.json')); s,w,n,e=g['bbox']
    ds=rasterio.open(f'detroit_{tag}_0.63m.tif')
    lats=np.linspace(s+half/dtmap.MLAT*1.02, n-half/dtmap.MLAT*1.02, nlat)
    lon=(w+e)/2
    def one(la, srch, ctr=(0.0,0.0)):
        S_,W_,N_,E_=la-half/dtmap.MLAT,lon-half/dtmap.MLON,la+half/dtmap.MLAT,lon+half/dtmap.MLON
        if not (W_>w and E_<e and S_>s and N_<n): return None
        cE=(lon-dtmap.LON0)*dtmap.MLON; cN=(la-dtmap.LAT0)*dtmap.MLAT
        div,nseg=bearing_diversity(S,cE,cN,half)
        if div<0.25: return None
        Wp=int((E_-W_)*dtmap.MLON/mpp); Hp=int((N_-S_)*dtmap.MLAT/mpp)
        win=rasterio.windows.Window((W_-w)/(e-w)*g['W'],(n-N_)/(n-s)*g['H'],
                                    (E_-W_)/(e-w)*g['W'],(N_-S_)/(n-s)*g['H'])
        I=ds.read(1,window=win,out_shape=(Hp,Wp)).astype(np.float32)
        if (I>0).mean()<0.45: return None
        minE=(W_-dtmap.LON0)*dtmap.MLON; maxN=(N_-dtmap.LAT0)*dtmap.MLAT
        r,valid=ridge(I,mpp)
        # shift the reference to the prior centre, search a small window around it
        Sc=S.copy(); Sc[:,[0,2]]+=ctr[0]; Sc[:,[1,3]]+=ctr[1]
        res=solve(r,valid,Sc,minE,maxN,mpp,mpp,search_m=srch)
        if res:
            res['dE']+=ctr[0]; res['dN']+=ctr[1]
            res.update(lat=float(la),div=div,nseg=nseg)
        return res
    A=[one(la,search) for la in lats]
    A=[x for x in A if x]
    if len(A)<3: return None,None
    P=np.array([[x['dE'],x['dN']] for x in A]); Wt=np.array([x['ratio'] for x in A])
    c=np.median(P,axis=0)
    for _ in range(6):
        d=np.hypot(P[:,0]-c[0],P[:,1]-c[1])
        sg=max(np.median(d)*1.4826,25.0)
        wt=Wt/(1+(d/(1.5*sg))**2)
        c=(P*wt[:,None]).sum(0)/wt.sum()
    B=[one(x['lat'],130.0,tuple(c)) for x in A]
    B=[x for x in B if x]
    return c, B

def solve_similarity(r, valid, S, minE, maxN, mppx, mppy, cE, cN,
                     search_m=400.0, rots=None, scales=(1.0,), width=3):
    """Full similarity absolute orientation: rotate/scale the arterial reference
    about the window centre, then FFT-correlate for translation. A rotated
    reference smears the correlation peak, which is why translation-only solves
    produced broad, ambiguous maxima."""
    H,W = r.shape
    A = r.copy(); A[~valid]=0
    A = A - A[valid].mean(); A[~valid]=0
    FA = np.fft.rfft2(A)
    if rots is None: rots = np.arange(-4.0, 4.01, 0.5)
    best=None
    for sc in scales:
        for rot in rots:
            th=math.radians(rot); c,s = math.cos(th), math.sin(th)
            dx0=S[:,0]-cE; dy0=S[:,1]-cN; dx1=S[:,2]-cE; dy1=S[:,3]-cN
            T=np.empty_like(S)
            T[:,0]=cE+sc*( c*dx0 - s*dy0); T[:,1]=cN+sc*( s*dx0 + c*dy0)
            T[:,2]=cE+sc*( c*dx1 - s*dy1); T[:,3]=cN+sc*( s*dx1 + c*dy1)
            M=mask_from(T,W,H,minE,maxN,mppx,mppy,width=width)
            if M.sum()<200: continue
            Mz=M-M.mean()
            C=np.fft.fftshift(np.fft.irfft2(FA*np.conj(np.fft.rfft2(Mz)), s=A.shape))
            lx=int(search_m/mppx); ly=int(search_m/mppy)
            cy,cx=H//2,W//2
            sub=C[max(0,cy-ly):cy+ly, max(0,cx-lx):cx+lx]
            pk=np.unravel_index(np.argmax(sub),sub.shape); peak=sub[pk]
            mm=sub.copy()
            exm=min(180.0,search_m*0.40)
            ey=max(1,int(exm/mppy)); ex=max(1,int(exm/mppx))
            mm[max(0,pk[0]-ey):pk[0]+ey, max(0,pk[1]-ex):pk[1]+ex]=-9e9
            second=mm.max()
            score=peak/math.sqrt(M.sum())
            if best is None or score>best['score']:
                best=dict(score=float(score), rot=float(rot), scale=float(sc),
                          dE=float((pk[1]-lx)*mppx), dN=float(-(pk[0]-ly)*mppy),
                          ratio=float(peak/max(second,1e-9)),
                          sigma=float((peak-sub.mean())/(sub.std()+1e-9)))
    return best

def solve_vs_image(rH, vH, rM, vM, mppx, mppy, search_m=380.0, rots=None):
    """Correlate the historical ridge response against the MODERN ridge response.
    Modern imagery carries continuous road structure (not 1.8% vector coverage),
    which is far more signal per station."""
    import numpy as _np
    H,W = rH.shape
    A=rH.copy(); A[~vH]=0
    A=A-A[vH].mean(); A[~vH]=0
    FA=_np.fft.rfft2(A)
    if rots is None: rots=_np.arange(-3.0,3.01,0.5)
    best=None
    for rot in rots:
        B = rM if abs(rot)<1e-9 else _rotate(rM, rot)
        vB = vM if abs(rot)<1e-9 else (_rotate(vM.astype(_np.float32), rot) > 0.5)
        if vB.mean()<0.25: continue
        Bz=B.copy(); Bz[~vB]=0
        Bz=Bz-Bz[vB].mean(); Bz[~vB]=0
        C=_np.fft.fftshift(_np.fft.irfft2(FA*_np.conj(_np.fft.rfft2(Bz)), s=A.shape))
        lx=int(search_m/mppx); ly=int(search_m/mppy)
        cy,cx=H//2,W//2
        sub=C[max(0,cy-ly):cy+ly, max(0,cx-lx):cx+lx]
        pk=_np.unravel_index(_np.argmax(sub),sub.shape); peak=sub[pk]
        mm=sub.copy()
        ey=max(1,int(150/mppy)); ex=max(1,int(150/mppx))
        mm[max(0,pk[0]-ey):pk[0]+ey, max(0,pk[1]-ex):pk[1]+ex]=-9e9
        second=mm.max()
        if best is None or peak>best['_p']:
            best=dict(_p=float(peak), rot=float(rot),
                      dE=float((pk[1]-lx)*mppx), dN=float(-(pk[0]-ly)*mppy),
                      ratio=float(peak/max(second,1e-9)),
                      sigma=float((peak-sub.mean())/(sub.std()+1e-9)))
    if best: best.pop('_p',None)
    return best

def _rotate(I, deg):
    import numpy as _np, math as _m
    h,w=I.shape; th=_m.radians(deg); c,s=_m.cos(th),_m.sin(th)
    Y,X=_np.mgrid[0:h,0:w]; dy=Y-h/2.0; dx=X-w/2.0
    sy=(-s*dx+c*dy)+h/2.0; sx=( c*dx+s*dy)+w/2.0
    yi=_np.clip(sy,0,h-1).astype(_np.int32); xi=_np.clip(sx,0,w-1).astype(_np.int32)
    o=I[yi,xi].astype(_np.float32)
    o[(sy<0)|(sy>=h)|(sx<0)|(sx>=w)]=0
    return o
