"""Self-calibrating mission geometry: derive px-per-metre and the flight frame's
rotation relative to north from the imagery + catalogue coords alone."""
import numpy as np, math, itertools
from PIL import Image

def mdeg(phi):
    p=math.radians(phi)
    return (111132.92-559.82*math.cos(2*p)+1.175*math.cos(4*p),
            111412.84*math.cos(p)-93.5*math.cos(3*p))

def pc_full(a,b):
    h=min(a.shape[0],b.shape[0]); w=min(a.shape[1],b.shape[1])
    ay=(a.shape[0]-h)//2; ax=(a.shape[1]-w)//2
    by=(b.shape[0]-h)//2; bx=(b.shape[1]-w)//2
    A=a[ay:ay+h,ax:ax+w]; B=b[by:by+h,bx:bx+w]
    A=A-A.mean(); B=B-B.mean()
    wn=np.outer(np.hanning(h),np.hanning(w))
    F=np.fft.rfft2(A*wn)*np.conj(np.fft.rfft2(B*wn))
    F/=(np.abs(F)+1e-9)
    c=np.fft.fftshift(np.fft.irfft2(F,s=(h,w)))
    pk=np.unravel_index(np.argmax(c),c.shape); peak=c[pk]
    m=c.copy(); m[max(0,pk[0]-8):pk[0]+8,max(0,pk[1]-8):pk[1]+8]=-9
    return pk[1]-w//2, pk[0]-h//2, float(peak/(c.std()+1e-12)), float(peak/(m.max()+1e-12))

def calibrate(frames, imgdir, max_sep=1700, min_sep=600, sharp=8, ratio=1.25, limit=400):
    """frames: list of {rec,lat,lon}.  Returns (px_per_m, rot_deg, n_used, diagnostics)."""
    lat0=np.median([f['lat'] for f in frames]); lon0=np.median([f['lon'] for f in frames])
    MLAT,MLON=mdeg(lat0)
    P={f['rec']:(((f['lon']-lon0)*MLON),((f['lat']-lat0)*MLAT)) for f in frames}
    recs=list(P)
    cand=[(a,b) for a,b in itertools.combinations(recs,2)
          if min_sep<math.hypot(P[a][0]-P[b][0],P[a][1]-P[b][1])<max_sep]
    cache={}
    def L(r):
        if r not in cache:
            if len(cache)>60: cache.clear()
            cache[r]=np.asarray(Image.open(f"{imgdir}/{r}.jpg").convert('L'),dtype=np.float32)
        return cache[r]
    S=[];A=[]
    for a,b in cand[:limit]:
        try: IA,IB=L(a),L(b)
        except Exception: continue
        dx,dy,sh,rt=pc_full(IA,IB)
        if sh<sharp or rt<ratio: continue
        dE=P[b][0]-P[a][0]; dN=P[b][1]-P[a][1]
        gm=math.hypot(dE,dN); pm=math.hypot(dx,dy)
        if gm<min_sep or pm<20: continue
        S.append(pm/gm)
        # bearing of the same baseline in each frame; image y is DOWN
        A.append(((math.atan2(-dy,dx)-math.atan2(dN,dE)+math.pi)%(2*math.pi))-math.pi)
    if len(S)<5: return None
    S=np.array(S); A=np.degrees(np.unwrap(np.array(A)))
    # robust
    s_med=float(np.median(S)); a_med=float(np.median(A))
    keep=(np.abs(S-s_med)<0.06)&(np.abs(A-a_med)<8)
    return dict(px_per_m=float(np.median(S[keep])), rot_deg=float(np.median(A[keep])),
                n=int(keep.sum()), n_try=len(S),
                s_mad=float(np.median(np.abs(S[keep]-np.median(S[keep])))),
                a_mad=float(np.median(np.abs(A[keep]-np.median(A[keep])))),
                lat0=lat0, lon0=lon0, MLAT=MLAT, MLON=MLON)
