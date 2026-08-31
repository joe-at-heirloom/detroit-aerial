"""Measure residual misalignment between two rasters that already share a grid.
Dense grid of windowed phase correlations -> per-window (dx,dy) in pixels."""
import numpy as np, math
from PIL import Image
Image.MAX_IMAGE_PIXELS=None

def highpass(I,k=16):
    """remove tonal/film-response differences so cross-epoch texture can correlate"""
    h,w=I.shape
    im=Image.fromarray(np.clip(I,0,255).astype(np.uint8))
    bg=np.asarray(im.resize((max(1,w//k),max(1,h//k)),Image.BILINEAR).resize((w,h),Image.BILINEAR),dtype=np.float32)
    out=I-bg
    m=I>0
    if m.any():
        s=out[m].std()
        if s>1e-6: out=out*(30.0/s)
    out[~m]=0
    return out

def _pc(A,B):
    h,w=A.shape
    A=A-A.mean(); B=B-B.mean()
    if A.std()<3 or B.std()<3: return None
    wn=np.outer(np.hanning(h),np.hanning(w))
    F=np.fft.rfft2(A*wn)*np.conj(np.fft.rfft2(B*wn))
    F/=(np.abs(F)+1e-9)
    c=np.fft.fftshift(np.fft.irfft2(F,s=(h,w)))
    pk=np.unravel_index(np.argmax(c),c.shape); peak=c[pk]
    m=c.copy(); m[max(0,pk[0]-6):pk[0]+6, max(0,pk[1]-6):pk[1]+6]=-9
    ratio=peak/(m.max()+1e-12); sharp=peak/(c.std()+1e-12)
    # parabolic subpixel
    def sub(i,axis):
        if i<=0 or i>=c.shape[axis]-1: return 0.0
        if axis==0: a,b,cc=c[i-1,pk[1]],c[i,pk[1]],c[i+1,pk[1]]
        else:       a,b,cc=c[pk[0],i-1],c[pk[0],i],c[pk[0],i+1]
        d=(a-2*b+cc)
        return 0.5*(a-cc)/d if abs(d)>1e-9 else 0.0
    dy=pk[0]-h//2+sub(pk[0],0); dx=pk[1]-w//2+sub(pk[1],1)
    return dx,dy,float(sharp),float(ratio)

def grid_offsets(imgA,imgB,win=320,step=None,sharp_min=6.0,ratio_min=1.25,edge=40):
    """imgA,imgB: same-size float arrays. Returns list of (cx,cy,dx,dy,sharp)."""
    H,W=imgA.shape
    step=step or win//2
    out=[]
    for cy in range(edge+win//2, H-edge-win//2, step):
        for cx in range(edge+win//2, W-edge-win//2, step):
            y0=cy-win//2; x0=cx-win//2
            A=imgA[y0:y0+win, x0:x0+win]; B=imgB[y0:y0+win, x0:x0+win]
            if (imgA[y0:y0+win,x0:x0+win]==0).mean()>0.10 or (imgB[y0:y0+win,x0:x0+win]==0).mean()>0.10: continue
            r=_pc(A,B)
            if not r: continue
            dx,dy,s,rt=r
            if s<sharp_min or rt<ratio_min: continue
            if abs(dx)>win*0.35 or abs(dy)>win*0.35: continue
            out.append((cx,cy,dx,dy,s))
    return out

def summarise(offs,mpp,label=''):
    if len(offs)<4: return f"  {label}: only {len(offs)} windows matched"
    a=np.array(offs); dx=a[:,2]; dy=a[:,3]
    d=np.hypot(dx,dy)
    return (f"  {label}: n={len(offs):3d}  median shift ({np.median(dx):+6.2f},{np.median(dy):+6.2f}) px "
            f"= {np.median(d)*mpp:5.2f} m | spread(MAD) {np.median(np.abs(d-np.median(d)))*mpp:4.2f} m "
            f"| max {d.max()*mpp:5.1f} m")

def fit_affine(offs, W, H):
    """Fit dx,dy = a + b*x + c*y  (robust, IRLS). Returns (coef_x, coef_y, resid_px)."""
    a=np.array(offs,dtype=float)
    x=a[:,0]/W - 0.5; y=a[:,1]/H - 0.5
    M=np.column_stack([np.ones_like(x), x, y])
    w=np.ones(len(x))
    for _ in range(6):
        Wm=np.sqrt(w)[:,None]
        cx=np.linalg.lstsq(M*Wm, a[:,2]*np.sqrt(w), rcond=None)[0]
        cy=np.linalg.lstsq(M*Wm, a[:,3]*np.sqrt(w), rcond=None)[0]
        rx=a[:,2]-M@cx; ry=a[:,3]-M@cy; r=np.hypot(rx,ry)
        s=max(np.median(r)*1.4826,0.4); w=1/(1+(r/(3*s))**2)
    return cx, cy, np.hypot(a[:,2]-M@cx, a[:,3]-M@cy)
