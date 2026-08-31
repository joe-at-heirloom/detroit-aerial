"""Relative orientation WITH rotation.

Frame-to-frame is where rotation is actually observable: ~60% overlap, identical
sensor, six-second interval. We rotate one window against the other over a small
angular range and keep the angle that maximises phase-correlation sharpness."""
import numpy as np, math
from PIL import Image
Image.MAX_IMAGE_PIXELS=None

def _pc(A,B):
    h,w=A.shape
    A=A-A.mean(); B=B-B.mean()
    wn=np.outer(np.hanning(h),np.hanning(w))
    F=np.fft.rfft2(A*wn)*np.conj(np.fft.rfft2(B*wn)); F/=(np.abs(F)+1e-9)
    c=np.fft.fftshift(np.fft.irfft2(F,s=(h,w)))
    pk=np.unravel_index(np.argmax(c),c.shape); peak=c[pk]
    m=c.copy(); m[max(0,pk[0]-8):pk[0]+8, max(0,pk[1]-8):pk[1]+8]=-9
    return pk[1]-w//2, pk[0]-h//2, float(peak/(c.std()+1e-12)), float(peak/(m.max()+1e-12))

def _rot(I,deg):
    if abs(deg)<1e-6: return I
    h,w=I.shape; th=math.radians(deg); c,s=math.cos(th),math.sin(th)
    Y,X=np.mgrid[0:h,0:w]; dy=Y-h/2.0; dx=X-w/2.0
    sy=(-s*dx + c*dy)+h/2.0; sx=( c*dx + s*dy)+w/2.0
    yi=np.clip(sy,0,h-1).astype(np.int32); xi=np.clip(sx,0,w-1).astype(np.int32)
    o=I[yi,xi].astype(np.float32)
    o[(sy<0)|(sy>=h)|(sx<0)|(sx>=w)]=0
    return o

def match(load,a,b,pa,pb,ppm,win=384,rots=np.arange(-3.0,3.01,0.5)):
    """Returns dict(dE,dN,drot,sharp,ratio) mapping b -> a, or None."""
    A=load(a); B=load(b)
    dE=pb[0]-pa[0]; dN=pb[1]-pa[1]
    ox,oy=dE/2,dN/2
    ax=A.shape[1]/2+ox*ppm; ay=A.shape[0]/2-oy*ppm
    bx=B.shape[1]/2-ox*ppm; by=B.shape[0]/2+oy*ppm
    def cut(I,cx,cy,pad=1.45):
        s=int(win*pad/2)
        x0=int(cx-s); y0=int(cy-s)
        if x0<0 or y0<0 or x0+2*s>I.shape[1] or y0+2*s>I.shape[0]: return None
        return I[y0:y0+2*s, x0:x0+2*s]
    PA=cut(A,ax,ay); PB=cut(B,bx,by)
    if PA is None or PB is None or PA.std()<7 or PB.std()<7: return None
    def centre(I):
        h,w=I.shape; s=win//2
        return I[h//2-s:h//2+s, w//2-s:w//2+s]
    CA=centre(PA)
    if CA.shape[0]<win or CA.shape[1]<win: return None
    best=None
    for rot in rots:
        RB=centre(_rot(PB,rot))
        if RB.shape!=CA.shape or RB.std()<5: continue
        dx,dy,sh,rt=_pc(CA,RB)
        if best is None or sh>best[2]:
            best=(dx,dy,sh,rt,rot)
    if best is None: return None
    dx,dy,sh,rt,rot=best
    return dict(dE=dE+dx/ppm, dN=dN-dy/ppm, drot=float(rot), sharp=float(sh), ratio=float(rt))
