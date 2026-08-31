"""Epoch-agnostic frame-to-frame registration.
Recovers (scale, rotation, translation) between two aerial frames with NO prior
assumption about mission geometry.  Fourier-Mellin: |FFT| is translation
invariant, so log-polar phase correlation of the spectra yields scale+rotation;
translation then follows from a second phase correlation."""
import numpy as np, math

def _win(h,w): return np.outer(np.hanning(h), np.hanning(w))

def _logpolar(M, nr=256, nt=360):
    h,w=M.shape; cy,cx=h/2.0,w/2.0
    rmax=min(cy,cx)
    # log-radial sampling, skip the DC core
    r=np.exp(np.linspace(math.log(3.0), math.log(rmax), nr))
    t=np.linspace(0, math.pi, nt, endpoint=False)     # spectra are symmetric -> pi
    R,T=np.meshgrid(r,t,indexing='ij')
    ys=(cy+R*np.sin(T)).astype(np.int32); xs=(cx+R*np.cos(T)).astype(np.int32)
    np.clip(ys,0,h-1,out=ys); np.clip(xs,0,w-1,out=xs)
    return M[ys,xs], r, t

def _pc(a,b):
    """phase correlation; returns (dy,dx,peak,ratio)"""
    A=np.fft.rfft2(a); B=np.fft.rfft2(b)
    R=A*np.conj(B); R/=(np.abs(R)+1e-9)
    c=np.fft.fftshift(np.fft.irfft2(R,s=a.shape))
    pk=np.unravel_index(np.argmax(c),c.shape); peak=c[pk]
    m=c.copy(); m[max(0,pk[0]-6):pk[0]+6, max(0,pk[1]-6):pk[1]+6]=-9
    ratio=peak/(m.max()+1e-12)
    return pk[0]-c.shape[0]//2, pk[1]-c.shape[1]//2, float(peak/(c.std()+1e-12)), float(ratio)

def spectrum(I):
    h,w=I.shape
    F=np.fft.fftshift(np.abs(np.fft.fft2((I-I.mean())*_win(h,w))))
    # high-emphasis filter suppresses the low-frequency dominance
    cy,cx=h/2,w/2
    Y,X=np.ogrid[:h,:w]
    rad=np.sqrt(((Y-cy)/cy)**2+((X-cx)/cx)**2)
    return np.log1p(F)*(0.5+rad)

def scale_rot(A,B,nr=256,nt=360):
    """relative scale and rotation taking B -> A"""
    SA,_ ,_=_logpolar(spectrum(A),nr,nt)+(None,) if False else (*_logpolar(spectrum(A),nr,nt),)
    SB,r,t=_logpolar(spectrum(B),nr,nt)
    SA=SA-SA.mean(); SB=SB-SB.mean()
    dr,dt,sharp,ratio=_pc(SA,SB)
    logstep=(math.log(r[-1])-math.log(r[0]))/(nr-1)
    scale=math.exp(dr*logstep)
    rot=math.degrees(dt*(math.pi/nt))
    return scale, rot, sharp, ratio

def warp(I, scale, rot_deg, out_shape=None):
    """rotate+scale I about its centre (nearest-neighbour, adequate at these sizes)"""
    h,w=I.shape
    oh,ow=out_shape or (h,w)
    th=math.radians(rot_deg); c,s=math.cos(th),math.sin(th)
    Y,X=np.mgrid[0:oh,0:ow]
    yc,xc=oh/2,ow/2
    dy=(Y-yc); dx=(X-xc)
    sy=( c*dy + s*dx)/scale + h/2
    sx=(-s*dy + c*dx)/scale + w/2
    yi=np.clip(sy.astype(np.int32),0,h-1); xi=np.clip(sx.astype(np.int32),0,w-1)
    out=I[yi,xi].astype(np.float32)
    out[(sy<0)|(sy>=h)|(sx<0)|(sx>=w)]=0
    return out

def register(A,B, try_180=True):
    """Full pairwise solve. Returns dict or None."""
    best=None
    for flip in ([0,180] if try_180 else [0]):
        Bf = B if flip==0 else np.rot90(B,2)
        sc,rot,sh,rt = scale_rot(A,Bf)
        if not (0.7<sc<1.4): continue
        for cand_rot in (rot, rot-180.0):
            Bw=warp(Bf, sc, -cand_rot, out_shape=A.shape)
            m=(Bw>0)
            if m.mean()<0.35: continue
            a=A-A.mean(); b=np.where(m,Bw-Bw[m].mean(),0)
            dy,dx,sharp,ratio=_pc(a*_win(*A.shape), b*_win(*A.shape))
            score=sharp
            if best is None or score>best['sharp']:
                best={'scale':sc,'rot':cand_rot+flip,'dx':dx,'dy':dy,'sharp':sharp,'ratio':ratio,'flip':flip}
    return best
