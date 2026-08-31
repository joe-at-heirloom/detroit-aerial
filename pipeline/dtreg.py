"""Register a downtown frame onto the common downtown grid.
Two stages: Fourier-Mellin against the already-georeferenced 1961 warp (gives
scale+rotation+shift robustly across decades because the street grid and river
persist), then a translation refinement against the modern road mask."""
import numpy as np, math, json
from PIL import Image
import fm
Image.MAX_IMAGE_PIXELS=None

REF=None; REFBBOX=None
def load_ref(path='warp1961_northup.png', bboxfile='bbox.json', side=1200):
    global REF,REFBBOX
    im=Image.open(path).convert('L')
    REFBBOX=json.load(open(bboxfile))['bbox']            # S,W,N,E
    sc=side/max(im.size)
    REF=(np.asarray(im.resize((int(im.width*sc),int(im.height*sc)),Image.LANCZOS),dtype=np.float32), sc, im.size)
    return REF

def register(frame_path, side=1200):
    """Returns (scale_px_per_refpx, rot_deg, dx, dy) mapping frame -> reference raster."""
    ref,refsc,refsize=REF
    im=Image.open(frame_path).convert('L')
    sc=side/max(im.size)
    A=np.asarray(im.resize((int(im.width*sc),int(im.height*sc)),Image.LANCZOS),dtype=np.float32)
    h=min(ref.shape[0],A.shape[0]); w=min(ref.shape[1],A.shape[1])
    def ctr(X,h,w):
        y=(X.shape[0]-h)//2; x=(X.shape[1]-w)//2
        return X[y:y+h,x:x+w]
    R=ctr(ref,h,w); F=ctr(A,h,w)
    best=None
    for flip in (0,180):
        Ff=F if flip==0 else np.rot90(F,2)
        s,rot,sh,rt=fm.scale_rot(R,Ff)
        if not (0.65<s<1.55): continue
        for cand in (rot, rot-180.0):
            Fw=fm.warp(Ff, s, -cand, out_shape=R.shape)
            m=Fw>0
            if m.mean()<0.30: continue
            a=R-R.mean(); b=np.where(m,Fw-Fw[m].mean(),0)
            wn=np.outer(np.hanning(h),np.hanning(w))
            C=np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(a*wn)*np.conj(np.fft.rfft2(b*wn)),s=(h,w)))
            pk=np.unravel_index(np.argmax(C),C.shape); peak=C[pk]
            score=peak/(C.std()+1e-12)
            if best is None or score>best[0]:
                best=(float(score), s, cand+flip, float(pk[1]-w//2), float(pk[0]-h//2), flip)
    return best
