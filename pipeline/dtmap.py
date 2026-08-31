"""Register one downtown frame straight to the map by correlating a rasterised
road network against the image. This is the method that scored 15 sigma on the
1961 downtown frame; here it also searches scale and rotation."""
import numpy as np, math, glob, json
from PIL import Image, ImageDraw
import os as _os
_DATA = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))), 'data')
def _d(p):
    return p if _os.path.isabs(p) or _os.path.exists(p) else _os.path.join(_DATA, _os.path.basename(p))

Image.MAX_IMAGE_PIXELS=None
LAT0,LON0=42.3340,-83.0450
def _md(phi):
    p=math.radians(phi)
    return (111132.92-559.82*math.cos(2*p)+1.175*math.cos(4*p),
            111412.84*math.cos(p)-93.5*math.cos(3*p))
MLAT,MLON=_md(LAT0)
def to_m(lat,lon): return ((lon-LON0)*MLON,(lat-LAT0)*MLAT)
KEEP={'primary','secondary','tertiary','residential','unclassified'}
_S=None
def roads():
    global _S
    if _S is None:
        segs=[]
        for p in sorted(glob.glob(_os.path.join(_DATA, 'osmtiles/r_*.json')))+['roads.json','streets.json']:
            try: els=json.load(open(p)).get('elements',[])
            except Exception: continue
            for w in els:
                if w.get('tags',{}).get('highway') not in KEEP: continue
                g=w.get('geometry')
                if not g or len(g)<2: continue
                pts=[to_m(q['lat'],q['lon']) for q in g]
                for a,b in zip(pts,pts[1:]): segs.append((a[0],a[1],b[0],b[1]))
        _S=np.array(segs,dtype=np.float32)
    return _S
def mask(cex,cey,w,h,ppm,rot_deg,width=2):
    S=roads()
    d2=(S[:,0]-cex)**2+(S[:,1]-cey)**2
    S=S[d2<(3200**2)]
    th=math.radians(rot_deg); ct,st=math.cos(th),math.sin(th)
    m=Image.new('L',(w,h),0); d=ImageDraw.Draw(m)
    def tf(x,y):
        dx,dy=x-cex,y-cey
        return ppm*(ct*dx-st*dy)+w/2, -ppm*(st*dx+ct*dy)+h/2
    X0,Y0=tf(S[:,0],S[:,1]); X1,Y1=tf(S[:,2],S[:,3])
    k=~(((X0<-40)&(X1<-40))|((X0>w+40)&(X1>w+40))|((Y0<-40)&(Y1<-40))|((Y0>h+40)&(Y1>h+40)))
    for a,b,c,e in zip(X0[k],Y0[k],X1[k],Y1[k]): d.line([(a,b),(c,e)],fill=255,width=width)
    return np.asarray(m,dtype=np.float32)/255.0
def highpass(I,k=25):
    h,w=I.shape
    im=Image.fromarray(np.clip(I,0,255).astype(np.uint8))
    return I-np.asarray(im.resize((max(1,w//k),max(1,h//k)),Image.BILINEAR).resize((w,h),Image.BILINEAR),dtype=np.float32)
def register(img, lat, lon, scales, rots):
    h,w=img.shape
    I=highpass(img); I=(I-I.mean())/(I.std()+1e-6)
    cex,cey=to_m(lat,lon)
    best=None
    for ppm in scales:
        for rot in rots:
            M=mask(cex,cey,w,h,ppm,rot)
            if M.sum()<400: continue
            Mz=M-M.mean()
            C=np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(I)*np.conj(np.fft.rfft2(Mz)),s=I.shape))
            pk=np.unravel_index(np.argmax(C),C.shape)
            sc=(C[pk]-C.mean())/(C.std()+1e-9)
            norm=sc/math.sqrt(M.sum())
            if best is None or norm>best[0]:
                dy=pk[0]-h//2; dx=pk[1]-w//2
                th=math.radians(rot); ct,st=math.cos(th),math.sin(th)
                # invert the placement to recover the true centre
                ex=cex-(ct*dx/ppm - st*(-dy)/ppm); ey=cey-(st*dx/ppm + ct*(-dy)/ppm)
                best=(norm,float(sc),ppm,rot,float(ex),float(ey))
    return best
