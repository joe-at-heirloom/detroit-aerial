import json, math, os, glob, numpy as np, urllib.request
from PIL import Image, ImageDraw
SCALE_FULL=1.59943; ROT_DEG=-0.1719; DS=4
BASE_SCALE=SCALE_FULL/DS
LAT0,LON0=42.3500,-83.1000
def _md(phi):
    p=math.radians(phi)
    return (111132.92-559.82*math.cos(2*p)+1.175*math.cos(4*p),
            111412.84*math.cos(p)-93.5*math.cos(3*p))
MLAT,MLON=_md(LAT0)
def to_m(lat,lon): return ((lon-LON0)*MLON,(lat-LAT0)*MLAT)
def to_ll(ex,ey): return (LAT0+ey/MLAT, LON0+ex/MLON)

_SEG=None;_GRID=None;CELL=2000.0
# roads that plausibly existed in 1961: exclude motorway/trunk (freeway era)
KEEP={'primary','secondary','tertiary','residential','unclassified'}
def load_roads():
    global _SEG,_GRID
    if _SEG is not None: return
    segs=[]
    for p in sorted(glob.glob('osmtiles/r_*.json')):
        for w in json.load(open(p)).get('elements',[]):
            if w.get('tags',{}).get('highway') not in KEEP: continue
            g=w.get('geometry')
            if not g or len(g)<2: continue
            pts=[to_m(q['lat'],q['lon']) for q in g]
            for a,b in zip(pts,pts[1:]): segs.append((a[0],a[1],b[0],b[1]))
    _SEG=np.array(segs,dtype=np.float32); grid={}
    for i,(x0,y0,x1,y1) in enumerate(_SEG):
        for cx in range(int(min(x0,x1)//CELL),int(max(x0,x1)//CELL)+1):
            for cy in range(int(min(y0,y1)//CELL),int(max(y0,y1)//CELL)+1):
                grid.setdefault((cx,cy),[]).append(i)
    _GRID=grid
def near(cx,cy,rad=3000):
    load_roads(); s=set()
    for gx in range(int((cx-rad)//CELL),int((cx+rad)//CELL)+1):
        for gy in range(int((cy-rad)//CELL),int((cy+rad)//CELL)+1):
            s.update(_GRID.get((gx,gy),()))
    return _SEG[sorted(s)] if s else np.empty((0,4),np.float32)
TH=math.radians(ROT_DEG);CT,ST=math.cos(TH),math.sin(TH)
def road_mask(cex,cey,w,h,sc,width=2):
    S=near(cex,cey)
    m=Image.new('L',(w,h),0);d=ImageDraw.Draw(m)
    if len(S):
        dx0,dy0=S[:,0]-cex,S[:,1]-cey; dx1,dy1=S[:,2]-cex,S[:,3]-cey
        X0= sc*(CT*dx0-ST*dy0)+w/2; Y0=-sc*(ST*dx0+CT*dy0)+h/2
        X1= sc*(CT*dx1-ST*dy1)+w/2; Y1=-sc*(ST*dx1+CT*dy1)+h/2
        for a,b,c,e in zip(X0,Y0,X1,Y1):
            if (a<-60 and c<-60) or (a>w+60 and c>w+60) or (b<-60 and e<-60) or (b>h+60 and e>h+60): continue
            d.line([(a,b),(c,e)],fill=255,width=width)
    return np.asarray(m,dtype=np.float32)/255.0
def highpass(I,k=25):
    h,w=I.shape
    im=Image.fromarray(np.clip(I,0,255).astype(np.uint8))
    bg=np.asarray(im.resize((max(1,w//k),max(1,h//k)),Image.BILINEAR).resize((w,h),Image.BILINEAR),dtype=np.float32)
    return I-bg
def _corr(I,M):
    Mz=M-M.mean()
    C=np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(I)*np.conj(np.fft.rfft2(Mz)),s=I.shape))
    pk=np.unravel_index(np.argmax(C),C.shape); peak=C[pk]
    # second peak outside an exclusion radius -> ambiguity ratio
    Y,X=np.ogrid[:C.shape[0],:C.shape[1]]
    far=((Y-pk[0])**2+(X-pk[1])**2) > 40**2
    second=C[far].max() if far.any() else 0.0
    ratio=peak/(second+1e-9)
    return pk,peak,ratio,C
def register(img,cex,cey,scales=None):
    h,w=img.shape
    I=highpass(img); I=(I-I.mean())/(I.std()+1e-6)
    best=None
    for f in (scales or [0.94,0.96,0.98,1.00,1.02,1.04,1.06]):
        sc=BASE_SCALE*f
        M=road_mask(cex,cey,w,h,sc)
        if M.sum()<400: continue
        pk,peak,ratio,_=_corr(I,M)
        norm=peak/math.sqrt(M.sum())
        if best is None or norm>best[0]: best=(norm,ratio,f,sc,pk,h,w)
    if best is None: return None
    norm,ratio,f,sc,pk,h,w=best
    dy=pk[0]-h//2; dx=pk[1]-w//2
    ex=cex-(CT*dx/sc - ST*(-dy)/sc); ey=cey-(ST*dx/sc + CT*(-dy)/sc)
    return {'norm':float(norm),'ratio':float(ratio),'sfac':f,'scale':sc,
            'dx':int(dx),'dy':int(dy),'shift_m':math.hypot(dx,dy)/sc,
            'cex':float(ex),'cey':float(ey)}
