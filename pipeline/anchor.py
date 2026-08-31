"""Anchor a solved block to the map. At block scale (10s of km) the street-grid
periodicity that defeats single-frame matching is no longer ambiguous."""
import json, glob, math, numpy as np
from PIL import Image, ImageDraw

def load_roads(lat0,lon0,MLAT,MLON,keep=None):
    segs=[]
    for p in sorted(glob.glob('osmtiles/r_*.json')):
        for w in json.load(open(p)).get('elements',[]):
            if keep and w.get('tags',{}).get('highway') not in keep: continue
            g=w.get('geometry')
            if not g or len(g)<2: continue
            pts=[(((q['lon']-lon0)*MLON),((q['lat']-lat0)*MLAT)) for q in g]
            for a,b in zip(pts,pts[1:]): segs.append((a[0],a[1],b[0],b[1]))
    return np.array(segs,dtype=np.float32)

def mosaic(pos, imgdir, ppm, rot_deg=0.0, pad=1500, maxpx=6000):
    E=[v[0] for v in pos.values()]; N=[v[1] for v in pos.values()]
    minE,maxE=min(E)-pad,max(E)+pad; minN,maxN=min(N)-pad,max(N)+pad
    W=int((maxE-minE)*ppm); H=int((maxN-minN)*ppm)
    sc=min(1.0, maxpx/max(W,H)); W=int(W*sc); H=int(H*sc); p=ppm*sc
    acc=np.zeros((H,W),np.float32); cnt=np.zeros((H,W),np.float32)
    for r,(e,n) in pos.items():
        try: im=Image.open(f"{imgdir}/{r}.jpg").convert('L')
        except Exception: continue
        if sc<1.0: im=im.resize((max(1,int(im.width*sc)),max(1,int(im.height*sc))),Image.LANCZOS)
        I=np.asarray(im,dtype=np.float32); h,w=I.shape
        x0=int((e-minE)*p-w/2); y0=int((maxN-n)*p-h/2)
        fy=np.minimum(np.arange(h),h-1-np.arange(h)); fx=np.minimum(np.arange(w),w-1-np.arange(w))
        fw=np.clip(np.minimum(fy[:,None],fx[None,:])/50.0,0,1).astype(np.float32)
        xs0,ys0=max(0,x0),max(0,y0); xs1,ys1=min(W,x0+w),min(H,y0+h)
        if xs1<=xs0 or ys1<=ys0: continue
        acc[ys0:ys1,xs0:xs1]+=I[ys0-y0:ys1-y0,xs0-x0:xs1-x0]*fw[ys0-y0:ys1-y0,xs0-x0:xs1-x0]
        cnt[ys0:ys1,xs0:xs1]+=fw[ys0-y0:ys1-y0,xs0-x0:xs1-x0]
    out=np.where(cnt>0.01,acc/np.maximum(cnt,1e-6),0)
    return out,(minE,maxN,p),(cnt>0.01)

def road_raster(S,minE,maxN,p,W,H,dE=0.0,dN=0.0,rot=0.0,scale=1.0,width=1):
    th=math.radians(rot); c,s=math.cos(th),math.sin(th)
    def tf(x,y):
        x=x*scale; y=y*scale
        X=c*x-s*y+dE; Y=s*x+c*y+dN
        return (X-minE)*p, (maxN-Y)*p
    m=Image.new('L',(W,H),0); d=ImageDraw.Draw(m)
    X0,Y0=tf(S[:,0],S[:,1]); X1,Y1=tf(S[:,2],S[:,3])
    keep=~(((X0<0)&(X1<0))|((X0>W)&(X1>W))|((Y0<0)&(Y1<0))|((Y0>H)&(Y1>H)))
    for a,b,cc,e in zip(X0[keep],Y0[keep],X1[keep],Y1[keep]): d.line([(a,b),(cc,e)],fill=255,width=width)
    return np.asarray(m,dtype=np.float32)/255.0

def hp(I,k=30):
    h,w=I.shape
    im=Image.fromarray(np.clip(I,0,255).astype(np.uint8))
    return I-np.asarray(im.resize((max(1,w//k),max(1,h//k)),Image.BILINEAR).resize((w,h),Image.BILINEAR),dtype=np.float32)

def refine(M,valid,S,geo,search_m=1200,rot_range=5.0,nrot=21,scales=(1.0,)):
    """FFT translation search x coarse rotation search. Returns (dE,dN,rot,score)."""
    minE,maxN,p=geo; H,W=M.shape
    I=hp(M)*valid; I=(I-I[valid].mean())/(I[valid].std()+1e-6); I*=valid
    best=None
    for scl in scales:
      for rot in np.linspace(-rot_range,rot_range,nrot):
        R=road_raster(S,minE,maxN,p,W,H,rot=rot,scale=scl)
        if R.sum()<500: continue
        Rz=R-R.mean()
        C=np.fft.fftshift(np.fft.irfft2(np.fft.rfft2(I)*np.conj(np.fft.rfft2(Rz)),s=I.shape))
        lim=int(search_m*p)
        cy,cx=H//2,W//2
        sub=C[max(0,cy-lim):cy+lim, max(0,cx-lim):cx+lim]
        pk=np.unravel_index(np.argmax(sub),sub.shape)
        peak=sub[pk]; sc=peak/ (C.std()+1e-9)
        dy=pk[0]+max(0,cy-lim)-cy; dx=pk[1]+max(0,cx-lim)-cx
        if best is None or sc>best[3]: best=(-dx/p, dy/p, float(rot), float(sc), float(scl))
    return best
