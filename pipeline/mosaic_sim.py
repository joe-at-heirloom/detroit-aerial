"""Composite a block from per-frame SIMILARITY placements (rotation + translation).

The previous compositor only translated frames, so any per-frame crab angle showed
up as crooked stitching. Here each output pixel is mapped back through its frame's
own rotation, and frames are still selected nearest-centre (Voronoi) so seams fall
midway between exposures."""
import numpy as np, math, os
from PIL import Image
Image.MAX_IMAGE_PIXELS=None

def regularise(sol, min_ratio=1.10, radius=4000.0):
    """Frames with a weak arterial lock inherit rotation from confident neighbours."""
    recs=[r for r,v in sol.items() if v.get('ok')]
    conf=[r for r in recs if sol[r]['ratio']>=min_ratio]
    out={}
    for r in recs:
        v=dict(sol[r])
        if v['ratio']<min_ratio and conf:
            d=[(math.hypot(sol[c]['e']-v['e'], sol[c]['n']-v['n']), c) for c in conf]
            d.sort()
            near=[c for dist,c in d if dist<radius][:6] or [d[0][1]]
            w=np.array([sol[c]['ratio'] for c in near])
            v['rot']=float(np.average([sol[c]['rot'] for c in near], weights=w))
            v['dE']=float(np.average([sol[c]['dE'] for c in near], weights=w))
            v['dN']=float(np.average([sol[c]['dN'] for c in near], weights=w))
            v['borrowed']=True
        else:
            v['borrowed']=False
        out[r]=v
    return out

def build(sol, imgdir, ppm_full, out_mpp, pad=900.0, chunk=1024, log=print, fill=True):
    """sol: {rec: {e,n,rot,dE,dN,gw,gh}}. Applies (dE,dN) and rot to each frame."""
    P={}
    for r,v in sol.items():
        if not v.get('ok'): continue
        P[r]=dict(e=v['e']-v['dE'], n=v['n']-v['dN'], rot=v['rot'], gw=v['gw'], gh=v['gh'])
    if not P: return None
    E=[p['e'] for p in P.values()]; N=[p['n'] for p in P.values()]
    gwm=max(p['gw'] for p in P.values()); ghm=max(p['gh'] for p in P.values())
    minE=min(E)-gwm/2-pad; maxE=max(E)+gwm/2+pad
    minN=min(N)-ghm/2-pad; maxN=max(N)+ghm/2+pad
    W=int((maxE-minE)/out_mpp); H=int((maxN-minN)/out_mpp)
    log(f"  canvas {W} x {H} px @ {out_mpp} m/px  ({(maxE-minE)/1000:.1f} x {(maxN-minN)/1000:.1f} km)")
    out=np.memmap(os.environ.get('MOSAIC_RAW','mosaic_sim.raw'), dtype=np.uint8, mode='w+', shape=(H,W)); out[:]=0
    C=np.array([[P[r]['e'],P[r]['n']] for r in P]); keys=list(P)
    done=0
    for r in keys:
        p=P[r]
        try: im=Image.open(f"{imgdir}/{r}.jpg").convert('L')
        except Exception: continue
        I=np.asarray(im); ih,iw=I.shape
        mpp_f=p['gw']/iw
        rad=math.hypot(p['gw'],p['gh'])/2
        x0=int((p['e']-rad-minE)/out_mpp); x1=int((p['e']+rad-minE)/out_mpp)
        y0=int((maxN-(p['n']+rad))/out_mpp); y1=int((maxN-(p['n']-rad))/out_mpp)
        x0=max(0,x0); y0=max(0,y0); x1=min(W,x1); y1=min(H,y1)
        if x1<=x0 or y1<=y0: continue
        th=math.radians(p['rot']); ct,st=math.cos(th),math.sin(th)
        others=np.array([[P[k]['e'],P[k]['n']] for k in keys if k!=r])
        xs=(minE+(np.arange(x0,x1)+0.5)*out_mpp).astype(np.float32)
        for ya in range(y0,y1,chunk):
            yb=min(y1,ya+chunk)
            ns=(maxN-(np.arange(ya,yb)+0.5)*out_mpp).astype(np.float32)[:,None]
            dE=xs[None,:]-p['e']; dN=ns-p['n']
            lx= ct*dE + st*dN
            ly=-st*dE + ct*dN
            px=lx/mpp_f + iw/2.0
            py=-ly/mpp_f + ih/2.0
            m=(px>=0)&(px<iw-1)&(py>=0)&(py<ih-1)
            if not m.any(): continue
            d0=dE*dE+dN*dN
            if len(others):
                best=np.full(d0.shape, np.inf, np.float32)
                for (ge,gn) in others:
                    ddx=xs[None,:]-ge; ddy=ns-gn
                    np.minimum(best, ddx*ddx+ddy*ddy, out=best)
                m &= (d0<=best)
            if not m.any(): continue
            xi=np.clip(px,0,iw-2); yi=np.clip(py,0,ih-2)
            xa=xi.astype(np.int32); yb_=yi.astype(np.int32)
            fx=(xi-xa).astype(np.float32); fy=(yi-yb_).astype(np.float32)
            v=(I[yb_,xa]*(1-fx)*(1-fy)+I[yb_,xa+1]*fx*(1-fy)
               +I[yb_+1,xa]*(1-fx)*fy+I[yb_+1,xa+1]*fx*fy)
            m &= (v>0)
            tgt=out[ya:yb, x0:x1]
            np.copyto(tgt, v.astype(np.uint8), where=m)
            out[ya:yb, x0:x1]=tgt
        done+=1
        if done%15==0: log(f"    {done}/{len(keys)} frames")
    if fill:
        # second pass: rotation moves footprints, so a pixel's nearest frame may not
        # actually cover it. Fill any remaining hole from whichever frame does.
        filled=0
        for r in keys:
            p=P[r]
            try: im=Image.open(f"{imgdir}/{r}.jpg").convert('L')
            except Exception: continue
            I=np.asarray(im); ih,iw=I.shape
            mpp_f=p['gw']/iw
            rad=math.hypot(p['gw'],p['gh'])/2
            x0=max(0,int((p['e']-rad-minE)/out_mpp)); x1=min(W,int((p['e']+rad-minE)/out_mpp))
            y0=max(0,int((maxN-(p['n']+rad))/out_mpp)); y1=min(H,int((maxN-(p['n']-rad))/out_mpp))
            if x1<=x0 or y1<=y0: continue
            th=math.radians(p['rot']); ct,st=math.cos(th),math.sin(th)
            xs=(minE+(np.arange(x0,x1)+0.5)*out_mpp).astype(np.float32)
            for ya in range(y0,y1,chunk):
                yb=min(y1,ya+chunk)
                tgt=out[ya:yb, x0:x1]
                hole=(tgt==0)
                if not hole.any(): continue
                ns=(maxN-(np.arange(ya,yb)+0.5)*out_mpp).astype(np.float32)[:,None]
                dE=xs[None,:]-p['e']; dN=ns-p['n']
                lx= ct*dE + st*dN
                ly=-st*dE + ct*dN
                px=lx/mpp_f + iw/2.0
                py=-ly/mpp_f + ih/2.0
                m=hole&(px>=0)&(px<iw-1)&(py>=0)&(py<ih-1)
                if not m.any(): continue
                xi=np.clip(px,0,iw-2); yi=np.clip(py,0,ih-2)
                xa=xi.astype(np.int32); yb_=yi.astype(np.int32)
                fx=(xi-xa).astype(np.float32); fy=(yi-yb_).astype(np.float32)
                v=(I[yb_,xa]*(1-fx)*(1-fy)+I[yb_,xa+1]*fx*(1-fy)
                   +I[yb_+1,xa]*(1-fx)*fy+I[yb_+1,xa+1]*fx*fy)
                m &= (v>0)
                if not m.any(): continue
                np.copyto(tgt, v.astype(np.uint8), where=m)
                out[ya:yb, x0:x1]=tgt
                filled+=int(m.sum())
        log(f"    hole-fill pass wrote {filled/1e6:.1f} Mpx")
    out.flush()
    return dict(arr=out, W=W, H=H, minE=minE, maxN=maxN, mpp=out_mpp)
