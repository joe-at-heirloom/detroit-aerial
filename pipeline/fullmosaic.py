"""Full-resolution nearest-frame mosaic.

Memory plan: the output is ~590 MP, so it lives in a uint8 memmap on disk.
We iterate over FRAMES (each loaded exactly once) rather than over output tiles,
and for each frame write only the pixels whose nearest frame-centre is that frame
(an exact Voronoi test against the frame's local competitors, evaluated in
row-chunks so no intermediate exceeds ~100 MB)."""
import numpy as np, math, os
from PIL import Image
Image.MAX_IMAGE_PIXELS=None

def plan(pos, ppm_full, pad_m=900):
    E=[v[0] for v in pos.values()]; N=[v[1] for v in pos.values()]
    minE,maxE=min(E)-pad_m,max(E)+pad_m; minN,maxN=min(N)-pad_m,max(N)+pad_m
    W=int(round((maxE-minE)*ppm_full)); H=int(round((maxN-minN)*ppm_full))
    return dict(minE=minE,maxN=maxN,W=W,H=H,p=ppm_full)

def competitors(pos, r, radius_m):
    e0,n0=pos[r]
    return [(k,)+tuple(pos[k]) for k in pos
            if k!=r and (pos[k][0]-e0)**2+(pos[k][1]-n0)**2 < radius_m**2]

def build(pos, imgdir, geo, out_path, chunk_rows=768, log=print):
    W,H,p,minE,maxN=geo['W'],geo['H'],geo['p'],geo['minE'],geo['maxN']
    out=np.memmap(out_path, dtype=np.uint8, mode='w+', shape=(H,W))
    out[:]=0
    # frame centres in output pixel space
    C={r:(( e-minE)*p, (maxN-n)*p) for r,(e,n) in pos.items()}
    done=0
    for r,(e,n) in pos.items():
        fp=f"{imgdir}/{r}.jpg"
        if not os.path.exists(fp): continue
        im=Image.open(fp).convert('L'); I=np.asarray(im); h,w=I.shape
        cx,cy=C[r]
        x0=int(round(cx-w/2)); y0=int(round(cy-h/2))
        xs0,ys0=max(0,x0),max(0,y0); xs1,ys1=min(W,x0+w),min(H,y0+h)
        if xs1<=xs0 or ys1<=ys0: continue
        comp=[(C[k][0],C[k][1]) for k in pos if k!=r
              and (C[k][0]-cx)**2+(C[k][1]-cy)**2 < (2.2*max(w,h))**2]
        xs=np.arange(xs0,xs1,dtype=np.float32)
        for ya in range(ys0,ys1,chunk_rows):
            yb=min(ys1,ya+chunk_rows)
            ys=np.arange(ya,yb,dtype=np.float32)[:,None]
            dx=xs[None,:]-cx; dy=ys-cy
            dmine=dx*dx+dy*dy
            best=np.full(dmine.shape, np.inf, np.float32)
            for (gx,gy) in comp:
                ddx=xs[None,:]-gx; ddy=ys-gy
                np.minimum(best, ddx*ddx+ddy*ddy, out=best)
            sub=I[ya-y0:yb-y0, xs0-x0:xs1-x0]
            m=(dmine<=best)&(sub>0)
            tgt=out[ya:yb, xs0:xs1]
            np.copyto(tgt, sub, where=m)
            out[ya:yb, xs0:xs1]=tgt
        done+=1
        if done%10==0: log(f"    composited {done}/{len(pos)} frames")
    out.flush()
    return out
