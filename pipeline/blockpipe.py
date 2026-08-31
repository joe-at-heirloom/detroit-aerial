"""Epoch-agnostic aerial block pipeline.
   calibrate -> iterative windowed tie matching -> robust bundle adjustment
   -> Procrustes anchoring to catalogue coords -> road-correlation refinement.
   No mission constants are hardcoded anywhere."""
import json, math, itertools, collections, numpy as np
from PIL import Image, ImageDraw

# ---------- image cache ----------
class Cache:
    def __init__(self,d,n=70): self.d=d; self.n=n; self.c={}
    def __call__(self,r):
        if r not in self.c:
            if len(self.c)>self.n: self.c.clear()
            self.c[r]=np.asarray(Image.open(f"{self.d}/{r}.jpg").convert('L'),dtype=np.float32)
        return self.c[r]

def _pc(A,B):
    h,w=A.shape
    A=A-A.mean(); B=B-B.mean()
    wn=np.outer(np.hanning(h),np.hanning(w))
    F=np.fft.rfft2(A*wn)*np.conj(np.fft.rfft2(B*wn)); F/=(np.abs(F)+1e-9)
    c=np.fft.fftshift(np.fft.irfft2(F,s=(h,w)))
    pk=np.unravel_index(np.argmax(c),c.shape); peak=c[pk]
    m=c.copy(); m[max(0,pk[0]-10):pk[0]+10,max(0,pk[1]-10):pk[1]+10]=-9
    return pk[1]-w//2, pk[0]-h//2, float(peak/(c.std()+1e-12)), float(peak/(m.max()+1e-12))

def match_pair(load,a,b,pa,pb,ppm,win=384):
    """Correlate the predicted overlap window. Returns (dE,dN,sharp,ratio) or None."""
    A=load(a); B=load(b)
    dE=pb[0]-pa[0]; dN=pb[1]-pa[1]
    ox,oy=dE/2,dN/2
    ax=A.shape[1]/2+ox*ppm; ay=A.shape[0]/2-oy*ppm
    bx=B.shape[1]/2-ox*ppm; by=B.shape[0]/2+oy*ppm
    def cut(I,cx,cy):
        x0=int(cx-win/2); y0=int(cy-win/2)
        if x0<0 or y0<0 or x0+win>I.shape[1] or y0+win>I.shape[0]: return None
        return I[y0:y0+win,x0:x0+win]
    PA=cut(A,ax,ay); PB=cut(B,bx,by)
    if PA is None or PB is None or PA.std()<7 or PB.std()<7: return None
    dx,dy,sh,rt=_pc(PA,PB)
    return dE+dx/ppm, dN-dy/ppm, sh, rt

def build_ties(frames,imgdir,calib,rounds=3,sharp=9,ratio=1.30,lo=150,hi=2800):
    load=Cache(imgdir)
    ppm=calib['px_per_m']; MLAT,MLON=calib['MLAT'],calib['MLON']
    lat0,lon0=calib['lat0'],calib['lon0']
    pos={f['rec']:(((f['lon']-lon0)*MLON),((f['lat']-lat0)*MLAT)) for f in frames}
    stated=dict(pos)
    edges={}; 
    for rd in range(rounds):
        cand=[(a,b) for a,b in itertools.combinations(pos,2)
              if (min(a,b),max(a,b)) not in edges
              and lo<math.hypot(pos[a][0]-pos[b][0],pos[a][1]-pos[b][1])<hi]
        added=0
        for a,b in cand:
            m=match_pair(load,a,b,pos[a],pos[b],ppm)
            if not m: continue
            dE,dN,sh,rt=m
            if sh>=sharp and rt>=ratio:
                edges[(min(a,b),max(a,b))]=(a,b,dE,dN,sh,rt); added+=1
        E=list(edges.values())
        pos=solve_blocks(list(pos), E, stated)[0]
        yield rd, added, len(edges), pos
    return

def components(recs,E):
    adj=collections.defaultdict(set)
    for a,b,*_ in E: adj[a].add(b); adj[b].add(a)
    seen=set(); out=[]
    for r in recs:
        if r in seen: continue
        st=[r]; c=[]
        while st:
            u=st.pop()
            if u in seen: continue
            seen.add(u); c.append(u); st+=[v for v in adj[u] if v not in seen]
        out.append(c)
    return sorted(out,key=len,reverse=True)

def bundle(recs,E,stated,lam=0.02,iters=6):
    idx={r:i for i,r in enumerate(recs)}; n=len(recs); m=len(E)
    st=np.array([stated[r] for r in recs])
    if m<1: return st[:,0],st[:,1],np.array([])
    w=np.ones(m)
    for _ in range(iters):
        A=np.zeros((m+n,n)); bE=np.zeros(m+n); bN=np.zeros(m+n); W=np.zeros(m+n)
        for k,(a,b,dE,dN,sh,rt) in enumerate(E):
            A[k,idx[a]]=-1; A[k,idx[b]]=1; bE[k]=dE; bN[k]=dN; W[k]=w[k]*min(math.sqrt(sh),8)
        for i in range(n): A[m+i,i]=1; bE[m+i]=st[i,0]; bN[m+i]=st[i,1]; W[m+i]=lam
        sw=np.sqrt(W)[:,None]
        sE=np.linalg.lstsq(A*sw,bE*np.sqrt(W),rcond=None)[0]
        sN=np.linalg.lstsq(A*sw,bN*np.sqrt(W),rcond=None)[0]
        r=np.array([math.hypot((sE[idx[b]]-sE[idx[a]])-dE,(sN[idx[b]]-sN[idx[a]])-dN) for a,b,dE,dN,_,_ in E])
        sg=max(np.median(r)*1.4826,5.0); w=1/(1+(r/(3*sg))**2)
    return sE,sN,r

def solve_blocks(recs,E,stated):
    comps=components(recs,E)
    pos={}; info=[]
    for c in comps:
        cs=set(c); Ec=[e for e in E if e[0] in cs and e[1] in cs]
        sE,sN,r=bundle(sorted(c),Ec,stated)
        for i,rr in enumerate(sorted(c)): pos[rr]=(float(sE[i]),float(sN[i]))
        info.append((len(c),len(Ec),float(np.median(r)) if len(r) else 0.0))
    for r in recs:
        if r not in pos: pos[r]=stated[r]
    return pos, comps, info

def run(frames,imgdir,calib,rounds=3,**kw):
    """returns (edges, pos, comps, info, stated)"""
    load=Cache(imgdir); ppm=calib['px_per_m']
    lat0,lon0,MLAT,MLON=calib['lat0'],calib['lon0'],calib['MLAT'],calib['MLON']
    stated={f['rec']:(((f['lon']-lon0)*MLON),((f['lat']-lat0)*MLAT)) for f in frames}
    pos=dict(stated); edges={}
    sharp=kw.get('sharp',9); ratio=kw.get('ratio',1.30); lo=kw.get('lo',150); hi=kw.get('hi',2800)
    for rd in range(rounds):
        cand=[(a,b) for a,b in itertools.combinations(pos,2)
              if (min(a,b),max(a,b)) not in edges
              and lo<math.hypot(pos[a][0]-pos[b][0],pos[a][1]-pos[b][1])<hi]
        for a,b in cand:
            m=match_pair(load,a,b,pos[a],pos[b],ppm)
            if not m: continue
            dE,dN,sh,rt=m
            if sh>=sharp and rt>=ratio: edges[(min(a,b),max(a,b))]=(a,b,dE,dN,sh,rt)
        E=list(edges.values())
        pos,comps,info=solve_blocks(list(stated),E,stated)
    return list(edges.values()),pos,comps,info,stated

def procrustes(local,geo):
    """similarity fit local->geo (both Nx2 metre arrays). Returns s,theta,t,resid"""
    A=np.asarray(local,float); B=np.asarray(geo,float)
    ma,mb=A.mean(0),B.mean(0); A0,B0=A-ma,B-mb
    za=A0[:,0]+1j*A0[:,1]; zb=B0[:,0]+1j*B0[:,1]
    k=(zb*np.conj(za)).sum()/((np.abs(za)**2).sum()+1e-12)
    s=abs(k); th=math.atan2(k.imag,k.real)
    pred=k*za
    resid=np.abs(pred-zb)
    return s,math.degrees(th),(mb[0]-(k*(ma[0]+1j*ma[1])).real, mb[1]-(k*(ma[0]+1j*ma[1])).imag),resid
