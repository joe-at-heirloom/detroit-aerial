"""Invert the solved anchor onto the frame placements.

solve_similarity finds T such that T(reference) matches the image; correcting the
image therefore means applying T^-1 to every frame's world placement."""
import math, json, dtmap
def corrected(sol, anchor, cE, cN):
    rot=anchor['rot']; dE=anchor['dE']; dN=anchor['dN']
    th=math.radians(-rot); c,s=math.cos(th),math.sin(th)
    out={}
    for r,v in sol.items():
        if not v.get('ok'): continue
        qx=v['e']-dE-cE; qy=v['n']-dN-cN
        out[r]=dict(v)
        out[r]['e']=cE+c*qx-s*qy
        out[r]['n']=cN+s*qx+c*qy
        out[r]['rot']=v['rot']-rot
    return out
