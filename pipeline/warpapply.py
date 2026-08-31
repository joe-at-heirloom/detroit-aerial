import numpy as np, json, math
def design(x,y,deg):
    cols=[np.ones_like(np.asarray(x,dtype=float))]
    for d in range(1,deg+1):
        for k in range(d+1): cols.append((np.asarray(x,float)**(d-k))*(np.asarray(y,float)**k))
    return np.column_stack(cols)
class Warp:
    def __init__(self,path):
        w=json.load(open(path)); self.__dict__.update(w)
    def at(self,cE,cN):
        X=(np.atleast_1d(cE)-self.mx)/self.sx; Y=(np.atleast_1d(cN)-self.my)/self.sy
        A=design(X,Y,self.deg)
        return (float((A@np.array(self.cE))[0]), float((A@np.array(self.cN))[0]), float((A@np.array(self.cR))[0]))
