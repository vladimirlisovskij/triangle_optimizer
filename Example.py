import numpy as np
import sympy as sp
from triangleOptimizer import TriOpt
    
func = lambda x,y: (x+1)**2 + y**2 + np.exp(x) + np.exp(y+1)
func_sp = (lambda x,y: (x+1)**2 + y**2 + sp.exp(x) + sp.exp(y+1))
koor = np.array([
    [1,1],
    [-1,-1],
    [-1,1]
])

print(TriOpt.opt(func=func,func_sp= func_sp,koor=koor,info = True,plot = (True,15)))
