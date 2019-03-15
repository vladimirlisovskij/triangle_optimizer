import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import seaborn
import array as arr
from triangleOptimizer import TriOpt


func = lambda x,y: (x+1)**2 + y**2 + np.exp(x) + np.exp(y+1)
func_sp = (lambda x,y: (x+1)**2 + y**2 + sp.exp(x) + sp.exp(y+1))
koor = np.array([
    [1,1],
    [-1,-1],
    [-1,1]
])

res = TriOpt.opt(
    func=func, # функция для минимизации
    func_sp = func_sp, # функция для SymPy
    koor=koor, # координаты вершин
    info = True, # дополнительная информация
    plot = (True,50), # построение графиков
    )
print (res)