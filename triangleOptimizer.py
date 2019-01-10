if __name__ == '__main__': 
    import matplotlib.pyplot as plt
    import matplotlib.colors as clr
    import matplotlib.cm as cm
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    import sympy as sp
    from math import sqrt, ceil, log
    import time


class triOpt():
    def __GS(self, a,b,e,f):
        if a[0] == b[0]: 
            if a[1] > b[1]:
                a,b = b,a
            A = a[1]
            B = b[1]
            F = (1 + sqrt(5))/2    
            while (B-A)/2 >= e:            
                y1 = B - (B-A)/F
                y2 = A + (B-A)/F
                if f(a[0],y1) >= f(a[0],y2):
                    A = y1
                else:
                    B = y2
            return [a[0],(A+B)/2]
        else:
            if a[0] > b[0]:
                a,b = b,a
            k = (a[1] - b[1])/(a[0] - b[0])
            bc = a[1] - a[0]*k
            A = a[0]
            B = b[0]
            F = (1 + sqrt(5))/2    
            while (B-A)/2 >= e:
                x1 = B - (B-A)/F
                y1 = x1*k + bc
                x2 = A + (B-A)/F
                y2 = x2*k + bc
                if f(x1, y1) >= f(x2, y2):
                    A = x1
                else:
                    B = x2
            return [(A+B)/2, ((A+B)/2)*k + bc]
    
    def __PLOT(self, kord, dot, func, color, n, dx, dy):
    
        a = np.array(kord[0])
        b = np.array(kord[1])
        c = np.array(kord[2])
        
        xe = []
        ye = []
        ze = []

        ans = {
        'AB' : True,
        'AC' : True,
        'BC' : True
        }

        if a[0] != b[0]:
            kAB = (a[1] - b[1])/(a[0] - b[0])
            bAB = a[1] - kAB*a[0]
            if kAB*c[0] + bAB < c[1]:
                ans['AB'] = "Yfalse"
            else:
                ans['AB'] = "Ytrue"
        else:        
            kAB = False
            bAB = False
            if a[0] < c[0]:
                ans['AB'] = "xfalse"
            else:
                ans['AB'] = "xtrue"

        if a[0] != c[0]:
            kAC = (a[1] - c[1])/(a[0] - c[0])
            bAC = a[1] - kAC*a[0]
            if kAC*b[0] + bAC < b[1]:
                ans['AC'] = "Yfalse"
            else:
                ans['AC'] = "Ytrue"
        else:
            kAC = False
            bAC = False
            if a[0] < b[0]:
                ans['AC'] = "xfalse"
            else:
                ans['AC'] = "xtrue"

        if c[0] != b[0]: 
            kBC = (c[1] - b[1])/(c[0] - b[0])
            bBC = c[1] - kBC*c[0]        
            if kBC*a[0] + bBC < a[1]:
                ans['BC'] = "Yfalse"
            else:
                ans['BC'] = "Ytrue"
        else:
            kBC = False
            bBC = False        
            if c[0] < a[0]:
                ans['BC'] = "xfalse"
            else:
                ans['BC'] = "xtrue"


        yMAX = np.max(kord[...,1])
        yMIN = np.min(kord[...,1])
        xMAX = np.max(kord[...,0])
        xMIN = np.min(kord[...,0])


        for X in np.linspace (xMIN, xMAX,n):
            for Y in np.linspace (yMIN, yMAX, n):
                if ans["AB"] == 'Ytrue' or ans["AB"] == 'Yfalse':
                    ytst = kAB*X + bAB
                    if ans["AB"] == "Ytrue":
                        if ytst < Y:
                            continue
                    elif ans["AB"] == "Yfalse":
                        if ytst > Y:
                            continue
                else:
                    if ans["AB"] == "xtrue":
                        if a[0] < X:
                            continue
                    elif ans["AB"] == "xfalse":
                        if a[0] > X:
                            continue
                            
                            
                if ans["BC"] == 'Ytrue' or ans["BC"] == 'Yfalse':
                    ytst = kBC*X + bBC
                    if ans["BC"] == "Ytrue":
                        if ytst < Y:
                            continue
                    elif ans["BC"] == "Yfalse":
                        if ytst > Y:
                            continue
                else:
                    if ans["BC"] == "xtrue":
                        if b[0] < X:
                            continue
                    elif ans["BC"] == "xfalse":
                        if b[0] > X:
                            continue  
                        
                if ans["AC"] == 'Ytrue' or ans["AC"] == 'Yfalse':
                    ytst = kAC*X + bAC

                    if ans["AC"] == "Ytrue":
                        if ytst < Y:
                            continue
                    elif ans["AC"] == "Yfalse":
                        if ytst > Y:
                            continue
                else:
                    if ans["AC"] == "xtrue":
                        if a[0] < X:
                            continue
                    elif ans["AC"] == "xfalse":
                        if a[0] > X:
                            continue 
                xe.append(X)
                ye.append(Y)
                ze.append(func(X,Y))

        fig = plt.figure()
        ax  = fig.add_subplot(111, projection = '3d')

        if color == True:
            cmap = cm.get_cmap('hot')
            normalize = clr.Normalize(vmin=min(ze), vmax=max(ze))
            colors = [cmap(normalize(value)) for value in ze]
        else:
            colors = None

        
        
        ax.scatter(xe, ye, ze, color = colors)
        ax.scatter(dot[0],dot[1], func(dot[0],dot[1]), color = 'b')
        plt.show()    
        
        plt.scatter(xe,ye, color = colors)
        plt.scatter(dot[0],dot[1], color = 'b')
        plt.plot(dx,dy, color = 'b')
        plt.show()

    def __grd (self, gr,kord): 
        return [float(gr[0].subs(sp.symbols('x'),kord[0])),float(gr[1].subs(sp.symbols('y'),kord[1])) ]

    def __norm(self, a):
        return sqrt(a[0]**2 + a[1]**2)

    def __step(self, a,b,c,func,opt,e,gr):
        AB = sqrt( (a[0] - b[0])**2 + (a[1] - b[1])**2 )
        AC = sqrt( (a[0] - c[0])**2 + (a[1] - c[1])**2 )
        BC = sqrt( (c[0] - b[0])**2 + (c[1] - b[1])**2 )
        g = 0
        if AC == max(AB,AC,BC):
            g = 1
            b,c = c,b
        if BC == max(AB,AC,BC):
            g = 2
            a,c = c,a
        AB = sqrt( (a[0] - b[0])**2 + (a[1] - b[1])**2 )
        AC = sqrt( (a[0] - c[0])**2 + (a[1] - c[1])**2 )
        BC = sqrt( (c[0] - b[0])**2 + (c[1] - b[1])**2 )
        
        n = [(a[0] + b[0])/2,(a[1]+b[1])/2]
        x = opt(n,c,e,func)
        grkord = self.__grd(gr,x)
        if grkord == [0,0]:
            return [None, None, n]
        gror = [x[0] + grkord[0],x[1] + grkord[1]]
        if c[0] != n[0]:
            k = (c[1] - n[1])/(c[0] - n[0])
            B = n[1] - n[0] * k
            if ((gror[0]*k + B)<gror[1] and (a[0] * k  + B)<a[1]) or ((gror[0]*k + B)>gror[1] and (a[0] * k  + B)>a[1]):
                return [b, c, n]
            else: 
                return [a, c, n]
        else:
            if ((gror[0]<c[0]) and (a[0] < c[0])) or ((gror[0]>c[0]) and (a[0] > c[0])) :
                return [b, c, n]
            else:
                return [a, c, n]

    def opt(self, kord, func, func_sp = None, m_const = None, l_const = None, E = 0.001, plot = False, info = False):

        """
        kord -- массив координат вершин треугольника;
        func -- функция, которую нужно минимизировать;
        func_sp -- func для SymPy, стоит использовать если заданая фунция содержит логарифмы и тригонометрические фунции;
        m_const -- константа Липшица для фунции;
        l_const -- константа Липшица для градиентов функции;
        E -- требуемая точность;
        plot -- построение графиков;
            первый элемент bool -- стоит ли делать заливку
            второй элемент int -- кучность точек
        info -- дополнительная информация.
        """

    
        a = kord[0]
        c = kord[1]
        b = kord[2]
        
        if info == True:
            start_time = time.time()
            t_meab = []
        
        if func_sp == None:
            func_sp = func(sp.symbols('x'),sp.symbols('y'))
            
        gr = [func_sp.diff(sp.symbols('x')),func_sp.diff(sp.symbols('y'))]
        
        if m_const == None:                
            gr_a = self.__grd(gr,kord[0])
            gr_b = self.__grd(gr,kord[1])
            gr_c = self.__grd(gr,kord[2])
            m_const = max(map(self.__norm,[gr_a,gr_b,gr_c]))
        if l_const == None:
            yMAX = np.max(kord[...,1])
            yMIN = np.min(kord[...,1])
            yl_max = max(gr[1].subs(sp.symbols('y'),yMAX),gr[1].subs(sp.symbols('y'),yMIN))
            xMAX = np.max(kord[...,0])
            xMIN = np.min(kord[...,0])
            xl_max = max(gr[0].subs(sp.symbols('x'),xMAX),gr[0].subs(sp.symbols('x'),xMIN))
            l_const = max([xl_max,yl_max])
        print (l_const,m_const)

        AB = sqrt( (a[0] - b[0])**2 + (a[1] - b[1])**2 )
        AC = sqrt( (a[0] - c[0])**2 + (a[1] - c[1])**2 )
        BC = sqrt( (c[0] - b[0])**2 + (c[1] - b[1])**2 )
        c_const = max(AB,AC,BC)

        n = ceil(log( (2*c_const*l_const)/E, 2/sqrt(3)))
        e = (E*l_const)/(4*(m_const*c_const + l_const)*n)
        
        resx,resy = [],[]
        if info == True:
            t1 = time.time()
        for i in range(n):
            t = self.__step(a,b,c, func,self.__GS, e, gr)
            a = t[0]
            b = t[1]
            c = t[2]
            resx.append(c[0])
            resy.append(c[1])
        if info == True:
            T = time.time()
            t2 = float(T - t1)/n               
            time_res = float(T - start_time)
        res  = {}
        res['кординаты точки'] = c
        res['значение функции'] = func(c[0],c[1])
        if info == True:   
            res['суммарное время'] = time_res
            res['среднее время итерации'] = t2
            res['константа L'] = l_const
            res['константа M'] = m_const
            res['колличество шагов'] = n
            res['точность одномерных задач'] = e
            
        if plot != False:
            self.__PLOT(kord, c, func, plot[0], plot[1], resx, resy)
        return res

