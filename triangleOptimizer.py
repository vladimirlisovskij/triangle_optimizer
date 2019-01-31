import array as arr
import time
from math import sqrt, ceil, log
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import mpl_toolkits.mplot3d
from matplotlib.cm import hot
import numpy as np
import sympy as sp


class TriOpt:
    @staticmethod
    def __GS(a, b, e, fun):
        """
        Золотое сечение

        a,b -- границы отрезка;
        e -- точность;
        fun -- функция.
        """
        if a[0] == b[0]:

            if a[1] > b[1]:
                a, b = b, a

            A = a[1]
            B = b[1]

            F = (1 + sqrt(5))/2

            while (B-A)/2 >= e:
                y1 = B - (B-A)/F
                y2 = A + (B-A)/F

                if fun(a[0], y1) >= fun(a[0], y2):
                    A = y1
                else:
                    B = y2

            return [a[0], (A+B)/2]
        else:

            if a[0] > b[0]:
                a, b = b, a

            k_coef = (a[1] - b[1])/(a[0] - b[0])
            b_coef = a[1] - a[0]*k_coef

            A = a[0]
            B = b[0]

            F = (1 + sqrt(5))/2

            while (B-A)/2 >= e:
                x1 = B - (B-A)/F
                y1 = x1*k_coef + b_coef
                x2 = A + (B-A)/F
                y2 = x2*k_coef + b_coef

                if fun(x1, y1) >= fun(x2, y2):
                    A = x1
                else:
                    B = x2

            return [(A+B)/2, ((A+B)/2)*k_coef + b_coef]

    @staticmethod
    def __PLOT(koor, dot, fun, color, n, dx, dy):
        """
        Построение треугольников

        koor -- координаты вершин;
        fun -- функция;
        color -- стоит ли делать заливку;
        n -- плотность точек;
        x_dots,y_dots -- результаты итераций.
        """

        a = koor[0]
        b = koor[1]
        c = koor[2]

        # точки треугольника
        xe = arr.array('f')
        ye = arr.array('f')
        ze = arr.array('f')

        # система неравенств
        ans = arr.array('i', [0, 0, 0])

        """
        Точки лежат ниже прямой -- 1;
        точки лежат выше прямой  -- 2;
        точки лежат левее прямой  -- 3;
        точки лежат правее прямой -- 4.
        
        AB -- ans[0];
        AC -- ans[1];
        BC -- ans[2].
        """

        if a[0] != b[0]:
            kAB = (a[1] - b[1])/(a[0] - b[0])
            bAB = a[1] - kAB*a[0]

            if kAB*c[0] + bAB < c[1]:
                ans[0] = 2
            else:
                ans[0] = 1
        else:
            if a[0] < c[0]:
                ans[0] = 4
            else:
                ans[0] = 3

        if a[0] != c[0]:
            kAC = (a[1] - c[1])/(a[0] - c[0])
            bAC = a[1] - kAC*a[0]

            if kAC*b[0] + bAC < b[1]:
                ans[1] = 2
            else:
                ans[1] = 1
        else:
            if a[0] < b[0]:
                ans[1] = 4
            else:
                ans[1] = 3

        if c[0] != b[0]:
            kBC = (c[1] - b[1])/(c[0] - b[0])
            bBC = c[1] - kBC*c[0]

            if kBC*a[0] + bBC < a[1]:
                ans[2] = 2
            else:
                ans[2] = 1
        else:
            if c[0] < a[0]:
                ans[2] = 4
            else:
                ans[2] = 3

        yMAX = np.max(koor[..., 1])
        yMIN = np.min(koor[..., 1])
        xMAX = np.max(koor[..., 0])
        xMIN = np.min(koor[..., 0])

        for X in np.linspace(xMIN, xMAX, n):
            for Y in np.linspace(yMIN, yMAX, n):
                if ans[0] == 1 or ans[0] == 2:
                    ytst = kAB*X + bAB

                    if ans[0] == 1:
                        if ytst < Y:
                            continue
                    elif ans[0] == 2:
                        if ytst > Y:
                            continue
                else:
                    if ans[0] == 3:
                        if a[0] < X:
                            continue
                    elif ans[0] == 4:
                        if a[0] > X:
                            continue

                if ans[2] == 1 or ans[2] == 2:
                    ytst = kBC*X + bBC

                    if ans[2] == 1:
                        if ytst < Y:
                            continue
                    elif ans[2] == 2:
                        if ytst > Y:
                            continue
                else:
                    if ans[2] == 3:
                        if b[0] < X:
                            continue
                    elif ans[2] == 4:
                        if b[0] > X:
                            continue

                if ans[1] == 1 or ans[1] == 2:
                    ytst = kAC*X + bAC

                    if ans[1] == 1:
                        if ytst < Y:
                            continue
                    elif ans[1] == 2:
                        if ytst > Y:
                            continue
                else:
                    if ans[1] == 3:
                        if a[0] < X:
                            continue
                    elif ans[1] == 4:
                        if a[0] > X:
                            continue

                xe.append(X)
                ye.append(Y)
                ze.append(fun(X, Y))

        # построение
        fig = plt.figure()
        ax = fig.add_subplot(121, projection='3d')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('f(X, Y)')
        ax.plot_trisurf(xe, ye, ze, cmap=hot) if color else ax.plot_trisurf(xe, ye, ze, cmap=hot)
        ax.scatter(dot[0], dot[1], fun(dot[0], dot[1]), color='b')

        x = tri.Triangulation(xe, ye)
        ax2 = fig.add_subplot(122)
        ax2.tricontourf(x, ze, cmap=hot) if color else ax2.fill(koor[..., 0], koor[..., 1])
        ax2.scatter(dot[0], dot[1], color='b')
        ax2.plot(dx, dy, color='b')

        plt.xlabel("x")
        plt.ylabel("y")

        plt.show()

    @staticmethod
    def __grd(gr, koor):
        """
        Градиент в точке

        gr -- формула градиента;
        koor -- координата точкию
        """
        return float(gr[0].subs(sp.symbols('x'), koor[0])), float(gr[1].subs(sp.symbols('y'), koor[1]))

    @staticmethod
    def __norm(a):
        """
        Норма

        a -- вектор.
        """
        return sqrt(a[0]**2 + a[1]**2)

    @staticmethod
    def __step(a, b, c, fun, e, gr):
        """
        Шаг алгоритма

        a,b,c -- вершины треугольника;
        fun -- функция;
        e -- точность;
        gr -- формула градиента.
        """
        AB = sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
        AC = sqrt((a[0] - c[0])**2 + (a[1] - c[1])**2)
        BC = sqrt((c[0] - b[0])**2 + (c[1] - b[1])**2)

        # приводим к удобному виду
        if AC == max(AB, AC, BC):
            b, c = c, b
        if BC == max(AB, AC, BC):
            a, c = c, a

        # вычислем точку для градиента
        n = ((a[0] + b[0])/2, (a[1]+b[1])/2)
        x = TriOpt.__GS(n, c, e, fun)

        # вычисляем градиент
        grkord = TriOpt.__grd(gr, x)

        if grkord == (0, 0):
            return n,

        # определяем направление
        gror = (x[0] + grkord[0], x[1] + grkord[1])

        if c[0] != n[0]:
            k = (c[1] - n[1])/(c[0] - n[0])
            B = n[1] - n[0] * k

            if ((gror[0]*k + B) < gror[1] and (a[0] * k + B) < a[1]) or ((gror[0]*k + B) > gror[1] and
                                                                         (a[0] * k + B) > a[1]):
                return b, c, n
            else:
                return a, c, n
        else:
            if ((gror[0] < c[0]) and (a[0] < c[0])) or ((gror[0] > c[0]) and (a[0] > c[0])):
                return b, c, n
            else:
                return a, c, n

    @staticmethod
    def opt(koor: np.ndarray, func: 'function for minimize', func_sp: 'function to sympy' = None, m_const: float = None,
            l_const: float = None, E: float = 0.001, plot: tuple = False, info: bool = False, n: int = None) -> dict:
        """
        Метод для минимизации функции

        in:
            koor [np.ndarray 3x2] -- массив координат вершин треугольника;
            func -- функция, которую нужно минимизировать;
            func_sp -- func для SymPy;
            m_const [float] -- константа Липшица для фунции;
            l_const [float]  -- константа Липшица для градиентов функции;
            E [float] -- требуемая точность;
            plot [tuple 2x1] -- построение графиков;
                первый элемент [bool] -- стоит ли делать заливку;
                второй элемент [int] -- кучность точек;
            n [int] -- колличество шагов;
            info [float] -- дополнительная информация.
        out:
            [dict] -- кооординаты точки, значение функции.
            Если был передан параметр info = True:
            время работы, среднее время итерации, значение констант Липшица, колличество шагов,
            точность решения одномерных задач.
        """
        a = koor[0]
        c = koor[1]
        b = koor[2]

        if info:
            start_time = time.time()  # засекаем время работы

        if func_sp is None:
            func_sp = func(sp.symbols('x'), sp.symbols('y'))  # создаем func для SymPy
        else:
            func_sp = func_sp(sp.symbols('x'), sp.symbols('y'))

        gr = [func_sp.diff(sp.symbols('x')), func_sp.diff(sp.symbols('y'))]  # формула для градиента

        # вычисляем константу m
        if m_const is None:
            gr_a = TriOpt.__grd(gr, koor[0])
            gr_b = TriOpt.__grd(gr, koor[1])
            gr_c = TriOpt.__grd(gr, koor[2])
            m_const = max(map(TriOpt.__norm, [gr_a, gr_b, gr_c]))
        # вычисляем константу l
        if l_const is None and n is None:
            yMAX = max(koor[..., 1])
            yMIN = min(koor[..., 1])
            yl_max = max(float(gr[1].subs(sp.symbols('y'), yMAX)), float(gr[1].subs(sp.symbols('y'), yMIN)))
            xMAX = max(koor[..., 0])
            xMIN = min(koor[..., 0])
            xl_max = max(float(gr[0].subs(sp.symbols('x'), xMAX)), float(gr[0].subs(sp.symbols('x'), xMIN)))
            l_const = max([xl_max, yl_max])

        # вычисляем диаметр
        AB = sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)
        AC = sqrt((a[0] - c[0])**2 + (a[1] - c[1])**2)
        BC = sqrt((c[0] - b[0])**2 + (c[1] - b[1])**2)
        c_const = max(AB, AC, BC)

        # вычисляем колличество шагов и точность
        if n is None:
            n = ceil(log((2*c_const*l_const)/E, 2/sqrt(3)))
        e = E/(m_const*c_const*n)

        if plot is not False:
            resx, resy = arr.array('f', []), arr.array('f', [])
        # засекаем время итераций
        if info:
            t1 = time.time()
        # итерации
        for _ in range(n):
            t = TriOpt.__step(a, b, c, func, e, gr)
            if len(t) != 1:
                a = t[0]
                b = t[1]
                c = t[2]
                if plot is not False:
                    resx.append(c[0])
                    resy.append(c[1])
            else:
                if plot is not False:
                    resx.append(t[0])
                    resy.append(t[1])
                break
        # вычисление времени работы и итераций
        if info:
            T = time.time()
            t2 = float(T - t1)/n
            time_res = float(T - start_time)

        # результат
        res = {}
        res['point coordinates'] = c
        res['function value'] = func(c[0], c[1])

        if info:
            res['total time'] = time_res
            res['average iteration time'] = t2
            if l_const:
                res['constant L'] = l_const
            res['constant M'] = m_const
            res['number of steps'] = n
            res['accuracy of one-dimensional tasks'] = e

        # построение графика
        if plot is not False:
            TriOpt.__PLOT(koor, c, func, plot[0], plot[1], resx, resy)
        return res
