# -*- coding: utf-8 -*-

import numpy as np

def make_shekel(lo, hi, m):
    """
    To look at different shekel functions:
        >>> f, a = plt.subplots()
        >>> lo, hi = 0, 10
        >>> step = (hi-lo)/100
        >>> xm, ym = np.mgrid[lo:hi:step, lo:hi:step]
        >>> def gen(lines=20, maxnum=10):
        ...:    while 1:
        ...:        s = make_shekel(lo, hi, maxnum)
        ...:        a.clear()
        ...:        a.contour(xm, ym, s(xm, ym), lines)
        ...:        yield
        ...:
        >>> g = gen()
        >>> next(g)
    """
    a = np.random.uniform(lo, hi, (m, 2))
    c = np.random.uniform(lo, hi, (m,))
    def shekel(x, y):
        return -np.sum( np.array([c[i]+(x-a[i,0])**2+(y-a[i,1])**2
                                  for i in range(m)])**(-1), axis=0 )
    return shekel


# global minimum: f(0, 0) = 0
# domain: -5.12 <= x(i) <= 5.12
def rastrigin(x, y):
    arr = np.array([x, y])
    return 20 + np.sum(arr**2 - 10*np.cos(2*np.pi*arr), axis=0)

# global minimum: f(0, 0) = 0
# domain: -5 <= x(i) <= 5
def ackley(x, y):
    arr = np.array([x, y])
    return -20 * np.exp(-0.2 * np.sqrt(np.sum(arr**2, axis=0)/2) ) \
            -np.exp(np.sum(np.cos(2*np.pi*arr), axis=0)/2 ) + np.exp(1) + 20

# global minimum: f(1, 1) = 0
# domain: -INF <= x(i) <= INF
def rosenbrock(x, y):
    return (1 - x)**2 + 100*(y - x**2)**2

# global minimum: f(0, -1) = 3
# domain: -2 <= x(i) <= 2
def goldstein_prince(x, y):
    return (1 + (x + y + 1)**2 * (19 - 14*x + 3*x**2 - 14*y + 6*x*y + 3*y**2))\
        *(30 + (2*x - 3*y)**2 * (18 - 32*x + 12*x**2 + 48*y - 36*x*y + 27*y**2))

# global minimum: f(-10, 1) = 0
# domain: -15 <= x <= -5, -3 <= y <= 3
def bukin_n6(x, y):
    return 100*np.sqrt(np.abs(y-0.01*x**2)) + 0.01*np.abs(x+10)

# global minimum: f(1, 1) = 0
# domain: -10 <= x(i) <= 10
def levi_n13(x, y):
    return np.sin(3*np.pi*x)**2 + (x-1)**2*(1+np.sin(3*np.pi*y)) \
            + (y-1)**2*(1+np.sin(2*np.pi*y))

# global minimum:
#   f(3, 2)                 \
#   f(-2.805118, 3.131312)  |   = 0
#   f(-3.779310, -3.283186) |
#   f(3.584428, -1.848126)  /
#
# domain: -5 <= x(i) <= 5
def himmelblau(x, y):
    return (x**2+y-11)**2 + (x+y**2-7)**2

# global minimum: f(np.pi, np.pi) = -1
# domain: -100 <= x(i) <= 100
def easom(x, y):
    arr = np.array([x, y])
    return -np.prod(np.cos(arr), axis=0) * np.exp(-np.sum((arr-np.pi)**2, axis=0))

# global minimum:
#     f(1.34941, -1.34941)  \
#     f(1.34941, 1.34941)   |   = -2.06261
#     f(-1.34941, 1.34941)  |
#     f(-1.34941, -1.34941) /
#
# domain: -10 <= x(i) <= 10
def cross_in_tray(x, y):
    arr = np.array([x, y])
    return -0.0001*( np.abs(np.prod(np.sin(arr), axis=0) \
        * np.exp(np.abs(100-np.sqrt(np.sum(arr**2, axis=0))/np.pi)))+1 )**(0.1)

# global minimum: f(512, 404.2319) = -959.6407
# domain: -512 <= x(i) <= 512
def eggholder(x, y):
    return -(y+47) * np.sin(np.sqrt(np.abs(x/2+y+47)))   \
        - x * np.sin(np.sqrt(np.abs(x-y+47)))

# global minimum:
#     f(8.05502, 9.66459)   \
#     f(-8.05502, 9.66459)  |   = -19.2085
#     f(8.05502, -9.66459)  |
#     f(-8.05502, -9.66459) /
#
# domain: -10 <= x(i) <= 10
def hoelder_table(x, y):
    return -np.abs(np.sin(x)*np.cos(y)*np.exp(np.abs(1-np.sqrt(np.sum(x**2,
                                                                      axis=0))/np.pi)))

# global minimum: -78.33234 < f(-2.903534, -2.903534) < -78.33232
# domain: -5 <= x(i) <= 5
def styblinski_tang(x, y):
    arr = np.array([x, y])
    return np.sum(arr**4 - 16*arr**2 + 5*arr, axis=0) / 2

# TODO global minimum
# domain: -3 <= x <= 3, -2 <= y <= 2
def tal(x, y):
    return -np.exp(-x**2-10.*y*(y-1.))*(1.+.8*np.sin(x*y))*np.cos(1.3*y+2.5*x)


test_functions = {
    'ackley'            :   (ackley,            [-5, 5]),
    'bukin-n6'          :   (bukin_n6,          [-15, -5], [-3, 3]),
    'cross-in-tray'     :   (cross_in_tray,     [-10, 10]),
    'easom'             :   (easom,             [-100, 100]),
    'eggholder'         :   (eggholder,         [-512, 512]),
    'goldstein-prince'  :   (goldstein_prince,  [-2, 2]),
    'himmelblau'        :   (himmelblau,        [-5, 5]),
    'hoelder-table'     :   (hoelder_table,     [-8, 8]),
    'levi-n13'          :   (levi_n13,          [-10, 10]),
    'rastrigin'         :   (rastrigin,         [-5.12, 5.12]),
    'rosenbrock'        :   (rosenbrock,        # [-100, 100]),
                             [-2, 2], [-3, 5]),
    'styblinski-tang'   :   (styblinski_tang,   [-5, 5]),
    'tal'               :   (tal,               [-3, 3], [-2, 2]),
    'shekel'            :   (make_shekel,       [-5, 5])
}
