# -*- coding: utf-8 -*-
'''function used in comfortmodels.py'''

STATIC_URL = '/static'

def bisect(a, b, fn, epsilon, target):
    '''bisect function'''

    while abs(b - a) > 2 * epsilon:
        midpoint = (b + a) / 2
        a_T = fn(a)
        b_T = fn(b)
        midpoint_T = fn(midpoint)
        if (a_T - target) * (midpoint_T - target) < 0:
            b = midpoint
        elif (b_T - target) * (midpoint_T - target) < 0:
            a = midpoint
        else:
            return -99
    return midpoint

def secant(a, b, fn, epsilon):
    '''secant # root-finding only'''

    f1 = fn(a)
    if abs(f1) <= epsilon:
        return a
    f2 = fn(b)

    if abs(f2) <= epsilon:
        return b
    for i in range(0, 100):
        slope = (f2 - f1) / (b - a)
        c = b - f2/slope
        f3 = fn(c)
        if abs(f3) < epsilon:
            return c
        a = b
        b = c
        f1 = f2
        f2 = f3
    return None

def getSensation(pmv):
    '''sensation from pmv'''
    if pmv < -2.5:
        return 'Cold'
    elif pmv < -1.5:
        return 'Cool'
    elif pmv < -0.5:
        return 'Slightly Cool'
    elif pmv < 0.5:
        return 'Neutral'
    elif pmv < 1.5:
        return 'Slightly Warm'
    elif pmv < 2.5:
        return 'Warm'
    return 'Hot'


def CtoF(x):
    '''convert celcius to farenheit'''
    return x * 9 / 5 + 32


def FtoC(x):
    '''fahrenheit to celcius'''
    return (x - 32) * 5 / 9


def secant_solve(f, x1, x2, ftol, xtol):
    '''secant solve'''
    f1 = f(x1)
    if abs(f1) <= ftol:
        return x1        # already effectively zero
    f2 = f(x2)
    if abs(f2) <= ftol:
        return x2        # already effectively zero
    while abs(x2 - x1) > xtol:
        slope = (f2 - f1)/(x2 - x1)
        if slope == 0:
            pass
            #sys.stderr.write("Division by 0 due to vanishing slope - exit!\n")
          #sys.exit(1)
        x3 = x2 - f2/slope               # the new approximate zero
        f3 = f(x3)                       # and its function value
        if abs(f3) <= ftol:
            break
        x1, f1 = x2, f2                    # copy x2,f2 to x1,f1
        x2, f2 = x3, f3                    # copy x3,f3 to x2,f2
    return x3
