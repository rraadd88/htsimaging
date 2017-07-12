import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import leastsq
# import numpy.random as npr
# import matplotlib.pyplot as plt
import matplotlib.pyplot as plt


def power(x, A, B):
    """power law equation."""
    return B*x**A
def power_residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    A,B = p
    err = y-power(x, A, B)
    return err
def power_peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B = p
    return power(x, A, B)

def line(x, m, C):
    """power law equation."""
    return m*x+C
def line_residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    m,C = p
    err = y-line(x, m,C)
    return err
def line_peval(x, p):
    """Evaluated value at x with current parameters."""
    m,C = p
    return line(x, m,C)

def logistic4(x, A, B, C, D):
    return ((A-D)/(1.0+((x/C)**B))) + D
def logistic4_residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D = p
    err = y-logistic4(x, A, B, C, D)
    return err
def logistic4_peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C,D = p
    return logistic4(x, A, B, C, D)

def logistic5(x, A, B, C, D, E):
    return  D+(A-D)/((1+(x/C)**B)**E)
def logistic5_residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""
    A,B,C,D,E = p
    err = y-logistic5(x, A, B, C, D, E)
    return err
def logistic5_peval(x, p):
    """Evaluated value at x with current parameters."""
    A,B,C,D,E = p
    return logistic5(x, A, B, C, D, E)

def fit_power(x,y,p0=[0, 1],plot=False):
    plsq,cov,infodict,mesg,ier= leastsq(power_residuals, p0, args=(y, x),full_output=1)
    ss_err=(infodict['fvec']**2).sum()
    ss_tot=((y-y.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    # if rsquared<0.85:
        # print "poor fit %d" % rsquared
    if plot==True:
        plt.plot(x,y,'o')  
        plt.plot(x,power_peval(x,plsq)) 
    return plsq[0], plsq[1],rsquared,power_peval(np.max(x),plsq)

def fit_line(x,y,p0=[0, 1],plot=False):
    plsq,cov,infodict,mesg,ier= leastsq(line_residuals, p0, args=(y, x),full_output=1)
    ss_err=(infodict['fvec']**2).sum()
    ss_tot=((y-y.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    # if rsquared<0.85:
    #     print "poor fit %d" % rsquared
    if plot==True:
        plt.plot(x,y,'o')  
        plt.plot(x,line_peval(x,plsq)) 
    return plsq[0], plsq[1],rsquared
