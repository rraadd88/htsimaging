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

from roux.lib.io_dfs import set_index
from os.path import basename
import matplotlib.cm as cm
def plot_kin(traj,dparams=None,ctime='time (s)',
             fit_eqn=None,
             smp_ctrl=None,
             plot_parami=1,
             ylabel='',
             ymargin=0.2,
             label='',
             color='b',
             ax1=None,
             plot_fh=None):
    traj_mean=set_index(traj,col_index=ctime)
    traj_mean=traj.T.mean()
    traj_std=traj.T.std()
    if fit_eqn=="logistic4":
        params=["Minimum asymptote","Hill's slope","Inflection point","Maximum asymptote"]
    elif fit_eqn=="logistic5":
        params=["Minimum asymptote","Hill's slope","Inflection point","Maximum asymptote","Asymmetry factor"]
    elif fit_eqn=='power':
        params=["power law exponent","amplitude"]
        
    if ax1 is None:
        plt.figure(figsize=(8,3))
        ax1=plt.subplot(131)
    ax1.set_xlabel(ctime)
    ax1.set_ylabel(ylabel)
    ax1.plot(traj_mean.index, traj_mean,'o',zorder=2, 
             markerfacecolor='none',markeredgecolor=color,alpha=0.5,
#             label='sd',
            )
    ax1.legend(label,loc='center', 
               bbox_to_anchor=(1.63, 0.5))
    ax1.fill_between(traj_std.index, \
                     traj_mean-traj_std, \
                     traj_mean+traj_std, \
                     color=color, alpha=0.2,zorder=0)
    
    if not dparams is None:
        params_mean=dparams.mean()
        params_std=dparams.std()
        p0=params_mean
    else:
        p0=np.zeros(len(params))
    if fit_eqn=="logistic4":
        plsq = leastsq(logistic4_residuals, p0, args=(traj_mean, traj_mean.index))
        ax1.plot(traj_mean.index.tolist(),logistic4_peval(traj_mean.index,plsq[0]),label=label,color=color,zorder=2)        
    elif fit_eqn=="logistic5":
        plsq = leastsq(logistic5_residuals, p0, args=(traj_mean, traj_mean.index))
        ax1.plot(traj_mean.index.tolist(),logistic5_peval(traj_mean.index,plsq[0]),label=label,color=color,zorder=2)        
    elif fit_eqn=="power":
        plsq = leastsq(power_residuals, p0, args=(traj_mean, traj_mean.index))
        ax1.plot(traj_mean.index.tolist(),power_peval(traj_mean.index,plsq[0]),label=label,color=color,zorder=2,lw=2)        
                
    ax1.set_ylim([0,ax1.get_ylim()[1]])
#     ax2=plt.subplot(133)
#     plot_param=params[plot_parami]
#     X=params_mean.loc[smps_plot,plot_param]
#     Y=range(len(X))
#     X_err=params_std.loc[smps_plot,plot_param]
#     for i in range(len(X)):
#         ax2.errorbar(X[i],Y[i],xerr=X_err[i],fmt="none",ecolor='black', capthick=2,zorder=0)
#         ax2.plot(X[i],Y[i],marker='o', linestyle='', ms=10,zorder=1)    
    
# #     ax2.errorbar(X,Y,xerr=X_err,fmt="none",ecolor='black', capthick=2,zorder=0)
# #     ax2.plot(X,Y,marker='o', linestyle='', ms=10,c='blue',alpha=0.5,zorder=1)    
# #     ax2.set_ylim(ax2.get_ylim()[::-1])
#     ax2.set_yticklabels([])
#     ax2.set_ylim([len(Y)-1+len(Y)*ymargin,0-len(Y)*ymargin])
# #     ax2.set_ylim([np.max(X_err)-1+len(Y)*ymargin,0-len(Y)*ymargin])
#     ax2.set_xlabel(plot_param)    

    plt.tight_layout()
    if plot_fh!=None:
        plt.savefig(plot_fh, format='pdf')
    return ax1

#from htsimaging.lib.io_strs import get_dir
from os.path import dirname
from htsimaging.lib.fit_kin import plot_kin
import matplotlib.pyplot as plt
import matplotlib.cm as cm
def plot_kin_all(expt_dh,imsd_fhs):
    fig=plt.figure(figsize=(8,3))
    colors=cm.rainbow(np.linspace(0,1,len(imsd_fhs)))
    ylimmx=[]
    for i,imsd_fh in enumerate(imsd_fhs):
        imsd_flt_fh='%s.imsd_flt' % imsd_fh
        imsd_flt=pd.read_csv(imsd_flt_fh).set_index('lag time [s]')
        params_flt_fh='%s.params_flt' % imsd_fh
        params_flt=pd.read_csv(params_flt_fh)
        ax1=plt.subplot(131)
        ax1=plot_kin(imsd_flt,
                 params_flt.loc[:,['power law exponent','power law constant']],
                 ctime='lag $t$',
                 fit_eqn='power',
                 smp_ctrl=None,
                 ylabel='MSD ($\mu m^{2}$)',
                 color=colors[i],
                 ymargin=0.2,
                 ax1=ax1,
                 label=basename(imsd_flt_fh).split('_replicate', 1)[0].replace('_',' '),
                 plot_fh=None)
        ylimmx.append(np.max(imsd_flt.max()))
    ax1.set_ylim([0,np.max(ylimmx)])
    ax1.legend(bbox_to_anchor=[1,1,0,0],loc=2)
    fig.savefig('%s/plot_flt_%s.pdf' % (expt_dh,dirname(expt_dh)))
