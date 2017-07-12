from os.path import splitext, join, exists, isdir,basename,abspath,dirname
from os import makedirs
import trackpy as tp
# import nd2reader
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pims
import pims_nd2
from glob import glob
import logging
import pandas as pd

# %matplotlib inline

@pims.pipeline
def average_z(image):
    return image.mean(axis=0) 
    
def plot_msd(imsd,emsd,ax=None,scale="log",plot_fh=None):
    if ax is None:
        plt.figure(figsize=(3, 3))
        ax=plt.subplot(111)

    imsd.plot(x=imsd.index,legend=None,alpha=0.75,ax=ax)
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')

    emsd.plot(x=emsd.index,style='o',legend=None,ax=ax)
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    plt.tight_layout()
    if not plot_fh is None: 
        ax.figure.savefig(plot_fh,format='pdf');plt.clf();plt.close()
    return ax

def plot_emsd(expt_data,ax=None,color='k',scale="log",plot_fh=None):
    if ax is None:
        plt.figure(figsize=(6, 3))
        ax=plt.subplot(121)
    expt_data.plot(x=expt_data.index,style='o',c=color,ax=ax,markeredgecolor = 'none')
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
                xlabel='lag time $t$')
    ax.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5)) #bbox_to_anchor=(2.2, 1.0)
    plt.tight_layout()
    if not plot_fh is None: 
        ax.figure.savefig(plot_fh,format='pdf');plt.clf();plt.close()
    return ax 

def nd2msd(nd_fh):
    frames=pims.ND2_Reader(nd_fh)
    if len(np.shape(frames))==4:
        frames = average_z(frames)
    f_batch = tp.batch(frames,diameter=11,threshold=np.percentile(frames,75))

    t = tp.link_df(f_batch, search_range=11, memory=3)
    t_flt = tp.filter_stubs(t, 3*int(len(frames)/4))
    d = tp.compute_drift(t_flt)
    t_cor = tp.subtract_drift(t_flt, d)
    imsd=tp.imsd(t_cor,mpp=0.0645,fps=0.2, max_lagtime=100, statistic='msd')
    emsd=tp.emsd(t_cor,mpp=0.0645,fps=0.2, max_lagtime=100)
    return imsd,emsd

def expt_dh2expt_info(expt_dh):
    try:
        expt_info=pd.read_csv(expt_dh+"info")
    except:
        expt_info=pd.read_csv(expt_dh+"info.csv")
    if "Unnamed" in expt_info.columns.tolist():
        cols_del.append(col)

    expt_info2=expt_info.set_index("smpn",drop=True).T
    expt_info2.columns=expt_info.loc[:,"smpn"]
    row_2_keep=[]
    for row in list(expt_info2.index):
        if "replicate" in row:
            row_2_keep.append(row)

    expt_info2=expt_info2.loc[row_2_keep,:]
    if len(expt_info2.index) != len(set(expt_info2.index)):
        print "duplicate values in index"
    for col in expt_info2.columns.tolist():
        for i in range(len(expt_info2)):
            if not pd.isnull(expt_info2.loc["replicate %d" % (i+1),col]):
                expt_info2.loc[("replicate %d" % (i+1)),col]=['replicate %d' % (i+1),expt_info2.loc[("replicate %d" % (i+1)),col]]
            else:
                expt_info2.loc[("replicate %d" % (i+1)),col]=[np.nan,np.nan]
    expt_info2=expt_info2.reset_index()
    del expt_info2["index"]
    return expt_info2

def expt2plots(expt_info,expt_dh):
    if not exists(expt_dh):
        makedirs(expt_dh)
    expt_data=pd.DataFrame()
    for test in expt_info:
        test_info=expt_info.loc[:,test]
        test_data=pd.DataFrame()    
        for rep in test_info:
            if not pd.isnull(rep[0]):
                repn="%s %s" % (test,rep[0])
    #             print repn
                nd_fh=rep[1]
                if exists(nd_fh):
                    out_fh="%s/%s" % (expt_dh,repn.replace(" ","_"))
    #                 print out_fh
                    plot_fh=out_fh+".imsd.pdf"
                    if not exists(plot_fh):
                        try:
                            imsd,emsd=nd2msd(nd_fh)
                        except:
                            continue
                        emsd=pd.DataFrame(emsd)
                        print repn
                        emsd.columns=[repn]
                        imsd.index=emsd.index
                        print imsd.head()
                        print out_fh+".imsd"
                        imsd.to_csv(out_fh+".imsd")
                        emsd.to_csv(out_fh+".emsd")
                        plot_msd(imsd,emsd,ax_imsd,scale='',plot_fh=plot_fh)
                    else:
                        emsd=pd.read_csv(out_fh+".emsd")
                        if "Unnamed: 0" in emsd.columns.tolist(): 
                            del emsd["Unnamed: 0"]
                else:
                    print "can not find"
                if len(test_data)==0:
                    test_data=emsd
                else:
        #                 test_data[repn]=emsd[repn]
                    test_data=pd.concat([test_data,emsd[repn]],axis=1)
            test_data.to_csv(expt_dh+test+".emsd")
        if len(expt_data)==0:
            expt_data=test_data
        else:
    #             expt_data.loc[:,test_data.columns.tolist()]=test_data.loc[:,test_data.columns.tolist()]
            expt_data=pd.concat([expt_data,test_data.loc[:,[col for col in test_data.columns.tolist() if col != "lagt"]]],axis=1)
    #     expt_data=expt_data.drop_duplicates("lagt",keep='first')
    expt_data.to_csv(expt_dh+"expt.emsd")
    plot_fh=expt_dh+"emsd.pdf"
    plot_emsd(expt_data, plot_fh=plot_fh)
    expt_data_log10=np.log10(expt_data)
    expt_data_log10.to_csv(expt_dh+"expt.emsdlog10")
    expt_data=expt_data.set_index('lagt',drop=True)
    return expt_data


def createinfo(expt_dh):
    try:
        info=pd.read_csv(expt_dh+"info")
    except:
        info=pd.read_csv(expt_dh+"info.csv")
    for i in range(len(info)):
        reps=glob("%s/%s*" % (info.loc[i,"dh"],info.loc[i,"fn_lead"]))
        for repi,rep in enumerate(reps):
            info.loc[i,"replicate %d" % (repi+1)]=rep
    #     break
    info.to_csv(expt_dh+"info",index=False)
    
    
from scipy.optimize import leastsq

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
def fit_power(x,y,p0=[0, 1],plot=False):
    plsq,cov,infodict,mesg,ier= leastsq(power_residuals, p0, args=(y, x),full_output=1)
    ss_err=(infodict['fvec']**2).sum()
    ss_tot=((y-y.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    if rsquared<0.85:
        print "poor fit %d" % rsquared
    if plot==True:
        plt.plot(x,y,'o')  
        plt.plot(x,power_peval(x,plsq)) 
    return plsq[0], plsq[1],rsquared,power_peval(np.max(x),plsq)
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
def fit_line(x,y,p0=[0, 1],plot=False):
    plsq,cov,infodict,mesg,ier= leastsq(line_residuals, p0, args=(y, x),full_output=1)
    ss_err=(infodict['fvec']**2).sum()
    ss_tot=((y-y.mean())**2).sum()
    rsquared=1-(ss_err/ss_tot)
    if rsquared<0.85:
        print "poor fit %d" % rsquared
    if plot==True:
        plt.plot(x,y,'o')  
        plt.plot(x,line_peval(x,plsq)) 
    return plsq[0], plsq[1],rsquared
