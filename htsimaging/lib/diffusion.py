#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

from os.path import splitext, join, exists, isdir,basename,abspath,dirname
from os import makedirs
import trackpy as tp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pims
# import pims_nd2
# import nd2reader
from glob import glob
from htsimaging.lib.io_nd_files import average_z
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # 

def plot_msd(imsd,emsd,ax):
    imsd.reset_index().plot(x="lag time [s]",legend=None,alpha=0.75,ax=ax)
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    ax.set_xscale('log')
    ax.set_yscale('log')

    emsd.plot(x="lagt",style='o',legend=None,ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    plt.tight_layout()
    return ax

def plot_emsd(expt_data,ax_emsd):
    # print expt_data.head()
    expt_data.plot(x="lagt",style='o',ax=ax_emsd)
    ax_emsd.set_xscale('log')
    ax_emsd.set_yscale('log')
    ax_emsd.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
                xlabel='lag time $t$')
    ax_emsd.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5)) #bbox_to_anchor=(2.2, 1.0)
    plt.tight_layout()

def nd2msd(nd_fh):
    # print nd_fh
    frames=pims.ND2_Reader(nd_fh)
    logging.info('number of frames = %d' % len(np.shape(frames)))
    if len(np.shape(frames))==4:
        frames = average_z(frames)
    threshold=np.percentile(frames,75)
    f_batch = tp.batch(frames,diameter=11,threshold=threshold)

    t = tp.link_df(f_batch, search_range=11, memory=3)
    t_flt = tp.filter_stubs(t, 3*int(len(frames)/4))
    try:
        d = tp.compute_drift(t_flt)
        t_cor = tp.subtract_drift(t_flt, d)
    except:
        t_cor=t_flt
        logging.info("drift correction excepted")    
    # plt.figure()
    # tp.plot_traj(t_flt)
    # plt.figure()
    # d.plot()
    imsd=tp.imsd(t_cor,0.1,0.2, max_lagtime=100, statistic='msd')
    emsd=tp.emsd(t_cor,0.1,0.2, max_lagtime=100)
    return imsd,emsd

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
                print repn
                nd_fh=rep[1]
                out_fn=splitext(basename(nd_fh))[0]
                out_fh=expt_dh+out_fn
                if not exists(out_fh+".imsd.png"):
                    # try:
                    print nd_fh
                    imsd,emsd=nd2msd(nd_fh)
                    # except:
                    #     continue
                    emsd=pd.DataFrame(emsd)
                    emsd.columns=[repn]
                    emsd=emsd.reset_index()
                    imsd.to_csv(out_fh+".imsd")
                    emsd.to_csv(out_fh+".emsd")

                    fig = plt.figure(figsize=(3, 3))
                    ax_imsd=plt.subplot(111)
                    ax_imsd=plot_msd(imsd,emsd,ax_imsd)
                    ax_imsd.figure.savefig(out_fh+".imsd.png");plt.clf();plt.close()
                else:
                    emsd=pd.read_csv(out_fh+".emsd")
                    del emsd["Unnamed: 0"]
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
    if len(expt_data)!=0:
        expt_data.to_csv(expt_dh+"expt.emsd")

        fig = plt.figure(figsize=(6, 3))
        ax_emsd=plt.subplot(121)
        ax_emsd=plot_emsd(expt_data,ax_emsd)
        plt.savefig(expt_dh+"emsd.png");plt.clf();plt.close()
        return expt_data
    else:
        logging.error('no data analyzed')