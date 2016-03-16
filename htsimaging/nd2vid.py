#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
import os
from os.path import exists,splitext,basename
import glob
import nd2reader
import pandas as pd
import string
import numpy as np
from scipy import stats,ndimage
from multiprocessing import Pool
import cv2
from skimage.segmentation import random_walker
# from skimage.data import binary_blobs
from skimage import io,exposure,restoration,filters,morphology,measure
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import subprocess
from nd2kins import nd2arr_list,raw2phasecorr,arr_list2regions

def arr_list2vid(arr_list,regions,kins_mean,vid_fh,xpixels, ypixels):
    dpi = 100
    png_dh=os.path.splitext(vid_fh)[0]
    if not os.path.exists(png_dh):
        try:
            os.makedirs(png_dh)
        except :
            print ">>> WARNING : race error data_lbl"
    fig = plt.figure(figsize=(ypixels/dpi*2, xpixels/dpi), dpi=dpi)
    ax_img = plt.subplot(1,2,1)
    ax_kin = plt.subplot(1,2,2)
    ax_img.set_axis_off()
    ax_img.set_aspect('equal')
    for i in range(len(arr_list)):
        ax_img.imshow(arr_list[i],cmap='gray',animated=True)
        ax_img.contour(regions, [0.25], linewidths=1.2, colors='r',animated=False)
        if len(kins_mean.columns)>1:
            kins_mean.plot(x='time',legend=False,ax=ax_kin)
            ax_kin.plot(kins_mean['time'],kins_mean.drop(['time'], axis=1).mean(axis=1),lw=6,color='k')
            ax_kin.set_ylabel("Fluorescence Intensity (FU)")
        ax_kin.set_xlim([kins_mean.loc[0,'time'],kins_mean.loc[len(kins_mean)-1,'time']])
        ax_kin.axvline(kins_mean.loc[i,'time'], color='r', linestyle='--',lw=2)
        fig.subplots_adjust(wspace=.4)
        fig.savefig(png_dh+'/%02d.png' % i)
        ax_kin.clear()
    bash_command=("ffmpeg -f image2 -r 4 -i "+png_dh+"/%02d.png -vcodec mpeg4 -y "+vid_fh)
    subprocess.Popen(bash_command, shell=True, executable='/bin/bash')

def main(fh_xls,well):
    info=pd.read_excel(fh_xls,'info')
    info=info.set_index('varname')
    for var in info.iterrows() :
        val=info['input'][var[0]]
        if not pd.isnull(val):
            exec("%s=info['input']['%s']" % (var[0],var[0]),locals(), globals())
        else:
            exec("%s=info['default']['%s']" % (var[0],var[0]),locals(), globals())

    data_job=pd.read_excel(fh_xls,'JobView')
    data_fns=pd.pivot_table(data_job,values='File Name',index='Loop_bleach Index',columns='Well Name', aggfunc=lambda x: x.iloc[0])
    data_fns_P =pd.pivot_table(data_job,values='File Name',index='TimeLapse1 Index',columns='Well Name', aggfunc=lambda x: x.iloc[0])
    data_fns   =pd.concat([data_fns,data_fns_P],axis=0)

    wells=[str(x) for x in list(data_job['Well Name'].unique())]
    if not any(x in well for x in wells):
        print "### ERROR : Could not find '%s'!" % well
        sys.exit(1)
    if not exists("%s.%sstb1.mp4" % (fh_xls,well)):
        print ">>> STATUS  : nd2vid : %s" % well 
        nd_fns=data_fns[well].dropna().unique()
        arr_list=nd2arr_list(nd_dh,nd_fns)
        arr_list_stb=raw2phasecorr(arr_list)
        plt.imshow(arr_list_stb[0])
        plt.savefig('%s.%s_ref_frame.png' % (fh_xls,well))
        regions,kins_mean=arr_list2regions(arr_list_stb,time_increment)
        kins_mean.to_csv('%s.%s_kins_mean.csv' % (fh_xls,well))
        arr_list2vid(arr_list_stb,regions,kins_mean,('%s.%sstb.mp4' % (fh_xls,well)),384, 384)
        # arr_list2vid(arr_list    ,regions,kins_mean,('%s.%sraw.mp4' % (fh_xls,well)),384, 384)
    else:
        print ">>> STATUS  : nd2vid :already done"

if __name__ == '__main__':
    # # GET INPTS
    if len(sys.argv) < 1: 
        print "### ERROR : check number of nput args required"
        sys.exit()
    if not exists(sys.argv[1]) :
        print >> sys.stderr, "### ERROR : Could not find '%s'!\n" % sys.argv[1]
        sys.exit(1)
    main(sys.argv[1],sys.argv[2])