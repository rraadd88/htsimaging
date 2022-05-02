#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

```
This converts the imaging data from microscopes to the usable images for analysis.
```

import sys
import os
from os import makedirs
from os.path import exists,splitext,basename
import glob
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


def arr2vid(arr_list,regions,kins_mean,vid_fh,xpixels, ypixels):
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

    
def makevid(gfp_list_stb,brf_list_stb,cmap_gfp,cmap_brf,vid_fh,conditionn=None,interval=None,dpi = 300):
    png_dh=splitext(vid_fh)[0]
    if not exists(png_dh):
        makedirs(png_dh)
    for i in range(len(gfp_list_stb)):
        fig = plt.figure(dpi=dpi)
        ax=plt.subplot(111)    
        ax.imshow(gfp_list_stb[i],cmap=cmap_gfp,zorder=1,alpha=0.5)
        ax.imshow(brf_list_stb[i],cmap=cmap_brf,zorder=0,alpha=0.5)
        ax.text(2,2,conditionn,color='r',fontsize=20,)
        ax.text(2,np.shape(gfp_list_stb[i])[1]*0.95,"time: %03d s" % ((i+1)*interval),color='b',fontsize=20,)
#         fig.subplots_adjust(wspace=.4)
        plt.tight_layout()
        plt.axis('off')
        fig.savefig('%s/%03d.png' % (png_dh,i))
        plt.close(fig)
    bash_command=("ffmpeg -f image2 -r 20 -i "+png_dh+"/%03d.png -vcodec mpeg4 -y "+vid_fh)
    print vid_fh
    subprocess.Popen(bash_command, shell=True)

