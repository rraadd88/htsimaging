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

import pims
import cv2
from skimage.segmentation import random_walker
# from skimage.data import binary_blobs
from skimage import io,exposure,restoration,filters,morphology,measure

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import subprocess



def nd2arr_list(nd_dh=None,nd_fns=[],nd_fh=None):
    arr_list=[]
    if nd_fh is None:
        for nd_fn in nd_fns:
            nd = nd2reader.Nd2("%s/%s" % (nd_dh,nd_fn))
            for ndi in nd:
                arr_list.append(np.array(ndi))
            del nd
        return arr_list
    elif not nd_fh is None:
        nd = nd2reader.Nd2("%s" % (nd_fh))
        for ndi in nd:
            arr_list.append(np.array(ndi))
        return arr_list
    
def raw2phasecorr(arr_list,clip=0): #cv
    cx = 0.0
    cy = 0.0
    stb_arr_list=[]
    prev_frame = arr_list[0]
    prev_image = np.float32(restoration.denoise_tv_chambolle(prev_frame.astype('uint16'), weight=0.1, multichannel=True)) #ref
    for frame in arr_list:           
        image = np.float32(restoration.denoise_tv_chambolle(frame.astype('uint16'), weight=0.1, multichannel=True))
        # TODO: set window around phase correlation
        dp = cv2.phaseCorrelate(prev_image, image)
        cx = cx - dp[0]
        cy = cy - dp[1]
        xform = np.float32([[1, 0, cx], [0, 1, cy]])
        stable_image = cv2.warpAffine(frame.astype('float32'), xform, dsize=(image.shape[1], image.shape[0]))
        prev_image = image
        #clip sides
        ht,wd=np.shape(stable_image)
#         clip=0.125 #0.25
        lt=int(wd*clip)
        rt=int(wd-wd*clip)
        up=int(ht*clip)
        dw=int(ht-ht*clip)
        stable_image_clipped=stable_image[up:dw,lt:rt]
        stb_arr_list.append(stable_image_clipped)
    return stb_arr_list

def arr_list2regions(arr_list, time_increment):
    pre_bleach=arr_list[0]
    denoised=restoration.denoise_bilateral(pre_bleach.astype('uint16'), sigma_range=0.01, sigma_spatial=15)
    smoothened = filters.median(denoised,np.ones((4,4)))
    markers = np.zeros(smoothened.shape, dtype=np.uint)
    markers[smoothened < filters.threshold_otsu(smoothened)] = 1
    markers[smoothened > filters.threshold_otsu(smoothened)] = 2
    labels = random_walker(smoothened, markers, beta=10, mode='bf')
    regions= measure.label(labels)
    props = measure.regionprops(regions,arr_list[0]) #ref
    regions_areas=np.array([prop.area  for prop in props])
    regions_lbls =np.array([prop.label for prop in props])
    regions_means=np.array([prop.mean_intensity  for prop in props])
    regions_index_large=np.where((regions_areas<10000) & (regions_areas>200) & (regions_means>2000))[0]
    regions_lbls_large=regions_lbls[regions_index_large] 
    regions_large=np.zeros((regions.shape), dtype=bool)
    for i in regions_lbls_large:
        booli = (regions == i) 
        regions_large=np.logical_or(regions_large,booli)
    kins_mean=pd.DataFrame(columns=regions_lbls_large, index=range(len(arr_list)))
    for i in range(len(arr_list)):
        props = measure.regionprops(regions,arr_list[i])
        means=np.array([prop.mean_intensity for prop in props])
        kins_mean.loc[i,:]=means[regions_index_large]
        del props
    kins_mean=kins_mean.loc[:, ~(kins_mean < 3000).any(axis=0)] #stitch
    kins_mean['time']=np.array(range(len(arr_list)))*time_increment #stitch
    return regions_large,kins_mean

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

@pims.pipeline
def average_z(frames):
    if len(np.shape(frames))==4:
        return frames.mean(axis=0)  # the same as image[0] + ... + image[4]
    else: 
        return frames

