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

def main(fh_xls,well):
    def nd2arr_list(nd_dh,nd_fns):
        arr_list=[]
        for nd_fn in nd_fns:
            nd = nd2reader.Nd2("%s/%s" % (nd_dh,nd_fn))
            for ndi in nd:
                arr_list.append(np.array(ndi))
            del nd
        return arr_list

    # def raw2phasecorr(arr_list): #skimage
    #     stb_arr_list=[]
    def threshold_otsu(smoothened):
        markers = np.zeros(smoothened.shape, dtype=np.uint)
        markers[smoothened < filters.threshold_otsu(smoothened)] = 1
        markers[smoothened > filters.threshold_otsu(smoothened)] = 2
        return markers

    def raw2phasecorr(arr_list): #cv
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
            stable_image = cv2.warpAffine(image, xform, dsize=(image.shape[1], image.shape[0]))
            prev_image = image
            #clip sides
            ht,wd=np.shape(stable_image)
            clip=0.125 #0.25
            lt=int(wd*clip)
            rt=int(wd-wd*clip)
            up=int(ht*clip)
            dw=int(ht-ht*clip)
            stable_image_clipped=stable_image[up:dw,lt:rt]
            stb_arr_list.append(stable_image_clipped)
        return stb_arr_list

    def arr_list2regions(arr_list, time_increment):
        pre_bleach=arr_list_stb[0]
        smoothened = filters.median(pre_bleach.astype('uint16'),np.ones((4,4)))
        markers = threshold_otsu(smoothened)
        labels = random_walker(smoothened, markers, beta=10, mode='bf')
        regions= measure.label(labels)
        label_objects, nb_labels = ndimage.label(regions)
        sizes = np.bincount(label_objects.ravel())
        mask_sizes = ((sizes > 200) & (sizes < 5000))
        mask_sizes[0] = 0
        regions_cleaned = mask_sizes[label_objects]
        props = measure.regionprops(regions,arr_list[0]) #ref
        regions_areas=np.array([prop.area for prop in props])
        regions_index_large=np.where((regions_areas<5000) & (regions_areas>200))[0]
        kins_mean=pd.DataFrame(columns=regions_index_large, index=range(len(arr_list)))
        for i in range(len(arr_list)):
            props = measure.regionprops(regions,arr_list[i])
            means=np.array([prop.mean_intensity for prop in props])
            kins_mean.loc[i,:]=means[regions_index_large]
            del props
        kins_mean=kins_mean.loc[:, ~(kins_mean < 2000).any(axis=0)] #stitch
        kins_mean['time']=np.array(range(len(arr_list)))*time_increment #stitch
        return regions_cleaned,kins_mean

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
            ax_img.contour(regions, [0.5], linewidths=1.2, colors='r',animated=False)
            if len(kins_mean.columns)>1:
                kins_mean.plot(x='time',legend=False,ax=ax_kin)
            ax_kin.set_xlim([kins_mean.loc[0,'time'],kins_mean.loc[len(kins_mean)-1,'time']])
            ax_kin.axvline(kins_mean.loc[i,'time'], color='r', linestyle='--',lw=2)
            plt.savefig(png_dh+'/%02d.png' % i)
            ax_kin.clear()
        bash_command=("ffmpeg -f image2 -r 4 -i "+png_dh+"/%02d.png -vcodec mpeg4 -y "+vid_fh)
        subprocess.Popen(bash_command, shell=True, executable='/bin/bash')

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
    if not exists("%s.%sstb.mp4" % (fh_xls,well)):
        print ">>> STATUS  : nd2vid : %s" % well 
        nd_fns=data_fns[well].dropna().unique()
        arr_list=nd2arr_list(nd_dh,nd_fns)
        arr_list_stb=raw2phasecorr(arr_list)
        regions,kins_mean=arr_list2regions(arr_list_stb,time_increment)
        arr_list2vid(arr_list_stb,regions,kins_mean,('%s.%sstb.mp4' % (fh_xls,well)),384, 384)
        arr_list2vid(arr_list    ,regions,kins_mean,('%s.%sraw.mp4' % (fh_xls,well)),384, 384)
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