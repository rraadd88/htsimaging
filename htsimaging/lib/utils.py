import sys
import os
from os import makedirs
from os.path import exists,splitext,basename
from glob import glob
import pims
# import nd2reader
import pandas as pd
import string
import numpy as np
from scipy import stats,ndimage
from multiprocessing import Pool
import cv2
from skimage.segmentation import random_walker
# from skimage.data import binary_blobs
from skimage import io,exposure,restoration,filters,morphology,measure
from scipy import ndimage as ndi
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib import colors
import subprocess

import trackpy as tp
from htsimaging.lib import spt

def filterframes(frames,cutoff=0):
    frames_cleaned=[]
    for framei in frames:
        frames_cleaned.append(filterframe(framei,cutoff))
    return frames_cleaned

def filterframe(frame,cutoff=0):
    frame_cleand=np.array(frame)
    frame_cleand[frame_cleandi<cutoff]=0
    return frame_cleaned

def filter_regions(img,regions,prop_type='area',mn=0,mx=0,check=False):
    regions_props       = measure.regionprops(regions.astype(int),img)
    regions_lbls        = np.array([prop.label for prop in regions_props])
    if prop_type=='area':
        regions_props_selected = np.array([prop.area  for prop in regions_props])
    elif prop_type=='mean_intensity':
        regions_props_selected = np.array([prop.mean_intensity  for prop in regions_props])
    elif prop_type=='eccentricity':
        regions_props_selected = np.array([prop.eccentricity  for prop in regions_props])        
    elif prop_type=='extent':
        regions_props_selected = np.array([prop.extent  for prop in regions_props])        

    if not check:
        regions_filtered_lbls = regions_lbls[np.where((regions_props_selected<mx) & (regions_props_selected>mn))[0]]
        regions_filtered=np.zeros(regions.shape)
        for lbli in regions_filtered_lbls:
            regions_filtered[np.where(regions==lbli)]=lbli
        return regions_filtered
    elif check:
        plt.hist(regions_props_selected)

def smoothen(img):
    denoised=restoration.denoise_bilateral(img.astype('uint16'), sigma_range=0.01, sigma_spatial=15)
    smoothened = filters.median(denoised,np.ones((4,4)))
    return smoothened

def smoothenframes(frames):
    frames_cleaned=[]
    for framei in frames:
        frames_cleaned.append(smoothen(np.array(framei)))
    return frames_cleaned
        
def get_regions(img):
    denoised=restoration.denoise_bilateral(img.astype('uint16'), sigma_range=0.01, sigma_spatial=15)
    smoothened = filters.median(denoised,np.ones((4,4)))
    markers = np.zeros(smoothened.shape, dtype=np.uint)
    markers[smoothened < filters.threshold_otsu(smoothened)] = 1
    markers[smoothened > filters.threshold_otsu(smoothened)] = 2
    labels = random_walker(smoothened, markers, beta=10, mode='bf')
    regions= measure.label(labels)
    return regions, denoised, smoothened,markers

def raw2phasecorr(arr_list,clip=0): #cv
    import cv2
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

def imlistcropper(imlist,loci):
    rowini,rowend,colini,colend=loci
    imlist_region=[]
    for im in imlist:
        imlist_region.append(im[rowini:rowend,colini:colend])
    return imlist_region


def imclipper(im_stb,clip):
    ht,wd=np.shape(im_stb)
#         clip=0.125 #0.25
    lt=int(wd*clip)
    rt=int(wd-wd*clip)
    up=int(ht*clip)
    dw=int(ht-ht*clip)
    im_stb_clipped=im_stb[up:dw,lt:rt]
    return im_stb_clipped

def phasecorr(imlist,imlist2=None,clip=0): #cv  [rowini,rowend,colini,colend]
    import cv2
    cx = 0.0
    cy = 0.0            
    imlist_stb=[]
    if imlist2!=None:
        imlist2_stb=[]

    imi=0
    im_prev = imlist[0]
    im_denoised_prev = np.float32(restoration.denoise_tv_chambolle(im_prev.astype('uint16'), weight=0.1, multichannel=True)) #ref
    for im in imlist:           
        im_denoised = np.float32(restoration.denoise_tv_chambolle(im.astype('uint16'), weight=0.1, multichannel=True))
        # TODO: set window around phase correlation
        dp = cv2.phaseCorrelate(im_denoised_prev, im_denoised)
        cx = cx - dp[0]
        cy = cy - dp[1]
        xform = np.float32([[1, 0, cx], [0, 1, cy]])
        im_stb = cv2.warpAffine(im.astype('float32'), xform, dsize=(im_denoised.shape[1], im_denoised.shape[0]))
        imlist_stb.append(imclipper(im_stb,clip))

        if imlist2!=None:
            im2=imlist2[imi]
            im2_stb=cv2.warpAffine(im2.astype('float32'), xform, dsize=(im_denoised.shape[1], im_denoised.shape[0]))
            imlist2_stb.append(imclipper(im2_stb,clip))

        im_denoised_prev = im_denoised
        imi+=1
    if imlist2!=None:
        return imlist_stb,imlist2_stb
    else:
        return imlist_stb