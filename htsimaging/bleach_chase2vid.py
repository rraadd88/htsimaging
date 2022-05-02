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
# from nd2kins import nd2arr_list,raw2phasecorr,arr_list2regions

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
        arr_list_stb=raw2phasecorr(arr_list,clip=0.125)
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