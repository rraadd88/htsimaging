#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
from os.path import exists
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # filename=cfg_xls_fh+'.log'
import nd2reader
import pandas as pd
import string
import numpy as np
from scipy import stats,ndimage
import os
from multiprocessing import Pool
from os.path import splitext,basename
import cv2
from skimage.segmentation import random_walker
# from skimage.data import binary_blobs
from skimage import io,exposure,restoration,filters,morphology,measure
import matplotlib.animation as animation
import matplotlib.pyplot as plt



def main(fh_xls):
    
    # fh_xls=sys.argv[1]
    info=pd.read_excel(fh_xls,'info')
    info=info.set_index('varname')
    for var in info.iterrows() :
        val=info['input'][var[0]]
        if not pd.isnull(val):
            exec("%s=info['input']['%s']" % (var[0],var[0]),locals(), globals())
        else:
            exec("%s=info['default']['%s']" % (var[0],var[0]),locals(), globals())

    #../tests/test.xlsx
    # fh_xls='../test/test.xlsx'
    data_job=pd.read_excel(fh_xls,'JobView')
    # nd_dh="/media/Transcend/20160219_000356_267"
    data_fns=pd.pivot_table(data_job,values='File Name',index='Loop_bleach Index',columns='Well Name', aggfunc=lambda x: x.iloc[0])
    data_fns_P =pd.pivot_table(data_job,values='File Name',index='TimeLapse1 Index',columns='Well Name', aggfunc=lambda x: x.iloc[0])
    data_fns   =pd.concat([data_fns,data_fns_P],axis=0)

    wells=[str(x) for x in list(data_job['Well Name'].unique())]
    wells.sort()
    wells_b, wells_u =wells[::2],wells[1::2]


    pt00s=range(24)
    for rowi in range(2,16):
        for coli in range(12):
            pt00s.append(coli) 
    pt00s_all=[None] * 384        
    pt00s_all[::2]=pt00s
    pt00s_all[1::2]=pt00s
    info_pt00s=pd.DataFrame({'well' : wells, \
                            'pt00': pt00s_all})
    info_pt00s=info_pt00s.set_index('well')
    time=np.array(range(12))*2
    # data_num=pd.DataFrame(columns=wells)
    # data_num.loc[:,'time']=time

    if not exists("%s.data_num_kin" % (fh_xls)):
        data_num_kin=pd.DataFrame(columns=wells)
        data_num_kin.loc[:,'time']=time
        for well in wells:
            logging.info("processing: %s" % well) 
            nd_fns=data_fns[well].dropna().unique()
            well_kin=nd2kins(nd_fns,nd_dh,time_increment)
            if not pd.isnull(info_pt00s.loc[well])[0]:
                pt00=int(info_pt00s.loc[well])
                data_num_kin[well]=well_kin[pt00:pt00+12].values
            else :   
                data_num_kin[well]=well_kin[0:12].values
        data_num_kin.to_csv("%s.data_num_kin" % (fh_xls))
    else:
        logging.info("already processed: %s.data_num_kin" % (fh_xls))
    data_num_kin=pd.read_csv("%s.data_num_kin" % (fh_xls))
    data_num_kin2ana(data_num_kin,wells,wells_b,wells_u)

if __name__ == '__main__':
    # # GET INPTS
    if len(sys.argv) < 1: 
        print "### ERROR : check number of nput args required"
        sys.exit()
    if not exists(sys.argv[1]) :
        print >> sys.stderr, "### ERROR : Could not find '%s'!\n" % sys.argv[1]
        sys.exit(1)
    # if not exists(sys.argv[2]) :
    #     print >> sys.stderr, "### ERROR : Could not find '%s'!\n" % sys.argv[2]
    #     sys.exit(1)

    main(sys.argv[1])
