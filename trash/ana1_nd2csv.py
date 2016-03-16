#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import nd2reader
import pandas as pd
import string
import numpy as np
from scipy import stats
import os
from multiprocessing import Pool
from os.path import splitext,basename

fh_xls='../test/test.xlsx'
#otpt_dh='/media/Transcend/20160219_000356_267_csvs'
nd2_dh="/media/Transcend/20160219_000356_267"
cores=16

data_job=pd.read_excel(fh_xls,'JobView')
nd_fhs=[ nd2_dh+'/'+str(fn) for fn in data_job['File Name'].unique()]

# convert t0 csvs
def nd2csv(nd_fh):
    nd = nd2reader.Nd2(nd_fh)
    frame_cnt=0
    for frame in nd:
        np.savetxt("%s/%s%02d.csv" % (nd2_dh,basename(nd_fh),frame_cnt), frame, delimiter=",")
        frame_cnt+=1
    nd.close
pool_data_lbl2data_fit=Pool(processes=int(cores)) 
pool_data_lbl2data_fit.map(nd2csv,nd_fhs)
pool_data_lbl2data_fit.close(); pool_data_lbl2data_fit.join()

#for wellid in wellids:
#    nd_fns = [filename for filename in os.listdir(nd2_dh) if filename.startswith("Well%s_Well%s" % (wellid,wellid))]
#    nd_fns = np.sort(nd_fns)
#    for nd_fn in nd_fns:
#        nd_fh  = "%s/%s" % (nd2_dh,nd_fn)
#            
#rowids=string.uppercase[:16]
#colids=range(1,25)
#wellids=[]
#for rowid in rowids:
#    for colid in colids: 
#        wellid="%s%02d" % (rowid,colid)
#        wellids.append(wellid)        
#data_mode=pd.DataFrame(columns=wellids)
#        data_mode_nd=[]
#
#            nd2 = nd2reader.Nd2(nd_fh)
#            nd2_data=np.concatenate(nd2[0])
#            data_mode_nd.append(stats.mode(nd2_data)[0][0])
#        data_mode[wellid]=data_mode_nd
#       
##hist=plt.hist(nd2_data,100)
##import matplotlib.pyplot as plt
##hist=plt.hist(nd2_data,100)
#
#fig, ax = plt.subplots()
#heatmap = ax.pcolor(nd2[0])
#
#plt.imshow(data)