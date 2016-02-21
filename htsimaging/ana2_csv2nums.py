#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
import pandas as pd
import string
import numpy as np
from scipy import stats
import os
from multiprocessing import Pool
from os.path import splitext,basename,exists

data_xls_fh=sys.argv[1]
info=pd.read_excel(data_xls_fh,'info')
info=info.set_index('varname')
for var in info.iterrows():
    val=info['input'][var[0]]
    #print var[0]
    if not pd.isnull(val):
        exec("%s=info['input']['%s']" % (var[0],var[0]))
    else:
        exec("%s=info['default']['%s']" % (var[0],var[0]))

def csv2nums(nd_fns_frameis_tps):
    nd_fn=nd_fns_frameis_tps[0]
    framei=nd_fns_frameis_tps[1]
    fh="%s/%s%02d.csv" % (nd_dh,nd_fn,framei)
    if exists(fh):
        print ">>> STATUS  : processing : %s" % fh 
        arr=np.genfromtxt(fh,delimiter=',')
        vec=np.concatenate(arr)
        nums=pd.DataFrame()
        nums.loc[0,'mode']=stats.mode(vec)[0][0]
        nums.loc[0,'mean']=np.mean(vec)
        nums.loc[0,'median']=np.median(vec)
        nums.to_csv(fh+".nums")
    else:
        print ">>> ERROR   : file no not exist : %s" % fh
#     data_job.to_csv('data_xls_fh_nums.csv')

# import matplotlib.pyplot as plt
# hist=plt.hist(vec,80)

if __name__ == '__main__':
    # # GET INPTS
    if 1 < len(sys.argv) < 1:
        print "### ERROR : check number of nput args required"
        sys.exit()
    if not exists(sys.argv[1]):
        print >> sys.stderr, "### ERROR : Could not find '%s'!\n" % sys.argv[1]
        exit(1)
        
    data_job=pd.read_excel(data_xls_fh,'JobView')
    data_job_del_cols=[str(s) for s in data_job.columns if 'Intensity' in str(s) ]
    for col in data_job_del_cols:
        exec("del data_job['%s']" % col,locals(), globals())
    nd_fns =[str(s) for s in data_job['File Name'].tolist()]
    frameis=[int(s) for s in data_job['File Frame Index'].tolist()]
    nd_fns_frameis_tps=[]
    for i in range(len(nd_fns)):
        nd_fns_frameis_tps.append([nd_fns[i],frameis[i]])
#     csv2nums(nd_fns_frameis_tps[0])
    pool=Pool(processes=int(cores)) 
    pool.map(csv2nums,nd_fns_frameis_tps)
    pool.close(); pool.join()
