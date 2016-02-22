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

def csv2nums(csv_fh):
    nums_fh=csv_fh+".nums"    
    if exists(csv_fh) and (not exists(nums_fh)) :
#         print ">>> STATUS  : processing : %s" % csv_fh 
        arr=np.genfromtxt(csv_fh,delimiter=',')
        vec=np.concatenate(arr)
        nums=pd.DataFrame()
        nums.loc[0,'mode']=stats.mode(vec)[0][0]
        nums.loc[0,'mean']=np.mean(vec)
        nums.loc[0,'median']=np.median(vec)
        nums.to_csv(nums_fh)
    else:
        print ">>> WARNING   : skipping : %s" % csv_fh
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
        
    data_job=pd.read_excel(data_xls_fh,'JobView')
    data_job_del_cols=[str(s) for s in data_job.columns if 'Intensity' in str(s) ]
    for col in data_job_del_cols:
        exec("del data_job['%s']" % col,locals(), globals())
    
    nd_fns =[str(s) for s in data_job['File Name'].tolist()]
    frameis=[int(s) for s in data_job['File Frame Index'].tolist()]
    
    csv_fhs =[]
    nums_fhs=[]
    nd_fn_framei_tps=[]
    for i in range(len(nd_fns)):
        csv_fhs.append("%s/%s%02d.csv" % (nd_dh,nd_fns[i],frameis[i]))
        nums_fhs.append("%s/%s%02d.csv.nums" % (nd_dh,nd_fns[i],frameis[i]))
        nd_fn_framei_tps.append([nd_fns[i],frameis[i]])

        #       csv2nums(nd_fns_frameis_tps[0])
    pool=Pool(processes=int(cores)) 
    pool.map(csv2nums,csv_fhs)
    pool.close(); pool.join()
    
    for i in range(len(nums_fhs)):  
        nums_fh=nums_fhs[i]
        nums_df=pd.read_csv(nums_fh)
        nd_fn=nd_fn_framei_tps[i][0]
        framei=nd_fn_framei_tps[i][1]
        data_job.loc[(nd_fn,framei),'mean']  =nums_df.loc[0,'mean']
        data_job.loc[(nd_fn,framei),'mode']  =nums_df.loc[0,'mode']
        data_job.loc[(nd_fn,framei),'median']=nums_df.loc[0,'median']
       
    data_job=data_job.reset_index()
    data_job.to_csv(data_xls_fh+".nums")
        
        
        
        
        
        
        
        
        