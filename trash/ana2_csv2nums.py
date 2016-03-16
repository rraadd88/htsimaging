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
from scipy.optimize import curve_fit
# from scipy.stats import gumbel_r

def pdf2cdfestimate(bincenters,pdf,percentile):
    pdf=pdf/pdf.sum()
    cdf=np.cumsum(pdf)
    for i in range(len(cdf)):
        if (percentile-cdf[i])<0:
            break
#     cdf.ix[cdfi[0]-50]
#     cdfi[0]-50
    Y=np.array([bincenters[i-1],bincenters[i]])
    X=np.array([cdf[i-1], cdf[i]])
    return np.interp(percentile, X, Y)

def csv2nums(csv_fh):
    nums_fh=csv_fh+".nums"    
    if exists(csv_fh):# and (not exists(nums_fh)) :
        print ">>> STATUS  : processing : %s" % csv_fh 
        arr=np.genfromtxt(csv_fh,delimiter=',')
        vec=np.concatenate(arr)
        nums=pd.DataFrame()
        
        nums.loc[0,'mode']=stats.mode(vec)[0][0]
        nums.loc[0,'mean']=np.mean(vec)
        nums.loc[0,'median']=np.median(vec)
        nums.loc[0,'Q25']=np.percentile(vec, 25)
        nums.loc[0,'Q50']=np.percentile(vec, 50)
        nums.loc[0,'Q75']=np.percentile(vec, 75)
        nums.loc[0,'Q90']=np.percentile(vec, 90)
        
        bins = np.arange(0, 5000, 50)#np.arange(5000)
        probs, binedges = np.histogram(vec, bins=bins, normed=True)
        bincenters = 0.5*(binedges[1:]+binedges[:-1])
        nums.loc[0,'cdf25']=pdf2cdfestimate(bincenters,probs,0.25)        
        nums.loc[0,'cdf50']=pdf2cdfestimate(bincenters,probs,0.5)        
        nums.loc[0,'cdf75']=pdf2cdfestimate(bincenters,probs,0.75)        
        nums.loc[0,'cdf90']=pdf2cdfestimate(bincenters,probs,0.9)            
        nums.to_csv(nums_fh)
    else:
        print ">>> WARNING   : skipping : %s" % csv_fh
# data_job.to_csv('data_xls_fh_nums.csv')

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

    # csv2nums(csv_fhs[0])
    pool=Pool(processes=int(cores)) 
    pool.map(csv2nums,csv_fhs)
    pool.close(); pool.join()
    
    data_job.set_index(['File Name', 'File Frame Index'], inplace=True)
    for i in range(len(nums_fhs)):  
        nums_fh=nums_fhs[i]
        nums_df=pd.read_csv(nums_fh)
        nd_fn=nd_fn_framei_tps[i][0]
        framei=nd_fn_framei_tps[i][1]
        data_job.loc[(nd_fn,framei),'mean']  =nums_df.loc[0,'mean']
        data_job.loc[(nd_fn,framei),'mode']  =nums_df.loc[0,'mode']
        data_job.loc[(nd_fn,framei),'median']=nums_df.loc[0,'median']
        data_job.loc[(nd_fn,framei),'Q25']=nums_df.loc[0,'Q25']
        data_job.loc[(nd_fn,framei),'Q50']=nums_df.loc[0,'Q50']
        data_job.loc[(nd_fn,framei),'Q75']=nums_df.loc[0,'Q75']
        data_job.loc[(nd_fn,framei),'Q90']=nums_df.loc[0,'Q90']

        data_job.loc[(nd_fn,framei),'cdf25']=nums_df.loc[0,'cdf25']
        data_job.loc[(nd_fn,framei),'cdf50']=nums_df.loc[0,'cdf50']
        data_job.loc[(nd_fn,framei),'cdf75']=nums_df.loc[0,'cdf75']
        data_job.loc[(nd_fn,framei),'cdf90']=nums_df.loc[0,'cdf90']

    data_job=data_job.reset_index()
    data_job.to_csv(data_xls_fh+".nums")      