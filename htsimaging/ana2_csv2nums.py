#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
import pandas as pd
import string
import numpy as np
from scipy import stats
import os
#from multiprocessing import Pool
from os.path import splitext,basename,exists

if len(sys.argv)>2:
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

def main(data_xls_fh):
    #data_xls_fh='../test/test.xlsx'
    data_job=pd.read_excel(data_xls_fh,'JobView')
    data_job_del_cols=[str(s) for s in data_job.columns if 'Intensity' in str(s) ]
    for col in data_job_del_cols:
        exec("del data_job['%s']" % col,locals(), globals())
    data_job
    for row in data_job.iterrows():
        print row[0]
        row_df=pd.DataFrame(row[1])
        fn_nd=str(row_df.loc['File Name',row[0]])
        framei=row_df.loc['File Frame Index',row[0]]
        fh="%s/%s%02d.csv" % (nd_dh,fn_nd,framei)
        if exists(fh):
            arr=np.genfromtxt(fh,delimiter=',')
            vec=np.concatenate(arr)
            data_job.loc[row[0],'mode']=stats.mode(vec)[0][0]
            data_job.loc[row[0],'mean']=np.mean(vec)
            data_job.loc[row[0],'median']=np.median(vec)
        else:
            print ">>> ERROR   : file no not exist : %s" % fh
        del row_df
    data_job.to_csv('data_xls_fh_nums.csv')

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
    main(sys.argv[1])
