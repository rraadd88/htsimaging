#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
from os.path import exists
import pandas as pd
import numpy as np
import seaborn as sns
#import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit

def main(data_xls_fh):
    def data2diff(data_num,b_well,u_well,info_pt00):    
        #def norm2max(pdseries):
        #    otpt_norm=pd.Series(index=pdseries.index)
        #    mx=pdseries.max()
        #    for i in range(len(pdseries)):
        #        print i
        #        otpt_norm.ix[i]=pdseries.ix[i]-mx+1
        #    return otpt_norm            
        #time   =data_num.loc[:,'time']
        b_data =data_num.loc[:,b_well]
        u_data =data_num.loc[:,u_well]
        b_data=b_data.loc[~b_data.isnull()]
        u_data=u_data.loc[~u_data.isnull()]
        b_data=b_data.reset_index(drop=True)
        u_data=u_data.reset_index(drop=True)    
        if not pd.isnull(info_pt00.loc[b_well])[0]:
            pt00=int(info_pt00.loc[b_well])
            otpt_b_data=b_data[pt00:pt00+12]
            #return 
            otpt_u_data=u_data[pt00:pt00+12]
            #return otpt_u_data
            #otpt_diff=np.log(otpt_u_data-otpt_b_data+10000)
            # otpt_diff=otpt_u_data-otpt_b_data
            otpt_diff=otpt_b_data-otpt_u_data
            #otpt_diff=norm2max(otpt_diff)
            # otpt_diff=otpt_diff - otpt_diff.min() # normalised to min
            #otpt_diff=(otpt_diff - otpt_diff.min())/(otpt_diff.max()-otpt_diff.min()) # normalised to max
            return otpt_b_data.reset_index(),otpt_u_data.reset_index(),otpt_diff.reset_index(drop=True)      
        else :   
            otpt_b_data=b_data[0:pt00+12]
            #return otpt_b_data
            otpt_u_data=u_data[0:pt00+12]
            #return otpt_u_data
            otpt_diff=np.log(otpt_u_data-otpt_b_data+10000)
            return otpt_b_data.reset_index(),otpt_u_data.reset_index(),otpt_diff.reset_index(drop=True)
    
    def exp1_growth(x, a, b, c):
        return a * np.exp( b * x) + c  

    def exp1_decay(x, a, b, c):
        return a * np.exp(-b * x) + c  
                                                                                                                                                      
    def getrateofrecov(y):
        y=(y - y.min())/(y.max()-y.min()) # norm to min0 max1
        y=y[~y.isnull()]#remove nans
        x=np.array(range(len(y)))*1.5
        otpt_diff_fitted=pd.Series(index=range(len(x)))
        #for totpts in range(12,4,-1):
            #rate, intercept, r_value, p_value, rate_std_error = stats.linregress(x[:totpts], y[:totpts])
        try:
            #popt, pcov = curve_fit(exp1, x, y,p0=[50,0,-50])
            popt, pcov = curve_fit(exp1_growth, x, y,p0=[0,0,1])
            #plt.figure()
            #plt.plot(x, y, 'ko', label="Original Noised Data")
            #plt.plot(x, exp1(x, *popt), 'r-', label="Fitted Curve")
            #plt.legend()
            #plt.show()
            perr=np.sqrt(np.diag(pcov)) #standard deviation errors
            #print r_value
            if perr[1] < 0.1:
                #predict_y = intercept + rate * x[:totpts]
                ## Plotting
                #plt.plot(x, y, 'o')
                #plt.plot(x[:totpts], predict_y, 'k-')
                #plt.show()
                rate=popt[1]
                otpt_diff_fitted=pd.Series(exp1_growth(x, *popt))
            else:
                rate=np.nan
                print ">>> WARNING : getrateofrecov : standard deviation error (%.4f) > 0.1 " % perr[1]
        except:
            rate=np.nan
            print ">>> WARNING : getrateofrecov : can not fit"
        #alpha=stats.linregress(x,y)
        return rate,otpt_diff_fitted

    def data_job2kin(data_job,num_type,info_pt00s,time,wells): #num_type mode, mean, meadian
        data_num   =pd.pivot_table(data_job,values=num_type,index='Loop_bleach Index',columns='Well Name')
        data_num_P =pd.pivot_table(data_job,values=num_type,index='TimeLapse1 Index',columns='Well Name')
        data_num   =pd.concat([data_num,data_num_P],axis=0)
        data_num   =data_num.reset_index()

        wells_b, wells_u =wells[::2],wells[1::2]

        data_num_kin=pd.DataFrame(columns=wells)
        data_num_kin.loc[:,'time']=time
        diff_df=pd.DataFrame(columns=wells_b)
        diff_df.loc[:,'time']=time
        diff_fitted_df=pd.DataFrame(columns=wells_b)
        diff_fitted_df.loc[:,'time']=time
        rateofrecov_df=pd.DataFrame(index=wells_b,columns=['rateofrecov'])
        rateofrecov_df.index.name='smp_well'

        for welli in range(len(wells_b)):   
            print wells_b[welli]
            data_num_kin.loc[:,wells_b[welli]]     =data2diff(data_num,wells_b[welli],wells_u[welli],info_pt00s)[0]
            data_num_kin.loc[:,wells_u[welli]]     =data2diff(data_num,wells_b[welli],wells_u[welli],info_pt00s)[1]
            diff_df.loc[:,wells_b[welli]]           =data2diff(data_num,wells_b[welli],wells_u[welli],info_pt00s)[2]
        #background correction
        diff_df_blanks=diff_df[['O13','O15','O17','O19','O21','O23','P13','P15','P17','P19','P21','P23']]
        diff_df_blanks_average=diff_df_blanks.mean(axis=1)
        diff_df_blanks.to_csv("%s.%s.diff_df_blanks" % (data_xls_fh,num_type))     
        for welli in range(len(wells_b)):            
            diff_df.loc[:,wells_b[welli]]=diff_df.loc[:,wells_b[welli]]-diff_df_blanks_average
            diff_df_mn=diff_df.loc[:,wells_b[welli]].min()
            diff_df_mx=diff_df.loc[:,wells_b[welli]].max()
            diff_df.loc[:,wells_b[welli]]=(diff_df.loc[:,wells_b[welli]]-diff_df_mn)/(diff_df_mx-diff_df_mn)
            rateofrecov_df.loc[wells_b[welli],'rateofrecov'] =getrateofrecov(diff_df.loc[:,wells_b[welli]])[0]        
            diff_fitted_df.loc[:,wells_b[welli]]             =getrateofrecov(diff_df.loc[:,wells_b[welli]])[1]
            del diff_df_mn,diff_df_mx
        wells_rows=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
        wells_rows_ctrl, wells_rows_test =wells_rows[::2],wells_rows[1::2] # odd,even
        rateofrecov_df.to_csv("%s.%s.rateofrecov" % (data_xls_fh,num_type))	    
        rateofrecov_df=rateofrecov_df.reset_index()    
        lbls=pd.read_excel(data_xls_fh,'lbls')
        for s in wells_rows_ctrl:
            rateofrecov_df.loc[(rateofrecov_df.smp_well.astype(str).str.contains(s)),'smp_type']=str(lbls.ix[0,0])
        for s in wells_rows_test:
            rateofrecov_df.loc[(rateofrecov_df.smp_well.astype(str).str.contains(s)),'smp_type']=str(lbls.ix[1,0])
        smpi=[]
        for rowi in range(1,9):
            for i in list(np.array(range(1,13))+12*(rowi-1)):
                smpi.append(i)
            for i in list(np.array(range(1,13))+12*(rowi-1)):
                smpi.append(i)

        len(smpi)
        len(rateofrecov_df)
        rateofrecov_df['smpi']=smpi
        rateofrecov_df['rateofrecov']=rateofrecov_df['rateofrecov'].astype(float)
        rateofrecov_per_smpi_df=pd.pivot_table(rateofrecov_df,values='rateofrecov',index='smpi',columns='smp_type')
        rateofrecov_per_smpi_df=rateofrecov_per_smpi_df.reset_index()
        rateofrecov_per_smpi_df=rateofrecov_per_smpi_df.reset_index()
        del rateofrecov_per_smpi_df.columns.name

        info_genes_df=pd.read_excel('/home/kclabws1/Documents/propro/writ/prjs/2_chem_chap_screens/data/data_sampling_yeast_gfp_genes/data/160208_corrected_smpi/151217_final_list90_data.xls')
        combo=pd.concat([info_genes_df,rateofrecov_per_smpi_df],axis=1)

        data_num.to_csv("%s.%s.data_num" % (data_xls_fh,num_type)) 
        data_num_kin.to_csv("%s.%s.data_num_kin" % (data_xls_fh,num_type))
        diff_df.to_csv("%s.%s.diff_df" % (data_xls_fh,num_type))  
        diff_fitted_df.to_csv("%s.%s.diff_fitted_df" % (data_xls_fh,num_type))
        rateofrecov_per_smpi_df.to_csv("%s.%s.rateofrecov_per_smpi_df" % (data_xls_fh,num_type))
        combo.to_csv("%s.%s.combo" % (data_xls_fh,num_type))

    #data_xls_fh='../../../data/yeast_gfp_half_life/data/160211_bleach_chase/160211_yeast_bleach_chase_data.xlsx'	    
    data_job=pd.read_csv(data_xls_fh+".nums")
    wells=[str(x) for x in list(data_job['Well Name'].unique())]
    wells.sort()
    wells_b, wells_u =wells[::2],wells[1::2]
    pt00s=range(24)
    for rowi in range(2,16):
        for coli in range(12):
            pt00s.append(coli)    
    info_pt00s=pd.DataFrame({'well' : wells_b, \
                            'pt00': pt00s})
    info_pt00s=info_pt00s.set_index('well')
    time=np.array(range(24))*2

    for num_type in ['sum_thr','pixels_thr']:#['mean','mode','median','peak','mean_thr','mode_thr','median_thr','peak_thr']: 
    	print num_type
        data_job2kin(data_job,num_type,info_pt00s,time,wells)

	    # data_mean  =data_job2data_num(data_job,'mean')
	    # data_mode  =data_job2data_num(data_job,'mode')
	    # data_median=data_job2data_num(data_job,'median')
	    
	    #data_num=pd.read_csv(data_xls_fh+'data_num.csv') 
	    #diff_df=pd.read_csv(data_xls_fh+'diff.csv')
	    #rateofrecov_per_smpi_df=pd.read_csv(data_xls_fh+'rateofrecov.csv')
	    #combo=pd.read_csv(data_xls_fh+'rateofrecov_info.csv')

	    #plot_diff_df=diff_df
	    ##plot_diff_df=plot_diff_df.set_index('smp_well')
	    #plot_diff_df=plot_diff_df.drop('smp_well',axis=1)
	    #plot_diff_df=plot_diff_df.ix[:][:11]
	    #plot_diff_df=plot_diff_df.set_index('time')
	    #plot_diff_df=plot_diff_df.stack()
	    #plot_diff_df=plot_diff_df.reset_index()
	    ##plot_diff_df=plot_diff_df.ix[1:][:]
	    ##plot_diff_df=plot_diff_df.reset_index(drop=True)
	    #plot_diff_df.columns=['time','well','diff']
	    #plot_diff_df=plot_diff_df.sort(columns='well')
	    #
	    #test=plot_diff_df.ix[:15]
	    #sns.set(style="ticks")
	    ## Initialize a grid of plots with an Axes for each walk
	    #grid = sns.FacetGrid(test, col="well", hue="well", col_wrap=4, size=1.5)    
	    ## Draw a horizontal line to show the starting point
	    #grid.map(plt.axhline, y=9, ls=":", c=".5")    
	    ## Draw a line plot to show the trajectory of each random walk
	    #grid.map(plt.plot, "time", "diff", marker="o", ms=4)    
	    ##fitted lines
	    #
	    ## Adjust the tick positions and labels
	    #grid.set(xticks=np.arange(0,20,5), yticks=np.arange(9.1,9.5,0.1),
	    #        xlim=(-0.5, 18.5), ylim=(9.2, 9.3))    
	    ## Adjust the arrangement of the plots
	    #grid.fig.tight_layout(w_pad=1)
	    #grid.fig.savefig(data_xls_fh+'diff.png')       
	    

	    
	    
	    #rateofrecov_df.to_csv('test3.csv')
	    #rateofrecov_ctrl=rateofrecov_df.loc[(rateofrecov_df.smp_type.astype(str).str.contains('ctrl')),'rateofrecov']
	    #rateofrecov_test=rateofrecov_df.loc[(rateofrecov_df.smp_type.astype(str).str.contains('test')),'rateofrecov']
	    #rateofrecov_ctrl=rateofrecov_ctrl.loc[~rateofrecov_ctrl.isnull()]
	    #rateofrecov_test=rateofrecov_test.loc[~rateofrecov_test.isnull()]
	    #
	    #plt.scatter( list(rateofrecov_ctrl),list(rateofrecov_test) )
	    
	    
	    #rateofrecov_df=pd.DataFrame({'smp_well': wells_b, \
	    #                        'rateofrecov' : rateofrecov_list}) 
	    #%matplotlib inline    
	    #data_num.plot(x="time", y=["C01","C02"])
    
if __name__ == '__main__':
    # # GET INPTS
    if 1 < len(sys.argv) < 1: 
	print "### ERROR : check number of nput args required"
	sys.exit()
    if not exists(sys.argv[1]) :
        print >> sys.stderr, "### ERROR : Could not find '%s'!\n" % sys.argv[1]
        exit(1)
    main(sys.argv[1])
