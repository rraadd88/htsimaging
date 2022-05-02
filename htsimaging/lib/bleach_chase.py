#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    



def arr_list_stb2kin(arr_list_stb,cut_off):
    ## cut_off=3000
    pre_bleach=arr_list_stb[0]
    index=np.where(pre_bleach>cut_off)
    kin=pd.DataFrame(columns=["time","mean"], index=range(len(arr_list_stb)))
    for arri in range(len(arr_list_stb)):
        arr=arr_list_stb[arri]
        kin.loc[arri,"mean"]=np.mean(arr[index])
        del arr
    kin['time']=np.array(range(len(arr_list_stb)))*time_increment
    return kin

def nd2kins(nd_fns,nd_dh,time_increment):
    arr_list=nd2arr_list(nd_dh,nd_fns)
    arr_list_stb=raw2phasecorr(arr_list)
    kin=arr_list_stb2kin(arr_list_stb,3000) #stitch
    kins_mean=kin.drop(['time'], axis=1) 
    return kins_mean

def data_num_kin2diff(data_num,b_well,u_well):    
    b_data =data_num.loc[:,b_well]
    u_data =data_num.loc[:,u_well]
    otpt_diff=np.log(u_data/b_data)
    return otpt_diff

def getrateofrecov(x,y):
    otpt_diff_fitted=pd.Series(index=range(len(x)))
    for totpts in range(12,4,-1):
        slope, intercept, r_value, p_value, slope_std_error = stats.linregress(x[:totpts], y[:totpts])
        #print r_value
        if r_value < -0.85:
            predict_y = intercept + slope * x[:totpts]
            otpt_diff_fitted=predict_y
            slope=slope
            break
        else:
            slope=np.nan
    #alpha=stats.linregress(x,y)
    return slope,otpt_diff_fitted

def data_num_kin2ana(data_num_kin,wells,wells_b,wells_u): #num_type mode, mean, meadian
    diff_df=pd.DataFrame(columns=wells_b)
    diff_df.loc[:,'time']=time
    for welli in range(len(wells_b)):
        logging.info("processing: %s" % wells_b[welli])
        diff_df.loc[:,wells_b[welli]]  =data_num_kin2diff(data_num_kin,wells_b[welli],wells_u[welli])
    diff_df.to_csv("%s.diff_df" % (fh_xls))  

    rateofrecov_df=pd.DataFrame(index=wells_b,columns=['rateofrecov'])
    rateofrecov_df.index.name='smp_well'
    diff_fitted_df=pd.DataFrame(columns=wells_b)
    diff_fitted_df.loc[:,'time']=time
    for welli in range(len(wells_b)): # get rate of recov
        rateofrecov_df.loc[wells_b[welli],'rateofrecov'] =getrateofrecov(diff_df.loc[:,'time'],diff_df.loc[:,wells_b[welli]])[0]        
        diff_fitted_df.loc[:,wells_b[welli]]             =getrateofrecov(diff_df.loc[:,'time'],diff_df.loc[:,wells_b[welli]])[1]
    diff_fitted_df.to_csv("%s.diff_fitted_df" % (fh_xls))


    wells_rows=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    wells_rows_ctrl, wells_rows_test =wells_rows[::2],wells_rows[1::2] # odd,even
    rateofrecov_df.to_csv("%s.rateofrecov" % (fh_xls))      

    rateofrecov_df=rateofrecov_df.reset_index()    
    lbls=pd.read_excel(fh_xls,'lbls')
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
    rateofrecov_per_smpi_df.to_csv("%s.rateofrecov_per_smpi_df" % (fh_xls))
    del rateofrecov_per_smpi_df.columns.name

    info_genes_df=pd.read_excel('/home/kclabws1/Documents/propro/writ/prjs/2_chem_chap_screens/data/data_sampling_yeast_gfp_genes/data/160208_corrected_smpi/151217_final_list90_data.xls')
    combo=pd.concat([info_genes_df,rateofrecov_per_smpi_df],axis=1)
    combo.to_csv("%s.combo" % (fh_xls))
