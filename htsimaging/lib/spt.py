from os.path import splitext, join, exists, isdir,basename,abspath,dirname
from rohan.global_imports import *
import trackpy as tp
tp.ignore_logging()
# import nd2reader
import pims
import pims_nd2
import logging
from htsimaging.lib.plot import *
from htsimaging.lib.stat import *

# %matplotlib inline

@pims.pipeline
def average_z(image):
    return image.mean(axis=0) 
    
def plot_msd(imsd,emsd,ax=None,scale="log",plot_fh=None,
    params_msd={"mpp":0.0645,
                "fps":0.2,
                "max_lagtime":100
                }):
    if ax is None:
        plt.figure(figsize=(3, 3))
        ax=plt.subplot(111)

    head=params_msd["fps"]*params_msd["max_lagtime"]
    head=int(head)

    imsd=imsd.head(head)
    emsd=emsd.head(head)

    imsd.plot(x=imsd.index,legend=None,alpha=0.75,ax=ax)
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')

    emsd.plot(x=emsd.index,style='o',legend=None,ax=ax)
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
           xlabel='lag time $t$')
    plt.tight_layout()
    if not plot_fh is None: 
        ax.figure.savefig(plot_fh);
#         plt.clf();plt.close()
    return ax

def plot_emsd(expt_data,ax=None,color='k',scale="log",plot_fh=None):
    if ax is None:
        plt.figure(figsize=(6, 3))
        ax=plt.subplot(121)
    expt_data.plot(x=expt_data.index,style='o',c=color,ax=ax,markeredgecolor = 'none')
    if scale=="log":
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set(ylabel=r'$\langle \Delta r^2 \rangle$ [$\mu$m$^2$]',
                xlabel='lag time $t$')
    ax.legend(loc = 'center left', bbox_to_anchor = (1.0, 0.5)) #bbox_to_anchor=(2.2, 1.0)
    plt.tight_layout()
    if not plot_fh is None: 
        ax.figure.savefig(plot_fh);
#         plt.clf();plt.close()
    return ax 

def nd2frames(nd_fh):
    if nd_fh.endswith('nd2'):
        frames=pims.ND2_Reader(nd_fh)
    elif nd_fh.endswith('mp4'):
        frames=pims.Video(nd_fh)
        from pimsviewer.utils import wrap_frames_sequence
        frames=wrap_frames_sequence(frames)

    if nd_fh.endswith('nd2'):
        if len(np.shape(frames))==4:
            frames = average_z(frames)
    else:
        frames.default_coords['c'] = 1
        frames.bundle_axes='yx'
        frames.iter_axes = 't'
    return frames

def get_params_locate(frame,diameter=15,minmass_percentile=92,out_fh=None,test=True,figsize=None,dbug=False):
    if dbug:
        out_fh=f"{out_fh}_diameter_{diameter:02d}"
    f = tp.locate(frame, diameter, invert=False)
    minmass=np.percentile(f['mass'],minmass_percentile)
    logging.info('feature count= %s, %spercentile= %s'  % (len(f),minmass_percentile,minmass))
                 
    f = tp.locate(frame, diameter, invert=False, minmass=np.percentile(f['mass'],minmass_percentile))
    logging.info('feature count= %s, %spercentile= %s'  % (len(f),minmass_percentile,
                                                           np.percentile(f['mass'],minmass_percentile)))
    
    if test:
        logging.info('getting plots annotate')
#         plt.clf()
        if not out_fh is None:
            fig=plt.figure(figsize=figsize)
            ax=plt.subplot(111)
            ax=tp.annotate(f, frame,ax=ax)
#             savefig('%s.annotate.pdf' % out_fh,format='pdf')
            savefig(f'{out_fh}.annotate.png')
#             plt.clf()

        if not out_fh is None:
            logging.info('getting plots hist')
            cols=['mass','size','ecc','signal','raw_mass','ep']
            fig=plt.figure()
            ax=plt.subplot(111)
            _=f.loc[:,cols].hist(ax=ax)
            savefig(f'{out_fh}.feature_props.png')
#             plt.clf()

        if not out_fh is None:
            logging.info('getting plots bias')
            fig=plt.figure()
            tp.subpx_bias(f);
            savefig(f'{out_fh}.subpx_bias.png')
#             plt.clf()

    params_locate={'diameter':diameter,
                  'minmass':minmass}
    return params_locate

def plot_traj(frame,traj):
    fig=plt.figure(figsize=[frame.shape[0]*0.02,frame.shape[1]*0.02])        
    ax=plt.subplot(111)
    ax.imshow(frame,cmap='binary_r',alpha=0.8)
    ax = tp.plot_traj(traj,label=False,ax=ax,lw=2)
    plt.tight_layout()
    return ax

def frames2coords(frames,out_fh,
                  params_locate,params_msd,params_link_df={'search_range':20,},
                  mass_cutoff=0.5,size_cutoff=0.5,ecc_cutoff=0.5,
                    filter_stubs=True,flt_mass_size=True,flt_incomplete_trjs=True,
                    force=False,test=False):
    dns=['f_batch','t','t1','t2']
    dn2dp={dn:f'{out_fh}.{dn}.tsv' for dn in dns}
    dn2df={}
    if not exists(dn2dp['t2']) or force:
        if not exists(dn2dp['t']) or force:
            dn2df['f_batch']=tp.batch(frames,engine='numba',**params_locate)
            dn2df['t']=tp.link_df(dn2df['f_batch'], **params_link_df)
            logging.info(params_link_df)
            to_table(dn2df['f_batch'],dn2dp['f_batch'])
            to_table(dn2df['t'],dn2dp['t'])
        else:
            dn2df['t']=read_table(dn2dp['t'])
        max_lagtime_stubs=params_msd["max_lagtime"]*params_msd["fps"]
        if filter_stubs:
            dn2df['t1'] = tp.filter_stubs(dn2df['t'], max_lagtime_stubs*1.25)
            logging.info('filter_stubs: particle counts: %s to %s' % (dn2df['t']['particle'].nunique(),dn2df['t1']['particle'].nunique()))
            if t1['particle'].nunique()==0:
                logging.error('filter_stubs: particle counts =0; using less stringent conditions')
                dn2df['t1'] = tp.filter_stubs(dn2df['t'], max_lagtime_stubs*1)
        else:
            dn2df['t1'] = dn2df['t'].copy()

        if test:        
            fig=plt.figure()
            ax=plt.subplot(111)
            tp.mass_size(dn2df['t1'].groupby('particle').mean(),ax=ax);
            plt.tight_layout()
            savefig('%s.mass_size.png' % out_fh)        
        if flt_mass_size:
            dn2df['t2'] = dn2df['t1'][((dn2df['t1']['mass'] > dn2df['t1']['mass'].quantile(mass_cutoff)) & (dn2df['t1']['size'] < dn2df['t1']['size'].quantile(size_cutoff)) &
                     (dn2df['t1']['ecc'] < ecc_cutoff))]
            logging.info('filter_mass_size: particle counts: %s to %s' % (dn2df['t1']['particle'].nunique(),dn2df['t2']['particle'].nunique()))
            if len(t2)==0:
                dn2df['t2'] = dn2df['t1'].copy()
                logging.warning('filter_mass_size produced 0 particles; using t2=t1.copy()')
        else:
            dn2df['t2'] = dn2df['t1'].copy()
        if test:        
            fig=plt.figure()
            ax=plt.subplot(111)
            tp.mass_size(dn2df['t2'].groupby('particle').mean(),ax=ax);
            plt.tight_layout()
            savefig('%s.mass_size_post_filtering.png' % out_fh)        
        if flt_incomplete_trjs:
            dn2df['t2']=dn2df['t2'].reset_index()
            vals=pd.DataFrame(dn2df['t2']['particle'].value_counts())
            partis=[i for i in vals.index if vals.loc[i,'particle']>=int(vals.max())*0.95 ]
            dn2df['t2']=dn2df['t2'].loc[[i for i in dn2df['t2'].index if (dn2df['t2'].loc[i,'particle'] in partis)],:]
        _particles=dn2df['t2']['particle'].unique()
        df=dn2df['t2'].groupby('particle').agg({'frame':lambda x : len(unique_dropna(x))>1})
        particles=df.loc[df['frame'],:].index.tolist()
        dn2df['t2']=dn2df['t2'].loc[dn2df['t2']['particle'].isin(particles),:]
        logging.info(f"removed single frame particles: {len(_particles)} to {len(particles)}")
        if len(dn2df['t2'])==0:
            return None
        to_table(dn2df['t2'],dn2dp['t2'])
    else:
        dn2df['t2']=read_table(dn2dp['t2'])
    if test:
        for traj in ['t','t1','t2']:
            ax=plot_traj(frames[-1],traj=dn2df[traj])
        logging.info('getting plots hist')
        cols=['mass','size','ecc','signal','raw_mass','ep']
        fig=plt.figure()
        ax=plt.subplot(111)
        _=dn2df['t2'].loc[:,cols].hist(ax=ax)        
    return dn2df['t2']

def frames2coords_cor(frames,params_locate_start={'diameter':11,'minmass_percentile':92},
                      params_filter={},
                      params_link_df={},
                     out_fh=None,
                     params_msd={},
                      subtract_drift=False,
                      force=False):
    t_fltp=f'{out_fh}.t2.tsv'
    if not exists(t_fltp) or force:
        dbug=False
        if dbug:
            for diameter in np.arange(5,20,2):
                d=params_locate_start.copy()
                d['diameter']=diameter
                params_locate=get_params_locate(frames[0],out_fh=out_fh,**d,dbug=True)
                print(diameter, params_locate)
            brk

        params_locate=get_params_locate(frames[0],out_fh=out_fh,**params_locate_start)
        logging.info(params_locate)
        logging.info('getting coords')
        t_flt=frames2coords(frames=frames,out_fh=out_fh,
                            params_locate=params_locate,params_msd=params_msd,params_link_df=params_link_df,
                            force=force,**params_filter)        
        if t_flt is None:
            return None
    else:
        t_flt=pd.read_csv(t_fltp,sep='\t')
    if subtract_drift:
        d = tp.compute_drift(t_flt)
        t_cor = tp.subtract_drift(t_flt, d)
        return t_cor
    else:
        return t_flt        

from htsimaging.lib.io_data_files import read_pkl,to_pkl
def nd2msd(nd_fh,
           params_locate_start={'diameter':11,'minmass_percentile':92},
           params_msd={'mpp':0.0645,'fps':0.2, 'max_lagtime':100},
          get_coords=True,out_fh=None):
    frames=nd2frames(nd_fh)
    t_cor=frames2coords_cor(frames,params_locate_start=params_locate_start,params_msd=params_msd,out_fh=out_fh)
    # debug
    imsd=tp.imsd(t_cor,statistic='msd',**params_msd)
    emsd=tp.emsd(t_cor,**params_msd)
    emsd=pd.DataFrame(emsd)
    emsd.index.name='lag time [s]'
    # print emsd
    
    import statsmodels.tsa.stattools as st
    acf=pd.DataFrame(columns=emsd.columns)
    for c in emsd:
        acf[c]=st.acf(emsd.loc[:,c],nlags=len(emsd))
    acf.index=emsd.index
    
    if not out_fh is None:
        figure=plt.figure()
        ax=plt.subplot(111)
        acf.plot(ax=ax)
        savefig('%s.acf.pdf' % out_fh)
#         plt.clf()
    
    if not out_fh is None:
        imsd.to_csv(out_fh+".imsd")
        emsd.to_csv(out_fh+".emsd")
        acf.to_csv(out_fh+".acf")
        dpkl={}
        dpkl['imsd']=imsd
        dpkl['emsd']=emsd
        dpkl['acf']=acf
#         dpkl['t']=t
#         dpkl['t_flt']=t_flt
#         dpkl['f_batch']=f_batch
        dpkl['t_cor']=t_cor
        to_pkl(dpkl,out_fh+".pkl")
        if get_coords:
            t_cor.to_csv(out_fh+".coords")
    return imsd,emsd

def expt_dh2expt_info(expt_dh):
    try:
        expt_info=pd.read_csv(expt_dh+"/info")
    except:
        expt_info=pd.read_csv(expt_dh+"/info.csv")
    if "Unnamed" in expt_info.columns.tolist():
        cols_del.append(col)

    expt_info2=expt_info.set_index("smpn",drop=True).T
    expt_info2.columns=expt_info.loc[:,"smpn"]
    row_2_keep=[]
    for row in list(expt_info2.index):
        if "replicate" in row:
            row_2_keep.append(row)

    expt_info2=expt_info2.loc[row_2_keep,:]
    if len(expt_info2.index) != len(set(expt_info2.index)):
        print("duplicate values in index")
    for col in expt_info2.columns.tolist():
        for i in range(len(expt_info2)):
            # print expt_info2.loc["replicate %d" % (i+1),col]
            if not pd.isnull(expt_info2.loc["replicate %d" % (i+1),col]):
                expt_info2.loc[("replicate %d" % (i+1)),col]=['replicate %d' % (i+1),expt_info2.loc[("replicate %d" % (i+1)),col]]
            else:
                expt_info2.loc[("replicate %d" % (i+1)),col]=[np.nan,np.nan]
    expt_info2=expt_info2.reset_index()
    del expt_info2["index"]
    return expt_info2

def expt2plots(expt_info,expt_dh,_cfg={},
               test=False,force=False):
    if not exists(expt_dh):
        makedirs(expt_dh)
    expt_data=pd.DataFrame()
    for smpl in expt_info:
        smpl_info=expt_info.loc[:,smpl]
        smpl_data=pd.DataFrame()    
        for rep in smpl_info:
            if not pd.isnull(rep[0]):
                repn="%s %s" % (smpl,rep[0])
    #             print repn
                nd_fh=rep[1]
                if not exists(nd_fh):    
                    if nd_fh.count('/')<2:
                        nd_fh='%s/%s' % (expt_dh,nd_fh)
                        print(nd_fh)
                if exists(nd_fh):
                    out_fh="%s/%s" % (expt_dh,repn.replace(" ","_"))
    #                 print out_fh
                    plot_fh=out_fh+".imsd.pdf"
                    if not exists(plot_fh) or force:
                        # try:
                        imsd,emsd=nd2msd(nd_fh,out_fh=out_fh,params_msd=_cfg)
                        # except:
                        #     continue
                        emsd=pd.DataFrame(emsd)
                        print(repn)
                        emsd.columns=[repn]
                        imsd.index=emsd.index
                        plot_msd(imsd,emsd,scale='linear',plot_fh=plot_fh,params_msd=_cfg)
                    else:
                        logging.warning('not exists or already processed: %s' % basename(plot_fh))
                        # emsd=pd.read_csv(out_fh+".emsd")
                        # if "Unnamed: 0" in emsd.columns.tolist(): 
                        #     del emsd["Unnamed: 0"]
                    if test:
                        break
                else:
                    print("can not find")
#                 if len(smpl_data)==0:
#                     smpl_data=emsd
#                 else:
# #                 smpl_data[repn]=emsd[repn]
#                     smpl_data=pd.concat([smpl_data,emsd[repn]],axis=1)
            else:
                logging.error('null in info')
            # print smpl_data.columns
            # print smpl_data.index.name
            # smpl_data=set_index(smpl_data,col_index='lag time [s]')
            # smpl_data.to_csv(expt_dh+smpl+".emsd")
    #     if len(expt_data)==0:
    #         expt_data=smpl_data
    #     else:
    # #             expt_data.loc[:,smpl_data.columns.tolist()]=smpl_data.loc[:,smpl_data.columns.tolist()]
    #         expt_data=pd.concat([expt_data,
    #                              smpl_data.loc[:,[col for col in smpl_data.columns.tolist() if col != "lagt"]]],axis=1)
    # #     expt_data=expt_data.drop_duplicates("lagt",keep='first')
    # expt_data=set_index(expt_data,col_index='lagt')
    # expt_data.to_csv(expt_dh+"expt.emsd")
    # plot_fh=expt_dh+"emsd.pdf"
    # plot_emsd(expt_data,scale='linear',plot_fh=plot_fh)
    # expt_data_log10=np.log10(expt_data)
    # expt_data_log10.to_csv(expt_dh+"expt.emsdlog10")
    return expt_data


def createinfo(expt_dh):
    try:
        info=pd.read_csv(expt_dh+"info")
    except:
        info=pd.read_csv(expt_dh+"info.csv")
    for i in range(len(info)):
        reps=glob("%s/%s*" % (info.loc[i,"dh"],info.loc[i,"fn_lead"]))
        for repi,rep in enumerate(reps):
            info.loc[i,"replicate %d" % (repi+1)]=rep
    #     break
    info.to_csv(expt_dh+"info",index=False)
    
from htsimaging.lib.fit_kin import fit_power
def get_params(imsd,fit_type='power',out_fh=None):
    if fit_type=='power':
        parameters=pd.DataFrame(index=imsd.columns.tolist(),
                                columns=["power law exponent",
                                         "power law constant",
                                         "rsquared of power",
                                         "amplitude",
#                                          "slope",
#                                          "y intercept",
#                                          "rsquared of line",
                                        ])
        for col in imsd:
        #     parameters=tp.utils.fit_powerlaw(imsd.loc[:,col],plot=True)
            parameters.loc[col,"power law exponent"],\
            parameters.loc[col,"power law constant"], \
            parameters.loc[col,"rsquared of power"], \
            parameters.loc[col,"amplitude"]= \
            fit_power(np.array(list(imsd.index)),np.array(list(imsd.loc[:,col])))
        if not out_fh is None:
            parameters.to_csv(out_fh)
        return parameters

def msd2params(imsd,
            flt_amplitude=True,mn_traj=3,
            out_fh=None,
            mpp=0.0645,
            fps=0.2,
            max_lagtime=100,
            flt_traj=False):
    
    l=fps*max_lagtime
    # print l #debug
    l=int(l)
    # for l in range(10,110,10):
    if not out_fh is None:
        params_fh='%s.params' % out_fh
        if flt_traj:
            params_flt_fh='%s.params_flt' % out_fh
            imsd_flt_fh='%s.imsd_flt' % out_fh
    parameters=get_params(imsd.head(l),fit_type='power',out_fh=None)
    parameters=parameters.T
    parameters.index.name='params'
    # print parameters
    parameters.to_csv(params_fh)

    if flt_traj:
        traj_flt1=parameters.loc[(parameters.loc[:,'power law exponent']<1) & \
                                (parameters.loc[:,'rsquared of power']>0.95) \
                                ,:].index.tolist()
        parameters=parameters.loc[traj_flt1,:]
    #     print parameters.head()
        #debug
        imsd.loc[:,traj_flt1].head(l).to_csv('imsd_traj_flt1.csv')    
        if flt_amplitude:
            for i in np.arange(0,1.1,0.1)[::-1]:
                traj_flt=parameters.loc[(parameters.loc[:,'amplitude']<(parameters.loc[:,'amplitude'].mean()+i*parameters.loc[:,'amplitude'].std())) & \
                                        (parameters.loc[:,'amplitude']>(parameters.loc[:,'amplitude'].mean()-i*parameters.loc[:,'amplitude'].std())) \
                                        ,:].index.tolist()
                if len(traj_flt)<mn_traj:
                    if 'traj_flt_prev' in locals():
                        traj_flt=traj_flt_prev
                    else:
                        traj_flt=traj_flt1
                    print('filtered: %s : %s' % (len(traj_flt1),len(traj_flt)))
                    break
                traj_flt_prev=traj_flt                
        else:
            traj_flt=traj_flt1
        print('filtered: %s : %s' % (len(traj_flt1),len(traj_flt)))
        imsd_flt=imsd.loc[:,traj_flt].head(l)
    #     emsd_flt=pd.DataFrame(imsd_flt.T.mean())
        params_flt=parameters.loc[traj_flt,:]
        if not out_fh is None:
            imsd_flt.to_csv(imsd_flt_fh)
            params_flt.to_csv(params_flt_fh)
        return imsd_flt,params_flt
    
# from htsimaging.lib.spt import frames2coords_cor    
def cellframes2distances(cellframes,cellframesmasked,out_fh=None,test=False,force=False):
    makedirs(dirname(out_fh),exist_ok=True)
    params_msd={'mpp':0.0645,'fps':0.2, 'max_lagtime':100}
    # for 150x150 images
    params_locate_start={'diameter':7,'minmass_percentile':90} 
    # round to odd number
    # diameter for cellboxdth=100: 7
    # diameter for cellboxdth=150: 7    
    params_link_df={'search_range':5,'memory':0,'link_strategy':'drop',}
    params_filter={'mass_cutoff':0.5,'size_cutoff':1,'ecc_cutoff':1,
                  'filter_stubs':False,'flt_mass_size':False,'flt_incomplete_trjs':False,
                  'test':test}
    makedirs(dirname(out_fh),exist_ok=True)
    t_cor=frames2coords_cor(frames=cellframesmasked,out_fh=out_fh,
                            params_locate_start=params_locate_start,
                            params_msd=params_msd,params_link_df=params_link_df,
                            params_filter=params_filter,
                            subtract_drift=True,
                            force=force)
    if t_cor is None:
        return None
    ddist=get_distance_travelled(frames=cellframesmasked,t_cor=t_cor,out_fh=out_fh,test=test,force=force)
    if not (out_fh is None or ddist is None):
        make_gif(cellframes,ddist,f"{dirname(out_fh)}/vid",force=force)    
        
def apply_cellframes2distances(cellcfgp):
    """
    wrapper around cellframes2distances for multiprocessing
    """
    celli=dirname(cellcfgp)
    print(celli);logging.info(celli)
    cellcfg=yaml.load(open(cellcfgp,'r'))
    cellframes2distances([np.load(p) for p in sorted(cellcfg['cellframeps'])],
                         [np.load(p) for p in sorted(cellcfg['cellframesmaskedps'])],
                         out_fh=f"{cellcfg['outp']}/plot_check",
                         test=cellcfg['test'],force=cellcfg['force'])
        