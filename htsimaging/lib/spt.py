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

# def get_params_locate(frame,diameter=15,minmass_percentile=92,outp=None,test=True,figsize=None,dbug=False):
#     if dbug:
#         outp=f"{outp}_diameter_{diameter:02d}"
#     f = tp.locate(frame, diameter, invert=False)
#     minmass=np.percentile(f['mass'],minmass_percentile)
#     logging.info('feature count= %s, %spercentile= %s'  % (len(f),minmass_percentile,minmass))
                 
#     f = tp.locate(frame, diameter, invert=False, minmass=np.percentile(f['mass'],minmass_percentile))
#     logging.info('feature count= %s, %spercentile= %s'  % (len(f),minmass_percentile,
#                                                            np.percentile(f['mass'],minmass_percentile)))
    
#     if test:
#         if not outp is None:
#             fig=plt.figure(figsize=figsize)
#             ax=plt.subplot(111)
#             ax=tp.annotate(f, frame,ax=ax)
#             savefig(f'{outp}/image_get_params_locate.png')            
#     params_locate={'diameter':diameter,
#                   'minmass':minmass}
#     return params_locate

def plot_traj(frame,traj):
    fig=plt.figure(figsize=[frame.shape[0]*0.02,frame.shape[1]*0.02])        
    ax=plt.subplot(111)
    ax.imshow(frame,cmap='binary_r',alpha=0.8)
    ax = tp.plot_traj(traj,label=False,ax=ax,lw=2)
    plt.tight_layout()
    return ax

def frames2coords(frames,outp,
                  params_locate,params_msd,params_link_df={'search_range':20,},
                  mass_cutoff=0.5,size_cutoff=0.5,ecc_cutoff=0.5,
                      subtract_drift=False,
                    filter_stubs=True,flt_mass_size=True,flt_incomplete_trjs=True,
                    force=False,test=False):
    dns=['f_batch','t','t1','t2']
    dn2dp={dn:f'{outp}.{dn}.tsv' for dn in dns}
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
            savefig('%s.mass_size.png' % outp)        
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
            savefig('%s.mass_size_post_filtering.png' % outp)        
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
#     return dn2df['t2']

# def frames2coords_cor(frames,params_locate_start={'diameter':11,'minmass_percentile':92},
#                       params_filter={},
#                       params_link_df={},
#                      outp=None,
#                      params_msd={},
#                       force=False):
#     t_fltp=f'{outp}.t2.tsv'
#     if not exists(t_fltp) or force:
#         logging.info(params_locate)
#         logging.info('getting coords')
#         t_flt=frames2coords(frames=frames,outp=outp,
#                             params_locate=params_locate,
#                             params_msd=params_msd,
#                             params_link_df=params_link_df,
#                             force=force,**params_filter)        
#         if t_flt is None:
#             return None
#     else:
#         t_flt=pd.read_csv(t_fltp,sep='\t')
    if subtract_drift:
        d = tp.compute_drift(dn2df['t2'])
        dn2df['t2_corrected'] = tp.subtract_drift(dn2df['t2'], d)
        return dn2df['t2_corrected']
    else:
        return dn2df['t2']
    
def test_locate_particles(cellcfg,params_locate,force=False):
        # test locate
    if exists(f"{cellcfg['outp']}/dlocate.tsv") and not force:
        return
    from htsimaging.lib.plot import dist_signal
    frame =np.load(cellcfg['cellgfpmaxp'])
    df1 = tp.locate(frame, **params_locate)
    df1['particle']=df1.index
    to_table(df1,f"{cellcfg['outp']}/dlocate_test.tsv")
    
    print(f"particles detected ={len(df1)}")
    # plot dist signal of the detected particles
    fig=plt.figure()
    ax=plt.subplot()
    dist_signal(np.load(cellcfg['cellgfpminp']),
                params_hist={'bins':20,'label':'gfp min','density':True,'color':'k'},ax=ax)
    dist_signal(np.load(cellcfg['cellgfpmaxp']),
                params_hist={'bins':20,'label':'gfp max','density':True,'color':'green'},ax=ax)
    dist_signal(df1['signal'],
                threshold=cellcfg['signal_cytoplasm'],label_threshold='signal_cytoplasm',
                params_hist={'bins':20,'label':'particles','density':True,'color':'r'},ax=ax)
    savefig(f"{cellcfg['plotp']}/dist_signal_locate_particles.png")  
    
    # plot image of the detected particles
    from htsimaging.lib.plot import image_locate_particles
    image_locate_particles(df1,frame)
    savefig(f"{cellcfg['plotp']}/image_locate_particles.png")  
            
    from htsimaging.lib.plot import plot_properties_cell
    plot_properties_cell(cellcfg,df1,cols_colorby=df1.select_dtypes('float').columns.tolist())
    savefig(f"{cellcfg['plotp']}/plot_properties_cell_locate_particles.png")
    
    
# from htsimaging.lib.spt import frames2coords_cor    
def cellcfg2distances(cellcfg,
                    # for 150x150 images
                    params={'locate':{'diameter':15, # round to odd number
                                      'noise_size':1,
                                      'separation':7,
                                      'threshold':4000,
                                      'preprocess':True,
                                      'invert':False,
                                      'max_iterations':50,
                                      'percentile':0,
                                      'engine':'numba',
                                      },
                    'link_df':{'search_range':5,
                               'memory':0,
                               'link_strategy':'drop',},
                    'filter_stubs':{'threshold':3},
#                     'filter':{'mass_cutoff':0.5,
#                               'size_cutoff':1,
#                               'ecc_cutoff':1,
#                               'filter_stubs':False,
#                               'flt_mass_size':False,
#                               'flt_incomplete_trjs':False,
# #                                   'test':test
#                                   },
#                     'msd':{'mpp':0.0645,'fps':0.2, 'max_lagtime':100},
                           },
                    test=False,force=False):
    params['locate']['threshold']=cellcfg['signal_cytoplasm']
    params['filter']['test']=test
    params['link_df']['search_range']=params['locate']['separation']*0.5
            
    to_dict(params,f"{cellcfg['outp']}/params.yml")
    
    test_locate_particles(cellcfg,params['locate'],force=force)
    return
    # get trajectories
    steps=['locate','link','filter_stubs','subtract_drift','distance']
    dn2dp={s:f"{cellcfg['outp']}/d{si}{s}.tsv" for si,s in enumerate(steps)}
    dn2plotp_suffix={s:f"{si}{s}.png" for si,s in enumerate(steps)}
    dn2df={}
    dn2df['locate']=tp.batch([np.load(p) for p in sorted(cellcfg['cellframesmaskedps'])],
                             **params['locate'])
    dn2df['link']=tp.link_df(dn2df['locate'], **params['link_df'])

    img_gfp=np.load(cellcfg['cellgfpmaxp'])
    img_bright=np.load(cellcfg['cellbrightp'])
#     %run ../htsimaging/htsimaging/lib/plot.py
    from htsimaging.lib.plot import image_trajectories
    image_trajectories(dtraj=dn2df['link'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['link']}")

    dn2df['filter_stubs']=tp.filter_stubs(dn2df['link'], threshold=params['filter_stubs']['threshold'])
    dn2df['filter_stubs'].index.name='index'
    dn2df['filter_stubs'].index=range(len(dn2df['filter_stubs']))
    image_trajectories(dtraj=dn2df['filter_stubs'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['filter_stubs']}")

    d = tp.compute_drift(dn2df['filter_stubs'])
    dn2df['subtract_drift'] = tp.subtract_drift(dn2df['filter_stubs'], d)
    image_trajectories(dtraj=dn2df['subtract_drift'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['subtract_drift']}")

    from htsimaging.lib.stat import get_distance_travelled
    dn2df['distance']=get_distance_travelled(t_cor=dn2df['subtract_drift'])

    if not (outp is None or ddist is None):
        make_gif([np.load(p) for p in sorted(cellcfg['cellframeps'])],ddist,f"{dirname(outp)}/vid",force=force)    
        
def apply_cellcfg2distances(cellcfgp):
    """
    wrapper around cellcfg2distances for multiprocessing
    """
    celli=dirname(cellcfgp)
    print(celli);logging.info(celli)
    cellcfg=yaml.load(open(cellcfgp,'r'))
    cellcfg2distances(cellcfg,
                      test=cellcfg['test'],
                      force=cellcfg['force'])