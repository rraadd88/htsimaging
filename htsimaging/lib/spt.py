from os.path import splitext, join, exists, isdir,basename,abspath,dirname
from roux.global_imports import *
import trackpy as tp
tp.ignore_logging()
# import nd2reader
import pims
import logging
# from htsimaging.lib.plot import *
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
    
def test_locate_particles(cellcfg,params_locate,frame=None,force=False,test=False):
    dlocate_testp=f"{cellcfg['outp']}/dlocate_test.tsv"
    if exists(dlocate_testp) and not force and not test:
        return True
    cellgfpmax=np.load(cellcfg['cellgfpmaxp'])
    frame=cellgfpmax if frame is None else frame
    df1 = tp.locate(frame, 
                    **params_locate)
    df1['particle']=df1.index
    if not test:
        to_table(df1,dlocate_testp)
    logging.info(f"particles detected ={len(df1)}") if not test else print(f"particles detected ={len(df1)}")
    if len(df1)==0:
        return False
    # plot dist signal of the detected particles
    fig=plt.figure()
    cellgfpmin=np.load(cellcfg['cellgfpminp'])
    from htsimaging.lib.plot import dist_signal
    ax=plt.subplot()
#     dist_signal(np.unique(cellgfpmin)[2:],
#                 params_hist={'bins':20,'label':'gfp min',
#                              'density':True,'color':'k'},ax=ax)
#     dist_signal(np.unique(cellgfpmax)[2:],
#                 params_hist={'bins':20,'label':'gfp max',
#                              'density':True,'color':'green'},ax=ax)
    dist_signal(np.unique(frame)[2:],
                params_hist={'bins':20,'label':'frame',
                             'density':True,'color':'lime'},ax=ax)        
    dist_signal(df1['signal'],
                threshold=np.unique(frame)[1],label_threshold='signal_cytoplasm',
                params_hist={'bins':20,'label':f'particles\n(total ={len(df1)})',
                             'density':True,'color':'r'},ax=ax)
    if not test:
        savefig(f"{cellcfg['plotp']}/dist_signal_locate_particles.png")
    # plot image of the detected particles
    from htsimaging.lib.plot import image_locate_particles
    ax=image_locate_particles(df1,
                           frame=frame,
                           img_region=np.load(cellcfg['cellbrightp']),
                           annotate_particles=False)
    _=ax.text(0,1,'\n'.join([f"{k}:{params_locate[k]}" for k in params_locate]),
              va='top',color='lime')
    if not test:
        savefig(f"{cellcfg['plotp']}/image_locate_particles.png")  
#     if len(df1)>=5:
#         from htsimaging.lib.plot import plot_properties_cell
#         plot_properties_cell(cellcfg,df1,cols_colorby=df1.select_dtypes('float').columns.tolist()) 
#         if not test:
#             savefig(f"{cellcfg['plotp']}/plot_properties_cell_locate_particles.png")
    return True
        
def trim_returns(df1):
    from htsimaging.lib.stat import get_inflection_point
    from htsimaging.lib.plot import plot_trajectories_stats
    # get inflection point if any
    df1=df1.groupby('particle').apply(lambda x: get_inflection_point(x))
    df2=df1.groupby('particle',as_index=False).apply(lambda df: df.set_index('frame').sort_index(axis=0).loc[df['inflection point'].unique()[0]:,:].reset_index())


    fig=plt.figure(figsize=[20,10])
    # plot before filtering
    ax=plt.subplot()
    df1.loc[df1['inflection point'].isnull(),:].groupby('particle').apply(lambda x: plot_trajectories_stats(x,
                                                            coly='distance effective from centroid per frame',
                                                            colx='frame',
                                                           fig=fig,ax=ax,
                                                           rescalex=False,
                                                          params_plot={'color':'lime','alpha':0.5}))
    df1.loc[~df1['inflection point'].isnull(),:].groupby('particle').apply(lambda x: plot_trajectories_stats(x,
                                                             coly='distance effective from centroid per frame',
                                                             colx='frame',
                                                           fig=fig,ax=ax,
                                                           rescalex=False,
                                                            axvlinex=x['inflection point'].unique()[0],
                                                          params_plot={'color':'r','alpha':0.75}))

    df2.loc[~df2['inflection point'].isnull(),:].groupby('particle').apply(lambda x: plot_trajectories_stats(x,
                                                             coly='distance effective from centroid per frame',
                                                             colx='frame',
                                                           fig=fig,ax=ax,
                                                           rescalex=False,
                                                            axvlinex=x['inflection point'].unique()[0],
                                                          params_plot={'color':'g','alpha':1}))
    ax.get_legend().remove()
    ax.set_ylim(df1['distance effective from centroid per frame'].min(),df1['distance effective from centroid per frame'].max())               
    ax.set_xlim(df1['frame'].min(),df1['frame'].max())
    return df2
                
def fill_frame_jumps(df1,jump_length):
    df1=df1.sort_values(by=['particle','frame'])
    df1['frame delta']=df1['frame']-([df1['frame'].tolist()[0]-1]+df1['frame'].tolist()[:-1])
    df=df1.loc[(df1['frame delta']==jump_length),:]
    df.loc[df.index,'frame']=df['frame']-1
    df1=pd.concat({'not':df1,'fill':df},axis=0)
    df1.index.name='frame fill or not'
    df1=df1.reset_index()
    df1=df1.sort_values(by=['particle','frame'])
    df1['frame delta']=df1['frame']-([df1['frame'].tolist()[0]-1]+df1['frame'].tolist()[:-1])
    logging.warning(sum((df1['frame']-([df1['frame'].tolist()[0]-1]+df1['frame'].tolist()[:-1]))==jump_length)==0)
    df1.index=range(len(df1))
    return df1                
            
def cellcfg2distances(cellcfg,
                    # for 150x150 images
                    params={'locate':{'diameter':11, # round to odd number
                                      'noise_size':1,
                                      'separation':15,
                                      'threshold':4000,
                                      'preprocess':True,
                                      'invert':False,
                                      'max_iterations':50,
                                      'percentile':0,
                                      'engine':'numba',
                                      },
                    'link_df':{
                               'search_range':5,
                               'memory':1,
                               'link_strategy':'drop',},
                    'filter_stubs':{'threshold':4},
                    'get_distance_from_centroid':{'center':[75,75]},
#                     'msd':{'mpp':0.0645,'fps':0.2, 'max_lagtime':100},
                           },
                    test=False,force=False):
    params['locate']['separation']=params['locate']['diameter']*1
    params['locate']['threshold']=cellcfg['signal_cytoplasm']*0.5
    params['link_df']['search_range']=params['locate']['diameter']*0.33
            
    to_dict(params,f"{cellcfg['outp']}/params.yml")
    
    if not test_locate_particles(cellcfg,params['locate'],force=force,test=False):
        print(cellcfg['cfgp'])
        return 
    # get trajectories
    steps=['locate','link_df','filter_stubs','filter_returns','subtract_drift','distance']
    dn2dp={s:f"{cellcfg['outp']}/d{si}{s}.tsv" for si,s in enumerate(steps)}
    dn2plotp_suffix={s:f"{si}{s}.png" for si,s in enumerate(steps)}
    steps_done=[k for k in dn2dp if exists(dn2dp[k])]
    
    if ('distance' in steps_done) and not force:
       print(cellcfg['cfgp'])
       return
           
    from htsimaging.lib.plot import image_trajectories
    from htsimaging.lib.stat import get_distance_from_centroid
           
    img_gfp=np.load(cellcfg['cellgfpmaxp'])
    img_bright=np.load(cellcfg['cellbrightp'])
    
    dn2df={}
    dn2df['locate']=tp.batch([np.load(p) for p in sorted(cellcfg['cellframes_masked_substracted'])],
                             **params['locate'])
    if len(dn2df['locate'])==0:
        return
    dn2df['locate']['frame']=dn2df['locate']['frame'].astype(np.integer)
    dn2df['link_df']=tp.link_df(dn2df['locate'], **params['link_df'])
#     if params['link_df']['memory']!=0:
    dn2df['link_df']=fill_frame_jumps(dn2df['link_df'],
                      jump_length=2 if params['link_df']['memory']==0 else params['link_df']['memory']+1)
#     to_table(dn2df['link_df'],'test.tsv')
#     to_table(dn2df['link_df'],dn2dp['link_df'])
    image_trajectories(dtraj=dn2df['link_df'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['link_df']}")
#     to_table(dn2df['link_df'],dn2dp['link_df'])

    
    dn2df['filter_stubs']=tp.filter_stubs(dn2df['link_df'], threshold=params['filter_stubs']['threshold'])
    dn2df['filter_stubs'].index.name='index'
    dn2df['filter_stubs'].index=range(len(dn2df['filter_stubs']))
    if len(dn2df['filter_stubs'])==0:
        to_table(dn2df['filter_stubs'],dn2dp['distance'])
        print(cellcfg['cfgp'])
        return 
    image_trajectories(dtraj=dn2df['filter_stubs'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['filter_stubs']}")
    
    dn2df['filter_returns']=get_distance_from_centroid(dn2df['filter_stubs'],**params['get_distance_from_centroid'])
    dn2df['filter_returns']=trim_returns(dn2df['filter_returns'])
    savefig(f"{cellcfg['plotp']}/image_trajectories_stats_trimming_{dn2plotp_suffix['filter_returns']}")

    dn2df['filter_returns']=tp.filter_stubs(dn2df['filter_returns'], threshold=params['filter_stubs']['threshold'])
    dn2df['filter_returns'].index.name='index'
    dn2df['filter_returns'].index=range(len(dn2df['filter_returns']))
    if len(dn2df['filter_returns'])==0:
        to_table(dn2df['filter_returns'],dn2dp['distance'])
        print(cellcfg['cfgp'])
        return             
    image_trajectories(dtraj=dn2df['filter_stubs'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['filter_returns']}")

    d = tp.compute_drift(dn2df['filter_returns'])
    dn2df['subtract_drift'] = tp.subtract_drift(dn2df['filter_stubs'], d)
    image_trajectories(dtraj=dn2df['subtract_drift'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['subtract_drift']}")

    dn2df['distance']=get_distance_from_centroid(dn2df['subtract_drift'],**params['get_distance_from_centroid'])
    from htsimaging.lib.stat import get_distance_travelled
    dn2df['distance']=get_distance_travelled(t_cor=dn2df['distance'])
    
    for k in dn2df:
        to_table(dn2df[k],dn2dp[k])
#     make_gif(cellcfg,
#              force=force)    
        
def apply_cellcfgp2distances(cellcfgp):
    """
    wrapper around cellcfg2distances for multiprocessing
    """
    celli='/'.join(dirname(cellcfgp).split('/')[-3:])
    print(f"processing {celli}: ");logging.info(celli)
    cellcfg=yaml.load(open(cellcfgp,'r'))
    cellcfg2distances(cellcfg,
                      test=cellcfg['test'],
                      force=cellcfg['force'])
