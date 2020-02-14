from os.path import splitext, join, exists, isdir,basename,abspath,dirname
from rohan.global_imports import *
import trackpy as tp
tp.ignore_logging()
# import nd2reader
import pims
import pims_nd2
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
    
def test_locate_particles(cellcfg,params_locate,force=False):
        # test locate
    dlocate_testp=f"{cellcfg['outp']}/dlocate_test.tsv"
    if exists(dlocate_testp) and not force:
        return
    from htsimaging.lib.plot import dist_signal
    frame =np.load(cellcfg['cellgfpmaxp'])
    df1 = tp.locate(frame, **params_locate)
    df1['particle']=df1.index
    to_table(df1,dlocate_testp)
    logging.info(f"particles detected ={len(df1)}")
    if len(df1)==0:
        return False
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
    image_locate_particles(df1,
                           frame=frame,
                           img_region=np.load(cellcfg['cellbrightp']))
    savefig(f"{cellcfg['plotp']}/image_locate_particles.png")  
    if len(df1)>=5:
        from htsimaging.lib.plot import plot_properties_cell
        plot_properties_cell(cellcfg,df1,cols_colorby=df1.select_dtypes('float').columns.tolist())
        savefig(f"{cellcfg['plotp']}/plot_properties_cell_locate_particles.png")
    return True

# df0=pd.DataFrame({'step name':steps,
# 'step #':range(len(steps)),}).set_index('step #')
# df0['dfp']=df0.apply(lambda x:f"{cellcfg['outp']}/d{'_'.join(df0.loc[range(x.name+1),'step name'].values.tolist())}.tsv" ,axis=1)
# df0['plotp suffix']=df0.apply(lambda x:f"_{'__'.join(df0.loc[range(x.name+1),'step name'].values.tolist())}.png" ,axis=1)
# cellcfg['track particles']=df0.set_index('step name')['dfp'].apply(lambda x: f"{basenamenoext(x)}p").to_dict()            
# class track_particles():            
#     def locate(cellcfg,params,dfp,plotp_suffix):
            
#         df1=tp.batch([np.load(p) for p in sorted(cellcfg['cellframesmaskedps'])],
#                                  **params)
#         to_table(df1,dfp)
#     def link(cellcfg,params,dfp,plotp_suffix):
#         df1=tp.link_df(df, **params)
#         to_table(cellcfg['']df1,dfp)
            
#         image_trajectories(dtraj=dn2df[step], 
#                            img_gfp=img_gfp, 
#                            img_bright=img_bright, fig=None, ax=None)
#         savefig(f"{cellcfg['plotp']}/image_trajectories_{plotp_suffix}")
#     def filter_stubs(cellcfg,params,dfp,plotp_suffix):
#         df1=tp.filter_stubs(dn2df['link'], threshold=params['filter_stubs']['threshold'])
#         df1.index.name='index'
#         df1.index=range(len(df1))
#         to_table(df1,dfp)
                
#         image_trajectories(dtraj=df1, 
#                            img_gfp=img_gfp, 
#                            img_bright=img_bright, fig=None, ax=None)
#         savefig(f"{cellcfg['plotp']}/image_trajectories_{plotp_suffix}")

#     def subtract_drift(cellcfg,params,dfp,plotp_suffix):
#         df1 = tp.subtract_drift(df, tp.compute_drift(df))
#         to_table(df1,dfp)
                
#         image_trajectories(dtraj=df1, 
#                            img_gfp=img_gfp, 
#                            img_bright=img_bright, fig=None, ax=None)
#         savefig(f"{cellcfg['plotp']}/image_trajectories_{plotp_suffix}")

#     def distance(cellcfg,params,dfp,plotp_suffix):
#         from htsimaging.lib.stat import get_distance_travelled
#         df1=get_distance_travelled(t_cor=df)
#         to_table(df1,dfp)
        
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
            
def cellcfg2distances(cellcfg,
                    # for 150x150 images
                    params={'locate':{'diameter':15, # round to odd number
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
                               'memory':0,
                               'link_strategy':'drop',},
                    'filter_stubs':{'threshold':4},
                    'get_distance_from_centroid':{'center':[75,75]},
#                     'msd':{'mpp':0.0645,'fps':0.2, 'max_lagtime':100},
                           },
                    test=False,force=False):
    params['locate']['separation']=params['locate']['diameter']*1.25
    params['locate']['threshold']=cellcfg['signal_cytoplasm']
    params['link_df']['search_range']=params['locate']['diameter']*0.33
            
    to_dict(params,f"{cellcfg['outp']}/params.yml")
    
    if not test_locate_particles(cellcfg,params['locate'],force=force):
        print(cellcfg['cfgp'])
        return 
    # get trajectories
    steps=['locate','link_df','filter_stubs','subtract_drift','distance']
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
    dn2df['locate']=tp.batch([np.load(p) for p in sorted(cellcfg['cellframesmaskedps'])],
                             **params['locate'])
    to_table(dn2df['locate'],dn2dp['locate'])

    dn2df['link_df']=tp.link_df(dn2df['locate'], **params['link_df'])
    to_table(dn2df['link_df'],dn2dp['link_df'])
    image_trajectories(dtraj=dn2df['link_df'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None)
    savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['link_df']}")
    to_table(dn2df['link_df'],dn2dp['link_df'])

    
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
    savefig(f"{cellcfg['plotp']}/image_trajectories_pre_{dn2plotp_suffix['filter_stubs']}")
    
    dn2df['filter_stubs']=get_distance_from_centroid(dn2df['filter_stubs'],**params['get_distance_from_centroid'])
    dn2df['filter_stubs']=trim_returns(dn2df['filter_stubs'])
    savefig(f"{cellcfg['plotp']}/image_trajectories_stats_trimming_{dn2plotp_suffix['filter_stubs']}")

    dn2df['filter_stubs']=tp.filter_stubs(dn2df['filter_stubs'], threshold=params['filter_stubs']['threshold'])
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

    d = tp.compute_drift(dn2df['filter_stubs'])
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

    make_gif(cellcfg,
             force=force)    
        
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