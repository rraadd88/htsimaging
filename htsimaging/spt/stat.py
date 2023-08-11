#!/usr/bin/env python
"""Statistical analysis of the single particle tracking."""

import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.INFO) # 

from os import makedirs
from os.path import exists, dirname

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from roux.lib.df import flatten_columns
from roux.lib.io import read_dict, to_dict, to_table

import trackpy as tp
tp.ignore_logging()

from scipy.spatial import distance

def test_locate_particles(
    cellcfg: dict,
    params_locate: dict,
    frame=None,
    force: bool=False,
    test: bool=False,
    ) -> bool:
    """
    Test locating of the particles.

    Args:
        cellcfg (dict): the cell level configuration. 
        params_locate (dict): parameters provided for the location.
        frame (np.array, optional): image frame. Defaults to None.
        force (bool, optional): over-write the outputs. Defaults to False.
        test (bool, optional): test mode. Defaults to False.

    Returns:
        bool
    """
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
    from htsimaging.viz.stat import dist_signal
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
        ax.set(title='QC: Particle detection')
        plt.savefig(f"{cellcfg['plotp']}/dist_signal_locate_particles.png")
    # plot image of the detected particles
    from htsimaging.endocytosis.viz import image_locate_particles
    ax=image_locate_particles(df1,
                           frame=frame,
                           img_region=np.load(cellcfg['cellbrightp']),
                           annotate_particles=False)
    _=ax.text(0,1,'\n'.join([f"{k}:{params_locate[k]}" for k in params_locate]),
              va='top',color='lime')
    if not test:
        ax.set(title='QC: Particle detection')
        plt.savefig(f"{cellcfg['plotp']}/image_locate_particles.png")  
    return True

def to_msd(
    frames: list,
    coff_intesity_perc : float = 75,
    diameter=11,
    cores : int = 4,
    test : bool = False,
    ) -> tuple:
    """
    MSD from the nd file.

    Args:
        frames (str): 2D frames.
        
        
    Returns:
        tuple: outputs.
    """
    logging.info('number of frames = %d' % len(frames))
    assert len(np.shape(frames))<4, f"{np.shape(frames)} perhaps use average_z(frames)"
    threshold=np.percentile(np.hstack(tuple(frames)),coff_intesity_perc)
    f_batch = tp.batch(frames,diameter=diameter,threshold=threshold,processes=cores)
    if test:
        plt.figure()
        ax=tp.annotate(f_batch.query("`frame` == 0"), frames[0])
        ax.set(title="Test: Detection of the particles")
    t = tp.link_df(f_batch, search_range=diameter, memory=3)
    t_flt = tp.filter_stubs(t, 3*int(len(frames)/4))
    try:
        d = tp.compute_drift(t_flt)
        t_cor = tp.subtract_drift(t_flt, d)
    except:
        t_cor=t_flt
        logging.info("drift correction excepted")
    imsd=tp.imsd(t_cor,0.1,0.2, max_lagtime=100, statistic='msd')
    emsd=tp.emsd(t_cor,0.1,0.2, max_lagtime=100)
    return imsd,emsd    
        
def trim_returns(
    df1: pd.DataFrame,
    ) -> pd.DataFrame:
    """
    Trim images.

    Args:
        df1 (pd.DataFrame): input dataframe.

    Returns:
        pd.DataFrame: output dataframe.
    """
    from htsimaging.lib.stat import get_inflection_point
    from htsimaging.spt.viz import plot_trajectories_stats
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
    ax.set(title='Returns trimmed')
    return df2
                
def fill_frame_jumps(
    df1: pd.DataFrame,
    jump_length,
    ) -> pd.DataFrame:
    """
    Fill the frame jumps.

    Args:
        df1 (pd.DataFrame): input dataframe.
        jump_length (_type_): length of the jump.

    Returns:
        pd.DataFrame: output dataframe.
    """
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
            
def cellcfg2distances(
    cellcfg: dict,
    params: dict,
    subtract_drift: bool=False,
    test: bool=False,
    force: bool=False,
    ):
    """
    Calculate distances from cell configuration.

    Args:
        cellcfg (dict): configuration
        params (_type_, optional): parameters. Defaults to { 'locate':{'diameter':11, # round to odd number 'noise_size':1, 'separation':15, 'threshold':4000, 'preprocess':True, 'invert':False, 'max_iterations':50, 'percentile':0, 'engine':'numba', }, 'link_df':{ 'search_range':5, 'memory':1, 'link_strategy':'drop',}, 'filter_stubs':{'threshold':4}, 'get_distance_from_centroid':{'center':[75,75]}, }.
        force (bool, optional): over-write the outputs. Defaults to False.

    """
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
           
    from htsimaging.endocytosis.viz import image_trajectories
    from htsimaging.spt.stat import get_distance_from_centroid
           
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
    image_trajectories(
        dtraj=dn2df['link_df'], 
        img_gfp=img_gfp, 
        img_bright=img_bright,
        fig=None,
        ax=None,
        ).set(title='Particles linked')
    makedirs(cellcfg['plotp'],exist_ok=True)
    plt.savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['link_df']}")
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
                       img_bright=img_bright, fig=None, ax=None).set(title='Stubs removed')
    plt.savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['filter_stubs']}")
    
    dn2df['filter_returns']=get_distance_from_centroid(dn2df['filter_stubs'],**params['get_distance_from_centroid'])
    dn2df['filter_returns']=trim_returns(dn2df['filter_returns'])
    plt.savefig(f"{cellcfg['plotp']}/image_trajectories_stats_trimming_{dn2plotp_suffix['filter_returns']}")

    dn2df['filter_returns']=tp.filter_stubs(dn2df['filter_returns'], threshold=params['filter_stubs']['threshold'])
    dn2df['filter_returns'].index.name='index'
    dn2df['filter_returns'].index=range(len(dn2df['filter_returns']))
    if len(dn2df['filter_returns'])==0:
        to_table(dn2df['filter_returns'],dn2dp['distance'])
        print(cellcfg['cfgp'])
        return             
    image_trajectories(dtraj=dn2df['filter_stubs'], 
                       img_gfp=img_gfp, 
                       img_bright=img_bright, fig=None, ax=None).set(title='Filtered returns')
    plt.savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['filter_returns']}")
    if subtract_drift:
        d = tp.compute_drift(dn2df['filter_returns'])
        dn2df['subtract_drift'] = tp.subtract_drift(dn2df['filter_stubs'], d)
        image_trajectories(dtraj=dn2df['subtract_drift'], 
                           img_gfp=img_gfp, 
                           img_bright=img_bright, fig=None, ax=None).set(title='Drift substracted')
        plt.savefig(f"{cellcfg['plotp']}/image_trajectories_{dn2plotp_suffix['subtract_drift']}")

    dn2df['distance']=get_distance_from_centroid(dn2df['subtract_drift' if subtract_drift else 'filter_returns'],**params['get_distance_from_centroid'])
    dn2df['distance']=get_distance_travelled(t_cor=dn2df['distance'])
    
    for k in dn2df:
        to_table(dn2df[k].reset_index(drop=True),dn2dp[k])
        
def apply_cellcfgp2distances(
    cellcfgp: str,
    ):
    """
    Wrapper around cellcfg2distances for multiprocessing.

    Args:
        cellcfgp (str): path to the configuration file.
    """
    celli='/'.join(dirname(cellcfgp).split('/')[-3:])
    print(f"processing {celli}: ");logging.info(celli)
    cellcfg=read_dict(cellcfgp)
    cellcfg2distances(cellcfg,
                      test=cellcfg['test'],
                      force=cellcfg['force'])
    

def get_distance_from_centroid(
    df1: pd.DataFrame,
    center: list=[75,75],
    ) -> pd.DataFrame:
    """
    Get distance from the centroid.

    Args:
        df1 (pd.DataFrame): input dataframe.
        center (list, optional): center point. Defaults to [75,75].

    Returns:
        pd.DataFrame: output dataframe.
    """
    # get distance from center
#     from scipy.spatial import distance
    if 'distance effective from centroid' in df1:
        df1=df1.drop(df1.filter(like='distance effective from centroid',axis=1).columns.tolist(),axis=1)
    df1['distance effective from centroid per frame']=df1.apply(lambda x: distance.euclidean(center,[x['x'],x['y']]),axis=1)
    df=df1.groupby('particle',as_index=False).agg({'distance effective from centroid per frame':[np.min,np.max]})
    df=flatten_columns(df)
    df=df.rename(columns={c:c if not 'per frame a' in c else c.replace('per frame a','') for c in df})
#     print(df.columns)
    df['distance effective from centroid']=df['distance effective from centroid max']-df['distance effective from centroid min']
    return df1.merge(df,on='particle',how='left')

## distances
def distance_effective(
    particle,
    frame1,
    frame2,
    t_cor: pd.DataFrame,
    ) -> float:
    """
    Effective distance between frames.

    Args:
        particle : particle
        frame1 (np.array): a frame.
        frame2 (np.array): another frame.
        t_cor (pd.DataFrame): t_cor.

    Returns:
        float: distance
    """
    from scipy.spatial import distance
    a=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame1)),['x','y']]
    b=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame2)),['x','y']]
    return distance.euclidean(a.values, b.values)

def get_distance_travelled(
    t_cor: pd.DataFrame,
    ) -> pd.DataFrame:
    """
    Distance travelled.

    Args:
        t_cor (pd.DataFrame): input dataframe.

    Returns:
        pd.DataFrame: output dataframe.
    """
    for f1,f2 in zip(list(range(0,t_cor['frame'].max())),
                list(range(1,t_cor['frame'].max()+1))):
        for p in t_cor['particle'].unique():
            a=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f1)),['x','y']]
            b=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),['x','y']]
            if len(a)!=0 and len(b)!=0:
                t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),'distance']=distance.euclidean(a.values, b.values)
    t_cor_distance=t_cor.groupby('particle').agg({'distance':sum})

    t_cor_rangeframes=t_cor.groupby('particle').agg({'frame':[min,max]})
    t_cor_rangeframes=flatten_columns(t_cor_rangeframes)
    t_cor_rangeframes['distance effective']=t_cor_rangeframes.apply(lambda x : distance_effective(x.name,x['frame min'],x['frame max'],t_cor) ,axis=1)

    t_cor=t_cor.merge(t_cor_distance,
                left_on='particle',right_index=True,suffixes=[' delta',' total'])
    t_cor=t_cor.merge(t_cor_rangeframes,
                left_on='particle',right_index=True)
    t_cor['distance delta']=t_cor['distance delta'].fillna(0)
    def get_distance_total_per_frame(df):
        df['distance total per frame']=df.apply(lambda x : df.set_index('frame').loc[[i for i in range(df['frame'].min(),int(x['frame'])+1)],'distance delta'].sum() ,axis=1)
        return df

    t_cor=t_cor.groupby('particle').apply(get_distance_total_per_frame)
    t_cor['distance effective per frame']=t_cor.apply(lambda x : distance_effective(particle=x['particle'],frame1=x['frame min'],frame2=x['frame'],t_cor=t_cor) ,axis=1)
    t_cor['intensity']=t_cor['mass']        
    return t_cor