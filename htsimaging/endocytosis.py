import sys
import argh                                       
from skimage import io,exposure,restoration,filters,morphology,measure
from glob import glob,iglob
from rohan.global_imports import *
from rohan.dandage.io_sys import runbashcmd
from htsimaging.lib.global_vars import *
from os.path import isdir
import pims
import trackpy as tp
import logging
import yaml
from htsimaging.lib.utils import filter_regions

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)

def segmentation2cells(imp,imsegp,fiterby_border_thickness=100,magnification=100,plotp=None):
    """
    prop_type='area',mn=100,mx=8000
    at 1.5X
    prop_type='area',mn=1500,mx=12000
    """
    im=io.imread(imp,as_gray=True)
    imseg=io.imread(imsegp,as_gray=True)
    regions=measure.label(imseg)
#     fiterby_border_thickness,im.shape[1]+fiterby_border_thickness
    regions=filter_regions(regions,im,prop_type='area',
                           mn=10*magnification,
                           mx=100000*magnification,
                           check=True,plotp=plotp)
    regions=filter_regions(regions,im,prop_type='eccentricity',mn=0,mx=0.8,check=False)
    regions=filter_regions(regions,im,prop_type='centroid_x',
                           mn=fiterby_border_thickness,mx=im.shape[0]+fiterby_border_thickness,check=False)
    regions=filter_regions(regions,im,prop_type='centroid_y',
                           mn=fiterby_border_thickness,mx=im.shape[1]+fiterby_border_thickness,check=False)
    return regions

from scipy.spatial import distance
def distance_effective(particle,frame1,frame2,t_cor):
    a=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame1)),['x','y']]
    b=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame2)),['x','y']]
    return distance.euclidean(a.values, b.values)

def get_distance_travelled(frames,t_cor,out_fh,test=False,force=False):
    ddistancesp=f"{out_fh}_distances.tsv"
    if not exists(ddistancesp) or force:    
        t_cor=read_table(f"{out_fh}.t2.tsv")
        for f1,f2 in zip(list(range(0,t_cor['frame'].max())),
                    list(range(1,t_cor['frame'].max()+1))):
            for p in t_cor['particle'].unique():
                a=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f1)),['x','y']]
                b=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),['x','y']]
                if len(a)!=0 and len(b)!=0:
                    t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),'distance']=distance.euclidean(a.values, b.values)
        if test:        
            print(t_cor.head())
        if 'distance' in t_cor:
            t_cor_distance=t_cor.groupby('particle').agg({'distance':sum})

            t_cor_rangeframes=t_cor.groupby('particle').agg({'frame':[min,max]})
            t_cor_rangeframes.columns=coltuples2str(t_cor_rangeframes.columns)
            t_cor_rangeframes['distance effective']=t_cor_rangeframes.apply(lambda x : distance_effective(x.name,x['frame min'],x['frame max'],t_cor) ,axis=1)

            t_cor=t_cor.merge(t_cor_distance,
                        left_on='particle',right_index=True,suffixes=[' delta',' total'])
            t_cor=t_cor.merge(t_cor_rangeframes,
                        left_on='particle',right_index=True)
            t_cor['distance delta']=t_cor['distance delta'].fillna(0)
            t_cor['distance total per frame']=t_cor.apply(lambda x : t_cor.set_index(['frame','particle']).loc[[(i,x['particle']) for i in range(int(x['frame'])+1)],'distance delta'].sum() ,axis=1)
            t_cor['distance effective per frame']=t_cor.apply(lambda x : distance_effective(particle=x['particle'],frame1=x['frame min'],frame2=x['frame'],t_cor=t_cor) ,axis=1)
            t_cor['intensity']=t_cor['mass']            
            to_table(t_cor,ddistancesp)
            ## plot dist distances
            plotp=f"{out_fh}_hist_distances.png"
            plt.figure()
            ax=plt.subplot()
            t_cor[['distance delta','distance total','distance effective']].dropna().hist(ax=ax)
            plt.tight_layout()
            savefig(plotp)   
            ## plot image labeled particles
            fig=plt.figure(figsize=[10,10])
            ax=plt.subplot()
            ax=plot_trajectories(t_cor, image=frames[0],label=True,colorby='frame',cmap='hsv',
                        ax=ax,
#params_text={'size':5},
                        )
#plot_trajectories(img=frames[-1],dtraj=t_cor,params_plot_traj={'label':True,'colorby':'frame','cmap':'hsv'})
            plotp=f"{out_fh}_trajectories.png"
            plt.tight_layout()
            savefig(plotp)    
        else:
            t_cor=pd.DataFrame(columns=t_cor.columns)
            to_table(t_cor,ddistancesp)
        return t_cor

from htsimaging.lib.spt import frames2coords_cor
from skimage.measure import label, regionprops_table, regionprops
def get_cellboxes(regions,plotp=None):
    import matplotlib.patches as mpatches
    if not plotp is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.imshow(regions,cmap='binary')
    cellbox_width=150
    cellboxes=[]
#     celli2props={}
    df1=pd.DataFrame(regionprops_table(regions.astype(int)))
    for regioni,region in enumerate(regionprops(regions.astype(int))):
        box_xmnxmxymnymx=[region.centroid[1]-(cellbox_width*0.5)
                          ,region.centroid[1]+(cellbox_width*0.5)
                          ,region.centroid[0]-(cellbox_width*0.5)
                          ,region.centroid[0]+(cellbox_width*0.5)
                         ]
        cellboxes.append([int(i) for i in box_xmnxmxymnymx])
#         celli2props[regioni+1]=region.area
        if not plotp is None:
            rect = mpatches.Rectangle([box_xmnxmxymnymx[0],box_xmnxmxymnymx[2]], cellbox_width, cellbox_width,
                                      fill=False, edgecolor='red', linewidth=2)
            ax.text(region.centroid[1],region.centroid[0],f"{regioni+1:d}",color='g')
#             ax.text(region.centroid[1],region.centroid[0],f"{region.extent:.2f}",color='g')
            ax.add_patch(rect)
    if not plotp is None:
        savefig(plotp)
#     df1=pd.DataFrame(pd.Series(celli2props))
#     df1.index.name='cell#'
#     df1.columns=['area']
    return cellboxes,df1#.reset_index()
    
def _plot(ax, coords, pos_columns, **plot_style):
    """ This function wraps Axes.plot to make its call signature the same for
    2D and 3D plotting. The y axis is inverted for 2D plots, but not for 3D
    plots.

    Parameters
    ----------
    ax : Axes object
        The axes object on which the plot will be called
    coords : DataFrame
        DataFrame of coordinates that will be plotted
    pos_columns : list of strings
        List of column names in x, y(, z) order.
    plot_style : keyword arguments
        Keyword arguments passed through to the `Axes.plot(...)` method

    Returns
    -------
    Axes object
    """
    if len(pos_columns) == 3:
        return ax.plot(coords[pos_columns[0]], coords[pos_columns[1]],
                       zs=coords[pos_columns[2]], **plot_style)
    elif len(pos_columns) == 2:
        return ax.plot(coords[pos_columns[0]], coords[pos_columns[1]],
                       **plot_style)
def plot_trajectories(traj,image, colorby='particle', mpp=None, label=False,
                        cmap=None, ax=None, t_column=None,
                        pos_columns=None, plot_style={},
                        params_text={'ha':'center','va':'center'}, **kwargs):
    """Plot traces of trajectories for each particle.
    Optionally image it on a frame from the video.

    Parameters
    ----------
    traj : DataFrame
        The DataFrame should include time and spatial coordinate columns.
    colorby : {'particle', 'frame'}, optional
    mpp : float, optional
        Microns per pixel. If omitted, the labels will have units of pixels.
    label : boolean, optional
        Set to True to write particle ID numbers next to trajectories.
    image : ndarray, optional
        Background image, default None
    cmap : colormap, optional
        This is only used in colorby='frame' mode. Default = mpl.cm.winter
    ax : matplotlib axes object, optional
        Defaults to current axes
    t_column : string, optional
        DataFrame column name for time coordinate. Default is 'frame'.
    pos_columns : list of strings, optional
        Dataframe column names for spatial coordinates. Default is ['x', 'y'].
    plot_style : dictionary
        Keyword arguments passed through to the `Axes.plot(...)` command

    Returns
    -------
    Axes object
    
    See Also
    --------
    plot_traj3d : the 3D equivalent of `plot_traj`
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from rohan.dandage.plot.colors import get_cmap_subset
#     if cmap is None:
#         cmap = get_cmap_subset('binary',vmin=0.15,vmax=0.05)
    if t_column is None:
        t_column = 'frame'
    if pos_columns is None:
        pos_columns = ['x', 'y']
    if len(traj) == 0:
        raise ValueError("DataFrame of trajectories is empty.")
    _plot_style = dict(linewidth=1)
    # Background image
    ax=plt.subplot() if ax is None else ax
    ax.imshow(image,
       cmap=get_cmap_subset('binary',vmin=0.2,vmax=0))
    # Trajectories
    if colorby == 'particle':
        # Unstack particles into columns.
        unstacked = traj.set_index(['particle', t_column])[pos_columns].unstack()
        for i, trajectory in unstacked.iterrows():
            _plot(ax, mpp*trajectory, pos_columns, **_plot_style)
    if colorby == 'frame':
        # Read http://www.scipy.org/Cookbook/Matplotlib/MulticoloredLine
        x = traj.set_index([t_column, 'particle'])['x'].unstack()
        y = traj.set_index([t_column, 'particle'])['y'].unstack()
        color_numbers = traj[t_column].values/float(traj[t_column].max())
#         logger.info("Drawing multicolor lines takes awhile. "
#                     "Come back in a minute.")
        for particle in x:
            points = np.array(
                [x[particle].values, y[particle].values]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap=cmap)
            lc.set_array(color_numbers)
            ax.add_collection(lc)
#             ax.set_xlim(x.apply(np.min).min(), x.apply(np.max).max())
#             ax.set_ylim(y.apply(np.min).min(), y.apply(np.max).max())
    if label:
        unstacked = traj.set_index([t_column, 'particle'])[pos_columns].unstack()
        first_frame = int(traj[t_column].min())
        coords = unstacked.fillna(method='backfill').stack().loc[first_frame]
        for particle_id, coord in coords.iterrows():
            ax.text(*coord.tolist(), s="%d" % particle_id,
                    **params_text)
    ax.set_xlim(0,image.shape[0])
    ax.set_ylim(0,image.shape[1])
    ax.set_ylim(ax.get_ylim()[::-1])
    return ax    

def make_gif(frames,t_cor,outd=None,test=False,force=False):
    from rohan.dandage.plot.colors import get_cmap_subset
    if not outd is None:
        makedirs(outd,exist_ok=True)
        gifp=f"{dirname(outd)}/vid.gif"
    else:
        test=True
        gifp=''
    if not exists(gifp) or force:
        for framei,frame in enumerate(frames):
            plotp=f'{outd}/{framei:02d}.png'
            plt.figure(figsize=[10,10])
            ax=plt.subplot(111)
            ax.imshow(frame,cmap=get_cmap_subset('binary_r',vmin=0.2,vmax=1),alpha=0.8,
                     )
            ax.text(frame.shape[0],frame.shape[1],f"frame={framei:05d}",ha='right',va='bottom',size='20',color='y')
            if not framei==0:
                # traj of particle
                for p in t_cor.loc[(t_cor['frame']==framei),'particle'].unique():
                    # traj of particle
                    framei_min=t_cor['frame min'].unique()[0]
                    df=t_cor.loc[((t_cor['particle']==p) \
                                  & (t_cor['frame'].isin(range(framei_min,framei+1)))\
                                  & (t_cor['x'].between(0,frame.shape[0]))\
                                  & (t_cor['y'].between(0,frame.shape[1]))),:].sort_values(by=['frame','particle'])
                    if test:
                        print(df['frame'])
                    if len(df)!=0:
                        if not 'move' in df:
                            df.plot.line(x='x',y='y',lw=7,
                                c='limegreen',
                                legend=False,ax=ax)                        
                        else:
                            df.plot.line(x='x',y='y',lw=7,
                                c='limegreen' if df.loc[:,['particle','move']].drop_duplicates().set_index('particle').loc[p,'move']==1 else 'magenta',
                                legend=False,ax=ax)
                ax.set_xlim(0,frame.shape[1])
                ax.set_ylim(0,frame.shape[1])
                plt.axis('off')
                if test:
                    if framei==3:
                        return ax
                if not test:    
                    makedirs(dirname(plotp),exist_ok=True)
                    plt.savefig(plotp)
                plt.clf();plt.close();
            _framei=framei
        if not test:    
            plt.close('all')
            com=f"convert -delay 10 -loop 0 {outd}/*.png {gifp}"
            runbashcmd(com)
    return gifp

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

from multiprocessing import Pool
def multiprocess_cellframes2distances(cellcfgp):
    celli=basename(dirname(cellcfgp))
    print(celli);logging.info(celli)
    cellcfg=yaml.load(open(cellcfgp,'r'))
    cellframes2distances([np.load(p) for p in cellcfg['cellframeps']],
                         [np.load(p) for p in cellcfg['cellframesmaskedps']],
                         out_fh=f"{cellcfg['outp']}/plot_check",
                         test=cellcfg['test'],force=cellcfg['force'])
    

def run_trials(prjd,bright_fn_marker,test=False,force=False,cores=4):
    """
    runs the analysis.   
    
    :param prjd: path to (project) directory containing individual runs of images
        eg. if the images are stored in such a way
        images_190919
            WT-WT_001
                tif ..
            WT-WT_002
                tif ..
            WT-WT_003
                tif ..
            WT-WT_004
                tif ..
        images_190919 will be prjd
        so the correct command will be
        python endocytosis.py run-trials /path/to/images_190919 _T1C1

    :param bright_fn_marker: _t if inhouse microscope else if chul: _T1C1    
    
    """    
    prjd=abspath(prjd)
    cfgp=f"{prjd}/cfg.yml"
    if not exists(cfgp) or force:
        print('making cfg')
        cfg={'prjd':prjd}
        cfg['cfgp']=cfgp
        cfg['cores']=cores
        cfg['bright_fn_marker']=bright_fn_marker
        if 'T1C1' in cfg['bright_fn_marker']:
            cfg['magnification']=150
        elif cfg['bright_fn_marker']=='_t':
            cfg['magnification']=100
        else:
            logging.error("unknown bright_fn_marker {cfg['bright_fn_marker']}")
            sys.exit()
        cfg['trials']={basename(d):{'datad':d} for d in glob(f"{cfg['prjd']}/*") if (isdir(d) and basename(d).replace('/','')!='segmentation_cell' and not basename(d).startswith('_'))}
        trials_bad=[]
        for k in cfg['trials']:
            cfg['trials'][k]['bright']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if cfg['bright_fn_marker'] in p and not p.endswith('.segmented.tif')]
            cfg['trials'][k]['gfp']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if not (cfg['bright_fn_marker'] in p or 'segmented' in p) ]
            cfg['trials'][k]['plotd']=f"{cfg['trials'][k]['datad']}/plot"
            makedirs(cfg['trials'][k]['plotd'],exist_ok=True)
            if len(cfg['trials'][k]['bright'])==0 or len(cfg['trials'][k]['gfp'])==0:
                trials_bad.append(k)
        for trial in trials_bad:
            del cfg['trials'][trial]
        yaml.dump(cfg,open(cfgp,'w'))   
        # QC
        ## 1 CHECK BLEACHING
        from rohan.dandage.plot.line import plot_mean_std                                                        
        for k in cfg['trials']:
            frames = pims.ImageSequence(np.sort(cfg['trials'][k]['gfp']), as_grey=True)
            dframes=pd.DataFrame({i:np.ravel(np.array(frame)) for i,frame in enumerate(frames)})
            plotp=f"{cfg['trials'][k]['plotd']}/plot_check_bleaching.png"
            makedirs(dirname(plotp),exist_ok=True)

            df=dframes.describe().T
            ax=plot_mean_std(df,cols=['mean','min','50%'])
            ax.set_xlabel('time points');ax.set_ylabel('intensity')
            savefig(plotp)
#         yaml.dump(cfg,open(cfgp,'w'))   
    else:
        cfg=yaml.load(open(cfgp,'r'))
    ## get segments from brightfield images
    if not 'flag_segmentation_done' in cfg or force:
        from htsimaging.lib.segment import run_yeastspotter
        cfg['yeastspotter_srcd']=f"{dirname(realpath(__file__))}/../deps/yeast_segmentation"
        logging.info(cfg.keys())
        cfg=run_yeastspotter(cfg,test=test)
        yaml.dump(cfg,open(cfgp,'w'))
        cfg['flag_segmentation_done']=True
        print('flag_segmentation_done')

#     if not '' in cfg:
    ## get and filter cells from segments images
    if not 'flag_cells_done' in cfg or force:
        for trial in cfg['trials']:
            if len(cfg['trials'][trial]['bright'])!=0:
                cellsps=[]
                for imp,imsegp in zip(cfg['trials'][trial]['bright'],cfg['trials'][trial]['bright_segmented']):
                    cellsp=f'{imsegp}.npy'
                    regions=segmentation2cells(imp,imsegp,magnification=cfg['magnification'],
                       plotp=f"{cfg['trials'][trial]['plotd']}/image_segmentation2cells.png")
                    np.save(cellsp, regions)
                    cellsps.append(cellsp)
                cfg['trials'][trial]['bright_segmented_cells']=cellsps                                        
        cfg['flag_cells_done']=True
        yaml.dump(cfg,open(cfgp,'w'))
        print('flag_cells_done')

    if not 'flag_cellframes_done' in cfg or force:    
        cellcfgps=[]
        for trial in cfg['trials']:
            frames = pims.ImageSequence(np.sort(cfg['trials'][trial]['gfp']), as_grey=True)
            cellsp=np.sort(cfg['trials'][trial]['bright_segmented_cells'])[0] # only bright field at the start
            cells=np.load(cellsp)
            cellboxes,dcellprops=get_cellboxes(cells,plotp=f"{cfg['trials'][trial]['plotd']}/image_get_cellboxes.png")
            to_table(dcellprops,f"{cellsp}.cellprops.tsv")
            for celli,cellbox in enumerate(cellboxes):
                logging.info(f"{trial};cell{celli+1:08d}")
                cellcfg={}
                cellcfg['outp']=f"{cfg['trials'][trial]['datad']}/cells/cell{celli+1:08d}/"
                cellcfg['cfgp']=f"{cellcfg['outp']}/cfg.yml"
                if not exists(cellcfg['cfgp']) or force:
                    cellcfg['test']=test
                    cellcfg['force']=force
                    cellcfg['cellbrightp']=f"{cellcfg['outp']}/cellbright.npy"
                    cellbright=cells[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]

                    if not exists(dirname(cellcfg['cellbrightp'])): 
                        makedirs(dirname(cellcfg['cellbrightp']),exist_ok=True)
                    np.save(cellcfg['cellbrightp'], cellbright) 
                    # only one cell per box
                    cellcfg['cellbrightmaskp']=f"{cellcfg['outp']}/cellbrightmask.npy"
                    cellbrightmask=filter_regions(cellbright.astype(int),prop_type='centroid_x',mn=70,mx=80)==0
                    np.save(cellcfg['cellbrightmaskp'], cellbrightmask)

                    cellframeps=[]
                    cellframesmaskedps=[]
                    for framei,frame in enumerate(frames):
                        cellframe=frame[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
                        cellframep=f"{cellcfg['outp']}/cellframe/frame{framei:08d}.npy"
                        if not exists(dirname(cellframep)): 
                            makedirs(dirname(cellframep),exist_ok=True)
                        np.save(cellframep, cellframe)
                        cellframeps.append(cellframep)

                        cellframemasked=cellframe.copy()
                        cellframemasked[cellbrightmask]=0
                        cellframemaskedp=f"{cellcfg['outp']}/cellframesmasked/frame{framei:08d}.npy"
                        if not exists(dirname(cellframemaskedp)): 
                            makedirs(dirname(cellframemaskedp),exist_ok=True)
                        np.save(cellframemaskedp, cellframemasked)
                        cellframesmaskedps.append(cellframemaskedp)

                    cellcfg['cellframeps']=cellframeps
                    cellcfg['cellframesmaskedps']=cellframesmaskedps
                    yaml.dump(cellcfg,open(cellcfg['cfgp'],'w'))
    #                 multiprocess_cellframes2distances(cellcfgp)
                cellcfgps.append(cellcfg['cfgp'])
        cfg['cellcfgps']=cellcfgps        
        cfg['flag_cellframes_done']=True
        yaml.dump(cfg,open(cfgp,'w'))
        print('flag_cellframes_done')
        
        # parallel processing
    if not 'flag_distances_done' in cfg or force:    
        cellcfgps=np.sort(cfg['cellcfgps'])
        if len(cellcfgps)!=0:
            print(f"{get_datetime()}: processing: {len(cellcfgps)} cells.")
            if not test:
                pool=Pool(processes=cfg['cores']) 
                pool.map(multiprocess_cellframes2distances, cellcfgps)
                pool.close(); pool.join()         
            else:
                for cellcfgp in cellcfgps:
                    logging.info(f'processing {cellcfgp}')
                    multiprocess_cellframes2distances(cellcfgp)
        cfg['flag_distances_done']=True
        yaml.dump(cfg,open(cfgp,'w'))
        print('flag_distances_done')
    print('finished')

## begin    
import sys
is_interactive_notebook=any([basename(abspath('.')).startswith(f'{i:02d}_') for i in range(10)])
# from rohan.dandage.io_sys import is_interactive_notebook
if not is_interactive_notebook:
    # assembling:
    parser = argh.ArghParser()
    parser.add_commands([run_trials])
    from rohan.dandage.io_strs import get_logger,get_datetime
    level=logging.ERROR
    logp=get_logger(program='htsimaging',
               argv=[get_datetime()],
               level=level,
               dp=None)        
    logging.info(f"start. log file: {logp}")
    print(f"start. log file: {logp}")    
    if __name__ == '__main__':
        logging.info('start')
        parser.dispatch()
else:
    print("can not run from notebook")
