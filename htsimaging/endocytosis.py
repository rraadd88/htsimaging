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

def segmentation2cells(imp,imsegp,fiterby_border_thickness=100,plotp=None):
    im=io.imread(imp,as_gray=True)
    imseg=io.imread(imsegp,as_gray=True)
    regions=measure.label(imseg)
#     fiterby_border_thickness,im.shape[1]+fiterby_border_thickness
    regions=filter_regions(regions,im,prop_type='area',mn=1000,mx=8000,check=True,plotp=plotp)
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
            to_table(t_cor,ddistancesp)
            if test:
                plotp=f"{out_fh}_hist_distances.png"
                plt.figure()
                ax=plt.subplot()
                t_cor[['distance delta','distance total','distance effective']].dropna().hist(ax=ax)
                plt.tight_layout()
                plt.savefig(plotp)   
            if test:
                plot_trajectories(img=frames[-1],dtraj=t_cor,params_plot_traj={'label':False})
                plotp=f"{out_fh}_trajectories.png"    
                plt.savefig(plotp)    
        else:
            to_table(pd.DataFrame(columns=t_cor.columns),ddistancesp)

from htsimaging.lib.spt import frames2coords_cor
from skimage.measure import label, regionprops
def get_cellboxes(regions,test=False):
    import matplotlib.patches as mpatches
    if test:
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.imshow(regions)
    cellboxes=[]
    for region in regionprops(regions.astype(int)):
        box_xmnxmxymnymx=[region.centroid[1]-50,region.centroid[1]+50,region.centroid[0]-50,region.centroid[0]+50]
        cellboxes.append([int(i) for i in box_xmnxmxymnymx])
        if test:
            rect = mpatches.Rectangle([box_xmnxmxymnymx[0],box_xmnxmxymnymx[2]], 100, 100,
                                      fill=False, edgecolor='red', linewidth=2)
            ax.text(region.centroid[1],region.centroid[0],f"{region.extent:.2f}",color='g')
            ax.add_patch(rect)
    return cellboxes

def plot_trajectories(img,dtraj,params_plot_traj={'label':False}):
    plt.figure()
    ax=plt.subplot(111)
    ax.imshow(img,cmap='binary_r',alpha=0.8,zorder=-1)
    ax = tp.plot_traj(dtraj,ax=ax,**params_plot_traj)
    ax.set_xlim(0,img.shape[0])
    ax.set_ylim(0,img.shape[1])
    plt.tight_layout()

def make_gif(frames,t_cor,outd=None,test=False,force=False):
    if outd is None:
        test=True            
    makedirs(outd,exist_ok=True)
    gifp=f"{dirname(outd)}/vid.gif"
    if not exists(gifp) or force:
        for framei,frame in enumerate(frames):
            plotp=f'{outd}/{framei:02d}.png'
            plt.figure(figsize=[5,5])
            ax=plt.subplot(111)
            ax.imshow(frame,cmap='binary_r',alpha=0.8)
            ax.text(frame.shape[0],frame.shape[1],f"frame={framei:02d}",ha='right',va='top',size='10',color='y')
            for p in t_cor['particle'].unique():
                df=t_cor.loc[((t_cor['particle']==p) \
                              & (t_cor['x'].between(0,frame.shape[0]))\
                              & (t_cor['y'].between(0,frame.shape[1]))),:]
                if len(df)!=0:
                    df.plot.line(x='x',y='y',lw=1,
                                      c='limegreen' if t_cor.loc[:,['particle','move']].drop_duplicates().set_index('particle').loc[p,'move']==1 else 'magenta',
                                      legend=False,ax=ax)
            ax.set_xlim(0,frame.shape[1])
            ax.set_ylim(0,frame.shape[1])
            plt.axis('off')
            if test:
                return ax
            makedirs(dirname(plotp),exist_ok=True)
            plt.savefig(plotp)
        plt.close('all')
        com=f"convert -delay 10 -loop 0 {outd}/*.png {gifp}"
        runbashcmd(com)
    return gifp

def cellframes2distances(cellframes,cellframesmasked,out_fh=None,test=False,force=False):
    makedirs(dirname(out_fh),exist_ok=True)
    params_msd={'mpp':0.0645,'fps':0.2, 'max_lagtime':100}
    # for 170x170 images
    params_locate_start={'diameter':7,'minmass_percentile':90} # for larger images increase the diameter
    params_link_df={'search_range':5,'memory':0,'link_strategy':'drop',}
    params_filter={'mass_cutoff':0.5,'size_cutoff':1,'ecc_cutoff':1,
                  'filter_stubs':False,'flt_mass_size':False,'flt_incomplete_trjs':False,
                  'test':test}
    makedirs(dirname(out_fh),exist_ok=True)
    t_cor=frames2coords_cor(frames=cellframesmasked,out_fh=out_fh,
                            params_locate_start=params_locate_start,
                            params_msd=params_msd,params_link_df=params_link_df,
                            params_filter=params_filter,
                            subtract_drift=False,
                            force=force)
    get_distance_travelled(frames=cellframesmasked,t_cor=t_cor,out_fh=out_fh,test=test,force=force)
    if not out_fh is None:
        make_gif(cellframes,t_cor,f"{dirname(out_fh)}/vid",force=force)
    
def run_trials(prjd,test=False,force=False):
    prjd=abspath(prjd)
    cfgp=f"{prjd}/cfg.yml"
    if not exists(cfgp) or force:
        print('making cfg')
        cfg={'prjd':prjd}
        cfg['cfgp']=cfgp
        cfg['trials']={basename(d):{'datad':d} for d in glob(f"{cfg['prjd']}/*") if (isdir(d) and basename(d).replace('/','')!='segmentation_cell' and not basename(d).startswith('_'))}
        for k in cfg['trials']:
            cfg['trials'][k]['gfp']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if '_t' in p]
            cfg['trials'][k]['bright']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if not '_t' in p]
            cfg['trials'][k]['plotd']=f"{cfg['trials'][k]['datad']}/plot"
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
            plt.savefig(plotp)
        yaml.dump(cfg,open(cfgp,'w'))   
    else:
        cfg=yaml.load(open(cfgp,'r'))
    ## get segments from brightfield images
    if not 'flag_segmentation_done' in cfg:
        print('flag_segmentation_done')
        from htsimaging.lib.segment import run_yeastspotter
        cfg['yeastspotter_srcd']=f"{dirname(realpath(__file__))}/../deps/yeast_segmentation"
        print(cfg.keys())
        cfg=run_yeastspotter(cfg,test=True)
        yaml.dump(cfg,open(cfgp,'w'))
        cfg['flag_segmentation_done']=True

#     if not '' in cfg:
    ## get and filter cells from segments images
    if not 'flag_cells_done' in cfg:
        print('flag_cells_done')
        for trial in cfg['trials']:
            if len(cfg['trials'][trial]['bright'])!=0:
                cellsps=[]
                for imp,imsegp in zip(cfg['trials'][trial]['bright'],cfg['trials'][trial]['bright_segmented']):
                    cellsp=f'{imsegp}.npy'
                    regions=segmentation2cells(imp,imsegp,
                       plotp=f"{cfg['trials'][trial]['plotd']}/image_segmentation2cells.png")
                    np.save(cellsp, regions)
                    cellsps.append(cellsp)
                cfg['trials'][trial]['bright_segmented_cells']=cellsps                                        
        cfg['flag_cells_done']=True
        yaml.dump(cfg,open(cfgp,'w'))

    if not 'flag_distances_done' in cfg:    
        print('flag_distances_done')
        for trial in cfg['trials']:
            frames = pims.ImageSequence(np.sort(cfg['trials'][trial]['gfp']), as_grey=True)
            for cellsp in cfg['trials'][trial]['bright_segmented_cells']:
                cells=np.load(cellsp)
                cellboxes=get_cellboxes(cells,test=test)
                for celli,cellbox in enumerate(cellboxes):
                    print(f"{trial};cell{celli+1:08d}")
                    logging.info(f"{trial};cell{celli+1:08d}")
                    outp=f"{cfg['trials'][trial]['datad']}/cells/cell{celli+1:08d}/"

                    cellbrightp=f"{outp}/cellbright.npy"
                    cellbright=cells[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
                    if not exists(dirname(cellbrightp)): 
                        makedirs(dirname(cellbrightp),exist_ok=True)
                    np.save(cellbrightp, cellbright) 
                    # only one cell per box
                    cellbrightmaskp=f"{outp}/cellbrightmask.npy"
                    cellbrightmask=filter_regions(cellbright.astype(int),prop_type='centroid_x',mn=45,mx=55)==0
                    np.save(cellbrightmaskp, cellbrightmask)
                    
                    cellframes=[]
                    cellframesmasked=[]
                    for fi,f in enumerate(frames):
                        cellframe=f[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
                        cellframes.append(f)
                        cellframep=f"{outp}/cellframe/frame{fi:08d}.npy"
                        if not exists(dirname(cellframep)): 
                            makedirs(dirname(cellframep),exist_ok=True)
                        np.save(cellframep, cellframe)
                        cellframe[cellbrightmask]=0
                        cellframesmasked.append(f)
                        cellframesmaskedp=f"{outp}/cellframesmasked/frame{fi:08d}.npy"
                        if not exists(dirname(cellframesmaskedp)): 
                            makedirs(dirname(cellframesmaskedp),exist_ok=True)
                        np.save(cellframesmaskedp, cellframesmasked)
                        
                    cellframes2distances(cellframes,cellframesmasked,
                                         out_fh=f"{outp}/plot_check",
                                         test=test,force=force)
import sys
exfromnotebook=any([basename(abspath('.')).startswith(f'{i:02d}_') for i in range(10)])
if not exfromnotebook:
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
        logging.info('done')
