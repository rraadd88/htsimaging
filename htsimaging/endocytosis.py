import argh
from skimage import io 
from skimage import io,exposure,restoration,filters,morphology,measure
from glob import glob,iglob
from rohan.global_imports import *
from os.path import isdir
import pims
import trackpy as tp
import logging
import yaml
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

from htsimaging.lib.utils import filter_regions
def segmentation2cells(imp,imsegp,fiterby_border_thickness=100,plotp=None):
    im=io.imread(imp,as_gray=True)
    imseg=io.imread(imsegp,as_gray=True)
    regions=measure.label(imseg)
    fiterby_border_thickness,im.shape[1]+fiterby_border_thickness
    regions=filter_regions(im,regions,prop_type='area',mn=1000,mx=8000,check=True,plotp=plotp)
    regions=filter_regions(im,regions,prop_type='eccentricity',mn=0,mx=0.8,check=False)
    regions=filter_regions(im,regions,prop_type='centroid_x',
                           mn=fiterby_border_thickness,mx=im.shape[0]+fiterby_border_thickness,check=False)
    regions=filter_regions(im,regions,prop_type='centroid_y',
                           mn=fiterby_border_thickness,mx=im.shape[1]+fiterby_border_thickness,check=False)
    return regions

def get_distance_travelled(frames,t_cor,out_fh):
    t_cor=f_batch=read_table(f"{out_fh}.t2.tsv")
    from scipy.spatial import distance
    for f1,f2 in zip(list(range(0,t_cor['frame'].max())),
                list(range(1,t_cor['frame'].max()+1))):
        for p in t_cor['particle'].unique():
            a=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f1)),['x','y']]
            b=t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),['x','y']]
            if len(a)!=0 and len(b)!=0:
                t_cor.loc[((t_cor['particle']==p) & (t_cor['frame']==f2)),'distance']=distance.euclidean(a.values, b.values)
    print(t_cor.head())
    t_cor_distance=t_cor.groupby('particle').agg({'distance':sum})

    t_cor_rangeframes=t_cor.groupby('particle').agg({'frame':[min,max]})
    t_cor_rangeframes.columns=coltuples2str(t_cor_rangeframes.columns)
    def distance_effective(particle,frame1,frame2,t_cor):
        a=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame1)),['x','y']]
        b=t_cor.loc[((t_cor['particle']==particle) & (t_cor['frame']==frame2)),['x','y']]
        return distance.euclidean(a.values, b.values)
    t_cor_rangeframes['distance effective']=t_cor_rangeframes.apply(lambda x : distance_effective(x.name,x['frame min'],x['frame max'],t_cor) ,axis=1)

    t_cor=t_cor.merge(t_cor_distance,
                left_on='particle',right_index=True,suffixes=[' delta',' total'])
    t_cor=t_cor.merge(t_cor_rangeframes,
                left_on='particle',right_index=True)

    ax=plt.subplot()
    t_cor[['distance delta','distance total','distance effective']].hist(ax=ax)

    t_cor['move']=t_cor['distance effective'].apply(lambda x : 1 if x>t_cor['distance effective'].quantile(0.6) else 0 )
    # print([t_cor_distance['distance'].quantile(q) for q in np.arange(0,1,0.1)])

    plotp=f"{out_fh}_trajectories.png"
    # t_cor_distance['move']=t_cor_distance['distance'].apply(lambda x : 1 if x>t_cor_distance['distance'].quantile(0.8) else 0 )
    plt.figure(figsize=[20,20])
    ax=plt.subplot(111)
    ax.imshow(frames[-1],cmap='binary_r',alpha=0.8)
    ax = tp.plot_traj(t_cor,label=False,ax=ax,
                      colorby='move',
                     )
    for p in t_cor['particle'].unique():
        t_cor.loc[(t_cor['particle']==p),:].plot.line(x='x',y='y',lw=3,
                          c='limegreen' if t_cor.loc[:,['particle','move']].drop_duplicates().set_index('particle').loc[p,'move']==1 else 'magenta',
                          legend=False,ax=ax)
    plt.savefig(plotp)    

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

def cellframes2distances(cellframes,out_fh=None,test=False,force=False):
    makedir(dirname(out_fh),exist_ok=True)
    params_msd={'mpp':0.0645,'fps':0.2, 'max_lagtime':100}
    # for 170x170 images
    params_locate_start={'diameter':5,'minmass_percentile':90} # for larger images increase the diameter
    params_link_df={'search_range':5,'memory':0,'link_strategy':'drop',}
    params_filter={'mass_cutoff':0.6,'size_cutoff':1,'ecc_cutoff':1,
                  'filter_stubs':False,'flt_mass_size':False,'flt_incomplete_trjs':False,
                  'test':test}
    makedirs(dirname(out_fh),exist_ok=True)
    t_cor=frames2coords_cor(frames=cellframes,
                            params_locate_start=params_locate_start,
                            params_msd=params_msd,params_link_df=params_link_df,
                            params_filter=params_filter,
                            out_fh=out_fh,force=force)
    get_distance_travelled(frames=cellframes,t_cor=t_cor,out_fh=out_fh)
        
def run_trials(prjd,test=False,force=False):
    cfgp=f"{prjd}/cfg.yml"
    if not exists(cfgp):
        cfg={'prjd':prjd}
        cfg['cfgp']=cfgp
        cfg['trials']={basename(d):{'datad':d} for d in glob(f"{cfg['prjd']}/*") if (isdir(d) and basename(d).replace('/','')!='segmentation_cell')}
        for k in cfg['trials']:
            cfg['trials'][k]['gfp']=[p for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if '_t' in p]
            cfg['trials'][k]['bright']=[p for p in glob(f"{cfg['trials'][k]['datad']}/*tif") if not '_t' in p]
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
        from htsimaging.lib.segment import run_yeastspotter
        cfg['yeastspotter_srcd']=f"{dirname(realpath(__file__))}/../deps/yeast_segmentation"
        print(cfg.keys())
        cfg=run_yeastspotter(cfg,test=True)
        yaml.dump(cfg,open(cfgp,'w'))
        cfg['flag_segmentation_done']=True

#     if not '' in cfg:
    ## get and filter cells from segments images
    if not 'flag_cells_done' in cfg:
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
        for trial in cfg['trials']:
            frames = pims.ImageSequence(np.sort(cfg['trials'][trial]['gfp']), as_grey=True)
            for cellsp in cfg['trials'][trial]['bright_segmented_cells']:
                cells=np.load(cellsp)
                cellboxes=get_cellboxes(cells,test=test)
                for celli,cellbox in enumerate(cellboxes):
                    cellframes=[f[cellbox[0]:cellbox[1],cellbox[2]:cellbox[3]] for f in frames]
                    cellframes2distances(cellframes,out_fh=f"{cfg['trials'][trial]['plotd']}/cell{celli:08d}/plot_check",
                                         test=test,force=force)
                    break
                break                                               
            break

# assembling:
parser = argh.ArghParser()
parser.add_commands([run_trials])

if __name__ == '__main__':
    logging.info('start')
    parser.dispatch()
    logging.info('done')
                                                        