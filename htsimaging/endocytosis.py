from rohan.global_imports import *
from rohan.dandage.io_sys import runbashcmd
import sys
import argh                                       
from glob import glob,iglob
from os.path import isdir

from skimage import io,exposure,restoration,filters,morphology,measure
from skimage.external import tifffile
import pims
import trackpy as tp

from htsimaging.lib.global_vars import *
from htsimaging.lib.utils import filter_regions,get_cellprops
from htsimaging.lib.plot import *
from htsimaging.lib.stat import *

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

from htsimaging.lib.spt import frames2coords_cor
def get_cellboxes(regions,plotp=None):
    import matplotlib.patches as mpatches
    if not plotp is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.imshow(regions,cmap='binary')
    cellbox_width=150
    cellboxes=[]
    for regioni,region in enumerate(measure.regionprops(regions.astype(int))):
        box_xmnxmxymnymx=[region.centroid[1]-(cellbox_width*0.5),
                          region.centroid[1]+(cellbox_width*0.5),
                          region.centroid[0]-(cellbox_width*0.5),
                          region.centroid[0]+(cellbox_width*0.5),
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
    return cellboxes
    

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
    celli=dirname(cellcfgp)
    print(celli);logging.info(celli)
    cellcfg=yaml.load(open(cellcfgp,'r'))
    cellframes2distances([np.load(p) for p in sorted(cellcfg['cellframeps'])],
                         [np.load(p) for p in sorted(cellcfg['cellframesmaskedps'])],
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
            dcellprops=get_cellprops(cells,
                                     intensity_imgtype2img={
            "gfp mean":np.mean(frames,axis=0),
            "bright":tifffile.imread(cfg['trials'][trial]['bright'])})
            to_table(dcellprops,f"{cellsp}.cellprops.tsv")
            
            cellboxes=get_cellboxes(cells,plotp=f"{cfg['trials'][trial]['plotd']}/image_get_cellboxes.png")
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
                for cellcfgp in cellcfgps:
                    cellcfg_=read_dict(cellcfgp)
                    cellcfg_['force']=force
                    to_dict(cellcfg_,cellcfgp)
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
# is_interactive_notebook=any([basename(abspath('.')).startswith(f'{i:02d}_') for i in range(10)])
# # from rohan.dandage.io_sys import is_interactive_notebook
# if not is_interactive_notebook:
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
# else:
#     print("can not run from notebook")
