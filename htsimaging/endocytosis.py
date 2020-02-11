from rohan.global_imports import *
import sys
import argh

from skimage.external import tifffile
import pims
from htsimaging.lib.global_vars import *
from htsimaging.lib.spt import apply_cellframes2distances

import warnings
warnings.filterwarnings("ignore")
warnings.simplefilter(action='ignore', category=FutureWarning)
from multiprocessing import Pool    

def run_trials(prjd,bright_fn_marker,test=False,force=False,cores=4):
    """
    runs the analysis.   
    
    Nomenclature:
    `project` is a folder containing many timelapses,
    timelapses are called `trials`,
    each `trial`, after analysis, would contain `cells`.
    
    :param prjd: path to (`project`) folder containing `trials`
        eg. if the images are stored in such a way
        
        images_190919 (`project`)
            WT-WT_001 (`trial`)
                tifs ..
            WT-WT_002 (`trial`)
                tifs ..
            WT-WT_003 (`trial`)
                tif ..
            WT-WT_004 (`trial`)
                tifs ..
        
        Here, `images_190919` will be prjd
        so the correct command will be
        python endocytosis.py run-trials /path/to/images_190919 _T1C1

    :param bright_fn_marker: _t if inhouse microscope else if chul: _T1C1    
    
    """    
    # make cfg for the project
    from htsimaging.lib.io_cfg import make_project_cfg
    make_project_cfg(prjd,bright_fn_marker,test,force,cores)
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
        from htsimaging.lib.segment import segmentation2cells
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
        from htsimaging.lib.segment import get_cellboxes
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
                # make cg for cell
                cellcfg=make_cell_cfg(cfg,cells,cellbox)
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
                pool.map(apply_cellframes2distances, cellcfgps)
                pool.close(); pool.join()         
            else:
                for cellcfgp in cellcfgps:
                    logging.info(f'processing {cellcfgp}')
                    apply_cellframes2distances(cellcfgp)
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
