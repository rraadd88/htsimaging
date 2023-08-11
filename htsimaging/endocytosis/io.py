#!/usr/bin/env python
"""Processing of the paths of input images to create configurations and metadata."""

import logging
from glob import glob
from os import makedirs
from os.path import isdir,abspath,exists,basename,dirname
import numpy as np
import matplotlib.pyplot as plt # used for QC plots
from roux.lib.io import read_dict, to_dict

def make_project_cfg(
    prjd: str,
    output_dir_path: str,
    bright_fn_marker: str=None,
    segmented_fn_marker:str=None,
    magnification: int=None,
    image_ext: str='tif',
    cores: int=1,
    test: bool=False,
    force: bool=False,
    )-> dict:
    """Make the confguration for the analysis run.

    Args:
        prjd (str): path to the directory with the images.
        output_dir_path (str): output directory path.
        bright_fn_marker (_type_): marker in the path of the bright field images.
        segmented_fn_marker (_type_): marker in the path of the segmented images.
        cores (int, optional): number of cores. Defaults to 1.
        test (bool, optional): test-mode. Defaults to False.
        force (bool, optional): over-write theoutputs. Defaults to False.

    Returns:
        dict: metadata
        
    Notes:
        Infer the magnification from the filenames:
        if 'T1C1' in cfg['bright_fn_marker']:
            cfg['magnification']=150
        elif cfg['bright_fn_marker']=='_t':
            cfg['magnification']=100
    """
    prjd=abspath(prjd)
    cfgp=f"{prjd}/cfg.yml"
    if not exists(cfgp) or force:    
        logging.info('making cfg')
        cfg={'prjd':prjd}
        cfg['outd']=output_dir_path
        cfg['cfgp']=cfgp
        cfg['cores']=cores
        cfg['bright_fn_marker']=bright_fn_marker
        cfg['segmented_fn_marker']=segmented_fn_marker
        assert not magnification is None
        cfg['magnification']=magnification
        # else:
        #     logging.error("unknown bright_fn_marker {cfg['bright_fn_marker']}")
        cfg['trials']={basename(d):{'datad':d} for d in glob(f"{cfg['prjd']}/*") if (isdir(d) and basename(d).replace('/','')!='segmentation_cell' and not basename(d).startswith('_'))}
        trials_bad=[]
        for k in cfg['trials']:
            if not cfg['bright_fn_marker'] is None:
                cfg['trials'][k]['bright']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*{image_ext}") if cfg['bright_fn_marker'] in p and not p.endswith(f'.segmented.{image_ext}')]
            else:
                cfg['trials'][k]['bright']=[]
            if not cfg['segmented_fn_marker'] is None:            
                cfg['trials'][k]['segmented']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*{image_ext}") if cfg['segmented_fn_marker'] in p]
            else:
                cfg['trials'][k]['segmented']=[]
            cfg['trials'][k]['gfp']=[abspath(p) for p in glob(f"{cfg['trials'][k]['datad']}/*{image_ext}")]
            cfg['trials'][k]['gfp']=list(set(cfg['trials'][k]['gfp'])-set(cfg['trials'][k]['bright']+cfg['trials'][k]['segmented']))
            cfg['trials'][k]['outd']=f"{cfg['outd']}/{k}/"
            cfg['trials'][k]['plotd']=f"{cfg['trials'][k]['outd']}/plot"
            makedirs(cfg['trials'][k]['plotd'],exist_ok=True)
            if len(cfg['trials'][k]['gfp'])==0:
                logging.warning(f"no gfp images for trial:{k}")
                trials_bad.append(k)
                continue
            if len(cfg['trials'][k]['bright'])+len(cfg['trials'][k]['segmented'])==0:
                logging.warning(f"no bright-field/segmented images for trial:{k}")
                trials_bad.append(k)
                continue
        for trial in trials_bad:
            logging.warning(f"removing trial:{trial}")
            del cfg['trials'][trial]
        to_dict(cfg,cfgp)   
        # QC
        ## 1 CHECK BLEACHING
        for k in cfg['trials']:
            from htsimaging.viz.stat import plot_summary_stats
            ax=plot_summary_stats(
                input_paths=sorted(cfg['trials'][k]['gfp'])
            )
            ax.set(xlabel='time points',ylabel='intensity',title="QC: bleaching test")
            plotp=f"{cfg['trials'][k]['plotd']}/plot_check_bleaching.png"
            makedirs(dirname(plotp),exist_ok=True)
            plt.savefig(plotp)
    else:
        cfg=read_dict(cfgp)
    return cfg

def make_cell_cfg(
    cfg: dict,
    frames: list,
    cells: list,
    trial: str,
    celli: int,
    cellbox: list,
    params_get_signal_summary_by_roi: dict={
        'xy_center':None,
        'width':20,
        'fun_summary_frame':'min',
        'fun_summary_frames':'median'
        },
    filterby_centroid: bool=False,
    scale_signal_cytoplasm: float=1.5,
    test: bool=False,
    force: bool=False,
    ) -> dict:
    """Make the configuration for an individual cell.

    Args:
        cfg (dict): metadata.
        frames (list): list of frames.
        cells (list): list of cells.
        trial (str): trial name.
        celli (int): index of the cell.
        cellbox (list): bounding box of the cell
        params_get_signal_summary_by_roi (dict, optional): parameters for the aggregation of the values at the ROI. Defaults to {'xy_center':None,'width':20, 'fun_summary_frame':'min', 'fun_summary_frames':'median' }.
        test (bool, optional): test-mode. Defaults to False.
        force (bool, optional): over-write the output. Defaults to False.

    Returns:
        dict: metadata
    """
    outp=f"{cfg['trials'][trial]['outd']}/cells/cell{celli:08d}/"
    cfgp=f"{outp}/cfg.yml"
    
    if not exists(cfgp) or force:
        cellcfg={}
        cellcfg['outp']=outp
        cellcfg['plotp']=f"{outp}/plot"
        cellcfg['cfgp']=cfgp
        cellcfg['test']=test
        cellcfg['force']=force
        cellcfg['cellbrightp']=f"{cellcfg['outp']}/cellbright.npy"       
        cellbright=cells[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
                                                              
        if not exists(dirname(cellcfg['cellbrightp'])): 
            makedirs(dirname(cellcfg['cellbrightp']),exist_ok=True)
        np.save(cellcfg['cellbrightp'], cellbright) 
        makedirs(dirname(cellcfg['outp']),exist_ok=True)
        makedirs(dirname(cellcfg['plotp']),exist_ok=True)
        
        cellcfg['cellbrightmaskp']=f"{cellcfg['outp']}/cellbrightmask.npy"
        if filterby_centroid:
            # only one cell per box
            from htsimaging.lib.utils import filter_regions
            cellbrightmask=filter_regions(cellbright.astype(int),prop_type='centroid_x',mn=70,mx=80)==0 ## get the region 
        else:
            cellbrightmask=(cellbright.astype(int))==0
            
        np.save(cellcfg['cellbrightmaskp'], cellbrightmask)
        
        cellframe_type2frames={'cellframes':{'frames':[],'ps':[]},
                               'cellframes_masked':{'frames':[],'ps':[]},
                               'cellframes_masked_substracted':{'frames':[],'ps':[]},
                              }
        for dn in sorted(list(cellframe_type2frames.keys())):
            if not exists(f"{cellcfg['outp']}/{dn}"): 
                makedirs(f"{cellcfg['outp']}/{dn}",exist_ok=True)
        from htsimaging.lib.utils import get_signal_summary_by_roi
        for framei,frame in enumerate(frames):
            cellframe=frame[cellbox[2]:cellbox[3],cellbox[0]:cellbox[1]]
            if test and framei==0:
                from htsimaging.viz.image import image_background
                import matplotlib.pyplot as plt
                
                plt.figure()
                image_background(
                    img=cellframe,
                    cmap='gfp',
                ).set(title="Original image")
                
            cellframe_masked=cellframe.copy()
            cellframe_masked[cellbrightmask]=0 ## background to zero
            if test and framei==0:
                plt.figure()
                image_background(
                    img=cellframe_masked,
                    cmap='gfp',
                ).set(title="Remove the extracellular regions")
            
            ## make the cytoplasm intensity uniform
            cellframe_masked_substracted=cellframe_masked.copy()
            signal_cytoplasm=get_signal_summary_by_roi(
                [cellframe_masked_substracted],
                **params_get_signal_summary_by_roi,
                )
            cellframe_masked_substracted=np.where(
                cellframe_masked_substracted<signal_cytoplasm*scale_signal_cytoplasm,
                signal_cytoplasm,
                cellframe_masked_substracted,
                ) 
            cellframe_masked_substracted[cellbrightmask]=0
            if test and framei==0:
                plt.figure()
                image_background(
                    img=cellframe_masked_substracted,
                    cmap='gfp',
                ).set(title="Cytoplasm intensity smoothening")                      

            cellframe_type2frames['cellframes']['frames'].append(cellframe)
            cellframe_type2frames['cellframes_masked']['frames'].append(cellframe_masked)
            cellframe_type2frames['cellframes_masked_substracted']['frames'].append(cellframe_masked_substracted)
            for dn in sorted(list(cellframe_type2frames.keys())):
                cellframe_type2frames[dn]['ps'].append(f"{cellcfg['outp']}/{dn}/frame{framei:08d}.npy")
                np.save(cellframe_type2frames[dn]['ps'][framei],
                        cellframe_type2frames[dn]['frames'][framei])

        for dn in sorted(list(cellframe_type2frames.keys())):
            cellcfg[dn]=cellframe_type2frames[dn]['ps']
        #gfp min max
        for funn in ['min','max']:
            cellcfg[f'cellgfp{funn}p']=f"{cellcfg['outp']}/cellgfp{funn}.npy"
            frame=getattr(np,f"a{funn}")(cellframe_type2frames['cellframes_masked_substracted']['frames'],axis=0)
            np.save(cellcfg[f'cellgfp{funn}p'], frame)
        del frame

        cellcfg['signal_cytoplasm']=get_signal_summary_by_roi(cellframe_type2frames['cellframes_masked']['frames'],
                                                              **params_get_signal_summary_by_roi)

        to_dict(cellcfg,cellcfg['cfgp'])
    else:
        cellcfg=read_dict(cfgp)
    return cellcfg