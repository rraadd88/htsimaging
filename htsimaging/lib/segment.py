import os
from rohan.dandage.io_strs import get_datetime
from rohan.dandage.io_files import copy
from rohan.dandage.io_sys import runbashcmd
from rohan.global_imports import *

def set_opts(indp,outdp,yeast_segmentation_srcdp):
    cfg={'input_directory':abspath(indp),
    'output_directory':abspath(outdp),
    'rescale':False,
    'scale_factor':2,
    'save_preprocessed':False,
    'save_compressed':False,
    'save_masks':True,
    'verbose':True,}
    optspyp=f"{yeast_segmentation_srcdp}/opts.py"
    with open(optspyp,'w') as f:    
        for k in cfg:
            if isinstance(cfg[k],str):
                cfg[k]=f"'{cfg[k]}'"
            f.write(f"{k}={cfg[k]}\n")
            
from os import symlink
from rohan.global_imports import *
def run_yeastspotter(cfg,test=False):
    cfg['segmentation_cell']={}
    cfg['segmentation_cell']['parentd']=f"{cfg['prjd']}/segmentation_cell"
    cfg['segmentation_cell']['inputd']=f"{cfg['segmentation_cell']['parentd']}/yeastspotter_input"
    cfg['segmentation_cell']['outputd']=f"{cfg['segmentation_cell']['parentd']}/yeastspotter_output"
    makedirs(cfg['segmentation_cell']['inputd'],exist_ok=True)
    makedirs(cfg['segmentation_cell']['outputd'],exist_ok=True)
    
    brightp2copyfromto={}
    for triald in cfg['trials']:
        brightps=cfg['trials'][triald]['bright']
        for brightp in brightps:
            brightp2copyfromto[brightp]=[f"{cfg['segmentation_cell']['outputd']}/masks/{basename(brightp)}",f"{dirname(brightp)}/{basename(brightp)}.segmented.tif"]
        cfg['trials'][triald]['bright_segmented']=[f"{dirname(p)}/{basename(p)}.segmented.tif" for p in brightps]
    print(brightp2copyfromto)
        
    for p in brightp2copyfromto:
        top=f"{cfg['segmentation_cell']['inputd']}/{basename(p)}"
        if not exists(top):
            print(p,top)        
            symlink(p,top)        
        else:
            logging.warning(f"exists: {top}")
    set_opts(cfg['segmentation_cell']['inputd'],cfg['segmentation_cell']['outputd'],
             yeast_segmentation_srcdp=cfg['yeastspotter_srcd'])

    runbashcmd(f"source activate htsimaging;python {cfg['yeastspotter_srcd']}/segmentation.py")
    for p in brightp2copyfromto:
        copy(brightp2copyfromto[p][0],brightp2copyfromto[p][1]) 
    return cfg