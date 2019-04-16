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
def run_yeastspotter(cfg,yeast_segmentation_srcdp,indp,outdp,test=False):
    brightps=[]
    for trial in cfg['trials'].keys():
        if 'bright' in cfg['trials'][trial]:
            if len(cfg['trials'][trial]['bright'])!=0:
                segmenteds=[]
                segmented_srcs=[]
                for p in cfg['trials'][trial]['bright']:
                    segmenteds.append(f"{dirname(p)}/{basenamenoext(p)}_segmented.tif")
                    segmented_srcs.append(f"{outdp}/masks/{basename(p)}")
                cfg['trials'][trial]['segmented']=segmenteds
                cfg['trials'][trial]['segmented_src']=segmented_srcs
                brightps.append(cfg['trials'][trial]['bright'])            

    makedirs(indp,exist_ok=True)
    if exists(outdp):
        os.rename('/'+outdp.strip('/'),f"{outdp}_{get_datetime()}")
    set_opts(indp,outdp,
             yeast_segmentation_srcdp=yeast_segmentation_srcdp)
    if test:
        print(np.ravel(brightps))  
    for p in list(np.ravel(brightps)):
        copy(p,f"{indp}/{basename(p)}")

    runbashcmd(f"source activate htsimaging;python {yeast_segmentation_srcdp}/segmentation.py")
    #     break

    for trial in cfg['trials'].keys():
        if 'bright' in cfg['trials'][trial]:
            if len(cfg['trials'][trial]['bright'])!=0:
                for outp,inp in zip(cfg['trials'][trial]['segmented'],cfg['trials'][trial]['segmented_src']):
                    copy(inp,outp) 
    return cfg