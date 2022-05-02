import os
from roux.lib.io_files import copyfile
from roux.lib.io_sys import runbash
from roux.global_imports import *
from skimage import io,measure
from htsimaging.lib.utils import filter_regions

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
from roux.global_imports import *
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
#     print(brightp2copyfromto)
        
    for p in brightp2copyfromto:
        top=f"{cfg['segmentation_cell']['inputd']}/{basename(p)}"
        if not exists(top):
#             print(p,top)   
            copyfile(p,top)
#             symlink(p,top)        
        else:
            logging.warning(f"exists: {top}")
    set_opts(cfg['segmentation_cell']['inputd'],cfg['segmentation_cell']['outputd'],
             yeast_segmentation_srcdp=cfg['yeastspotter_srcd'])
    print('running yeastspotter')
    runbash(f"conda activate htsimaging;python {cfg['yeastspotter_srcd']}/segmentation.py",test=test)
#                 &> {cfg['segmentation_cell']['parentd']}/log_yeastspotter.txt
#     print(brightp2copyfromto)
    for p in brightp2copyfromto:
        if exists(brightp2copyfromto[p][0]):
            copyfile(brightp2copyfromto[p][0],brightp2copyfromto[p][1]) 
    return cfg

def segmentation2cells(imp,imsegp,fiterby_border_thickness=100,magnification=100,plotp=None,
                      **kws):
    """
    prop_type='area',mn=100,mx=8000
    at 1.5X
    prop_type='area',mn=1500,mx=12000
    """
    if isinstance(imp,(str)):
        im=io.imread(imp,as_gray=True)
    else:
        im=imp
    if isinstance(imsegp,(str)):        
        imseg=io.imread(imsegp,as_gray=True)
    else:
        imseg=imsegp
        
    regions=measure.label(imseg)
#     fiterby_border_thickness,im.shape[1]+fiterby_border_thickness
    regions=filter_regions(regions,im,prop_type='area',
                           mn=10*magnification,
                           mx=100000*magnification,
                           plotp=plotp,**kws)
    regions=filter_regions(regions,im,prop_type='eccentricity',mn=0,mx=0.8,check=False)
    regions=filter_regions(regions,im,prop_type='centroid_x',
                           mn=fiterby_border_thickness,mx=im.shape[0]+fiterby_border_thickness,check=False)
    regions=filter_regions(regions,im,prop_type='centroid_y',
                           mn=fiterby_border_thickness,mx=im.shape[1]+fiterby_border_thickness,check=False)
    return regions               
               
def get_cellboxes(regions,plotp=None,cellbox_width=150):
    import matplotlib.patches as mpatches
    if not plotp is None:
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.imshow(regions,cmap='binary')
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
               