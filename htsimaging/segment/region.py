#!/usr/bin/env python
"""Processing of the segmented regions."""

import logging
import numpy as np
from skimage import measure

from htsimaging.lib.utils import filter_regions
from htsimaging.lib.io import read_image

def segmentation2cells(
    # imp: str,
    imsegp: str,
    kind : str = 'yeast',
    fiterby_border_thickness: int=None,
    magnification: int=100,
    test: bool=False,
    **kws: dict,
    ) -> list:
    """
    Segment the image to the single cells.

    Args:
        imsegp (str): _description_
        fiterby_border_thickness (int, optional): _description_. Defaults to 100.
        magnification (int, optional): _description_. Defaults to 100.
        plotp (str, optional): _description_. Defaults to None.

    Returns:
        list: _description_

    Examples:
        1. Parameters:
            prop_type='area',mn=100,mx=8000
            at 1.5X
            prop_type='area',mn=1500,mx=12000
    """
    # im=read_image(imp)
    imseg=read_image(imsegp)
        
    regions=measure.label(imseg)
#     fiterby_border_thickness,im.shape[1]+fiterby_border_thickness
    if kind=='yeast':
        if test: logging.info(f"number of regions = {len(np.unique(regions))}")
        regions=filter_regions(regions,imseg,prop_type='area',
                               mn=10*magnification,
                               mx=100000*magnification,
                               plotp=plotp,
                               **kws,
                              )
        if test: logging.info(f"number of regions = {len(np.unique(regions))}")
        regions=filter_regions(regions,imseg,prop_type='eccentricity',mn=0,mx=0.8,test=test,
                               **kws,
                              )
    if test: logging.info(f"number of regions = {len(np.unique(regions))}")
    if not fiterby_border_thickness is None:
        regions=filter_regions(regions,imseg,prop_type='centroid_x',
                               mn=fiterby_border_thickness,mx=imseg.shape[0]+fiterby_border_thickness,test=test,
                               **kws,
                              )
        if test: logging.info(f"number of regions = {len(np.unique(regions))}")
        regions=filter_regions(regions,imseg,prop_type='centroid_y',
                               mn=fiterby_border_thickness,mx=imseg.shape[1]+fiterby_border_thickness,test=test,
                               **kws,
                              )
        if test: logging.info(f"number of regions = {len(np.unique(regions))}")
    if len(np.unique(regions))==1:
        logging.error("the output does not contain any regions")
    return regions               
               
def get_cellboxes(
    regions: list,
    cellbox_width: int=150,
    test: bool=False,
    ) -> list:
    """
    Get the bounding boxes of the cells.

    Args:
        regions (list): regions.
        cellbox_width (int, optional): width of the bounding box of the cell. Defaults to 150.
        test (bool, optional): test-mode. Defaults to False.
    Returns:
        list: list of the bounding boxes for cells. 
    """
    import matplotlib.patches as mpatches
    if test:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.imshow(regions,cmap='binary')
    cellboxes=[]
    for regioni,region in enumerate(measure.regionprops(regions.astype(int))):
        box_xmnxmxymnymx=[
            region.centroid[1]-(cellbox_width*0.5),
            region.centroid[1]+(cellbox_width*0.5),
            region.centroid[0]-(cellbox_width*0.5),
            region.centroid[0]+(cellbox_width*0.5),
                         ]
        cellboxes.append([int(i) for i in box_xmnxmxymnymx])
#         celli2props[regioni+1]=region.area
        if test:
            rect = mpatches.Rectangle([box_xmnxmxymnymx[0],box_xmnxmxymnymx[2]], cellbox_width, cellbox_width,
                                      fill=False, edgecolor='red', linewidth=2)
            ax.text(region.centroid[1],region.centroid[0],f"{regioni+1:d}",color='g')
#             ax.text(region.centroid[1],region.centroid[0],f"{region.extent:.2f}",color='g')
            ax.add_patch(rect)
    return cellboxes

def arr_list2regions(
    arr_list: list,
    time_increment: int,
    ) -> tuple:
    """Parameterized cell-segmentation for the time lapse images. 

    Args:
        arr_list (list): frames of images.
        time_increment (int): time interval.

    Returns:
        tuple: regions and table with intensities.
        
    """
    pre_bleach=arr_list[0]
    denoised=restoration.denoise_bilateral(pre_bleach.astype('uint16'), sigma_range=0.01, sigma_spatial=15)
    smoothened = filters.median(denoised,np.ones((4,4)))
    markers = np.zeros(smoothened.shape, dtype=np.uint)
    markers[smoothened < filters.threshold_otsu(smoothened)] = 1
    markers[smoothened > filters.threshold_otsu(smoothened)] = 2
    labels = random_walker(smoothened, markers, beta=10, mode='bf')
    regions= measure.label(labels)
    props = measure.regionprops(regions,arr_list[0]) #ref
    regions_areas=np.array([prop.area  for prop in props])
    regions_lbls =np.array([prop.label for prop in props])
    regions_means=np.array([prop.mean_intensity  for prop in props])
    regions_index_large=np.where((regions_areas<10000) & (regions_areas>200) & (regions_means>2000))[0]
    regions_lbls_large=regions_lbls[regions_index_large] 
    regions_large=np.zeros((regions.shape), dtype=bool)
    for i in regions_lbls_large:
        booli = (regions == i) 
        regions_large=np.logical_or(regions_large,booli)
    kins_mean=pd.DataFrame(columns=regions_lbls_large, index=range(len(arr_list)))
    for i in range(len(arr_list)):
        props = measure.regionprops(regions,arr_list[i])
        means=np.array([prop.mean_intensity for prop in props])
        kins_mean.loc[i,:]=means[regions_index_large]
        del props
    kins_mean=kins_mean.loc[:, ~(kins_mean < 3000).any(axis=0)] #stitch
    kins_mean['time']=np.array(range(len(arr_list)))*time_increment #stitch
    return regions_large,kins_mean
