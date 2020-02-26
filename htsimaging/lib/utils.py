from rohan.global_imports import *
import pims
# import nd2reader
import string
from scipy import stats,ndimage
from multiprocessing import Pool
from skimage.segmentation import random_walker
# from skimage.data import binary_blobs
from skimage import io,exposure,restoration,filters,morphology,measure
from scipy import ndimage as ndi
import matplotlib.animation as animation
from matplotlib import colors
import subprocess

import trackpy as tp
from htsimaging.lib import spt

def filterframes(frames,cutoff=0):
    frames_cleaned=[]
    for framei in frames:
        frames_cleaned.append(filterframe(framei,cutoff))
    return frames_cleaned

def filterframe(frame,cutoff=0):
    frame_cleand=np.array(frame)
    frame_cleand[frame_cleandi<cutoff]=0
    return frame_cleaned

def get_data_by_regions(regions,img=None,prop_type='area'):
    if prop_type=='mean_intensity':
        ValueError("arg img is required")
    regions_props= measure.regionprops(regions.astype(int),intensity_image=img)
    regions_lbls = np.array([prop.label for prop in regions_props])
    if prop_type=='area':
        regions_props_selected = np.array([getattr(prop,prop_type)  for prop in regions_props])
    df1=pd.DataFrame({'labels':regions_lbls,
                      prop_type:regions_props_selected,
                     })
    return df1

def filter_regions(regions,img=None,prop_type='area',mn=0,mx=0,
                   check=False,plotp=None):
    if prop_type=='mean_intensity':
        ValueError("arg img is required")
    regions_props= measure.regionprops(regions.astype(int),intensity_image=img)
    regions_lbls = np.array([prop.label for prop in regions_props])
    if prop_type=='area':
        regions_props_selected = np.array([prop.area  for prop in regions_props])
    elif prop_type=='mean_intensity':
        regions_props_selected = np.array([prop.mean_intensity  for prop in regions_props])
    elif prop_type=='eccentricity':
        regions_props_selected = np.array([prop.eccentricity  for prop in regions_props])        
    elif prop_type=='extent':
        regions_props_selected = np.array([prop.extent  for prop in regions_props])     
    elif prop_type=='centroid_x':
        regions_props_selected = np.array([prop.centroid[0] for prop in regions_props])     
    elif prop_type=='centroid_y':
        regions_props_selected = np.array([prop.centroid[1] for prop in regions_props])     

    if check:
        print(regions_props_selected)
        print(regions_lbls)
        
    regions_filtered_lbls = regions_lbls[np.where((regions_props_selected<mx) & (regions_props_selected>mn))[0]]
    regions_filtered=np.zeros(regions.shape)
    for lbli in regions_filtered_lbls:
        regions_filtered[np.where(regions==lbli)]=lbli
    if check:
        plt.figure()
        ax=plt.subplot(111)        
        ax.hist(regions_props_selected)
        plt.figure()
        ax=plt.subplot(111)
        if not img is None:
            ax.imshow(img, cmap=plt.cm.gray, interpolation='nearest')        
        ax.contour(regions, [0.5], linewidths=1.2, colors='r')
        ax.contour(regions_filtered, [0.5], linewidths=1.2, colors='g')
        plt.tight_layout()
        if not plotp is None:
            makedirs(dirname(plotp),exist_ok=True)
            plt.savefig(plotp)
    return regions_filtered

def smoothen(img):
    denoised=restoration.denoise_bilateral(img.astype('uint16'), sigma_range=0.01, sigma_spatial=15)
    smoothened = filters.median(denoised,np.ones((4,4)))
    return smoothened

def smoothenframes(frames):
    frames_cleaned=[]
    for framei in frames:
        frames_cleaned.append(smoothen(np.array(framei)))
    return frames_cleaned
        
def get_regions(img):
    denoised=restoration.denoise_bilateral(img.astype('uint16'), 
#                                            sigma_range=0.01, 
                                           sigma_spatial=15,
                                          multichannel=False)
    smoothened = filters.median(denoised,np.ones((4,4)))
    markers = np.zeros(smoothened.shape, dtype=np.uint)
# otsu works only for only for multi-channel images     
#     markers[smoothened < filters.threshold_otsu(smoothened)] = 1
#     markers[smoothened > filters.threshold_otsu(smoothened)] = 2
    markers[smoothened < filters.median(smoothened)] = 1
    markers[smoothened > filters.median(smoothened)] = 2

    labels = random_walker(smoothened, markers, beta=10, mode='bf')
    regions= measure.label(labels)
    return regions, denoised, smoothened,markers

def raw2phasecorr(arr_list,clip=0): #cv
    import cv2
    cx = 0.0
    cy = 0.0
    stb_arr_list=[]
    prev_frame = arr_list[0]
    prev_image = np.float32(restoration.denoise_tv_chambolle(prev_frame.astype('uint16'), weight=0.1, multichannel=True)) #ref
    for frame in arr_list:           
        image = np.float32(restoration.denoise_tv_chambolle(frame.astype('uint16'), weight=0.1, multichannel=True))
        # TODO: set window around phase correlation
        dp = cv2.phaseCorrelate(prev_image, image)
        cx = cx - dp[0]
        cy = cy - dp[1]
        xform = np.float32([[1, 0, cx], [0, 1, cy]])
        stable_image = cv2.warpAffine(frame.astype('float32'), xform, dsize=(image.shape[1], image.shape[0]))
        prev_image = image
        #clip sides
        ht,wd=np.shape(stable_image)
#         clip=0.125 #0.25
        lt=int(wd*clip)
        rt=int(wd-wd*clip)
        up=int(ht*clip)
        dw=int(ht-ht*clip)
        stable_image_clipped=stable_image[up:dw,lt:rt]
        stb_arr_list.append(stable_image_clipped)
    return stb_arr_list

def imlistcropper(imlist,loci):
    rowini,rowend,colini,colend=loci
    imlist_region=[]
    for im in imlist:
        imlist_region.append(im[rowini:rowend,colini:colend])
    return imlist_region


def imclipper(im_stb,clip):
    ht,wd=np.shape(im_stb)
#         clip=0.125 #0.25
    lt=int(wd*clip)
    rt=int(wd-wd*clip)
    up=int(ht*clip)
    dw=int(ht-ht*clip)
    im_stb_clipped=im_stb[up:dw,lt:rt]
    return im_stb_clipped

def phasecorr(imlist,imlist2=None,clip=0): #cv  [rowini,rowend,colini,colend]
    import cv2
    cx = 0.0
    cy = 0.0            
    imlist_stb=[]
    if imlist2!=None:
        imlist2_stb=[]

    imi=0
    im_prev = imlist[0]
    im_denoised_prev = np.float32(restoration.denoise_tv_chambolle(im_prev.astype('uint16'), weight=0.1, multichannel=True)) #ref
    for im in imlist:           
        im_denoised = np.float32(restoration.denoise_tv_chambolle(im.astype('uint16'), weight=0.1, multichannel=True))
        # TODO: set window around phase correlation
        dp = cv2.phaseCorrelate(im_denoised_prev, im_denoised)
        cx = cx - dp[0]
        cy = cy - dp[1]
        xform = np.float32([[1, 0, cx], [0, 1, cy]])
        im_stb = cv2.warpAffine(im.astype('float32'), xform, dsize=(im_denoised.shape[1], im_denoised.shape[0]))
        imlist_stb.append(imclipper(im_stb,clip))

        if imlist2!=None:
            im2=imlist2[imi]
            im2_stb=cv2.warpAffine(im2.astype('float32'), xform, dsize=(im_denoised.shape[1], im_denoised.shape[0]))
            imlist2_stb.append(imclipper(im2_stb,clip))

        im_denoised_prev = im_denoised
        imi+=1
    if imlist2!=None:
        return imlist_stb,imlist2_stb
    else:
        return imlist_stb

## features
def get_cellprops(regions,intensity_imgtype2img,properties=['area',
    'bbox_area',
    'convex_area',
    'eccentricity',
    'equivalent_diameter',
    'euler_number',
    'extent',
    'filled_area',
    'label',
    'major_axis_length',
    'max_intensity',
    'mean_intensity',
    'min_intensity',
    'minor_axis_length',
    'orientation',
    'perimeter',
    'solidity',
    'centroid']):
    from skimage.external import tifffile
    dn2df={}
    for imgtype in intensity_imgtype2img:
        df=pd.DataFrame(measure.regionprops_table(regions.astype(int),
                                           intensity_image=intensity_imgtype2img[imgtype],
                                           properties=properties))
        df.index=df.index+1
        dn2df[imgtype]=dmap2lin(df,coln='property',idxn='cell #')
    df=pd.concat(dn2df,axis=0,names=['image type']).droplevel(1)#.reset_index()
    return df

def get_signal_summary_by_roi(cellframes,xy_center=None,width=20,
                       fun_summary_frame='min',fun_summary_frames='median'):
    """
    place of the roi in the image is defined by
    :param xy_center:
    :param width:    
    """
    width_half=int(width*0.5)
    frame=cellframes[0]
    if xy_center is None:
        xy_center=[int(i/2) for i in cellframes[0].shape]
    signal_cytoplasm=getattr(np,fun_summary_frames)([getattr(np,fun_summary_frame)(frame[xy_center[1]-width_half:xy_center[1]+width_half,xy_center[0]-width_half:xy_center[0]+width_half]) for frame in cellframes])
    return int(signal_cytoplasm)