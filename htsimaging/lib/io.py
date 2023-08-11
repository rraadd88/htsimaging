#!/usr/bin/env python
"""I/O"""

import os
from os import makedirs
from os.path import exists,splitext,basename

import numpy as np
import matplotlib.pyplot as plt

import subprocess
from skimage import io

## image io
def read_image(
    imp: str
    ):
    """Read image.

    Args:
        imp (str): path to the image file.

    Returns:
        np.array
        
    TODOs: 
        For a tiff file:
        from skimage.external import tifffile
    """
    from skimage import io
    if isinstance(imp,(str)) and not isinstance(imp,(list)):
        if not imp.endswith('.npy'):
            im=io.imread(imp,as_gray=True)
        else:
            im=np.load(imp)
    else:
        im=imp
    return im

## make a video
def arr2vid(
    arr_list: list,
    regions: list,
    kins_mean: float,
    vid_fh: str,
    xpixels: list,
    ypixels: list,
    dpi: int = 100,
    ) -> str:
    """From array to video.

    Args:
        arr_list (list): list of frames.
        regions (list): regions 
        kins_mean (float): kinetics 
        vid_fh (str): video file path
        xpixels (list): pixels allong x-axis.
        ypixels (list): pixels allong y-axis.
        dpi (int, optional): DPI resolution. Defaults to 100.

    Returns:
        str: path of the video
    """
    png_dh=os.path.splitext(vid_fh)[0]
    if not os.path.exists(png_dh):
        try:
            os.makedirs(png_dh)
        except :
            logging.warning("error data_lbl")
    fig = plt.figure(figsize=(ypixels/dpi*2, xpixels/dpi), dpi=dpi)
    ax_img = plt.subplot(1,2,1)
    ax_kin = plt.subplot(1,2,2)
    ax_img.set_axis_off()
    ax_img.set_aspect('equal')
    for i in range(len(arr_list)):
        ax_img.imshow(arr_list[i],cmap='gray',animated=True)
        ax_img.contour(regions, [0.25], linewidths=1.2, colors='r',animated=False)
        if len(kins_mean.columns)>1:
            kins_mean.plot(x='time',legend=False,ax=ax_kin)
            ax_kin.plot(kins_mean['time'],kins_mean.drop(['time'], axis=1).mean(axis=1),lw=6,color='k')
            ax_kin.set_ylabel("Fluorescence Intensity (FU)")
        ax_kin.set_xlim([kins_mean.loc[0,'time'],kins_mean.loc[len(kins_mean)-1,'time']])
        ax_kin.axvline(kins_mean.loc[i,'time'], color='r', linestyle='--',lw=2)
        fig.subplots_adjust(wspace=.4)
        fig.savefig(png_dh+'/%02d.png' % i)
        ax_kin.clear()
    bash_command=("ffmpeg -f image2 -r 4 -i "+png_dh+"/%02d.png -vcodec mpeg4 -y "+vid_fh)
    subprocess.Popen(bash_command, shell=True, executable='/bin/bash')
    return vid_fh
    
def makevid(
    gfp_list_stb: list,
    brf_list_stb: list,
    cmap_gfp: str,
    cmap_brf: str,
    vid_fh: str,
    conditionn: int=None,
    interval=None,
    dpi : int= 300,
    ) -> str:
    """Convert to a video.

    Args:
        gfp_list_stb (list): channel (e.g. GFP) images.  
        brf_list_stb (list): segmented (e.g. bright-field) images. 
        cmap_gfp (str): colormap for the channel images.
        cmap_brf (str): colormap for the segmented images.
        vid_fh (str): path to the video file.
        conditionn (int, optional): title. Defaults to None.
        interval (_type_, optional): inerval of the frames. Defaults to None.
        dpi (int, optional): DPI resolution. Defaults to 300.

    Returns:
        str: path to the video file.
    """
    png_dh=splitext(vid_fh)[0]
    if not exists(png_dh):
        makedirs(png_dh)
    for i in range(len(gfp_list_stb)):
        fig = plt.figure(dpi=dpi)
        ax=plt.subplot(111)    
        ax.imshow(gfp_list_stb[i],cmap=cmap_gfp,zorder=1,alpha=0.5)
        ax.imshow(brf_list_stb[i],cmap=cmap_brf,zorder=0,alpha=0.5)
        ax.text(2,2,conditionn,color='r',fontsize=20,)
        ax.text(2,np.shape(gfp_list_stb[i])[1]*0.95,"time: %03d s" % ((i+1)*interval),color='b',fontsize=20,)
#         fig.subplots_adjust(wspace=.4)
        plt.tight_layout()
        plt.axis('off')
        fig.savefig('%s/%03d.png' % (png_dh,i))
        plt.close(fig)
    bash_command=("ffmpeg -f image2 -r 20 -i "+png_dh+"/%03d.png -vcodec mpeg4 -y "+vid_fh)
    print(vid_fh)
    subprocess.Popen(bash_command, shell=True)
    return vid_fh

def nd2arr_list(
    nd_dh: str=None,
    nd_fns: list=[],
    nd_fh: str=None,
    ) -> list:
    """Raw image to list of arrays.

    Args:
        nd_dh (str, optional): directory containing raw files e.g. nd2. Defaults to None.
        nd_fns (list, optional): file names. Defaults to [].
        nd_fh (str, optional): path to the files. Defaults to None.

    Returns:
        list: list of arrays
    """
    import nd2reader
    arr_list=[]
    if nd_fh is None:
        for nd_fn in nd_fns:
            nd = nd2reader.Nd2("%s/%s" % (nd_dh,nd_fn))
            for ndi in nd:
                arr_list.append(np.array(ndi))
            del nd
        return arr_list
    elif not nd_fh is None:
        nd = nd2reader.Nd2("%s" % (nd_fh))
        for ndi in nd:
            arr_list.append(np.array(ndi))
        return arr_list

def to_csv(
    fh_xls='../test/test.xlsx',
    nd2_dh="/media/Transcend/20160219_000356_267",
    cores=16,
    ):
    """
    Convert nd2 files to csv using parallel processing.

    Args:
        fh_xls (str, optional): metadata file. Defaults to '../test/test.xlsx'.
        nd2_dh (str, optional): path of the directory containing raw images. Defaults to "/media/Transcend/20160219_000356_267".
        cores (int, optional): number of cores. Defaults to 16.
    """
    import nd2reader
    import pandas as pd
    data_job=pd.read_excel(fh_xls,'JobView')
    nd_fhs=[ nd2_dh+'/'+str(fn) for fn in data_job['File Name'].unique()]

    # convert t0 csvs
    def nd2csv(
        nd_fh: str,
        ):
        nd = nd2reader.Nd2(nd_fh)
        frame_cnt=0
        for frame in nd:
            np.savetxt("%s/%s%02d.csv" % (nd2_dh,basename(nd_fh),frame_cnt), frame, delimiter=",")
            frame_cnt+=1
        nd.close
    from multiprocessing import Pool
    pool_data_lbl2data_fit=Pool(processes=int(cores)) 
    pool_data_lbl2data_fit.map(nd2csv,nd_fhs)
    pool_data_lbl2data_fit.close(); pool_data_lbl2data_fit.join()