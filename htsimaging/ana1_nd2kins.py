import nd2reader
import pandas as pd
import string
import numpy as np
from scipy import stats
import os
from multiprocessing import Pool
from os.path import splitext,basename
import cv2

from skimage.segmentation import random_walker
# from skimage.data import binary_blobs
# import skimage
from skimage import io,exposure,restoration,filters,morphology,measure

fh_xls='../../../data/yeast_gfp_half_life/data/160305_bleachto60ormore/160305_bleach10minB06.xlsx'
data_job=pd.read_excel(fh_xls,'JobView')

nd_dh='../../../data/yeast_gfp_half_life/data/160305_bleachto60ormore/160305_bleachchase_bleach10min/20160305_234331_875'
nd_fns=['WellB06_Seq0000.nd2',\
        'WellB06_Seq0002.nd2',\
        'WellB06_Seq0007.nd2']

def raw2phasecorr(arr_list):
    cx = 0.0
    cy = 0.0
    stb_arr_list=[]
    prev_image = numpy.float32(arr_list[0]) #ref
    for frame in arr_list:
        image = numpy.float32(frame)
        # TODO: set window around phase correlation
        dp = cv2.phaseCorrelate(prev_image, image)
        cx = cx - dp[0]
        cy = cy - dp[1]
        xform = numpy.float32([[1, 0, cx], [0, 1, cy]])
        stable_image = cv2.warpAffine(image, xform, dsize=(image.shape[1], image.shape[0]))
        prev_image = image
        stb_arr_list.append(stable_image)
    return stb_arr_list

def arrlist2vid(arrlist,vid_fh):
    dpi = 100
    xpixels, ypixels = 512, 512
    fig = plt.figure(figsize=(ypixels/dpi, xpixels/dpi), dpi=dpi)
    ax = plt.Axes(fig, [0., 0., 1, 1])
    fig.add_axes(ax)
    ax.set_axis_off()
    ax.set_aspect('equal')
    im = ax.imshow(arr_list[0],cmap='gray')
    # plt.savefig('test.png')

    def update_img(n):
        tmp = arr_list[1:][n]
        im.set_data(tmp)
        return im

    # #legend(loc=0)
    ani = animation.FuncAnimation(fig,update_img,np.arange(len(arr_list)-1),interval=60,blit=False)
    writer = animation.writers['ffmpeg'](fps=1)

    ani.save(vid_fh,writer=writer)#dpi=dpi

def arr_list2kins(arr_list)
	pre_bleach=arr_list[0]
	smoothened = filters.median(pre_bleach.astype(uint16),np.ones((10,10)))
	markers = np.zeros(smoothened.shape, dtype=np.uint)
	markers[smoothened < filters.threshold_otsu(edges)] = 1
	markers[smoothened > filters.threshold_otsu(edges)] = 2
	labels = random_walker(smoothened, markers, beta=10, mode='bf')
	regions=measure.label(labels)
	kins_area=pd.DataFrame(index=range(lbl.max()), columns=range(len(arr_list)))
	kins_mean=pd.DataFrame(index=range(lbl.max()), columns=range(len(arr_list)))

	for i in range(len(arr_list)):
		props = measure.regionprops(regions,arr_list[i])
		kins_area.loc[:,i]=np.array([prop.area for prop in props])
		kins_mean.loc[:,i]=np.array([prop.intensity_mean for prop in props])
	return kins_area,kins_mean

# for well in wells:
arr_list=[]
for nd_fn in nd_fns:
    nd = nd2reader.Nd2("%s/%s" % (nd_dh,nd_fn))
    for ndi in nd:
        arr_list.append(np.array(ndi))
    del nd
stb_arr_list=raw2phasecorr(arr_list)
kins_area,kins_mean=arr_list2kins(stb_arr_list)


