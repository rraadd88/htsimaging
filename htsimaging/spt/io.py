#!/usr/bin/env python
"""I/O"""

from os.path import exists
from glob import glob

import numpy as np
import pandas as pd 

def to_frames(
    input_path: str,
    channeli=None,
    ):
    """
    Convert to frames.

    Args:
        input_path (str): path to the raw data.

    Returns:
        list: list of frames.
    """
    if input_path.endswith('nd2'):
        frames=pims.ND2_Reader(input_path)
    elif input_path.endswith('mp4'):
        import imageio
        vid = imageio.get_reader(input_path,  'ffmpeg')
        frames=[frame for frame in vid.iter_data()]
        if not channeli is None:
            # get a channel 
            frames=[frame[:,:,1] for frame in vid.iter_data()]
    if input_path.endswith('nd2'):
        if len(np.shape(frames))==4:
            frames = map(lambda x: x.mean(axis=0), frames)
    # else:
    #     frames.default_coords['c'] = 1
    #     frames.bundle_axes='yx'
    #     frames.iter_axes = 't'
    return frames

def expt_dh2expt_info(
    expt_dh: str,
    ) -> pd.DataFrame:
    """
    Make configuration using the directory structure for an experiment.

    Args:
        expt_dh (str): str

    Returns:
        pd.DataFrame: output dataframe.
    """
    if exists(expt_dh+"info"):
        expt_info=pd.read_csv(expt_dh+"info")
    elif exists(expt_dh+"info.csv"):
        expt_info=pd.read_csv(expt_dh+"info.csv")
    cols_del=["dh","fn_lead"]
    for col in expt_info.columns.tolist():
        if "Unnamed" in col:
            cols_del.append(col)
    for col in cols_del:
        if col in expt_info.columns.tolist():
            del expt_info[col]    
    expt_info2=expt_info.drop(["smpn"],axis=1).T
    expt_info2.columns=expt_info.loc[:,"smpn"]
    
    # print expt_info2
    for col in expt_info2.columns.tolist():
        for i in range(len([i for i in expt_info2.index if 'replicate' in i])):
            if not pd.isnull(expt_info2.loc["replicate %d" % (i+1),col]):
                expt_info2.loc[("replicate %d" % (i+1)),col]=['replicate %d' % (i+1),expt_info2.loc[("replicate %d" % (i+1)),col]]
            else:
                expt_info2.loc[("replicate %d" % (i+1)),col]=[np.nan,np.nan]
    expt_info2=expt_info2.reset_index()
    del expt_info2["index"]
    return expt_info2

def createinfo(
    expt_dh: str,
    ) -> str:
    """
    Create information file.

    Args:
        expt_dh (str): path to the directory containing the raw data.

    Returns: 
        path to the file containing the metadata. 
    """
    info=pd.read_csv(expt_dh+"info")
    for i in range(len(info)):
        reps=glob("%s/%s*" % (info.loc[i,"dh"],info.loc[i,"fn_lead"]))
        for repi,rep in enumerate(reps):
            info.loc[i,"replicate %d" % (repi+1)]=rep
    #     break
    output_path=expt_dh+"info"
    info.to_csv(output_path)
    return output_path