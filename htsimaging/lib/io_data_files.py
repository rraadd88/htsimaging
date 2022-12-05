#!/usr/bin/env python

from os.path import splitext, join, exists, isdir,basename,abspath,dirname
from os import makedirs
import numpy as np
import pandas as pd
from glob import glob

import logging

import pickle

def to_pkl(
    data,
    fh,
    ):
    if not fh is None:
        with open(fh, 'wb') as f:
            try:
                pickle.dump(data, f, -1)
            except:
                pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
                
def read_pkl(
    fh,
    ):
    with open(fh,'rb') as f:
        return pickle.load(f) 


def expt_dh2expt_info(
    expt_dh: str,
    ):
    """
    Make configuration using the directory structure for an experiment.
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
    ):
    """
    Create information file.
    """
    info=pd.read_csv(expt_dh+"info")
    for i in range(len(info)):
        reps=glob("%s/%s*" % (info.loc[i,"dh"],info.loc[i,"fn_lead"]))
        for repi,rep in enumerate(reps):
            info.loc[i,"replicate %d" % (repi+1)]=rep
    #     break
    info.to_csv(expt_dh+"info")

