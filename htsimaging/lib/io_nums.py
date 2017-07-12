#!usr/bin/python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com,rohan@igib.in>
# This program is distributed under General Public License v. 3.  

"""
================================
``io_nums``
================================
"""
import numpy as np
import pandas as pd 

def is_numeric(obj):
    """
    This detects whether an input object is numeric or not.
    """
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        return True
    
def str2num(x):
    """
    This extracts numbers from strings. eg. 114 from M114R.
    """
    return int(''.join(ele for ele in x if ele.isdigit()))

def plog(x,p = 0.5):
    return np.log2(x+p)

def glog(x,l = 2):
    return np.log((x+np.sqrt(x**2+l**2))/2)/np.log(l)

def float2int(x):
    if not pd.isnull(x):
        if is_numeric(x):
            x=int(x)
    return x    
