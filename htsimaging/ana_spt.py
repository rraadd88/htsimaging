#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
import logging
import pandas as pd
from os.path import exists, basename, dirname, abspath
import numpy as np
import trackpy as tp
from glob import glob

logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # 
# from htsimaging.lib import spt #.expt2plots,spt.fit_power,spt.fit_line
from htsimaging.lib.spt import expt2plots,expt_dh2expt_info,flt_traj

def main(expt_dh):
    # expt_dh=expt_dh
    if not exists(expt_dh):    
        if expt_dh.count('/')<2:
            expt_dh='%s/../test_dataset/%s' % (dirname(abspath(__file__)),expt_dh) 
            print expt_dh     
    if exists(expt_dh):    
        # print expt_dh
        _cfg_fh='%s/_cfg.csv' % expt_dh
        if exists(_cfg_fh):
            _cfg=pd.read_csv(_cfg_fh).set_index('var')
            _cfg=_cfg['val'].to_dict()
            # print _cfg #debug
        else:
            _cfg={}

        expt2plots(expt_dh2expt_info(expt_dh),expt_dh,_cfg=_cfg)
        imsd_fhs=glob('%s/*.imsd' % expt_dh)
        for imsd_fh in imsd_fhs:  
            imsd=pd.read_csv(imsd_fh).set_index('lagt')
            imsd_flt,params_flt=flt_traj(imsd,flt_amplitude=True,out_fh=imsd_fh,**_cfg)
        from htsimaging.lib.fit_kin import plot_kin_all        
        plot_kin_all(expt_dh,imsd_fhs)
    else:
        logging.error('path does not exist %s' % expt_dh)

    
if __name__ == '__main__':
    main(sys.argv[1])