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
from htsimaging.lib.spt import expt2plots,expt_dh2expt_info,msd2params
from roux.lib.io_dfs import set_index

import argparse
def main():
    logging.info("start")
    parser = argparse.ArgumentParser(description='data2ml')
    parser.add_argument("prj_dh", help="path to project directory", 
                        action="store", default=False)    
    parser.add_argument("--test", help="Debug mode on", dest="test", 
                        action='store_true', default=False)    
    parser.add_argument("--force", help="Debug mode on", dest="force", 
                        action='store_true', default=False)    
    args = parser.parse_args()
    pipeline(args.prj_dh,
             test=args.test,
             force=args.force)
    logging.shutdown()

def pipeline(expt_dh,test=False,force=False):
    if not exists(expt_dh):    
        if expt_dh.count('/')<2:
            expt_dh='%s/../test_dataset/%s' % (dirname(abspath(__file__)),expt_dh) 
            print expt_dh
    if exists(expt_dh):    
        _cfg_fh='%s/_cfg.json' % expt_dh
        if exists(_cfg_fh):
            if _cfg_fh.endswith('json'):
                with open(_cfg_fh, 'r') as f:
                    import json
                    _cfg = json.load(f)
        elif exists('%s/_cfg.json' % expt_dh):
            _cfg_fh='%s/_cfg.json' % expt_dh
            if _cfg_fh.endswith('csv'):
                _cfg=pd.read_csv(_cfg_fh).set_index('var')
                _cfg=_cfg['val'].to_dict()
                # print _cfg #debug
        else:
            _cfg={}
            
        print _cfg
        expt2plots(expt_dh2expt_info(expt_dh),expt_dh,_cfg=_cfg,
                  test=test,force=force)
        emsd_fhs=glob('%s/*.emsd' % expt_dh)
        for fh in emsd_fhs:  
            emsd=pd.read_csv(fh)
            emsd=set_index(emsd,'lag time [s]')
            msd2params(emsd,out_fh=fh,**_cfg)
            # emsd_flt,params_flt=flt_traj(imsd,flt_amplitude=True,out_fh=imsd_fh,**_cfg)
        # from htsimaging.lib.fit_kin import plot_kin_all        
        # plot_kin_all(expt_dh,imsd_fhs)
    else:
        logging.error('path does not exist %s' % expt_dh)

    
if __name__ == '__main__':
    main()
