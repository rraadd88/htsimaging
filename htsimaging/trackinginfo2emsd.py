#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
import logging
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # 
from htsimaging.lib.io_data_files import expt_dh2expt_info
from htsimaging.lib.diffusion import expt2plots

def main(expt_dh):
    expt_data_chem_trial=expt2plots(expt_dh2expt_info(expt_dh),expt_dh)

if __name__ == '__main__':
    main(sys.argv[1])
