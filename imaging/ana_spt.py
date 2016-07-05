#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    

import sys
import logging
import pandas as pd
import numpy as np
import trackpy as tp
logging.basicConfig(format='[%(asctime)s] %(levelname)s\tfrom %(filename)s in %(funcName)s(..): %(message)s',level=logging.DEBUG) # 
from imaging.lib import spt #.expt2plots,spt.fit_power,spt.fit_line

def main(expt_dh):
    expt_dh=expt_dh+"/"
    expt_data_chem_trial=spt.expt2plots(spt.expt_dh2expt_info(expt_dh),expt_dh)

    expt_data=pd.read_csv(expt_dh+"expt.emsd")
    if expt_data.index.name!='lagt':
        expt_data=expt_data.set_index('lagt',drop=True)

    if 'Unnamed: 0' in expt_data.columns.tolist():
        del expt_data['Unnamed: 0']

    parameters=pd.DataFrame(index=expt_data.columns.tolist(), columns=["power law exponent","power law constant","rsquared of power","slope","y intercept","rsquared of line"])

    for col in expt_data.columns.tolist():
    #     print col
        values=tp.utils.fit_powerlaw(expt_data.loc[:,col],plot=False)
        values=values.reset_index()
        parameters.loc[col,"n"]= values.loc[0,"n"]
        parameters.loc[col,"A"]= values.loc[0,"A"]
        
        parameters.loc[col,"power law exponent"],\
        parameters.loc[col,"power law constant"], \
        parameters.loc[col,"rsquared of power"] = \
        spt.fit_power(np.array(list(expt_data.index))[:60],np.array(list(expt_data.loc[:,col]))[:60])

        parameters.loc[col,"slope"],\
        parameters.loc[col,"y intercept"], \
        parameters.loc[col,"rsquared of line"] = \
        spt.fit_line(np.log10(np.array(list(expt_data.index))[:60]),np.log10(np.array(list(expt_data.loc[:,col]))[:60]))
    parameters.to_csv("%s/parameters" % expt_dh)
    

    
if __name__ == '__main__':
    main(sys.argv[1])