#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rohan@igib.in,rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    


"""
========
setup.py
========

installs package

USAGE :
python setup.py install

Or for local installation:

python setup.py install --prefix=/your/local/dir

"""

import sys
try:
    from setuptools import setup, find_packages, Extension 
except ImportError:
    from distutils.core import setup, find_packages, Extension


if (sys.version_info[0], sys.version_info[1]) != (2, 7):
    raise RuntimeError('Python 2.7 required ')
               
# main setup command
setup(
name='htsimaging',
author='Rohan Dandage',
author_email='rraadd_8@hotmail.com,rohan@igib.in',
version='1',
url='https://github.com/rraadd88/..',
download_url='https://github.com/rraadd88/htsimaging/..',
description='..',
long_description='https://github.com/rraadd88/../README.md',
license='General Public License v. 3',
install_requires=[  'pandas >= 0.18.0',
                    'scipy >= 0.17.0',
                    'scikit_learn >= 0.17.1',
                    'matplotlib >= 1.5.1',
                    'seaborn >= 0.6.0',
                    'pims == 0.3.3',
                    # 'opencv',
                    'scikit-image == 0.11.3',
                    'nd2reader',
                    'trackpy',
                    ],
platforms='Tested on Linux (Ubuntu 12.04)',
keywords=['1','1','1'],
packages=find_packages(),
)