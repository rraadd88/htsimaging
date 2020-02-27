#!/usr/bin/env python

# Copyright 2019, Rohan Dandage <rohanadandage@gmail.com,rraadd_8@hotmail.com>
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


if (sys.version_info[0], sys.version_info[1]) != (3, 6):
    raise RuntimeError('Python 3.6 required ')
    
# update requirements.txt from environment.yml
# with open('requirements.txt','w') as f:
#     f.write('\n'.join(dropna(yaml.load(open('environment.yml','r'))['dependencies'][2]['pip'])))

with open('requirements.txt') as f:
    required = f.read().splitlines()

    
# main setup command
setup(
name='htsimaging',
author='Rohan Dandage',
author_email='rraadd_8@hotmail.com',
version='1.0.2',
url='https://github.com/rraadd88/..',
download_url='https://github.com/rraadd88/htsimaging',
description='For analysis of microscopy images',
long_description='https://github.com/rraadd88/htsimaging/README.md',
license='General Public License v. 3',
# install_requires=required,
platforms='Tested on Ubuntu 16.04',
keywords=['microscopy','imaging','trackpy'],
packages=find_packages(),    
entry_points={
    'console_scripts': ['endocytosis = htsimaging.endocytosis:parser.dispatch',],
    },
    
)

#for git 
# setup version change -> commit
# git tag -a v$(python setup.py --version) -m 'still beta'
# for pypi
# pandoc -o README.rst README.md
