# htsimaging

[![build status](
  http://img.shields.io/travis/rraadd88/htsimaging/master.svg?style=flat)](
 https://travis-ci.org/rraadd88/htsimaging)

[![PyPI version](https://badge.fury.io/py/htsimaging.svg)](https://pypi.python.org/pypi/htsimaging)

Python package for analysis of big bunch of images obtained from high throughput imaging.

## Installation

	cd htsimaging_master
	pip install .

## Contents of the package

	htsimaging/  
	|-- bleach_chase2kins.py (bleach-chase images to kinetics)  
	|-- bleach_chase2vid.py (bleach-chase images to video)  
	|-- trackinginfo2emsd.py (single particle tracking analysis)  
	|-- __init__.py  
	`-- lib (under-the-hood libraries)  
	    |-- bleach_chase.py   
	    |-- diffusion.py  
	    |-- __init__.py  
	    |-- io_data_files.py  
	    `-- io_nd_files.py  
