# `htsimaging`

[![build status](http://img.shields.io/travis/rraadd88/htsimaging/master.svg?style=flat)](https://travis-ci.org/rraadd88/htsimaging) [![PyPI version](https://badge.fury.io/py/htsimaging.svg)](https://pypi.python.org/pypi/htsimaging)

Python package for analysis of big bunch of images obtained from high throughput imaging.

Table of Contents
-----------------

1. [Installation](#installation)
2. [Configuration](#configuration-file)
3. [Input format](#input-format)
4. [Output format](#output-format)
5. [API](#api)

Installation
------------
Basic requirements: [`Anaconda package manager`](https://www.anaconda.com/download/#linux). See [requirements.md](https://github.com/rraadd88/test_beditor/blob/master/requirements.md) for set of bash commands that would install it.
0. Installing Ana(mini)conda:

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc
```

1.  Once all the requirements are satisfied, create a python 3.6 virtual environment.

``` {.sourceCode .text}
wget https://raw.githubusercontent.com/rraadd88/htsimaging/master/environment.yml
conda env create -f environment.yml
```

2.  Activate the virtual environment.

``` {.sourceCode .text}
source activate htsimaging
```

3. Install `htsimaging` python package from github.

``` {.sourceCode .text}
git clone https://github.com/rraadd88/htsimaging.git
cd htsimaging;pip install -e .
```

"Pull" updates
------------
`htsimaging` is constantly under development. In order to update local repo to the latest version, first change direcotry to (`cd`) `htsimaging` folder, and run following command:

```
git pull origin master
```

Often you may need to update the dependencies also. Dependencies such as (`rohan`).

Dependencies 
-------------

1. `yeastspotter`: for segmentation of bright field images

``` {.sourceCode .bash}
# go to htsimaging repo (`cd htsimaging`)
mkdir deps; cd deps
git clone https://github.com/rraadd88/yeast_segmentation.git
```
