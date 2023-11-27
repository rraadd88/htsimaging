#!/usr/bin/env python
"""
## Development

### release new version

    git commit -am "version bump";git push origin master
    python setup.py --version
    git tag -a v$(python setup.py --version) -m "upgrade";git push --tags

"""

import sys
if (sys.version_info[0]) != (3):
    raise RuntimeError('Python 3 required ')

import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()
    
requirements = {
'base':[
    'roux>=0.0.9', # helper functions
    'scikit-image', # image processing.
    'argh', # for command-line convinience
    'seaborn',
    'fastparquet',
    ],
'spt':[
    'pims==0.4.1', 
    'trackpy', # for tracking the particles
    'numba', # required by trackpy for faster processing
],
## development and maintenance
'dev':[
    'pytest',
    'jupyter','ipywidgets','ipykernel',
    # 'sphinx','recommonmark',
    'black',
    'coveralls == 3.*',
    'flake8',
    'isort',
    'pytest-cov == 2.*',
    'testbook',
    'papermill',
    'lazydocs', 'regex', ## docs
],
}

extras_require={k:l for k,l in requirements.items() if not k=='base'}
## all: extra except dev
extras_require['all']=[l for k,l in extras_require.items() if not k=='dev']
### flatten
extras_require['all']=[s for l in extras_require['all'] for s in l]
### unique
extras_require['all']=list(set(extras_require['all']))

# main setup command
setuptools.setup(
    name='htsimaging',
    version='1.0.5',
    description='High-Throughput Single-cell Imaging analysis.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='http://github.com/rraadd88/htsimaging',
    author='rraadd88',
    author_email='rohanadandage@gmail.com',
    license='General Public License v. 3',
    packages=setuptools.find_packages('.',exclude=['test','tests', 'unit','deps','data','examples']),
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=requirements['base'],
    extras_require=extras_require,
    # entry_points={
    # 'console_scripts': ['htsimaging = htsimaging.run:parser.dispatch',],
    # },
    python_requires='>=3.7, <4',
)