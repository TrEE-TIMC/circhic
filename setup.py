#! /usr/bin/env python
#
# License: 3-clause BSD

import os
import sys
from setuptools import setup
import circhic

DISTNAME = 'circhic'
DESCRIPTION = 'Visualizing HiC data of circular chromosomes'
with open('README.md') as f:
    LONG_DESCRIPTION = f.read()
MAINTAINER = 'Ivan Junier'
MAINTAINER_EMAIL = 'ivan.junier@univ-grenoble-alpes.fr'
URL = 'https://github.com/tree-timc/circHiC'
DOWNLOAD_URL = 'https://pypi.org/project/circHiC/#files'
LICENSE = 'new BSD'
PROJECT_URLS = {
    'Bug Tracker': 'https://github.com/tree-timc/circHiC/issues',
    'Documentation': 'https://tree-timc.github.io/circhic/',
    'Source Code': 'https://github.com/tree-timc/circHiC'
}

# Let's avoid having to change version number here as well as in the
# __init__.py
VERSION = circhic.__version__

SCIPY_MIN_VERSION = '0.19.1'
NUMPY_MIN_VERSION = '1.16'
MATPLOTLIB_MIN_VERSION = '3.0.0'

setup(
    name=DISTNAME,
    version=VERSION,
    author=MAINTAINER,
    author_email=MAINTAINER_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
    include_package_data=True,
)
