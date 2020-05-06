#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.
  
# License: 3-clause BSD
set -e

pip install --upgrade pip pytest pytest-cov
pip install --upgrade numpy scipy matplotlib

# For the documentation
pip install --upgrade sphinx sphinx-gallery numpydoc

python setup.py develop

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "import matplotlib; print('matplotlib %s' % matplotlib.__version__)"
