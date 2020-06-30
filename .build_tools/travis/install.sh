#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.
  
# License: 3-clause BSD
set -e

pip install --upgrade pip pytest pytest-cov

if [[ $NUMPY_VERSION != "*" ]]; then
    pip install --upgrade \
        numpy==$NUMPY_VERSION
else
    pip install numpy --upgrade
fi

if [[ $SPHINX_VERSION != "*" ]]; then
    pip install --upgrade \
        sphinx==$SPHINX_VERSION
else
    pip install sphinx --upgrade
fi



if [[ $MATPLOTLIB_VERSION != "*" ]]; then
    pip install matplotlib==$MATPLOTLIB_VERSION
else
    pip install matplotlib --upgrade
fi

if [[ $SCIPY_VERSION != "*" ]]; then
    pip install --upgrade scipy==$SCIPY_VERSION
else
    pip install scipy --upgrade
fi

pip install --upgrade pandas iced

# For the documentation
pip install --upgrade sphinx-gallery numpydoc
pip install --upgrade pillow

python setup.py develop

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "import matplotlib; print('matplotlib %s' % matplotlib.__version__)"
