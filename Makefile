PYTHON ?= python
CYTHON ?= cython
PYTEST ?= pytest
CTAGS ?= ctags

all: clean inplace test

inplace:
	$(PYTHON) setup.py build_ext -i

install:
	$(PYTHON) setup.py install

test: test-code

test-code: inplace
	$(PYTEST) --showlocals -v circhic --durations=20

test-coverage:
	rm -rf coverage .coverage
	$(PYTEST) circhic --showlocals -v --cov=circhic

clean-ctags:
	rm -f tags

clean: clean-ctags
	$(PYTHON) setup.py clean
	rm -rf dist
	rm -rf build

trailing-spaces:
	find circhic -name "*.py" -exec perl -pi -e 's/[ \t]*$$//' {} \;

