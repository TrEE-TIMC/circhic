#!/bin/bash
# This script is meant to be called by the "script" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

set -e
set -v

run_tests() {
    TEST_CMD="pytest --showlocals --pyargs"

    # Get into a temp directory to run test from the installed scikit learn
    # and
    # check if we do not leave artifacts
    mkdir -p $TEST_DIR
    pushd $TEST_DIR

    if [[ "$COVERAGE" == "true" ]]; then
        TEST_CMD="$TEST_CMD --cov=circhic --cov-report=xml"
    fi
    $TEST_CMD circhic
    popd
}

compile_documentation() {
    pushd docs/
    make html SPHINXOPTS="-W"
    popd
}

run_tests
compile_documentation
