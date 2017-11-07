#!/bin/bash

NEXTSV_ROOT=$PWD

echo "#### installation of PBHoney ####"
cd $NEXTSV_ROOT/aligners_and_callers
tar -xzf PBSuite_15.8.24.tgz
cd PBSuite_15.8.24
echo "
export SWEETPATH=$NEXTSV_ROOT/aligners_and_callers/PBSuite_15.8.24
export PYTHONPATH=\$PYTHONPATH:\$SWEETPATH
export PATH=\$PATH:\$SWEETPATH/bin
" > $NEXTSV_ROOT/setup-env.sh
cd $NEXTSV_ROOT

echo "
########################################################################

Please add the following environmental variables to your ~/.bashrc file:

export SWEETPATH=$NEXTSV_ROOT/aligners_and_callers/PBSuite_15.8.24
export PYTHONPATH=\$PYTHONPATH:\$SWEETPATH
export PATH=\$PATH:\$SWEETPATH/bin

########################################################################

"

