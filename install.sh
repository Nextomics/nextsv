#!/bin/bash

NEXTSV_ROOT=$PWD

mkdir -p bin

echo "#### installation of Sniffles ####"
cd $NEXTSV_ROOT/aligners_and_callers
sniffles_build_dir=$NEXTSV_ROOT/aligners_and_callers/Sniffles-1.0.5/build
if [ -d "$sniffles_build_dir" ]
then
    echo "cleaning failed files"
    rm -rf $sniffles_build_dir/*
else
    tar -xzf Sniffles-1.0.5.tar.gz
    cd Sniffles-1.0.5
    mkdir build 
fi

echo "building Sniffles"
cd $sniffles_build_dir
cmake ..
make
mv ../bin/sniffles-core-1.0.5/sniffles $NEXTSV_ROOT/bin/

echo "#### installation of PBHoney ####"
cd $NEXTSV_ROOT/aligners_and_callers
pip install --user argparse
pip install --user pysam
pip install --user networkx
pip install --user h5py
pip install --user PyIntervalTree

tar -xzf PBSuite_15.8.24.tgz
cd PBSuite_15.8.24
echo "
export SWEETPATH=$NEXTSV_ROOT/aligners_and_callers/PBSuite_15.8.24
export PYTHONPATH=\$PYTHONPATH:\$SWEETPATH
export PATH=\$PATH:\$SWEETPATH/bin
" > $NEXTSV_ROOT/setup-env.sh
cd $NEXTSV_ROOT

echo "#### installation of PBHoney ####"
mv $NEXTSV_ROOT/aligners_and_callers/blasr $NEXTSV_ROOT/bin/blasr

echo "
########################################################################

Please add the following environmental variables to your ~/.bashrc file:

export SWEETPATH=$NEXTSV_ROOT/aligners_and_callers/PBSuite_15.8.24
export PYTHONPATH=\$PYTHONPATH:\$SWEETPATH
export PATH=\$PATH:\$SWEETPATH/bin

########################################################################

"

echo "Installation finished"
