#!/bin/bash

NEXTSV_ROOT=$PWD

mkdir -p bin
mkdir -p src

echo "#### installation of samtools-1.3 and htslib ####"
cd $NEXTSV_ROOT/aligners_and_callers
tar -xjf samtools-1.3.1.tar.bz2
cd samtools-1.3.1/
./configure
make
cp samtools $NEXTSV_ROOT/bin/samtools1.3
mkdir -p $NEXTSV_ROOT/bamstat/lib/ $NEXTSV_ROOT/bamstat/include
cp htslib-1.3.1/libhts.a $NEXTSV_ROOT/bamstat/lib/libhts.a
cp htslib-1.3.1/htslib/* $NEXTSV_ROOT/bamstat/include/

echo "#### installation of bamstat ####"
cd $NEXTSV_ROOT/bamstat
make
cp bamstat $NEXTSV_ROOT/bin

echo "#### installation of Sniffles ####"
cd $NEXTSV_ROOT/aligners_and_callers
tar -xzf Sniffles-1.0.5.tar.gz
cd Sniffles-1.0.5
mkdir build 
cd build
cmake ..
make
mv ../bin/sniffles-core-1.0.5/sniffles $NEXTSV_ROOT/bin

echo "#### installation of ngmlr ####"
cd $NEXTSV_ROOT/aligners_and_callers
tar -xzf ngmlr-0.2.3.tar.gz
cd ngmlr-0.2.3
mkdir build
cd build
cmake ..
make
mv ../bin/ngmlr-0.2.3/ngmlr $NEXTSV_ROOT/bin


echo "#### installation of bwa ####"
cd $NEXTSV_ROOT/aligners_and_callers
tar -xzf bwa-0.7.15.tar.gz
cd bwa-0.7.15
make
mv bwa $NEXTSV_ROOT/bin


echo "#### installation of PBHoney ####"
cd $NEXTSV_ROOT/aligners_and_callers
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
cp $NEXTSV_ROOT/aligners_and_callers/blasr $NEXTSV_ROOT/bin/blasr

echo "
########################################################################

Please add the following environmental variables to your ~/.bashrc file:

export SWEETPATH=$NEXTSV_ROOT/aligners_and_callers/PBSuite_15.8.24
export PYTHONPATH=\$PYTHONPATH:\$SWEETPATH
export PATH=\$PATH:\$SWEETPATH/bin

########################################################################

"

echo "Installation finished"
