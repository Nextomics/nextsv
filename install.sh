#!/bin/bash

NEXTSV_ROOT=$PWD

mkdir -p bin
mkdir -p src

echo "#### installation of samtools-1.3 and htslib ####"
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar -xjf samtools-1.3.1.tar.bz2
mv samtools-1.3.1.tar.bz2 src
cd samtools-1.3.1/
./configure
make
cp samtools $NEXTSV_ROOT/bin/samtools1.3
mkdir -p ../bamstat/lib/ ../bamstat/include
cp htslib-1.3.1/libhts.a ../bamstat/lib/libhts.a
cp htslib-1.3.1/htslib/* ../bamstat/include/

echo "#### installation of bamstat ####"
cd $NEXTSV_ROOT/bamstat
make
cp bamstat ../bin
cd $NEXTSV_ROOT

echo "#### installation of Sniffles ####"
wget https://github.com/fritzsedlazeck/Sniffles/archive/v1.0.5.tar.gz
mv v1.0.5.tar.gz Sniffles-1.0.5.tar.gz
tar -xzf Sniffles-1.0.5.tar.gz
mv Sniffles-1.0.5.tar.gz src/
cd Sniffles-1.0.5
mkdir build 
cd build
cmake ..
make
mv ../bin/sniffles-core-1.0.5/sniffles $NEXTSV_ROOT/bin
cd $NEXTSV_ROOT

echo "#### installation of ngmlr ####"
wget https://github.com/philres/ngmlr/archive/v0.2.3.tar.gz
mv v0.2.3.tar.gz ngmlr-0.2.3.tar.gz
tar -xzf ngmlr-0.2.3.tar.gz
mv ngmlr-0.2.3.tar.gz src
cd ngmlr-0.2.3
mkdir build
cd build
cmake ..
make
mv ../bin/ngmlr-0.2.3/ngmlr $NEXTSV_ROOT/bin
cd $NEXTSV_ROOT


echo "#### installation of bwa ####"
wget https://github.com/lh3/bwa/archive/v0.7.15.tar.gz
mv v0.7.15.tar.gz bwa-0.7.15.tar.gz
tar -xzf bwa-0.7.15.tar.gz
mv bwa-0.7.15.tar.gz src/
cd bwa-0.7.15
make
mv bwa $NEXTSV_ROOT/bin
cd $NEXTSV_ROOT


echo "#### installation of PBHoney ####"
pip install --user pysam
pip install --user networkx
pip install --user h5py
pip install --user PyIntervalTree

wget https://pilotfiber.dl.sourceforge.net/project/pb-jelly/PBSuite_15.8.24.tgz
tar -xzf PBSuite_15.8.24.tgz
mv PBSuite_15.8.24.tgz src
cd PBSuite_15.8.24
echo "
export SWEETPATH=$NEXTSV_ROOT/PBSuite_15.8.24
export PYTHONPATH=\$PYTHONPATH:\$SWEETPATH
export PATH=\$PATH:\$SWEETPATH/bin
" > $NEXTSV_ROOT/setup-env.sh
cd $NEXTSV_ROOT


echo "#### installation of blasr ####"
git clone https://github.com/PacificBiosciences/pitchfork.git
cd pitchfork
echo PREFIX=$NEXTSV_ROOT/blasr > settings.mk
make init
make blasr
mv $NEXTSV_ROOT/blasr/bin/blasr $NEXTSV_ROOT/bin
echo "export LD_LIBRARY_PATH=$NEXTSV_ROOT/blasr/lib:\$LD_LIBRARY_PATH" >> $NEXTSV_ROOT/setup-env.sh

echo "
########################################################################

Please add the following environmental variables to your ~/.bashrc file:

export SWEETPATH=$NEXTSV_ROOT/PBSuite_15.8.24
export PYTHONPATH=\$PYTHONPATH:\$SWEETPATH
export PATH=\$PATH:\$SWEETPATH/bin
export LD_LIBRARY_PATH=$NEXTSV_ROOT/blasr/lib:\$LD_LIBRARY_PATH

########################################################################

"

echo "Installation finished"
