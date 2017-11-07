#!/bin/bash

NEXTSV_ROOT=$PWD

mkdir -p bin

echo "#### installation of Sniffles ####"
cd $NEXTSV_ROOT/aligners_and_callers
sniffles_build_dir=$NEXTSV_ROOT/aligners_and_callers/Sniffles-1.0.5/build/
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

if [ -f "$NEXTSV_ROOT/bin/sniffles" ]
then
	echo "Installation of Sniffles finished"
else
	echo "Installation of Sniffles failed"
fi

