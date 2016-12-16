#!/bin/bash
mkdir -p bin
mkdir -p src
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar -xjf samtools-1.3.1.tar.bz2
mv samtools-1.3.1.tar.bz2 src
cd samtools-1.3.1/
./configure
make
cp samtools ../bin/samtools1.3
mkdir -p ../bamstat/lib/ ../bamstat/include
cp htslib-1.3.1/libhts.a ../bamstat/lib/libhts.a
cp htslib-1.3.1/htslib/* ../bamstat/include/
cd ../bamstat
make
cp bamstat ../bin
