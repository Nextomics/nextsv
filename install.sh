#!/bin/bash
cd longreadqc
sh build.sh
cp longreadqc ../bin/
cd ../pigz
make
cp pigz ../bin/
cp unpigz ../bin/
