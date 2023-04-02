#!/bin/bash
cd longreadqc && make && cp longreadqc ../bin/ && make clean 
cd ..
cd pigz && make && cp pigz ../bin/ && cp unpigz ../bin/ && make clean 
