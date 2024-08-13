#!/bin/bash

cmake --preset release -DPASTA_BIT_VECTOR_BUILD_BENCHMARKS=On
cmake --build --preset release

[ -f result.txt ] && rm result.txt 

for FILLRATE in 10 20 30 40 50 60 70 80 90
do
    for BITVECTORSIZE in 10M 20M 30M
    do
        ./build/bit_vector_benchmark -b $BITVECTORSIZE -q 10M -f $FILLRATE | tee -a result.txt
    done
    
done

sqlplot-tools plots.tex
