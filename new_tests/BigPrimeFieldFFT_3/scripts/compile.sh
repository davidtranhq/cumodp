#!/bin/bash

# PROJECT_NAME="BigPrimeFieldFFT_3"
# PROJECT="$CUMODP_HOME/new_tests/$PROJECT_NAME"

#################################################
. ./projectProperties.conf

COMPUTE_CAPABILITY=compute_20
GPU_SM=sm_20

COMPUTE_CAPABILITY=compute_30
GPU_SM=sm_30

GENCODE_ARCH_COMPUTE="-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30"
#################################################

g++ -DLINUXINTEL64 -I $CUMODP_HOME/include -O2 -std=c++11 -w $PROJECT/src/generateData.cpp -o $PROJECT/bin/generateData.bin 

#################################################
nvcc -DLINUXINTEL64 --cudart shared --compiler-options '-fPIC -Wall -Wno-unused-but-set-variable -Wno-comment' $GENCODE_ARCH_COMPUTE -w -O2 -Xptxas -dlcm=ca $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin
#################################################