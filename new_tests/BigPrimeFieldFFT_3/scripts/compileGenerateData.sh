#!/bin/bash

# PROJECT_NAME="BigPrimeFieldFFT_3"
# PROJECT="$CUMODP_HOME/new_tests/$PROJECT_NAME"

. ./projectProperties.conf

# echo $CUMODP_HOME
# echo $PROJECT
g++ -DLINUXINTEL64 -I $CUMODP_HOME/include -O2 -std=c++11 -w $PROJECT/src/generateData.cpp -o $PROJECT/bin/generateData.bin 

# nvcc -DLINUXINTEL64 --cudart shared -w --ptxas-options=-v -O2 $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin 2> $PROJECT/test/compileMessage

# nvcc -DLINUXINTEL64 --cudart shared -w -O2 -Xptxas -dlcm=ca $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin 2> $PROJECT/test/compileMessage

# nvcc -DLINUXINTEL64 --compiler-options '-fPIC -Wall -Wno-unused-but-set-variable -Wno-comment' --cudart shared -w -O2 -Xptxas -dlcm=cg $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin 2> $PROJECT/test/compileMessage

# nvcc -DLINUXINTEL64 --cudart shared -w -O2 -Xptxas -dlcm=ca $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin 2> $PROJECT/test/compileMessage

# nvcc -DLINUXINTEL64 --cudart shared -w -O2 -Xptxas -dlcm=ca $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin 2> $PROJECT/test/compileMessage

COMPUTE_CAPABILITY=compute_20
GPU_SM=sm_20

COMPUTE_CAPABILITY=compute_30
GPU_SM=sm_30

GENCODE_ARCH_COMPUTE="-gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30"
# -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30
# -gencode $GENCODE_ARCH_COMPUTE
# printf "new compile script for $PROJECT_NAME"

# nvcc -DLINUXINTEL64 --cudart shared --compiler-options '-fPIC -Wall -Wno-unused-but-set-variable -Wno-comment' --gpu-architecture=$COMPUTE_CAPABILITY --gpu-code=$GPU_SM -w -O2 -Xptxas -dlcm=ca $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin

# nvcc -DLINUXINTEL64 --cudart shared --compiler-options '-fPIC -Wall -Wno-unused-but-set-variable -Wno-comment' $GENCODE_ARCH_COMPUTE -w -O2 -Xptxas -dlcm=ca $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin

#nvcc -DLINUXINTEL64 -w -O2 -Xptxas -dlcm=ca $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin
# nvcc -DLINUXINTEL64 $PROJECT/src/bigArithmetic_test.cu -I $CUMODP_HOME/include -L $CUMODP_HOME/src -lcumodp -o $PROJECT/bin/bigArithmetic_test.bin
# nvcc -DLINUXINTEL64 testCTPowers.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o ctpowers
# 2> $PROJECT/test/compileMessage

# nvcc -DLINUXINTEL64 -I/$CUMODP_HOME/src -I/$INC -std=c++11 -w --ptxas-options=-v -gencode arch=compute_30,code=sm_30 ../src/bigArithmetic_test.cu -o ../bin/bigArithmetic_test.bin 2> ../test/compileMessage

# less $PROJECT/test/compileMessage
