#!/bin/bash

# /* opAdd:=0
# * opMultiplication:=1
# * opSubtraction:=2
# * opCyclicShift:=3
# * opFFT_16:=4
# * opFFT_256:=5
# * opFFT_4k:=6
# * opFFT_16k:=7
# * opFFT_64k:=8
# */

# PROJECT_NAME="BigPrimeFieldFFT_3"
# PROJECT="$CUMODP_HOME/new_tests/$PROJECT_NAME"
. ./projectProperties.conf

rm -f $CUMODP_HOME/src/bigArithmeticKernels.cu
rm -f $PROJECT/src/bigArithmeticKernels.cuh
rm -f $CUMODP_HOME/include/bigArithmeticKernels.h
rm -f $CUMODP_HOME/include/bigPrimeField.h
rm -f $CUMODP_HOME/include/bigArithmetic_addition.h
rm -f $CUMODP_HOME/include/bigArithmetic_cyclicShift.h
rm -f $CUMODP_HOME/include/bigArithmetic_fft.h
rm -f $CUMODP_HOME/include/bigArithmetic_multiplication.h
rm -f $CUMODP_HOME/include/bigArithmetic_subtraction.h



# include_list=()
#for compiling without cumodp library (-lcumodp)
ln -sf $CUMODP_HOME/src/$PROJECT_NAME/bigArithmeticKernels_P3.cu $CUMODP_HOME/src/bigArithmeticKernels_P3.cu

ln -sf $CUMODP_HOME/src/bigArithmeticKernels_P3.cu $PROJECT/src/bigArithmeticKernels.cuh

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigArithmeticKernels.h $CUMODP_HOME/include/bigArithmeticKernels.h

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigPrimeField.h $CUMODP_HOME/include/bigPrimeField.h

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigArithmetic_addition.h $CUMODP_HOME/include/bigArithmetic_addition.h

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigArithmetic_cyclicShift.h $CUMODP_HOME/include/bigArithmetic_cyclicShift.h

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigArithmetic_fft.h $CUMODP_HOME/include/bigArithmetic_fft.h

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigArithmetic_multiplication.h $CUMODP_HOME/include/bigArithmetic_multiplication.h

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigArithmetic_subtraction.h $CUMODP_HOME/include/bigArithmetic_subtraction.h

# ln -sf $CUMODP_HOME/include/$PROJECT_NAME/bigFFT_reduction.h $CUMODP_HOME/include/bigFFT_reduction.h

cd $PROJECT/scripts/mapleScripts/
	ln -sf maple_operation_add.mpl maple_operation_0.mpl
	ln -sf maple_operation_mult.mpl maple_operation_1.mpl
	ln -sf maple_operation_mult_revised.mpl maple_operation_1_revised.mpl
	ln -sf maple_operation_sub.mpl maple_operation_2.mpl
	ln -sf maple_operation_cyclicShift.mpl maple_operation_3.mpl
	ln -sf maple_operation_fft_16.mpl maple_operation_4.mpl
	ln -sf maple_operation_fft_256.mpl maple_operation_5.mpl
	ln -sf maple_operation_fft_4k.mpl maple_operation_6.mpl
	ln -sf maple_operation_fft_16k.mpl maple_operation_7.mpl
	ln -sf maple_operation_fft_64k.mpl maple_operation_8.mpl