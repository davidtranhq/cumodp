##################################################
# Compiler variables and flags
##################################################
#For Linux
ARCH=LINUXINTEL64

##################################################
# Generating code for a specific sm and compute 
# capability
# GENCODE_ARCH_COMPUTE=-gencode arch=compute_20,code=sm_20
# GENCODE_ARCH_COMPUTE=-gencode arch=compute_30,code=sm_30
# GENCODE_ARCH_COMPUTE=-gencode arch=compute_50,code=sm_50

##################################################
ifeq ($(ARCH), LINUXINTEL64)
	CC=gcc
	GXX=g++ -std=c++14
	NVCC=nvcc -std=c++14 --extended-lambda -ccbin $(GXX) -O2 -lineinfo
	
	CFLAGS= -c -O2 -DLINUXINTEL64 -I$(CUMODP_HOME)/include -w
	
	CUFLAGS=--compiler-options '-fPIC -Wall -Wno-unused-but-set-variable -Wno-comment' -Xptxas -dlcm=cg --cudart shared $(GENCODE_ARCH_COMPUTE)

	LIBFLAGS= -shared -o 
	LDFLAGS=-lm -lcudart -lgmp -lgmpxx
	CUMODP=cumodp
endif

ifeq ($(ARCH), WINDOWS)
	CC=cl
	LINK=link
	NVCC=nvcc
	CFLAGS=
	CUFLAGS=
	LIBFLAGS=
	LDFLAGS=
	CUMODP=cumodp
endif

# Flags for the C/C++ compiler
# Useful flags:
# 	-O2   --  level of optimization
# 	-m64  --  needed to get 64-bit longs
# 	-g    --  debugging
#
# Flags for the CUDA compiler
# Useful flags:
#   --ptxas-options -v     

##################################################

 SOURCES=fft_aux.cu\
	list_stockham.cu\
	list_inv_stockham.cu\
	stockham_mont.cu\
	stockham.cu\
 	ct_fft_mont.cu\
	rdr_poly.cu\
	subres_p.cu\
	subres.cu\
	scube.cu\
	modpn.cu\
	condensationReal.cu\
	condensationFiniteField.cu\
	naive_poly_mul.cu\
	opt_plain_mul.cu\
	fastPolyEvaluation.cu\
	fastPolyInterpolation.cu\
	taylor_shift_kernel.cu\
	taylor_shift_cpu.cu\
	taylor_shift.cu\
	cudaDiv.cu\
	cudaGcd.cu\
	list_pointwise_mul.cu\
	taylor_shift_fft.cu\
	ts_chinese_remainder.cu\
	taylor_comb.cu\
	ts_VCA_count.cu\
	ts_VCA_all.cu\
	new_prefix_mul.cu\
	new_taylor_shift.cu\
	new_pkc.cu\
	new_mixed_radix_rep.cu\
	mixed_radix_reverse.cu\
	BigPrimeFieldFFT_3/bigArithmetic_addition_P3.cu\
	BigPrimeFieldFFT_3/bigArithmetic_subtraction_P3.cu\
	BigPrimeFieldFFT_3/bigArithmetic_cyclicShift_P3.cu\
	BigPrimeFieldFFT_3/bigArithmetic_multiplication_P3.cu\
	BigPrimeFieldFFT_3/bigArithmetic_fft_P3.cu\
	BigPrimeFieldFFT_3/bigArithmetic_aux_P3.cu\
	BigPrimeFieldFFT_3/bigFFT_reduction_P3.cu\
	BigPrimeFieldFFT_4/bigArithmetic_addition_P4.cu\
	BigPrimeFieldFFT_4/bigArithmetic_subtraction_P4.cu\
	BigPrimeFieldFFT_4/bigArithmetic_cyclicShift_P4.cu\
	BigPrimeFieldFFT_4/bigArithmetic_multiplication_P4.cu\
	BigPrimeFieldFFT_4/bigArithmetic_fft_P4.cu\
	BigPrimeFieldFFT_4/bigArithmetic_aux_P4.cu\
	BigPrimeFieldFFT_4/bigFFT_reduction_P4.cu\
	two_convolution_poly_mul.cu\
	two_convolution_constants_check.cpp\
	
##################################################
#  Building 
##################################################

ifeq ($(ARCH), LINUXINTEL64) 

OBJECTS := $(SOURCES:.cu=.o)
OBJECTS := $(OBJECTS:.cpp=.o)

##################################################
all: so

##################################################
so: $(OBJECTS)
	$(NVCC) $(LIBFLAGS) lib$(CUMODP).so $(OBJECTS) $(LDFLAGS)

##################################################
%.o : %.cu
	$(NVCC) $(CFLAGS) $(CUFLAGS) $*.cu -o $*.o

##################################################

%.o : %.cpp
	$(GXX) -fPIC --std=c++20 $(CFLAGS) $*.cpp -o $*.o

else

echo "To build for Windows platforms is not supported!"

endif

clean:	
	rm -f *.o *~ *.linkinfo *.out *.ptx *.log *.prof *.obj
	rm -f lib$(CUMODP).*
	rm -f main exe
	rm -f BigPrimeFieldFFT_3/*.o
	rm -f BigPrimeFieldFFT_4/*.o

####################END OF THE FILE#########################"