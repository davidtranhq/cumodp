## cumodp/new_tests/FFTTest/Makefile

LDFLAGS=-L/usr/local/32-bit/lib
all:
	nvcc -DLINUXINTEL64 bigAdditionWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigAdditionWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigAdditionWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigAdditionWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalAdditionWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalAdditionWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalAdditionWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalAdditionWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalDoubleMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalDoubleMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalDoubleMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalDoubleMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalMontMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalMontMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalMontMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalMontMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 unsignedMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o unsignedMultiplicationWithSM -arch=sm_20
check: clean
	nvcc -DLINUXINTEL64 bigAdditionWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigAdditionWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigAdditionWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigAdditionWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalAdditionWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalAdditionWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalAdditionWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalAdditionWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalDoubleMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalDoubleMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalDoubleMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalDoubleMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalMontMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalMontMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalMontMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalMontMultiplicationWithoutSM -arch=sm_20
  nvcc -DLINUXINTEL64 unsignedMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o unsignedMultiplicationWithSM -arch=sm_20
test: dat
	nvcc -DLINUXINTEL64 bigAdditionWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigAdditionWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigAdditionWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigAdditionWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalAdditionWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalAdditionWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalAdditionWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalAdditionWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 bigMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o bigMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalDoubleMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalDoubleMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalDoubleMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalDoubleMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalMontMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalMontMultiplicationWithSM -arch=sm_20
	nvcc -DLINUXINTEL64 normalMontMultiplicationWithoutSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o normalMontMultiplicationWithoutSM -arch=sm_20
	nvcc -DLINUXINTEL64 unsignedMultiplicationWithSM.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o unsignedMultiplicationWithSM -arch=sm_20
	./bigAdditionWithSM
	./bigAdditionWithoutSM
	./normalAdditionWithSM
	./normalAdditionWithoutSM
	./bigMultiplicationWithSM
	./bigMultiplicationWithoutSM
	./normalDoubleMultiplicationWithSM
	./normalDoubleMultiplicationWithoutSM
	./normalMontMultiplicationWithSM
	./normalMontMultiplicationWithoutSM
	./unsignedMultiplicationWithSM

clean:
	rm -f *.o *~ *.linkinfo a.out *.log *.dat bigAdditionWithSM bigAdditionWithoutSM normalAdditionWithSM normalAdditionWithoutSM bigMultiplicationWithSM bigMultiplicationWithoutSM normalDoubleMultiplicationWithSM normalDoubleMultiplicationWithoutSM normalMontMultiplicationWithSM normalMontMultiplicationWithoutSM unsignedMultiplicationWithSM

dat:
	rm -f W_4096.dat W_16777216.dat FFT_Data_16.dat FFT_4096_input.dat FFT_256_input.dat  FFT_1048576_input.dat addition_benchmark_data.dat W_1048576.dat FFT_Data_Block.dat FFT_65536_input.dat FFT_16777216_input.dat small_number_benchmark.dat big_number_benchmark.dat W_65536.dat
	ln -s /share/d1/cumodp/data/W_65536.dat
	ln -s /share/d1/cumodp/data/W_4096.dat 
	ln -s /share/d1/cumodp/data/W_16777216.dat
	ln -s /share/d1/cumodp/data/FFT_Data_16.dat
	ln -s /share/d1/cumodp/data/FFT_4096_input.dat
	ln -s /share/d1/cumodp/data/FFT_256_input.dat 
	ln -s /share/d1/cumodp/data/FFT_1048576_input.dat
	ln -s /share/d1/cumodp/data/addition_benchmark_data.dat
	ln -s /share/d1/cumodp/data/W_1048576.dat          
	ln -s /share/d1/cumodp/data/FFT_Data_Block.dat
	ln -s /share/d1/cumodp/data/FFT_65536_input.dat
	ln -s /share/d1/cumodp/data/FFT_16777216_input.dat
	ln -s /share/d1/cumodp/data/small_number_benchmark.dat
	ln -s /share/d1/cumodp/data/big_number_benchmark.dat
