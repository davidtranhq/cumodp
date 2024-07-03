## cumodp/new_tests/FastEvaluation/Makefile
SRC=test_fastEvaluation
##################################################
all:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin

##################################################
check: clean
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 4 469762049 
	./$(SRC).bin 5 469762049 
	./$(SRC).bin 6 469762049 

##################################################
test:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 4 469762049 1	
	./$(SRC).bin 5 469762049 1
	./$(SRC).bin 6 469762049 1
	./$(SRC).bin 7 469762049 1
	./$(SRC).bin 8 469762049 1
	./$(SRC).bin 9 469762049 1
	./$(SRC).bin 10 469762049 1
	./$(SRC).bin 11 469762049 1
	./$(SRC).bin 12 469762049 1
	./$(SRC).bin 13 469762049 1
#	./$(SRC).bin 14 469762049 1	
#	./$(SRC).bin 15 469762049 1

##################################################
benchmark:
	nvcc -DLINUXINTEL64 testFastEvaluation.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o FastEvaluation
	./$(SRC).bin 4 469762049 
	./$(SRC).bin 5 469762049 
	./$(SRC).bin 6 469762049 
	./$(SRC).bin 7 469762049 
	./$(SRC).bin 8 469762049 
	./$(SRC).bin 9 469762049 
	./$(SRC).bin 10 469762049
	./$(SRC).bin 11 469762049
	./$(SRC).bin 12 469762049
	./$(SRC).bin 13 469762049
	./$(SRC).bin 14 469762049	
	./$(SRC).bin 15 469762049
	./$(SRC).bin 16 469762049 
	./$(SRC).bin 17 469762049 
	./$(SRC).bin 18 469762049 
	./$(SRC).bin 19 469762049 
	./$(SRC).bin 20 469762049
	./$(SRC).bin 21 469762049
	./$(SRC).bin 22 469762049
	./$(SRC).bin 23 469762049

##################################################
clean:
	rm -f *.o *~ *.linkinfo a.out *.log test $(SRC).bin
	rm -f KP.dat  Points.dat  PolyF.dat  PolyFgpu.dat  PolyM.dat  PolyMgpu.dat  PolyMinv.dat  PolyMinvMgpu.dat  value.dat





