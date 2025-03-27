## cumodp/new_tests/Multiplication/Makefile

SRC=test_PolyMult
##################################################
all:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -lm -o $(SRC).bin

interactive:
	nvcc -DLINUXINTEL64 interactive_test_PolyMult.cu -g -G -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -lm -o interactive_test_PolyMult.bin

##################################################
check: clean
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 500 500 469762049 	
	./$(SRC).bin 1000 500 469762049 
	./$(SRC).bin 2000 1000 469762049 
	./$(SRC).bin 2000 2000 469762049

##################################################
test:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 500 500 469762049 1	
	./$(SRC).bin 1000 500 469762049 1
	./$(SRC).bin 2000 1000 469762049 1
	./$(SRC).bin 2000 2000 469762049 1
	./$(SRC).bin 3000 1000 469762049 1
	./$(SRC).bin 3000 2000 469762049 1
	./$(SRC).bin 3000 3000 469762049 1
	./$(SRC).bin 4000 2000 469762049 1
	./$(SRC).bin 4000 4000 469762049 1
#	./$(SRC).bin 5000 1000 469762049 1
#	./$(SRC).bin 5000 5000 469762049 1	
#	./$(SRC).bin 6000 1000 469762049 1
#	./$(SRC).bin 6000 6000 469762049 1
#	./$(SRC).bin 7000 1000 469762049 1
#	./$(SRC).bin 7000 7000 469762049 1
#	./$(SRC).bin 8000 1000 469762049 1
	./$(SRC).bin 8000 8000 469762049 1
##################################################
benchmark:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 500 500 469762049 	
	./$(SRC).bin 1000 500 469762049 
	./$(SRC).bin 2000 1000 469762049 
	./$(SRC).bin 2000 2000 469762049
	./$(SRC).bin 3000 1000 469762049 
	./$(SRC).bin 3000 2000 469762049 
	./$(SRC).bin 3000 3000 469762049 
	./$(SRC).bin 4000 2000 469762049 
	./$(SRC).bin 4000 4000 469762049 
	./$(SRC).bin 5000 1000 469762049 
	./$(SRC).bin 5000 5000 469762049 	
	./$(SRC).bin 6000 1000 469762049 
	./$(SRC).bin 6000 6000 469762049 
	./$(SRC).bin 7000 1000 469762049 
	./$(SRC).bin 7000 7000 469762049 
	./$(SRC).bin 8000 1000 469762049 
	./$(SRC).bin 8000 8000 469762049 

##################################################
clean:
	rm -f *.o *~ *.linkinfo a.out *.log testPolyMul *.dat
	rm -f $(SRC).bin




