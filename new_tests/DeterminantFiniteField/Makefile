## cumodp/new_tests/DeterminantFiniteField/Makefile

LDFLAGS=-L/usr/local/32-bit/lib
SRC=test_determinantFF
SRC_NoPivot=test_determinantFFNoPivotGauss
##################################################
all:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	nvcc -DLINUXINTEL64 $(SRC_NoPivot).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC_NoPivot)
##################################################
check: clean
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	nvcc -DLINUXINTEL64 $(SRC_NoPivot).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC_NoPivot)
	
	./$(SRC).bin 100 469762049 4 128 8 0
	./$(SRC).bin 200 469762049 4 128 8 0
	./$(SRC).bin 500 469762049 4 128 8 0
	./$(SRC_NoPivot) 128 469762049
	./$(SRC_NoPivot) 256 469762049
	./$(SRC_NoPivot) 512 469762049

##################################################
test:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	nvcc -DLINUXINTEL64 $(SRC_NoPivot).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC_NoPivot)

	./$(SRC).bin 50 469762049 4 128 8 1
	./$(SRC).bin 100 469762049 4 128 8 1
	./$(SRC).bin 200 469762049 4 128 8 1
	./$(SRC).bin 300 469762049 4 128 8 1
	./$(SRC_NoPivot) 128 469762049
	./$(SRC_NoPivot) 256 469762049
	./$(SRC_NoPivot) 512 469762049
benchmark:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	nvcc -DLINUXINTEL64 $(SRC_NoPivot).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC_NoPivot)

	./$(SRC).bin 100 469762049 4 128 8 0
	./$(SRC).bin 200 469762049 4 128 8 0
	./$(SRC).bin 500 469762049 4 128 8 0
	./$(SRC).bin 1000 469762049 4 128 8 0
	./$(SRC).bin 2000 469762049 4 128 8 0
	./$(SRC).bin 3000 469762049 4 128 8 0
	./$(SRC_NoPivot) 128 469762049
	./$(SRC_NoPivot) 256 469762049
	./$(SRC_NoPivot) 512 469762049

##################################################
clean:
	rm -f *.o *~ *.linkinfo a.out *.log *.dat $(SRC).bin $(SRC_NoPivot)




