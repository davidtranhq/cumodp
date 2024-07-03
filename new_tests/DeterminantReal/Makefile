## cumodp/new_tests/DeterminantReal/Makefile

# LDFLAGS=-L/usr/local/32-bit/lib

SRC=test_determinantReal

##################################################
all:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin

##################################################
check: clean
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 100 
	./$(SRC).bin 200

##################################################	
test:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 50 4 3 1
	./$(SRC).bin 100 4 3 1
	./$(SRC).bin 200 4 3 1
	./$(SRC).bin 300 4 3 1

##################################################
benchmark:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 100 
	./$(SRC).bin 200
	./$(SRC).bin 500
	./$(SRC).bin 1000 
	./$(SRC).bin 2000 
	./$(SRC).bin 3000

##################################################
clean:
	rm -f *.o *~ *.linkinfo a.out *.log *.dat $(SRC).bin