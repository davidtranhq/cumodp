
SRC=test_fastInterpolation

##################################################
all:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin 

##################################################
check: test
##################################################
test:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin 
	./$(SRC).bin 4  469762049 1 > results.dat
	./$(SRC).bin 5  469762049 1 >> results.dat
	./$(SRC).bin 6  469762049 1 >> results.dat
	./$(SRC).bin 7  469762049 1 >> results.dat
	./$(SRC).bin 8  469762049 1 >> results.dat
	./$(SRC).bin 9  469762049 1 >> results.dat
	./$(SRC).bin 10 469762049 1 >> results.dat
	./$(SRC).bin 11 469762049 1 >> results.dat
	./$(SRC).bin 12 469762049 1 >> results.dat
	./$(SRC).bin 13 469762049 0 >> results.dat
	./$(SRC).bin 14 469762049 0 >> results.dat
	./$(SRC).bin 15 469762049 0 >> results.dat
	./$(SRC).bin 16 469762049 0 >> results.dat
	./$(SRC).bin 17 469762049 0 >> results.dat
	./$(SRC).bin 18 469762049 0 >> results.dat
	./$(SRC).bin 19 469762049 0 >> results.dat
	./$(SRC).bin 20 469762049 0 >> results.dat
	./$(SRC).bin 21 469762049 0 >> results.dat
	./$(SRC).bin 22 469762049 0 >> results.dat
	./$(SRC).bin 23 469762049 0 >> results.dat
##################################################
clean:
	rm -f *.o *~ *.linkinfo a.out *.log test *.dat $(SRC).bin



