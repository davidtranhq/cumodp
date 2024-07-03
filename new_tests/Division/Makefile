## cumodp/new_tests/Division/Makefile

LDFLAGS=-L/usr/local/32-bit/lib
SRC=test_division
##################################################
all:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
##################################################
check: clean
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 500 100 469762049 
	./$(SRC).bin 500 200 469762049
	./$(SRC).bin 500 300 469762049
	./$(SRC).bin 500 400 469762049

##################################################
test:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 500 100 469762049 1
	./$(SRC).bin 500 200 469762049 1
	./$(SRC).bin 500 300 469762049 1
	./$(SRC).bin 500 400 469762049 1
	./$(SRC).bin 1000 500 469762049 1
	./$(SRC).bin 1000 800 469762049 1
	./$(SRC).bin 2000 1000 469762049 1
	./$(SRC).bin 2000 1500 469762049 1
	./$(SRC).bin 3000 800 469762049 1
	./$(SRC).bin 3000 1000 469762049 1
	./$(SRC).bin 3000 2500 469762049 1
	./$(SRC).bin 4000 3000 469762049 1
	./$(SRC).bin 5000 3000 469762049 1
	./$(SRC).bin 6000 3000 469762049 1
	./$(SRC).bin 7000 3000 469762049 1
	./$(SRC).bin 8000 3000 469762049 1
	./$(SRC).bin 9000 3000 469762049 1
	./$(SRC).bin 10000 3000 469762049 1

##################################################
benchmark:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o $(SRC).bin
	./$(SRC).bin 500 100 469762049 
	./$(SRC).bin 500 200 469762049
	./$(SRC).bin 500 300 469762049
	./$(SRC).bin 500 400 469762049
	./$(SRC).bin 1000 500 469762049
	./$(SRC).bin 1000 800 469762049
	./$(SRC).bin 2000 1000 469762049
	./$(SRC).bin 2000 1500 469762049
	./$(SRC).bin 3000 800 469762049  
	./$(SRC).bin 3000 1000 469762049
	./$(SRC).bin 3000 2500 469762049
	./$(SRC).bin 4000 3000 469762049

##################################################
## The target naive compiles a stand alone division
## which does not use the "common+provate parts"
## techniques.
naive:
	nvcc -O2 test_naiveDivision.cu -o test_naiveDivision.bin
	./test_naiveDivision.bin 1000 100 469762049
	./test_naiveDivision.bin 2000 100 469762049
	./test_naiveDivision.bin 2000 1000 469762049
	./test_naiveDivision.bin 3000 1000 469762049
	./test_naiveDivision.bin 3000 2000 469762049
	./test_naiveDivision.bin 4000 3000 469762049

##################################################
ntl:
	g++ -m32 -O2 $(LDFLAGS) test_divisionNTL.c -o $(SRC)NTL.bin  -lntl  -lgmp -lm
	# ./$(SRC)NTL.bin 500 200 469762049
	# ./$(SRC)NTL.bin 500 300 469762049
	# ./$(SRC)NTL.bin 500 400 469762049
	# ./$(SRC)NTL.bin 1000 500 469762049
	# ./$(SRC)NTL.bin 1000 800 469762049
	# ./$(SRC)NTL.bin 2000 1000 469762049
	# ./$(SRC)NTL.bin 2000 1500 469762049
	# ./$(SRC)NTL.bin 3000 800 469762049 
	# ./$(SRC)NTL.bin 3000 1000 469762049
	# ./$(SRC)NTL.bin 3000 2500 469762049
	# ./$(SRC)NTL.bin 4000 3000 469762049

##################################################
clean:
	rm -f *.o *~ *.linkinfo a.out *.log *.dat $(SRC).bin $(SRC)NTL.bin test_naiveDivision.bin




