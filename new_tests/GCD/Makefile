## cumodp/new_tests/GCD/Makefile

# LDFLAGS=-L/usr/local/32-bit/lib
SRC=test_gcd
##################################################
all:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o ./$(SRC).bin
	@# nvcc -DLINUXINTEL64 gcdCuda.cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o naiveGCD
	@#g++ -m32 -O2 $(LDFLAGS) testGcdNTL.c -o testGcdNTL  -lntl  -lgmp -lm
##################################################
check: clean
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o ./$(SRC).bin
	./$(SRC).bin 469762049 500 100 
	./$(SRC).bin 469762049 500 200 
	./$(SRC).bin 469762049 500 300 
##################################################
test:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o ./$(SRC).bin
	./$(SRC).bin 469762049 500 100 1 1 1
	./$(SRC).bin 469762049 500 200 1 1 1 
	./$(SRC).bin 469762049 500 300 100 1 1
	./$(SRC).bin 469762049 500 400 10 1 1
	./$(SRC).bin 469762049 1000 500 100 1 1
	./$(SRC).bin 469762049 1000 800 200 1 1
	./$(SRC).bin 469762049 2000 1000 10 1 1
	./$(SRC).bin 469762049 2000 1500 1 1 1
	./$(SRC).bin 469762049 3000 800 50 1 1
	./$(SRC).bin 469762049 3000 1000 25 1 1 1
	./$(SRC).bin 469762049 3000 2500 100 1 1
	./$(SRC).bin 469762049 5000 5000 100 1 1
	./$(SRC).bin 469762049 6000 6000 100 1 1
	./$(SRC).bin 469762049 7000 7000 100 1 1
	./$(SRC).bin 469762049 8000 8000 100 1 1
	./$(SRC).bin 469762049 9000 9000 100 1 1
	./$(SRC).bin 469762049 10000 10000 100 1 1
	./$(SRC).bin 469762049 11000 11000 100 1 1
	./$(SRC).bin 469762049 12000 12000 100 1 1
	./$(SRC).bin 469762049 13000 13000 100 1 1
## ./$(SRC).bin <prime> <deg A - deg Gcd> <deg B - deg Gcd> <deg Gcd> <check?> 

##################################################
benchmark:
	nvcc -DLINUXINTEL64 $(SRC).cu -I$(CUMODP_HOME)/include -L $(CUMODP_HOME)/src -lcumodp -o ./$(SRC).bin
	./$(SRC).bin 469762049 500 100 
	./$(SRC).bin 469762049 500 200 
	./$(SRC).bin 469762049 500 300 
	./$(SRC).bin 469762049 500 400 
	./$(SRC).bin 469762049 1000 500
	./$(SRC).bin 469762049 1000 800
	./$(SRC).bin 469762049 2000 1000 
	./$(SRC).bin 469762049 2000 1500 
	./$(SRC).bin 469762049 3000 800 
	./$(SRC).bin 469762049 3000 1000 
	./$(SRC).bin 469762049 3000 2500 
	./$(SRC).bin 469762049 4000 3000 

##################################################
## The target old compiles a stand alone Gcd
## which does not use the "common+private parts"
## techniques.
##################################################
naiveGcd:
	nvcc -O2 test_naiveGcd.cu -o test_naiveGcd.bin
	./naiveGcd.bin 1000 100 469762049
	./naiveGcd.bin 2000 100 469762049
	./naiveGcd.bin 2000 1000 469762049
	./naiveGcd.bin 3000 1000 469762049
	./naiveGcd.bin 3000 2000 469762049
	./naiveGcd.bin 4000 3000 469762049
##################################################
ntl:
	g++ -O2 $(LDFLAGS) test_gcdNTL.c -o test_gcdNTL.bin -lntl  -lgmp -lm
	./test_gcdNTL.bin 500 200 469762049
	./test_gcdNTL.bin 500 300 469762049
	./test_gcdNTL.bin 500 400 469762049
	./test_gcdNTL.bin 1000 500 469762049
	./test_gcdNTL.bin 1000 800 469762049
	./test_gcdNTL.bin 2000 1000 469762049
	./test_gcdNTL.bin 2000 1500 469762049
	./test_gcdNTL.bin 3000 800 469762049 
	./test_gcdNTL.bin 3000 1000 469762049
	./test_gcdNTL.bin 3000 2500 469762049
	./test_gcdNTL.bin 4000 3000 469762049

##################################################
clean:
	@rm -f *.o *~ *.linkinfo *.log *.dat 
	@rm -f test_gcdNTL.bin $(SRC).bin test_naiveGcd.bin


