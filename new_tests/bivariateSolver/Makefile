## cumodp/new_tests/bivariateSolver/Makefile

all: bivariateSolver

bivariateSolver: libmodpn
	gcc -w -DLINUXINTEL64 bivariateSolver.c -I$(MODPN_HOME)/src -I$(CUMODP_HOME)/include -L$(MODPN_HOME)/src -L$(CUMODP_HOME)/src -o bivariateSolver -lcumodp -lmodpn

bivariateSolverTester: libmodpn
	gcc -w -DLINUXINTEL64 bivariateSolverTester.c -I$(MODPN_HOME)/src -I$(CUMODP_HOME)/include -L$(MODPN_HOME)/src -L$(CUMODP_HOME)/src -o bivariateSolverTester -lcumodp -lmodpn

libmodpn:
	@cd $(MODPN_HOME)/src && if [ ! -f libmodpn.so ]; then make; fi;

check: bivariateSolver
	./bivariateSolver Sys86 469762049 
	./bivariateSolver Sys86 469762049 1
	./bivariateSolver Sys86 469762049 2
	./bivariateSolver Sys93 469762049 
	./bivariateSolver Sys93 469762049 1
	./bivariateSolver Sys93 469762049 2
	./bivariateSolver Sys203 469762049 
	./bivariateSolver Sys203 469762049 1
	./bivariateSolver Sys203 469762049 2

test: bivariateSolverTester
	./bivariateSolverTester Sys86
	./bivariateSolverTester Sys93
	./bivariateSolverTester Sys203

clean:
	rm -f *.o *~ *.linkinfo a.out *.log bivariateSolver *.dat *.txt bivariateSolverTester
	