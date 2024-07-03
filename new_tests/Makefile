SHELL = bash

# Files location
TOPDIR = $(PWD)
SRCS = Multiplication Division GCD FastEvaluation FastInterpolation SubproductTree DeterminantFiniteField DeterminantReal BigPrimeFieldFFT_3 BigPrimeFieldFFT_4 bivariateSolver 
## Work in progress:
## Taylor_shift  BigPrimeFieldFFT  FFT newGCD bivariateSolver2 

CHECKARGS = $(SRCS:%=%.check)
TESTTARGS = $(SRCS:%=%.test)
CLEANTARGS = $(SRCS:%=%.clean)


## Sanity check
check: $(CHECKARGS)

## Longer tests
test: $(TESTTARGS)

clean: $(CLEANTARGS)

%.check:
	(cd $(TOPDIR)/$* ; make check)

%.test:
	(cd $(TOPDIR)/$* ; make test)

%.clean:
	(cd $(TOPDIR)/$* ; make clean)
