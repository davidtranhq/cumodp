#For AMD 64-bit (default)
ARCH=LINUXINTEL64

#For Intel Linux 32-bit
#ARCH=LINUXINTEL32

#For MAC 64-bit
#ARCH=MAC64

#For Intel Mac 32-bit
#ARCH=MAC32

#For Sparc Solaris 64-bit
#ARCH=SOLARIS64

#For PPC 64-bit
#ARCH=PPC64

#For Windows (32/64-bit)
#ARCH=WINDOWS

########################################
convertMapleConverter_h:= \
	cp MapleCConverter.h MapleCConverter.h.bak && \
	awk '{gsub(/__declspec/,"//__declspec"$1); print;}' MapleCConverter.h > MapleCConverter.h.tmp && \
	mv MapleCConverter.h.tmp MapleCConverter.h 

convertMapleConverter_c:=\
	cp MapleCConverter.c MapleCConverter.c.bak && \
	awk '{gsub(/__declspec/,"//__declspec"$1); print;}' MapleCConverter.c > MapleCConverter.c.tmp &&\
	mv MapleCConverter.c.tmp MapleCConverter.c

convertBack:= if [ -e MapleCConverter.h.bak ]; then mv MapleCConverter.h.bak MapleCConverter.h; fi && \
			  if [ -e MapleCConverter.c.bak ]; then mv MapleCConverter.c.bak MapleCConverter.c ; fi 

##########################################
# Setting compiler flags based on machine architecture
########################################
ifeq ($(ARCH), LINUXINTEL64)
	CC = gcc 
	CFLAGS=-c -O3 -Wall -fPIC -DLINUXINTEL64=1 -I$(CUMODP_HOME)/include 
	SOFLAGS= -shared 
	LDFLAGS=-lm -L$(CUMODP_HOME)/src -lcumodp
	# should we also add -lcumodp to other targets?
endif
####################
ifeq ($(ARCH), LINUXINTEL32)
	CC = gcc
	CFLAGS=-c  -O3 -fPIC -DLINUXINTEL32=1 
	SOFLAGS= -shared 
	LDFLAGS=-lm 
endif
####################
ifeq ($(ARCH), MAC64)
	CC = gcc
	CFLAGS= -O3 -fomit-frame-pointer -fno-common -fPIC -arch ppc -arch i386 -arch ppc64 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -DMAC64=1
	SOFLAGS= -bundle -fPIC -arch ppc -arch ppc64
	LDFLAGS=-lm 
endif
####################
ifeq ($(ARCH), MAC32)
	CC = gcc
	CFLAGS= -O3 -fomit-frame-pointer -fno-common -fPIC -arch ppc -arch i386 -arch ppc64 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -DMAC32=1
	SOFLAGS= -bundle -fPIC -arch ppc -arch i386  -arch ppc64
	LDFLAGS=-l
endif
####################
ifeq ($(ARCH), SOLARIS64)
	CC = gcc 
	CFLAGS=-c  -O3 -fPIC -DSOLARIS64=1
	SOFLAGS= -shared
	LDFLAGS=-lm
endif
####################
ifeq ($(ARCH), PPC64)
	CC = gcc
	CFLAGS= -O3 -fomit-frame-pointer -fno-common -fPIC -arch ppc -arch i386 -arch ppc64 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -DPPC64=1
	SOFLAGS= -bundle -fPIC -arch ppc -arch ppc64
	LDFLAGS=-lm
endif
####################
ifeq ($(ARCH), WINDOWS)
	CC = cl
	LINK = link
	CFLAGS = /nologo /w /O2 /MD 
	LIBFLAGS = /DLL /nologo /Gz
endif
########################################
SRC=generalFuncs.c \
	utilityFunctions.c \
	FMUL.c \
	FDIV.c \
	GCD.c \
	HGCD.c \
	matrix.c \
	MultiDFFT.c \
	MPMMTS.c \
	HashTable.c \
	UniHensel.c \
	FINTERP.c \
	Factorization.c \
	SubResultantSeq.c \
	IteratedResultant.c \
	LinkedList.c \
	IsInvertible.c \
	RegularGcd.c \
	CONSTANTS.c \
	MapleCConverter.c \
	solve2.c
########################################
ifeq ($(ARCH), WINDOWS)
	OBJ=$(SRC:.c=.obj)
	OBJEXE=$(SRCEXE:.c=.obj)
	dllLib=modpn.dll
else
	OBJ=$(SRC:.c=.o)
	OBJEXE=$(SRCEXE:.c=.o)
	soLib=libmodpn.so
endif

SRCEXE=$(SRC) TEST.c
EXEFILE=TEST

########################################
ifeq ($(ARCH), WINDOWS)

all:maincc exe
####################
%.obj:  $(SRC)
	$(CC) $(CFLAGS) /c $*.c
####################
maincc: $(OBJ)
	$(LINK) $(LIBFLAGS) /out:$(dllLib) $(OBJ)
	mt /manifest modpn.dll.manifest /outputresource:modpn.dll\;
####################
%.obj: $(SRCEXE)
	$(CC) $(CFLAGS) /c $*.c

####################
exe: $(OBJEXE)
	$(LINK) /out:$(EXEFILE) $(OBJEXE)

########################################
else
########################################
all: libcumodp maincc
###################
libcumodp:
	cd $(CUMODP_HOME)/src && if [ ! -f libcumodp.so ]; then make; fi;
####################
maincc: $(OBJ)
	@ $(convertMapleConverter_h)
	@ $(convertMapleConverter_c)

	$(CC)  $(SOFLAGS) -o $(soLib) $(OBJ)  $(LDFLAGS)
####################
test: 
	@ $(convertMapleConverter_h)
	@ $(convertMapleConverter_c)
	
	$(CC)  $(CFLAGS) TEST.c 
	$(CC) $(OBJEXE) -o $(EXEFILE) $(LDFLAGS)
####################
run:
	./TEST 0 257 999 234  #Testing plain multiplication
	./TEST 1 257 64 63 #Testing FFT-based multiplication
	./TEST 2 257 64 63 #Testing TFT-based multiplication
	./TEST 3 257 100 95 #Testing univariate fast division code
	./TEST 6 469762049 2000 1000 #Testing univariate fast quotient only
	./TEST 3 257 100 5
	./TEST 3 257 132 5
	./TEST 4 257 132 76 0 #Testing EEA 
	./TEST 4 257 132 76 1
	./TEST 4 257 132 76 2
	./TEST 4 257 132 76 5
	./TEST 5 257 100 50 # Testing Half GCD
	./TEST 5 469762049 500 1000
	./TEST 8 469762049 63 1 1 63 1 1 # Testing resultant3
	./TEST 8 469762049 7 7 7 7 7 7  # Testing resultant3 and raising a BUG
	./TEST 8 469762049 63 5 1 63 5 5 # Testing resultant3 and raising a BUG
	./TEST 7 $(SYSPOOL_HOME)/Systems/Sys203 469762049 0  # bivariateSolver with CUDA
	./TEST 7 $(SYSPOOL_HOME)/Systems/Sys203 469762049 # bivariateSolver without CUDA
	./TEST 9 469762049 3 3 #normal form
####################
bivariateSolver:
	@ $(convertMapleConverter_h)
	@ $(convertMapleConverter_c)

	$(CC) $(CFLAGS) bivariateSolver.c
	$(CC) $(OBJEXE) -o  bivariateSolver $(LDFLAGS)
########################################
endif
########################################
# doc:
# 	perl gdoc -html AS.h > AS.html
# 	perl gdoc -html inlineFuncs.h > inlineFuncs.html
# 	perl gdoc -html generalFuncs.c > generalFuncs.html
# 	perl gdoc -html FMUL.c > FMUL.html
# 	perl gdoc -html FDIV.c > FDIV.html
# 	perl gdoc -html GCD.c > GCD.html
# 	perl gdoc -html HGCD.c > HGCD.html
# 	perl gdoc -html matrix.c > matrix.html
# 	perl gdoc -html MultiDFFT.c > MultiDFFT.html
# 	perl gdoc -html MPMMTS.c > MPMMTS.html
# 	perl gdoc -html HashTable.c > HashTable.html
# 	perl gdoc -html UniHensel.c > UniHensel.html
# 	perl gdoc -html MapleCConverter.c > MapleCConverter.html
# 	perl gdoc -html FINTERP.c > FINTERP.html
# 	perl gdoc -html Factorization.c > Factorization.html
# 	perl gdoc -html SubResultantSeq.c >  SubResultantSeq.html
# 	perl gdoc -html IteratedResultant.c > IteratedResultant.html
# 	perl gdoc -html LinkedList.c > LinkedList.html
# 	perl gdoc -html IsInvertible.c > IsInvertible.html
# 	perl gdoc -html RegularGcd.c > RegularGcd.html
# 	perl gdoc -html CONSTANTS.c > CONSTANTS.html
########################################
clean:
	@ $(convertBack)
	rm -f main *~ *.o *.op *.DAT  *.ini a.out mk test *.obj
	rm -f callgrind.out.* ex TEST  TEST-RESULT *.exp *.manifest
	rm -f \#*  TEST.lib
	rm -f *~
	rm -f libmodpn.*
