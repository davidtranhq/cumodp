
all: clean
	(cd $(CUMODP_HOME)/src; make)

##################################################	
install:
	./install_cumodp.sh
##################################################
clean:
	(cd $(CUMODP_HOME); rm -fr tmp)
	(cd $(CUMODP_HOME)/src; make clean)
	(cd $(CUMODP_HOME)/new_tests; make clean)
	(rm -rf help/html help/latex)

##################################################
check:
	(cd $(CUMODP_HOME)/new_tests; make check)
##################################################
doc:
	doxygen
