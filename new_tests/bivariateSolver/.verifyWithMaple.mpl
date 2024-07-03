# read in the original data file


fd := fopen("inputFile.txt", READ):
RRR := readdata(fd, [string, string]):
fclose(fd):

s1 := RRR[1][1]:
s2 := RRR[1][2]:

read s1:
read s2:
print(s1);


#data_file_name := "data_file";
#read data_file_name;




# read the output from the BiSolver
read "rc.txt":


with(RegularChains):

# substitue back the variables names
#print("regular chain from bivariate solver: "):
lrcs := subs([a=vars[1], b=vars[2]], rc mod p):

# p is the characteristic
R := PolynomialRing([vars[2], vars[1]], p):

erc := ChainTools:-Empty(R):

# make regular chain data structure
#print("regular chain from bivariate solver: "):
lrcs := map(t->ChainTools:-Chain(t, erc, R ), lrcs):

lrses := map(ConstructibleSetTools:-RegularSystem, lrcs, R):

# build the constructible set from output of BiSolve
cs1 := ConstructibleSetTools:-ConstructibleSet(lrses, R):


# build the constructible set via Triangularize

cs2 := Triangularize(eqs[1..2], [1], R, 'output'='lazard'):

#verify the result, both command are expecting true
#print("From Maple: "):
#Display(cs2,R);
#print("From bivariate solver: "):
#Display(cs1,R):
 ConstructibleSetTools:-IsContained(cs1, cs2, R);

 ConstructibleSetTools:-IsContained(cs2, cs1, R);





