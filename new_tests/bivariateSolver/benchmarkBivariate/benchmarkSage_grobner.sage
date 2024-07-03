import time;
import os;
a=[];
b=[];

if os.path.exists('bivariateSystem.dat'):
	print("file is there")
	execfile("bivariateSystem.py")
else:
	a_0=40;
	b_0=20;
	p = 469762049;
	nElements=10;
	stepSize=15	;

writeFileGrobner=open("results_sage_groebner.txt", 'w')

writeStr="a \t b \t elapsedTime \t algorithm \n------------------------------- \n"
print(writeStr)
writeFileGrobner.write(writeStr);

# defining polynomial ring
R.<x,y>=PolynomialRing(GF(p),2, order='lex');

for i in range(nElements):
	a.append(stepSize*(i+1) + a_0);
	b.append(stepSize*(i+1) + b_0);

for i in range(nElements):
	f1=x^a[i] + y -1;
	f2=x + y^b[i] -1;
	F=[f1,f2];
	
	I=Ideal(F);
	###############################################
	startTime=time.time();
	v= I.groebner_basis();
	elapsedTime=time.time()-startTime;
	#print (a[i], b[i], "%.6f" %elapsedTime, "(groebner_basis)");

	writeStr=str(a[i])+"\t"+str(b[i])+"\t"+str("%.6f" %elapsedTime)+ "\t(groebner_basis)\n"
	print(writeStr)
	writeFileGrobner.write(writeStr)
	writeFileGrobner.write("\n")
	###############################################

writeStr="\n------------------------------- \n"
print(writeStr)