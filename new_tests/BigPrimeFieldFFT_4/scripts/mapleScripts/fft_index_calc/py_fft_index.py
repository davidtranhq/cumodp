#!/usr/bin/python
#####################################
def ruler():
	print ("===========================");
#####################################
def prepare_for_threads(indexList,stride):
	tmpDict=dict()
	for i in range(baseCase):
		tmpDict[i]=indexList[i];
	print(tmpDict);

	for i in range(baseCase):
		tmpDict[i]=tmpDict[i]-i;
	print(tmpDict);
	print(tmpDict.values())
#####################################
baseCase=32;
def L(indexList,k,m):
	tmpList=baseCase*[0];
	steps=baseCase/(k*m);
	for s in range(steps):
		for i in range(m):
			for j in range(k):
				tmpList[k*s*m+m*j+i]=indexList[k*s*m+k*i+j];

	# print(tmpList);
	return tmpList;

#####################################
#computing twiddle T_{m}^{km}
def T(indexList,k,m):
	# tmpList=baseCase*[0];
	tmpList=dict();
	steps=baseCase/(k*m);
	for s in range(steps):
		for i in range(m):
			for j in range(k):
				key=indexList[k*s*m+k*i+j];
				# if (k*s*m+k*i+j)%(k*m)==((k*m)-1):
				# 	print(key);
				tmpList[key]=(i*j)*(baseCase/(k*m));
	# for x in sorted(tmpList):
	# 	print(tmpList[x])
	valueList=tmpList.values()
	return valueList;


#####################################
baseCase=32;
indexList=baseCase*[0];
for i in range(baseCase):
	indexList[i]=i;

print("indexList",indexList);
ruler();

indexList=L(indexList,2,baseCase/2)
indexList=L(indexList,2,baseCase/4)
indexList=L(indexList,2,baseCase/8)
indexList=L(indexList,2,baseCase/16)

print("step1", indexList)
ruler()
################D4
twiddle=T(indexList,2,2);
print("twiddle,D4",twiddle)
ruler()

indexList=L(indexList,2,2);
print("step2", indexList)
ruler()

indexList=L(indexList,2,2);

################D8
twiddle=T(indexList,2,4);
print("twiddle,D8",twiddle)
ruler()

indexList=L(indexList,4,2);
print("step3", indexList)
ruler()

indexList=L(indexList,2,4);
################D16
twiddle=T(indexList,2,8);
print("twiddle,D16",twiddle)
ruler()

indexList=L(indexList,8,2);
print("step4", indexList)
ruler()

indexList=L(indexList,2,8);
# print("step7", indexList)
# ruler()
################D32
twiddle=T(indexList,2,16);
print("twiddle,D32",twiddle)
ruler()	
indexList=L(indexList,16,2);

print("step5", indexList)
ruler()

indexList=L(indexList,2,16);
print("step6", indexList)
ruler()

prepare_for_threads(indexList,baseCase/2);


step6=indexList
for i in range(baseCase):
	indexList[i]=i;

indexList=L(indexList,2,16)
indexList=L(indexList,2,8)
indexList=L(indexList,2,4)
indexList=L(indexList,2,2)
print(indexList)
# for x in range(baseCase):
# 	if step6[x]!=indexList[x]:
# 		print("fail")
# 		break;
# 	if x==baseCase-1:
# 		print("equal")

# print(step6)	