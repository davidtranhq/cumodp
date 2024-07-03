#!/bin/bash

./compileGenerateData.sh
xSeed=11
ySeed=22
size=1024
size=1048576
for i in {0,1,2,3,4,5}; do
	./generateData.sh $i $xSeed $ySeed $size
done;