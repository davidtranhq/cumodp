#!/bin/bash

PROJECT_NAME="BigPrimeFieldFFT_3"
PROJECT=$CUMODP_HOME/new_tests/$PROJECT_NAME

reportName=$PROJECT/test/profilingReports/$1
if [ ! $# -eq 1 ]; then
	printf "\n--: Please specifiy a name for report and try again\n"
	printf "*********************************************\n"
	exit
fi

if [ -e $reportName ]; then
		rm $reportName
	fi;

for f in $PROJECT/test/profilingReports/*.log; 
do
#	cat $f | awk 'NR==3 || NR>5 {print}' >> $reportName
	cat $f >> $reportName
	rm $f
	echo "********************************************************************************">>$reportName
done; 
