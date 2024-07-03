#!/bin/bash

	# for j in $(seq 0 0); do
	# 	for i in $(seq 0 2); do ./bigArithmetic1 0 10 $i; done;
	# 	echo "###########################################################"
	# done;

	for i in $(seq 0 2); do ./profile.sh "./bigArithmetic1 0 10 $i" >> log; done;