#!/bin/bash

	for i in $(seq 0 0); do
		for j in $(seq 0 2); do ./bigArithmetic1 0 0 1 $i; done;
		echo "###########################################################"
	done;