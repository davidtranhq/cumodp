#!/bin/bash
. ./projectProperties.conf

radix_high_pow=63
radix_low_pow=34
K=16

##################################################
g++ -w $PROJECT/src/precompute_powers_of_omega.cpp -o $PROJECT/bin/precompute -I $CUMODP_HOME/include -O3 -lntl -lgmp -lm

##################################################
for e in $(seq 1 5); do
	$PROJECT/bin/precompute $K $e $radix_high_pow $radix_low_pow $PROJECT/src/"omega_K$e" $PROJECT/src/"powers_of_omega_K$e"
done;