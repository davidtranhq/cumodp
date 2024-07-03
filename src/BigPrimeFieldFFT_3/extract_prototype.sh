#!/bin/bash


# src=code.cu
src=$1
device_prototypes_only=${src/.h/_device.h}
cp $src $device_prototypes_only

kernel_prototypes_only=${src/.h/.cu}
cp $src $kernel_prototypes_only

kernel_prototypes_only=${src/.h/.h}
f=$(cat $src)


# parse all prototypes
sed ':again;$!N;$!b again; :b; s/{[^{}]*}/;/g; t b' $src  > all_prototypes

# sed -e ':again;$!N;$!b again; :b; s/__device__[^(]*([^)]*[^;]*;//g; t b;' all_prototypes > $kernel_prototypes_only 

sed -e ':again;$!N;$!b again; :b; s/__device__[^(]*([^)]*[^;]*;/__device_removed__/g; t b;' all_prototypes > tmp1

# remove unmatched }
sed -e ':again;$!N;$!b again; :b; s/__device_removed__[^}]*}/__device_removed__/g; t b;' tmp1 > tmp2

# remove __device_removed__ flag
sed -e ':again;$!N;$!b again; :b; s/__device_removed__//g; t b;' tmp2 > $kernel_prototypes_only 

# less $kernel_prototypes_only

# sed ':again;$!N;$!b again; :b; s/{[^{}]*}/;/g; t b' $src  > all_prototypes

rm -f tmp1 tmp2 all_prototypes