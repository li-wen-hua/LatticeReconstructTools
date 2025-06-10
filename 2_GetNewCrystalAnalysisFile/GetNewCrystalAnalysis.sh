#!/bin/bash


#========================
# Manually specify the location of the dislocation loops or dislocation lines to be removed
x=101.096
y=96.4252
z=-99.2811
#========================









# parameter changing
replacement="loops(1,1:3) = ($x, $y, $z)"
sed -i "62s#.*#${replacement}#" loopremoving.f90
# count the start-reading line
bash word_count.sh
# loop removing
gfortran loopremoving.f90
./a.out > log.dat
# counting vertex
n1=$(awk 'NR==5 {print $3}' log.dat)
n2=$(awk 'NR==10 {print $3}' log.dat)
echo "Number of Vertex in Old Crystal Analysis file : $n1"
echo "Number of Vertex in New Crystal Analysis file : $n2"
# modify nodes.dump file
sed -i "4s/.*/$n2/" nodes.dump
# outputs
mv newca.input datafiles_outputs/disline_new.input
mv nodes.dump datafiles_outputs/

echo 'Selection is DONE!'

