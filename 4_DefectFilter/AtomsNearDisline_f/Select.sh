#!/bin/bash

echo '-----------Word Counting-----------'
sh word_count.sh
gfortran  SelectingAtomsAwayDislines.f90 
./a.out

echo 'Selecting atoms away from dislines, DONE!'

