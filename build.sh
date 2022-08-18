#!/bin/bash

gfortran --free-form -c   ./src/utility/*.f90 ./src/*.f90 ./src/test/*.f90  -J ./output/mod
mv *.o ./output/obj

gfortran --free-form ./output/obj/test_dwt_2d.o ./output/obj/pgm.o ./output/obj/csv_file.o ./output/obj/basic.o ./output/obj/integral.o ./output/obj/wavelet_bank.o ./output/obj/wavelet_transform.o  -o test_dwt_2d 

mv test_*  ./bin