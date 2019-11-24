#!/bin/bash

# Ejemplo de parámetros
# Parametros de entrada
# numAsteroides = 5;
# numIteraciones = 200;
# numPlanetas = 8;
# semilla = 100;


if [ -f "nasteroids-seq-no-opti" ]; then
    rm -f nasteroids-seq-no-opti
fi 
if [ -f "stepbystep.txt" ]; then
    rm -f stepbystep.txt
fi  

#clear

printf "\nCompilando programa\n"  
g++ main.cpp -o nasteroids-seq-no-opti -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -g

#valgrind --tool=cachegrind --cachegrind-out-file=vgrind_seq_out --I1=16384,8,32 --LL=131072,8,64 --branch-sim=yes ./nasteroids-seq 5 2 8 100
#cg_annotate vgrind_seq_out --auto=yes

# Ejecutando programa
./nasteroids-seq-no-opti 5 2 5 2000
