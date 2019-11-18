#!/bin/bash

# Ejemplo de parámetros
# Parametros de entrada
# numAsteroides = 5;
# numIteraciones = 200;
# numPlanetas = 8;
# semilla = 100;

ulimit -s unlimited
#KMP_STACKSIZE = 1m ./nasteroids-par

if [ -f "stepbystep" ]; then
    rm -f stepbystep
fi  

clear

printf "\nCompilando programa\n"  
g++ main.cpp -o nasteroids-par -std=c++14 -O3 -DNDEBUG -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -g

#valgrind --tool=cachegrind --cachegrind-out-file=vgrind_seq_out --I1=16384,8,32 --LL=131072,8,64 --branch-sim=yes ./nasteroids-seq 5 2 8 100
#cg_annotate vgrind_seq_out --auto=yes

# Ejecutando programa
./nasteroids-par 10 250 5 2000
