#!/bin/bash

# Ejemplo de parámetros
# Parametros de entrada
# numAsteroides = 5;
# numIteraciones = 200;
# numPlanetas = 8;
# semilla = 100;

if [ -f "nasteroids-seq" ]; then
    rm -f nasteroids-seq
fi  
if [ -f "init_conf.txt" ]; then
    rm -f init_conf.txt
fi  
if [ -f "out.txt" ]; then
    rm -f out.txt
fi  
clear

# TODO: (ver enunciado) preparar tests, añadiendo los tiempos des del programa original propio
# En el programa de referencia de aula global debe escribir por separado los tiempos

printf "\nCompilando programa:\n"  
g++ launch_tests.cpp -o launch_tests -std=c++14 -O3 -DNDEBUG -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -g


./launch_tests 10 250 5 2000

# To debug: g++ main.cpp -o nasteroids-seq -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -g
g++ main.cpp -o nasteroids-seq -std=c++14 -Wall -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors -g -DNDEBUG
#./nasteroids2018-base_v20 5 2 8 100

#valgrind --tool=cachegrind --cachegrind-out-file=vgrind_seq_out --I1=16384,8,32 --LL=131072,8,64 --branch-sim=yes ./nasteroids-seq 5 2 8 100
#cg_annotate vgrind_seq_out --auto=yes

printf "\nEjecutando programa:\n"  
./nasteroids-seq 10 250 5 2000
