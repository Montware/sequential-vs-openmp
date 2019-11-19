#!/bin/bash

# Ejemplo de parámetros
# Parametros de entrada
# numAsteroides = 5;
# numIteraciones = 200;
# numPlanetas = 8;
# semilla = 100;

# Borrado del step_by_step.txt proporcionado
if [ -f "../seq2019-ref/step_by_step.txt" ]; then
    rm -f ../seq2019-ref/step_by_step.txt
fi  

# Reset de los tests.csv
if [ -f "../seq_no_opti/tests.csv" ]; then
    rm -f ../seq_no_opti/tests.csv
fi  

if [ -f "../seq/tests.csv" ]; then
    rm -f ../seq/tests.csv
fi  

if [ -f "../par/tests.csv" ]; then
    rm -f ../par/tests.csv/media/dev/DATA/CLOUD/0_MEGA_UC3M/WS/'VISUAL STUDIO'/UNI/ARCOS/sequential-vs-openmp
fi  

clear


printf "\nEjecutando tests de los programas:\n"  
#../seq2019-ref/build.sh 10 250 5 2000
../seq2019-ref/build.sh
../seq_no_opti/build.sh
../seq/build.sh
../par/build.sh
