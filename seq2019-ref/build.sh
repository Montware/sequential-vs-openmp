#!/bin/bash

# Ejemplo de parámetros
# Parametros de entrada
# numAsteroides = 5;
# numIteraciones = 200;
# numPlanetas = 8;
# semilla = 100;

ulimit -s unlimited
#KMP_STACKSIZE = 1m ./nasteroids-par

if [ -f "step_by_step.txt" ]; then
    rm -f step_by_step.txt
fi  


# Ejecutando programa
./nasteroids2019 100 50 50 2000
