#!/bin/bash

# Ejemplo de parámetros
# Parametros de entrada
# numAsteroides = 5;
# numIteraciones = 200;
# numPlanetas = 8;
# semilla = 100;

ulimit -s unlimited
#KMP_STACKSIZE = 1m ./nasteroids-par        // TODO: Ver si necesario

if [ -f "stepbystep" ]; then
    rm -f stepbystep
fi  


# Ejecutando programa
./nasteroids2019 5 2 1 2000
