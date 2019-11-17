# Parametros de entrada
# numAsteroides = 5;
# numIteraciones = 200;
# numPlanetas = 8;
# semilla = 100;
ulimit -s unlimited
#KMP_STACKSIZE = 1m ./nasteroids-par

if [ -f "nasteroids-seq" ]; then
    rm -f nasteroids-par
fi  
if [ -f "init_conf.txt" ]; then
    rm -f init_conf.txt
fi  
if [ -f "out.txt" ]; then
    rm -f out.txt
fi 
clear 
echo "Compilando archivo"  
g++ main.cpp -o nasteroids-par -std=c++14 -Wall -fopenmp -Wextra -Wno-deprecated -Werror -pedantic -pedantic-errors

#./nasteroids2018-base_v20 5 2 8 100

##valgrind --tool=cachegrind --cachegrind-out-file=vgrind_par_out --I1=16384,8,32 --LL=131072,8,64 --branch-sim=yes ./nasteroids-par 5 2 8 100
##cg_annotate vgrind_par_out --auto=yes
./nasteroids-par 10 250 5 2000