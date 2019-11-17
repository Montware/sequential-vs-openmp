#ifndef test
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <math.h>
#include <omp.h>
#endif

#include "metodos.h"
using namespace std;
// Definicion de constantes del proyecto
#define GRAVITY 6.674E-5
#define PERIODO 0.1
#define DISTMIN 2.0
#define ANCHURA 200
#define ALTURA 200
#define MEDIADISTRIBUCIONMASAS 1000
#define DESVIACIONSDM 50
#define SPACE " "

#define INITFILE "init_conf.txt"
#define STEPSFILE "stepbystep.txt"
#define OUTFILE "out.txt"

int main(int argc, char const *argv[])
{
    cout << "Ejecutando archivo" << endl;
    double t1 = omp_get_wtime();
    // Prepara el programa para procesamientos en paralelo
    int n_threads = omp_get_max_threads();
    cout << "Numero maximo de nucleos = " << n_threads << endl;

    //checkeo de los datos introducidos como parametros
    if (argc <= 4)
    {
        printf("nasteroids-seq: Wrong arguments.\n");
        printf("Correct use:\n");
        printf("nasteroids-seq num_asteroides num_iteraciones num_planetas semilla\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        if (argc != 5)
        {
            printf("nasteroids-seq: Wrong arguments.\n");
            printf("Correct use:\n");
            printf("nasteroids-seq num_asteroides num_iteraciones num_planetas semilla\n");
            exit(EXIT_FAILURE);
        }
        else
        {
            int num_asteroides = stoi(argv[1]);
            int num_iteraciones = stoi(argv[2]);
            int num_planetas = stoi(argv[3]);
            int semilla = stoi(argv[4]);
            if (num_asteroides < 0 || num_iteraciones < 0 || num_planetas < 0 || semilla < 0)
            {
                printf("nasteroids-seq: Wrong arguments.\n");
                printf("Correct use:\n");
                printf("nasteroids-seq num_asteroides num_iteraciones num_planetas semilla\n");
                exit(EXIT_FAILURE);
            }
            else
            {
                //precarga el programa y genera datos -- paso 1
                vector<Asteriode> gast = initAsteriodes(num_asteroides, semilla);
                vector<Planeta> planetas = initPlanetas(num_planetas, semilla);
                generateFile(INITFILE, gast, planetas, num_asteroides, num_iteraciones, num_planetas, semilla);
                vector<double> velocidades_finales_x;
                vector<double> velocidades_finales_y;
                //calcula fuerzas y genera datos -- paso 2
                for (int ia = 0; ia <= num_iteraciones - 1; ia++)
                {    
                    #pragma omp parallel for ordered num_threads(n_threads) firstprivate(planetas)                     
                    for (size_t i = 0; i <= gast.size() - 1; ++i)
                    {
                        #pragma omp ordered
                        {
                            calcularDistancia(gast[i], gast, planetas, n_threads);
                            calcularMovimientoNormal(gast[i], gast, planetas, n_threads);
                            calcularFuerzaX(gast[i], gast, planetas, n_threads);
                            calcularFuerzaY(gast[i], gast, planetas, n_threads);
                            calcularMovimientoAsteriode(gast[i], n_threads);
                            calcularRebotePared(gast[i]);
                        }
                    }
                    calcularReboteAsteroides(gast, n_threads);
                    generateStepFile(STEPSFILE, gast, planetas, ia);
                }
                //acaba el prorama imprimiendo el archivo de salida -- paso 3
                generateOutFile(OUTFILE, gast);
            }
        }
    }
    double t2 = omp_get_wtime();
    double t_ejec = t2 - t1;
    cout << "Tiempo total de ejecucion = " << t_ejec << " segundos" << endl;
    return 0;
}