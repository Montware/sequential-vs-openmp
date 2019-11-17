#ifndef test
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <math.h>
#include <chrono>
#include <ctime>
#endif

#include "metodos.h"
//#include "estructuras.h"
using namespace std;
using namespace std::chrono;
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
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
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
                    
                    for (size_t i = 0; i <= gast.size() - 1; i++)
                    {
                        //printf("%f %f\n", gast.at(i).getCorx(), gast.at(i).getCory());
                        //Asteriode selast = gast[i];
                        calcularDistancia(gast[i], gast, planetas);
                        //cout << "distancia asteroide introducida de asteroide " << i << " con asteroide 4 = " << gast[i].getDistAsteroides()[4] << endl;

                        calcularMovimientoNormal(gast[i], gast, planetas);
                        calcularFuerzaX(gast[i], gast, planetas);
                        calcularFuerzaY(gast[i], gast, planetas);
                        calcularMovimientoAsteriode(gast[i]);
                        calcularRebotePared(gast[i]);
                    }
                    calcularReboteAsteroides(gast);
                    generateStepFile(STEPSFILE, gast, planetas, ia);
                }
                //acaba el prorama imprimiendo el archivo de salida -- paso 3
                generateOutFile(OUTFILE, gast);
            }
        }
    }

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << "Tiempo total de ejecucion = " << time_span.count() << " segundos" << endl;

    return 0;
}