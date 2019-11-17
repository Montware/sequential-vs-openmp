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

using namespace std;
using namespace std::chrono;

/* Definicion de constantes del proyecto */
#define GRAVITY 6.674E-5
#define PERIODO 0.1
#define ANCHURA 200
#define ALTURA 200
#define MEDIADISTRIBUCIONMASAS 1000
#define DESVIACIONSDM 50

/* Definicion de constantes con los nombres de los archivos de entrada y salida */
#define INITFILE "init_conf.txt"
#define STEPSFILE "stepbystep.txt"
#define OUTFILE "out.txt"
#define TESTFILE "tests.csv"



/* Ejecución principal */
int main(int argc, char const *argv[])
{
    cout << "Ejecutando nasteroides-seq no optimizado" << endl;
    int num_iteraciones;
    int num_asteroides;
    int num_planetas;


    /* Inicio del temportizador para métrica de tiempo del programa*/
    high_resolution_clock::time_point program_start_time = high_resolution_clock::now();
    high_resolution_clock::time_point loops_start_time;

    /* Comprobación de los parámetros de entrada para la ejecución del programa*/
    
    /* Comprobamos que hay 5 parámetros */
    if (argc != 5)
    {
        printf("nasteroids-seq: Wrong arguments.\n");
        printf("Correct use:\n");
        printf("nasteroids-seq num_asteroides num_iteraciones num_planetas semilla\n");
        exit(EXIT_FAILURE);
    } else {
        /* Obtenemos los parámetros introducidos */
        num_asteroides = stoi(argv[1]);
        num_iteraciones = stoi(argv[2]);
        num_planetas = stoi(argv[3]);
        int semilla = stoi(argv[4]);

        /* Comprobamos ningún parámetro es negativo o semilla no positiva */
        if (num_asteroides < 0 || num_iteraciones < 0 || num_planetas < 0 || semilla <= 0)
        {
            printf("nasteroids-seq: Wrong arguments.\n");
            printf("Correct use:\n");
            printf("nasteroids-seq num_asteroides num_iteraciones num_planetas semilla\n");
            exit(EXIT_FAILURE);
        } else {
            /* Paso 1 (inicial) */
            /* Preparación de los vectores de objetos para el programa y generaciónd e datos y archivo init_config.txt */
            vector<Asteroide> asteroides = init_asteroides(num_asteroides, semilla);
            vector<Planeta> planetas = init_planetas(num_planetas, semilla);
            gen_init_file(INITFILE, asteroides, planetas, num_asteroides, num_iteraciones, num_planetas, semilla);
            vector<double> velocidades_finales_x;
            vector<double> velocidades_finales_y;

            /* Inicio del temportizador para métrica de tiempo de los bucles */
            loops_start_time = high_resolution_clock::now();

            /* Paso 2 (cálculo del movimiento de asteroides) */
            /* Cálculo de fuerzas y movimientos y actualización de la info de asteroides en cada iteración */
            for (int i = 0; i <= num_iteraciones - 1; ++i)
            {   
                /* Cálculos de cada asteroide */
                for (size_t j = 0; j <= asteroides.size() - 1; ++j)
                {
                    calc_distancias(asteroides[j], asteroides, planetas);
                    calc_movs_normales(asteroides[j], asteroides, planetas);
                    calc_fuerzas_x(asteroides[j], asteroides, planetas);
                    calc_fuerzas_y(asteroides[j], asteroides, planetas);
                    calc_mov_asteroide(asteroides[j]);
                    calc_rebote_pared(asteroides[j]);
                }

                calc_rebote_asteroides(asteroides);
                gen_step_file(STEPSFILE, asteroides, planetas, i);
            }

            // Acaba el programa imprimiendo el archivo de salida -- paso 3
            gen_out_file(OUTFILE, asteroides);
        }
    }
    

    /* Fin del temportizador para métrica de tiempo */
    high_resolution_clock::time_point program_end_time = high_resolution_clock::now();

    /* Cálculo de tiempo medio de cada iteración */
    duration<double> duracion_ejecucion_loops = duration_cast<duration<double>>(program_end_time - loops_start_time);
    duration<double> duracion_media_iteracion = duracion_ejecucion_loops / num_iteraciones;
    cout << "Tiempo medio de cada iteración = " << duracion_media_iteracion.count() << " segundos" << endl;

    /* Cálculo de tiempo total del programa */
    duration<double> duracion_ejecucion = duration_cast<duration<double>>(program_end_time - program_start_time);
    cout << "Tiempo total de ejecucion = " << duracion_ejecucion.count() << " segundos" << endl;

    /* Almacenamiento de tiempos para los tests de evaluación con métrica */
    gen_test_file(TESTFILE, num_iteraciones, num_asteroides, num_planetas, duracion_ejecucion, duracion_media_iteracion);

    return 0;
}


