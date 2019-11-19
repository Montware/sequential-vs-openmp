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
#define GRAVITY 6.674e-5
#define PERIODO 0.1
#define DISTMIN 5.0
#define ANCHURA 200
#define ALTURA 200
#define MEDIADISTRIBUCIONMASAS 1000
#define DESVIACIONSDM 50

/* Definicion de constantes con los nombres de los archivos de entrada y salida */
#define INITFILE "init_conf.txt"
#define STEPSFILE "stepbystep.txt"
#define OUTFILE "out.txt"
#define TESTFILE "tests.csv"



/* Predeclaración de funciones */
void print_program_info(int num_asteroides, int num_iteraciones, string init_fpath,
                        string out_fpath, int num_planetas, int semilla, double gravity,
                        double delta, double min_dist, int anchura, int altura);


/* Ejecución principal */
int main(int argc, char const *argv[])
{
    //cout << "Ejecutando nasteroides-seq optimizado" << endl;
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
            /* Impresión por pantalla d ela información del programa */
            print_program_info(num_asteroides, num_iteraciones, INITFILE, OUTFILE, num_planetas,
                               semilla, GRAVITY, PERIODO, DISTMIN, ANCHURA, ALTURA);
            
            /*Generación de semilla aleatoria*/
            default_random_engine semilla_re{semilla};
                   
            /* Paso 1 (inicial) */          
            /* Preparación de los vectores de objetos para el programa y generaciónd e datos y archivo init_config.txt */
            vector<Asteroide> asteroides = init_asteroides(num_asteroides, semilla_re);
            vector<Planeta> planetas = init_planetas(num_planetas, semilla_re);
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
                    calc_fuerzas(asteroides[j], asteroides, planetas);
                    calc_mov_asteroide(asteroides[j]);
                    calc_rebote_pared(asteroides[j]);
                }

                calc_rebote_asteroides(asteroides);
                gen_step_file(STEPSFILE, asteroides, planetas);
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
    cout << "\nTiempo medio de cada iteración = " << duracion_media_iteracion.count() << " segundos" << endl;

    /* Cálculo de tiempo total del programa */
    duration<double> duracion_ejecucion = duration_cast<duration<double>>(program_end_time - program_start_time);
    cout << "Tiempo total de ejecucion = " << duracion_ejecucion.count() << " segundos" << endl;

    /* Almacenamiento de tiempos para los tests de evaluación con métrica */
    gen_test_file(TESTFILE, num_iteraciones, num_asteroides, num_planetas,
                  duracion_ejecucion.count(), duracion_media_iteracion.count());
    return 0;
}


/* Función de impresión de la información de la ejecución */
void print_program_info(int num_asteroides, int num_iteraciones, string init_fpath,
                        string out_fpath, int num_planetas, int semilla, double gravity, double delta, double min_dist, int anchura, int altura){
    cout << "Execution setup" << endl;
    cout << "\nNumber of bodies: " << num_asteroides << endl;
    cout << "Number of iterations: " << num_iteraciones << endl;
    cout << "Initial file: " << init_fpath << endl;
    cout << "Output file: " << out_fpath << endl;
    cout << "Number of planets: " << num_planetas << endl;
    cout << "Seed: " << semilla << endl;
    cout << "\nNumber of bodies: " << num_asteroides << endl;
    cout << "Gravity: " << gravity << endl;
    cout << "Delta time: " << delta << endl;
    cout << "Number of steps: " << num_iteraciones << endl;
    cout << "Min. distance: " << min_dist << endl;
    cout << "Width: " << anchura << endl;
    cout << "Height: " << altura << endl;
}

