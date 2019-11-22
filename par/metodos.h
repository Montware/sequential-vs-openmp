#include <random>
#include <vector>
#include <math.h>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <iomanip>


/* Constantes predefinidas de la simulación */
#define GRAVITY 6.674e-5
#define PERIODO 0.1
#define DISTMIN 5.0
#define ANCHURA 200
#define ALTURA 200
#define MEDIADISTRIBUCIONMASAS 1000
#define DESVIACIONSDM 50

#include "estructuras.h"

using namespace std;
using namespace std::chrono;

/* Predeclaración de funciones */
vector<Asteroide> init_asteroides(unsigned int num_asteroides, default_random_engine semilla_re);
vector<Planeta> init_planetas(unsigned int num_planetas, default_random_engine semilla_re);
void gen_init_file(string init_file_path, vector<Asteroide> asteroides, vector<Planeta> planetas,
                   unsigned int num_asteroides, unsigned int num_iteraciones,
                   unsigned int num_planetas, unsigned int semilla);
void gen_step_file(string step_file_path, vector<Asteroide> asteroides, vector<Planeta> planetas,
                   unsigned int iteration);
void gen_out_file(string out_file_path, vector<Asteroide> asteroides);
void gen_test_file(string out_file_path, int num_iteraciones, int num_asteroides, int num_planetas,
                   double duracion_ejecucion, double duracion_media_iteracion);
void calc_distancias(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas);
void calc_movs_normales(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas);
void calc_fuerzas(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas);
void calc_mov_asteroide(Asteroide& asteroide);
void calc_rebote_pared(Asteroide& asteroide);
void calc_rebote_asteroides(vector<Asteroide> asteroides);
Asteroide* clonar_asteroide(const Asteroide& orig);


/* Funciones para generacion automática de valores (en base a una distribución y desviación)
    de posicion (coordenadas X e Y) dentro de las dimensiones (200x200) del escenario de simulación.
*/
uniform_real_distribution<double> xdist(0.0, std::nextafter(ANCHURA,
                                        std::numeric_limits<double>::max()));
uniform_real_distribution<double> ydist(0.0, std::nextafter(ALTURA,
                                        std::numeric_limits<double>::max()));
normal_distribution<double> mdist(MEDIADISTRIBUCIONMASAS, DESVIACIONSDM);


/* Funciones de inicialización de los objetos */

/* Función de inicialización de asteroides y sus posiciones aleatorias.
    Recibe el número de asteroides introducidos y el valor de la semilla para la inicialización aleatoria.
    Devuelve un vector de asteroides con las posiciones X e Y y las masas ya definidas.
 */
vector<Asteroide> init_asteroides(unsigned int num_asteroides, default_random_engine semilla_re)
{
    vector<Asteroide> asteroides_vect;

    //#pragma omp parallel for ordered
    for(unsigned int i = 0; i <= num_asteroides - 1; ++i)
    {
        double pos_x = xdist(semilla_re);
        double pos_y = ydist(semilla_re);
        double masa = mdist(semilla_re);
        Asteroide asteroide(pos_x, pos_y, masa);
        //#pragma omp ordered
        asteroides_vect.push_back(asteroide); 
    }

    return asteroides_vect;
}


/* Función de inicialización de planetas y sus posiciones aleatorias a lo largo del marco del escenario.
    Recibe el número de planetas introducidos y el valor de la semilla para la inicialización aleatoria.
    Devuelve un vector de planetas con las posiciones X e Y y las masas ya definidas.
*/
vector<Planeta> init_planetas(unsigned int num_planetas, default_random_engine semilla_re)
{
    vector<Planeta> planetas_vect;
    double pos_x = 0.0, pos_y = 0.0, masa = 0.0;

    for(unsigned int i = 0; i <= num_planetas - 1; ++i)
    {
        /* Colocación de los planetas en los laterales del marco de forma repartida */
        if(i % 4 == 0)
        {
            pos_x = 0.0;
            pos_y = ydist(semilla_re);
        }

        if(i % 4 == 1)
        {
            pos_x = xdist(semilla_re);
            pos_y = 0.0;
        }

        if(i % 4 == 2)
        {
            pos_x = ANCHURA;
            pos_y = ydist(semilla_re);
        }

        if(i % 4 == 3)
        {
            pos_x = xdist(semilla_re);
            pos_y = ALTURA;
        }

        /* Mayoración x10 de la masa de los planetas */
        masa = mdist(semilla_re) * 10;

        Planeta planeta(pos_x, pos_y, masa);
        planetas_vect.push_back(planeta);
    }

    return planetas_vect;
}


/* Generación de archivos */

/* Generación del archivo init_conf.txt con las posiciones iniciales y los parámetros de ejecución
    del programa.
    Recibe el path para al archivo init_conf.txt, los vectores con la info de los asteroides y
    planetas, y los parámetros de ejecución del programa.
    No devuelve nada.
*/
void gen_init_file(string init_file_path, vector<Asteroide> asteroides, vector<Planeta> planetas,
                   unsigned int num_asteroides, unsigned int num_iteraciones,
                   unsigned int num_planetas, unsigned int semilla)
{
    /* Preparación para la escritura del archivo */
    ofstream initconf;
    initconf.open (init_file_path, ios::out | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);

    /* Escritura de parámetros de ejecución del programa */
    initconf << num_asteroides << " " <<  num_iteraciones << " " <<  num_planetas << " " <<  semilla << "\n";

    /* Escritura de info de los asteroides */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {
        Asteroide a = asteroides.at(i);
        initconf << a.get_pos_x() << " " <<  a.get_pos_y() << " " <<  a.get_masa() << "\n";
    }

    /* Escritura de info de los planetas */
    for(size_t j = 0; j <= planetas.size() - 1; ++j)
    {
        Planeta planeta = planetas.at(j);
        initconf << planeta.get_pos_x() << " " <<  planeta.get_pos_y() << " " <<  planeta.get_masa() << "\n";
    }

    initconf.close();
}


/* Generación del archivo stepbystep.txt con las fuerzas y el ángulo de influencia de éstas por
    las que se ve influenciado cada asteroide.
    Recibe el path para al archivo stepbystep.txt, los vectores con la info de los asteroides y
    planetas, y los parámetros de ejecución del programa.
    No devuelve nada.
*/
void gen_step_file(string step_file_path, vector<Asteroide> asteroides, vector<Planeta> planetas)
{
    /* Preparación para la escritura del archivo */
    ofstream initconf;
    initconf.open (step_file_path, ios::out | ios::app | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);

    /* Escritura de las fuerzas de cada asteroide en cada iteración */
    /* Escritura de asteroide vs asteroides */
    initconf << "--- asteroids vs asteroids ---\n";
    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {
        Asteroide asteroide = asteroides.at(i);
        
        for(size_t j = 0; j <= asteroides.size() - 1; ++j)
        {
            if (j > i)
            {
                initconf << i << " " << j << " " << std::fixed << std::setprecision(6) << asteroide.get_fuerzas_x()[j] << 
                " " << std::fixed << std::setprecision(6) << asteroide.get_ang_influencia()[j] << "\n";
            }
        }
    }
    
    /* Escritura de asteroide vs planetas */
    initconf << "--- asteroids vs planets ---\n";
    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {
        Asteroide asteroide = asteroides.at(i);

        for(size_t k = 0; k <= planetas.size() - 1; ++k)
        {
           initconf << i << " " << k << " " << std::fixed << std::setprecision(6) << asteroide.get_fuerzas_x()[k + asteroides.size()] <<
           " " << std::fixed << std::setprecision(6) << asteroide.get_ang_influencia()[k + asteroides.size()] << "\n";
        }
    }

    initconf << "\n******************** ITERATION ********************\n";

    initconf.close();
}


/* Generación del archivo out.txt con las posiciones finales de los asteroides, velocidades y masas.
    Recibe el path para al archivo out.txt y los vectores con la info de los ateroides.
    No devuelve nada.
*/
void gen_out_file(string out_file_path, vector<Asteroide> asteroides)
{
    /* Preparación para la escritura del archivo */
    ofstream initconf;
    initconf.open (out_file_path, ios::out | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);

    /* Escritura de la posición final de los asteroides */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {
        Asteroide asteroide = asteroides.at(i);
        initconf << asteroide.get_pos_x() << " " <<  asteroide.get_pos_y() << " " <<
        asteroide.get_vel_x() << " " << asteroide.get_vel_y()  << " " <<  asteroide.get_masa() << "\n";
    }

    initconf.close();
}

/* Generación del archivo tests.csv con las posiciones finales de los asteroides, velocidades y masas.
    Recibe el path para al archivo tests.txt y los vectores con la info de los ateroides.
    No devuelve nada.
*/
void gen_test_file(string out_file_path, int num_iteraciones, int num_asteroides, int num_planetas,
                   double duracion_ejecucion, double duracion_media_iteracion)
{
    /* Preparación para la escritura del archivo */
    ofstream initconf;
    initconf.open (out_file_path, ios::out | ios::app | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);

    initconf << num_iteraciones << ", " <<  num_asteroides << ", " << num_planetas << ", " <<
        duracion_ejecucion << ", " << duracion_media_iteracion << "\n";
}


/* Cálculo de distancia entre un asteroide con los demás asteroides y planetas y actualización de su info.
    Recibe el asteroide a evaluar, un vector con todos los asteroides y otro con los planetas.
    No devuelve nada.
*/
void calc_distancias(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas)
{
    /* Reset de las distancias respecto a asteroides y planetas*/
    asteroide.clear_dists_asteroides();
    asteroide.clear_dists_planetas();

    /* Cálculo de la distancia con los demás asteroides */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {
        Asteroide asteroide_tmp = asteroides[i];
        double dist = sqrt(pow((asteroide.get_pos_x() - asteroide_tmp.get_pos_x()), 2) +
                           pow((asteroide.get_pos_y() - asteroide_tmp.get_pos_y()), 2));
        asteroide.add_dist_asteroides(dist);
    }

    /* Cálculo de la distancia con los demás planetas */
    for(size_t j = 0; j <= planetas.size() - 1; ++j)
    {
        Planeta planeta_temp = planetas[j];
        double dist = sqrt(pow((asteroide.get_pos_x() - planeta_temp.get_pos_x()), 2) +
                           pow((asteroide.get_pos_y() - planeta_temp.get_pos_y()), 2));
        asteroide.add_dist_planetas(dist);
    }
}


/* Cálculo del movimiento normal provocado por las por los demás asteroides y planetas y actualización de su info.
    Recibe el asteroide a evaluar, un vector con todos los asteroides y otro con los planetas.
    No devuelve nada.
*/
void calc_movs_normales(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas)
{
    vector<double> pendiente = asteroide.get_pendiente();

    /* Reset de movimientos normales*/
    asteroide.clear_movs_normales();

    /* Cálculo del movimiento normal provocado en un asteroide por los demás */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {
        /* Si el asteroide es él mismo o la distancia entre ellos es menor que 5 la pendiente es ignorada */
        if(asteroide.get_dist_asteroides()[i] >= DISTMIN){
            /* Cálculo de pendiente */
            double pendiente_gen = (asteroide.get_pos_y() - asteroides[i].get_pos_y()) /
            (asteroide.get_pos_x() - asteroides[i].get_pos_x());
            
            /* Correción antes de almacenar */
            if (pendiente_gen > 1)
            {
                pendiente_gen = 1;
            } else if (pendiente_gen < -1)
            {
                pendiente_gen = -1;
            }

            /* Comprobamos que es un número correcto */
            if(isnan(pendiente_gen) > 0)
            {
                pendiente_gen = 0.0;
            }

            asteroide.add_pendiente(pendiente_gen);

            /* Cálculo del ángulo con la arcotangente */
            asteroide.add_ang_influencia(atan(pendiente_gen));
        }
    }

    /* Cálculo del movimiento normal provocado en un asteroide por los planetas */
    for(size_t j = 0; j <= planetas.size() - 1; ++j)
    {
        /* Si el asteroide es él mismo o la distancia entre ellos es menor que 5 la pendiente es ignorada */  
        if(asteroide.get_dist_planetas()[j] >= DISTMIN)
        {
            /* Cálculo de pendiente */
            double pendiente_gen2 = (asteroide.get_pos_y() - planetas[j].get_pos_y()) /
            (asteroide.get_pos_x() - planetas[j].get_pos_x());
            
            /* Correción antes de almacenar */
            if (pendiente_gen2 > 1)
            {
                pendiente_gen2 = 1;
            } else if (pendiente_gen2 < -1)
            {
                pendiente_gen2 = -1;
            }

            /* Comprobamos que es un número */
            if(isnan(pendiente_gen2) > 0)
            {
                pendiente_gen2 = 0.0;
            }

            asteroide.add_pendiente(pendiente_gen2);

            /* Cálculo del ángulo con la arcotangente */
            asteroide.add_ang_influencia(atan(pendiente_gen2));
        }
    }
}


/* Cálculo de fuerzas de atracción X e Y sobre un asteroide ejercicas por los demás asteroides y planetas
    y actualización de su info.
    Recibe el asteroide a evaluar, un vector con todos los asteroides y otro con los planetas.
    No devuelve nada.
*/
void calc_fuerzas(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas)
{
    /* Obtención de distancias */    
    vector<double> dists_asteroides = asteroide.get_dist_asteroides();
    vector<double> dists_planetas = asteroide.get_dist_planetas();
    vector<double> angs_influencia = asteroide.get_ang_influencia();
    
    /* Reset de fuerzas */
    asteroide.clear_fuerzas_x();
    asteroide.clear_fuerzas_y();

    /* Optimización A1: Fusión de bucles*/
    /* Cálculo de componentes X e Y de la fuerza de atracción sobre un asteroide ejercida por los demás */
    for(size_t i = 0; i <= dists_asteroides.size() - 1; ++i)
    {
        double fuerza_x;
        double fuerza_y;

        /* Si el asteroide es él mismo o la distancia entre ellos es menor que 5 la fuerza es ignorada (nula) */
        if(asteroide.get_dist_asteroides()[i] < DISTMIN)
        {
            fuerza_x = 0.0;
            fuerza_y = 0.0;
        } else
        {
            fuerza_x = ((GRAVITY * asteroide.get_masa() * asteroides[i].get_masa()) /
                         pow(asteroide.get_dist_asteroides()[i], 2)) * cos(angs_influencia[i]);
            fuerza_y = ((GRAVITY * asteroide.get_masa() * asteroides[i].get_masa()) /
                         pow(asteroide.get_dist_asteroides()[i], 2)) * sin(angs_influencia[i]);
        }

        /* Comprobación de fuerzas superiors a 100 */
        if (fuerza_x > 100.0)
        {
            fuerza_x = 100.0;
        }

        if (fuerza_y > 100.0)
        {
            fuerza_y = 100.0;
        }

        asteroide.add_fuerza_x(fuerza_x);
        asteroide.add_fuerza_y(fuerza_y);

    }

    /* Optimización A2: Fusión de bucles*/
    /* Cálculo de de componentes X e Y de la fuerza de atracción sobre un asteroide ejercida por los planetas */
    for(size_t j = 0; j <= dists_planetas.size() - 1; ++j)
    {
        double fuerza_x = ((GRAVITY * asteroide.get_masa() * planetas[j].get_masa()) /
                            pow(asteroide.get_dist_planetas()[j], 2)) *
                            cos(angs_influencia[j + asteroides.size()]);
        double fuerza_y = ((GRAVITY * asteroide.get_masa() * planetas[j].get_masa()) /
                            pow(asteroide.get_dist_planetas()[j], 2)) *
                             sin(angs_influencia[j + asteroides.size()]);
        
        /* Comprobación de fuerzas superiors a 100 */        
        if (fuerza_x > 100.0)
        {
            fuerza_x = 100.0;
        }

        if (fuerza_y > 100.0)
        {
            fuerza_y = 100.0;
        }

        asteroide.add_fuerza_x(fuerza_x);
        asteroide.add_fuerza_y(fuerza_y);
    }
}


/* Cálculo de movimiento final de un asteroide y actualización de su info.
    Recibe el asteroide a evaluar.
    No devuelve nada.
*/
void calc_mov_asteroide(Asteroide& asteroide)
{
    vector<double> fuerzas_x = asteroide.get_fuerzas_x();
    vector<double> fuerzas_y = asteroide.get_fuerzas_y();

    /* Reset de fuerzas totales */
    asteroide.set_fuerza_tot_x(0.0);
    asteroide.set_fuerza_tot_y(0.0);

    /* Optimización A3: Fusión de bucles*/
    /* Sumatorio de las fuerzas totales X e Y optimizado */
    for(size_t i = 0; i <= fuerzas_x.size() - 1 && i <= fuerzas_y.size() - 1; ++i)
    {
        asteroide.set_fuerza_tot_x(asteroide.get_fuerza_tot_x() + fuerzas_x[i]);
        asteroide.set_fuerza_tot_y(asteroide.get_fuerza_tot_y() + fuerzas_y[i]);
    }

    /* Cálculo de la aceleración */
    /* Cálculo de la aceleración X*/
    asteroide.set_acel_x(asteroide.get_fuerza_tot_x() / asteroide.get_masa());
    double acel_x = asteroide.get_acel_x();

    /* Cálculo de la aceleración Y*/
    asteroide.set_acel_y(asteroide.get_fuerza_tot_y() / asteroide.get_masa());
    double acel_y = asteroide.get_acel_y();

    /* Cálculo de la velocidad */
    /* Cálculo de la velocidad X */
    double velx = asteroide.get_vel_x();
    velx = velx + acel_x * PERIODO;

    /* Cálculo de la velocidad Y */
    double vely = asteroide.get_vel_y();
    vely = vely + acel_y * PERIODO;

    asteroide.set_vel_x(velx);
    asteroide.set_vel_y(vely);

    /* Cálculo de la posición */
    /* Cálculo de la posición X */
    double pos_x = asteroide.get_pos_x();
    pos_x = pos_x + velx * PERIODO;

    /* Cálculo de la posición Y */
    double pos_y = asteroide.get_pos_y();
    pos_y = pos_y + vely * PERIODO;

    asteroide.set_pos_x(pos_x);
    asteroide.set_pos_y(pos_y);
}


/* Cálculo de rebote de un asteroide con la pared.
    Recibe el asteroide a evaluar.
    No devuelve nada.
*/
void calc_rebote_pared(Asteroide& asteroide)
{
    double pos_x = asteroide.get_pos_x();
    double pos_y = asteroide.get_pos_y();

    /* Cuando un el asteroide está a menos de la distancia mínima de los bordes, sale rebotado
        cambiando el signo de su velocidad */
    if(pos_x <= 0)
    {
        asteroide.set_pos_x(DISTMIN);
        asteroide.set_vel_x(asteroide.get_vel_x() * -1);
    } else if (pos_y <= 0)
    {
        asteroide.set_pos_y(DISTMIN);
        asteroide.set_vel_y(asteroide.get_vel_y() * -1);

    } else if (pos_x >= ANCHURA)
    {
        asteroide.set_pos_x(ANCHURA - DISTMIN);
        asteroide.set_vel_x(asteroide.get_vel_x() * -1);

    } else if (pos_y >= ALTURA)
    {
        asteroide.set_pos_y(ALTURA - DISTMIN);
        asteroide.set_vel_y(asteroide.get_vel_y() * -1);
    }
}


/* Cálculo de rebote entre asteroides.
    Recibe el vector de asteroides.
    No devuelve nada.
*/
void calc_rebote_asteroides(vector<Asteroide> asteroides)
{
    /* Copia temporal de vector asteroides para no perder las velocidades antes de los cambios */
    vector<Asteroide> asteroides_temp_copy;  

    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {
        Asteroide *asteroide_orig = &asteroides[i];
        Asteroide *asteroide_temp = clonar_asteroide(*asteroide_orig);        
        asteroides_temp_copy.push_back(*asteroide_temp);
    }

    /* Cálculo de intercambio de velocidades de los asteroides si estos rebotan (dist <= DISTMIN) */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i)
    {

        for(size_t j = 0; j <= asteroides_temp_copy.size() - 1; ++j)
        {
            /* Comprueba que la distancia es menor a la mínima y que no está comparando el asteroide consigo mismo */
            if(asteroides_temp_copy[i].get_dist_asteroides().at(j) < DISTMIN &&
               asteroides[i].get_pos_x() != asteroides_temp_copy[j].get_pos_x() &&
               asteroides[i].get_pos_y() != asteroides_temp_copy[j].get_pos_y())
               {
                cout << "REBOTE ENTRE ASTEROIDES " << i << " con pos " << asteroides[i].get_pos_x() <<
                ", " << asteroides[i].get_pos_y() << " y " << j << " con pos " <<
                asteroides_temp_copy[j].get_pos_x() << ", " << asteroides_temp_copy[j].get_pos_y() << endl;

                /* Intercambio de velocidades*/
                asteroides[i].set_vel_x(-1 * asteroides_temp_copy[j].get_vel_x());
                asteroides[i].set_vel_y(-1 * asteroides_temp_copy[j].get_vel_y());
            }
        }
    }
}


/* Helpers */

/* Función para clonar un objeto Asteroide.
    Recibe el asteroide a clonar.
    Devuelve el asteroide nuevo.
*/
Asteroide* clonar_asteroide(const Asteroide& orig)
{
    Asteroide* new_asteroide = new Asteroide();
    *new_asteroide = orig;
    return new_asteroide;
}



