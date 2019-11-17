#include <random>
#include <vector>
#include <math.h>
#include <iomanip>

/* Constantes predefinidas de la simulación */
#define GRAVITY 6.674E-5
#define PERIODO 0.1
#define DISTMIN 5.0      // TODO: Borrar
#define ANCHURA 200
#define ALTURA 200
#define MEDIADISTRIBUCIONMASAS 1000
#define DESVIACIONSDM 50

#include "estructuras.h"

using namespace std;

/* TODO: Revisar saltos de línea para no salirse d elos márgenes con líneas muy largas */


/* Funciones para generacion automática de valores (en base a una distribución y desviación)
    de posicion (coordenadas X e Y) dentro de las dimensiones (200x200) del escenario de simulación.

*/
uniform_real_distribution<double> xdist(0.0, std::nextafter(ANCHURA, std::numeric_limits<double>::max()));
uniform_real_distribution<double> ydist(0.0, std::nextafter(ALTURA, std::numeric_limits<double>::max()));
normal_distribution<double> mdist(MEDIADISTRIBUCIONMASAS, DESVIACIONSDM);


/* Funciones de inicialización de los objetos */

/* Función de inicialización de asteroides y sus posiciones aleatorias.
    Recibe el número de asteroides introducidos y el valor de la semilla para la inicialización aleatoria.
    Devuelve un vector de asteroides con las posiciones X e Y y las masas ya definidas.
 */
vector<Asteroide> init_asteroides(unsigned int num_asteroides, unsigned int val_sem){
    vector<Asteroide> asteroides_vect;
    default_random_engine semilla{val_sem};

    for(unsigned int i = 0; i <= num_asteroides - 1; ++i){
        double pos_x = xdist(semilla);
        double pos_y = ydist(semilla);
        double masa = mdist(semilla);
        Asteroide asteroide(pos_x, pos_y, masa);
        asteroides_vect.push_back(asteroide); 
    }

    return asteroides_vect;
}


/* Función de inicialización de planetas y sus posiciones aleatorias a lo largo del marco del escenario.
    Recibe el número de planetas introducidos y el valor de la semilla para la inicialización aleatoria.
    Devuelve un vector de planetas con las posiciones X e Y y las masas ya definidas.
*/
vector<Planeta> init_planetas(unsigned int num_planetas, unsigned int val_sem){
    vector<Planeta> planetas_vect;
    default_random_engine semilla{val_sem};
    double pos_x = 0.0, pos_y = 0.0, masa = 0.0;

    for(unsigned int i = 0;i <= num_planetas - 1; ++i){

        /* Colocación de los planetas en los laterales del marco de forma uniformemente distribuida */
        if(i % 1 == 0){
            pos_x = xdist(semilla);
            pos_y = 0.0;
        }

        if(i % 2 == 0){
            pos_x = ANCHURA;
            pos_y = ydist(semilla);
        }

        if(i % 3 == 0){
            pos_x = xdist(semilla);
            pos_y = ALTURA;
        }

        if(i % 4 == 0){
            pos_x = 0.0;
            pos_y = ydist(semilla);
        }

        /* Mayoración de la masa de los planetas */
        masa = mdist(semilla) * 10;

        Planeta planeta(pos_x, pos_y, masa);
        planetas_vect.push_back(planeta);
    }
    return planetas_vect;
}


/* Generación de archivos */

/* Generación del archivo init_conf.txt con las posiciones iniciales y los parámetros de ejecución del programa.
    Recibe el path para al archivo init_conf.txt, los vectores con la info de los asteroides y planetas, y los parámetros de ejecución del programa.
    No devuelve nada.
*/
void gen_init_file(string init_file_path, vector<Asteroide> asteroides, vector<Planeta> planetas, unsigned int num_asteroides, unsigned int num_iteraciones, unsigned int num_planetas, unsigned int semilla){
    /* Preparación para la escritura del archivo */
    /* TODO: Revisar los parámetros de .open y los std:: */
    ofstream initconf;
    initconf.open (init_file_path, ios::out | ios::app | ios::binary);          // TODO: Revisar si ios:app es conveniente o mejor no hacer append
    initconf << std::fixed;
    initconf << std::setprecision(3);

    /* Escritura de parámetros de ejecución del programa */
    initconf << num_asteroides << " " <<  num_iteraciones << " " <<  num_planetas << " " <<  semilla << "\n";

    /* Escritura de info de los asteroides */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i){
        Asteroide a = asteroides.at(i);
        initconf << a.get_pos_x() << " " <<  a.get_pos_y() << " " <<  a.get_masa() << "\n";
    }

    /* Escritura de info de los planetas */
    for(size_t i = 0; i <= planetas.size() - 1; ++i){
        Planeta planeta = planetas.at(i);
        initconf << planeta.get_pos_x() << " " <<  planeta.get_pos_y() << " " <<  planeta.get_masa() << "\n";
    }

    initconf.close();
}


/* Generación del archivo stepbystep.txt con las fuerzas y el ángulo de influencia de éstas por las que se ve influenciado cada asteroide.
    Recibe el path para al archivo stepbystep.txt, los vectores con la info de los asteroides y planetas, y los parámetros de ejecución del programa.
    No devuelve nada.
*/
void gen_step_file(string step_file_path, vector<Asteroide> asteroides, vector<Planeta> planetas, unsigned int iteration){
    /* Preparación para la escritura del archivo */
    /* TODO: Revisar los parámetros de .open y los std:: */
    ofstream initconf;
    initconf.open (step_file_path, ios::out | ios::app | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);

    initconf << "*******ITERATION " << (iteration + 1) << "*******\n";

    /* Escritura de las fuerzas de cada asteroide en cada iteración */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i){
        Asteroide asteroide = asteroides.at(i);
        initconf << "--- asteroid " << i << " vs asteroids ---\n";
        
        /* Escritura de asteroide vs asteroides */
        for(size_t j = 0; j <= asteroides.size() - 1; ++j){
            initconf << i << " " << j << " " << asteroide.get_fuerzas_x()[j] << " " << asteroide.get_ang_influencia()[j] << "\n";
         }
        
        /* Escritura de asteroide vs planetas */
        initconf << "--- asteroid " << i << " vs planets ---\n";
        for(size_t j = 0; j <= planetas.size() - 1; ++j){
            // TODO: Revisar qué fuerza/s hay que guardar en el archivo step_by_Step
           initconf << i << " " << j << " " << asteroide.get_fuerzas_x()[j + asteroides.size()] << " " << asteroide.get_ang_influencia()[j + asteroides.size()] << "\n";
        }
    }
    initconf.close();
}


/* Generación del archivo out.txt con las posiciones finales de los asteroides, velocidades y masas.
    Recibe el path para al archivo out.txt y los vectores con la info de los ateroides.
    No devuelve nada.
*/
void gen_out_file(string out_file_path, vector<Asteroide> asteroides){
    /* Preparación para la escritura del archivo */
    /* TODO: Revisar los parámetros de .open y los std:: */
    ofstream initconf;
    initconf.open (out_file_path, ios::out | ios::app | ios::binary);          // TODO: Revisar si ios:app es conveniente o mejor no hacer append
    initconf << std::fixed;
    initconf << std::setprecision(3);

    /* Escritura de la posición final de los asteroides */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i){
        Asteroide asteroide = asteroides.at(i);
        initconf << asteroide.get_pos_x() << " " <<  asteroide.get_pos_y() << " " << asteroide.get_vel_x() << " " << asteroide.get_vel_y()  << " " <<  asteroide.get_masa() << "\n";
    }

    initconf.close();
}


/* Cálculos de fuerzas */

/* Cálculo de fuerzas de atracción X sobre un asteroide ejercicas por los demás asteroides y planetas y actualización de su info.
    Recibe el asteroide a evaluar, un vector con todos los asteroides y otro con los planetas.
    No devuelve nada.
*/
void calc_fuerzas_x(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas){
    /* Obtención de distancias */
    vector<double> dists_asteroides = asteroide.get_dist_asteroides();
    vector<double> dists_planetas = asteroide.get_dist_planetas();
    vector<double> angs_influencia = asteroide.get_ang_influencia();
    
    /* Reset de fuerzas X */
    asteroide.clear_fuerzas_x();

    /* Cálculo de componentes X de la fuerza de atracción sobre un asteroide ejercida por los demás */
    //for(size_t i = 0; i <= dist_asteroide.size() - 1; ++i){       // TODO: Borrar
    for(size_t i = 0; i <= dists_asteroides.size() - 1; ++i){

        /* Si el asteroide es él mismo o la distancia entre ellos es menor que 5 la fuerza es ignorada (nula) */
        if(asteroide.get_dist_asteroides()[i] < DISTMIN){
            asteroide.add_fuerza_x(0);
        } else {
            double fuerza_x = ((GRAVITY * asteroide.get_masa() * asteroides[i].get_masa()) / pow(asteroide.get_dist_asteroides()[i], 2)) * cos(angs_influencia[i]);
            asteroide.add_fuerza_x(fuerza_x);
        }
    }

    /* Cálculo de de componentes X de la fuerza de atracción sobre un asteroide ejercida por los planetas */
    for(size_t i = 0; i <= dists_planetas.size() - 1; ++i){
        double fuerza_x = ((GRAVITY * asteroide.get_masa() * planetas[i].get_masa()) / pow(asteroide.get_dist_planetas()[i], 2)) * cos(angs_influencia[i + asteroides.size()]);
        asteroide.add_fuerza_x(fuerza_x);
    }

    // TODO: Borrar impresión de prueba
    vector<double> fuerza_x_0 = asteroide.get_fuerzas_x();
    cout << fuerza_x_0[0] << "\n" << endl;
    cout << fuerza_x_0[1] << "\n" << endl;
}


/* Cálculo de fuerzas de atracción Y sobre un asteroide ejercicas por los demás asteroides y planetas y actualización de su info.
    Recibe el asteroide a evaluar, un vector con todos los asteroides y otro con los planetas.
    No devuelve nada.
*/
void calc_fuerzas_y(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas){
    vector<double> dists_asteroides = asteroide.get_dist_asteroides();
    vector<double> dists_planetas = asteroide.get_dist_planetas();
    vector<double> ang_influencia = asteroide.get_ang_influencia();

    /* Reset de fuerzas Y */
    asteroide.clear_fuerzas_x();

    /* Cálculo de componentes Y de la fuerza de atracción sobre un asteroide ejercida por los demás */
    //for(size_t i = 0; i <= dist_asteroide.size() - 1; ++i){       // TODO: Borrar
    for(size_t i = 0; i <= dists_asteroides.size() - 1; ++i){
        
        /* Si el asteroide es él mismo o la distancia entre ellos es menor que 5 la fuerza es ignorada (nula) */
        if(asteroide.get_dist_asteroides()[i] < DISTMIN){
            asteroide.add_fuerza_y(0);
        } else {
            double fuerza_y = ((GRAVITY * asteroide.get_masa() * asteroides[i].get_masa()) / pow(asteroide.get_dist_asteroides()[i], 2)) * sin(ang_influencia[i]);
            asteroide.add_fuerza_y(fuerza_y);
        }
    }

    /* Cálculo de de componentes Y de la fuerza de atracción sobre un asteroide ejercida por los planetas */
    for(size_t i = 0; i <= dists_planetas.size() - 1; ++i){
        double fuerza_y = ((GRAVITY * asteroide.get_masa() * planetas[i].get_masa()) / pow(asteroide.get_dist_planetas()[i], 2)) * sin(ang_influencia[i + asteroides.size()]);
        asteroide.add_fuerza_y(fuerza_y);
    }
}


/* Cálculo del movimiento normal provocado por las por los demás asteroides y planetas y actualización de su info.
    Recibe el asteroide a evaluar, un vector con todos los asteroides y otro con los planetas.
    No devuelve nada.
*/
void calc_movs_normales(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas){
    vector<double> pendiente = asteroide.get_pendiente();
    //vector<double> ang_influencia = asteroide.get_ang_influencia();       // TODO: Borrar
    //vector<double> dists_asteroides = asteroide.get_dist_asteroides();
    //vector<double> dists_planetas = asteroide.get_dist_planetas();

    /* Reset de movimientos normales*/
    asteroide.clear_movs_normales();

    /* Cálculo del movimiento normal provocado en un asteroide por los demás */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i){
        /* Cálculo de pendiente */
        double pendiente_gen = (asteroide.get_pos_y() - asteroides[i].get_pos_y()) / (asteroide.get_pos_x() - asteroides[i].get_pos_x());
        
        /* Correción antes de almacenar */
        if (pendiente_gen > 1 || pendiente_gen < -1){
            pendiente_gen = pendiente_gen - ((int) pendiente_gen / 1);
        }

        /* Comprobamos que es un número */
        if(isnan(pendiente_gen) > 0){         // TODO: Revisar qué es "isnan". Ver si es necesario este if
            pendiente_gen = 0.0;
        }
        asteroide.add_pendiente(pendiente_gen);
        asteroide.add_ang_influencia(atan(pendiente_gen));
    }

    /* Cálculo del movimiento normal provocado en un asteroide por los planetas */
    for(size_t i = 0; i <= planetas.size() - 1; ++i){
        /* Cálculo de pendiente */
        double pendiente_gen2 = (asteroide.get_pos_y() - planetas[i].get_pos_y()) / (asteroide.get_pos_x() - planetas[i].get_pos_x());
        
        /* Correción antes de almacenar */
        if (pendiente_gen2 > 1 || pendiente_gen2 < -1){
            pendiente_gen2 = pendiente_gen2 - ((int) pendiente_gen2 / 1);
        }

        /* Comprobamos que es un número */
        if(isnan(pendiente_gen2) > 0){         // TODO: Revisar qué es "isnan". Ver si es necesario este if
            pendiente_gen2 = 0.0;
        }

        asteroide.add_pendiente(pendiente_gen2);
        asteroide.add_ang_influencia(atan(pendiente_gen2));
    }
}


/* Cálculo de distancia entre un asteroide con los demás asteroides y planetas y actualización de su info.
    Recibe el asteroide a evaluar, un vector con todos los asteroides y otro con los planetas.
    No devuelve nada.
*/
void calc_distancias(Asteroide& asteroide, vector<Asteroide> asteroides, vector<Planeta> planetas){
    /* Reset de las distancias respecto a asteroides y planetas*/
    asteroide.clear_dists_asteroides();
    asteroide.clear_dists_planetas();

    /* Cálculo de la distancia con los demás asteroides */
    for(size_t i = 0; i <= asteroides.size() - 1; ++i){
        Asteroide asteroide_tmp = asteroides[i];
        double dist = sqrt(pow((asteroide.get_pos_x() - asteroide_tmp.get_pos_x()), 2) + pow((asteroide.get_pos_y() - asteroide_tmp.get_pos_y()), 2));
        asteroide.add_dist_asteroides(dist);
    }

    /* Cálculo de la distancia con los demás planetas */
    for(size_t j = 0; j <= planetas.size() - 1; ++j){
        Planeta planeta_temp = planetas[j];
        double dist = sqrt(pow((asteroide.get_pos_x() - planeta_temp.get_pos_x()), 2) + pow((asteroide.get_pos_y() - planeta_temp.get_pos_y()), 2));
        asteroide.add_dist_planetas(dist);
    }
}


/* Cálculo de movimiento final de un asteroide y actualización de su info.
    Recibe el asteroide a evaluar.
    No devuelve nada.
*/
void calc_mov_asteroide(Asteroide& asteroide){
    vector<double> fuerzas_x = asteroide.get_fuerzas_x();
    vector<double> fuerzas_y = asteroide.get_fuerzas_y();
    //vector<double> ang_influencia = asteroide.get_ang_influencia();       // TODO: Borrar

    /* Reset de fuerzas totales */
    asteroide.set_fuerza_tot_x(0.0);
    asteroide.set_fuerza_tot_y(0.0);

    /* Sumatorio de las fuerzas totales X */
    for(size_t i = 0; i <= fuerzas_x.size() - 1; ++i){
                            // TODO: borrar prints
                                                                                    printf("ESTOY AQUI\n");
                                                                                    cout << "Obeteniendo fuerza_x i=" << i << " del asteroide\n" << endl;
                                                                                    cout << fuerzas_x[i] << endl;
                                                                                 printf("ESTOY AQUI\n");


        asteroide.set_fuerza_tot_x(asteroide.get_fuerza_tot_x() + fuerzas_x[i]);

    }

    /* Sumatorio de las fuerzas totales Y */
    for(size_t i = 0; i <= fuerzas_y.size() - 1; ++i){
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
    /* Cálculo de la velocidad X*/
    double velx = asteroide.get_vel_x();
    velx = velx + acel_x * PERIODO;

    /* Cálculo de la velocidad Y*/
    double vely = asteroide.get_vel_y();
    vely = vely + acel_y * PERIODO;

    asteroide.set_vel_x(velx);
    asteroide.set_vel_y(vely);

    /* Cálculo de la posición */
    /* Cálculo de la posición X*/
    double pos_x = asteroide.get_pos_x();
    pos_x = pos_x + velx * PERIODO;

    /* Cálculo de la posición X*/
    double pos_y = asteroide.get_pos_y();
    pos_y = pos_y + vely * PERIODO;

    asteroide.set_pos_x(pos_x);
    asteroide.set_pos_y(pos_y);
}


// TODO: Dividir funciónde cálculo de movimiento para sacar fuera el cálculo de aceleración y de velocidad
/* Cálculo de aceleración */

/* Cálculo de velocidad */



/* Cálculo de rebote de un asteroide con la pared.
    Recibe el asteroide a evaluar.
    No devuelve nada.
*/
void calc_rebote_pared(Asteroide& asteroide){
    double pos_x = asteroide.get_pos_x();
    double pos_y = asteroide.get_pos_y();

    /* Cuando un el asteroide está a menos de la distancia mínima de los bordes, sale rebotado cambiando el signo de su velocidad */
    // TODO: Revisar si cuando posx <= 0 o posx <= DISTMIN o posx <= 2 y lo mismo en los demas laterales
    // TODO: En caso de ser <= DISTMIN, borrar sentencias asteroide.set_pos_x(DISTMIN);
    if(pos_x <= DISTMIN)
    {
        asteroide.set_pos_x(DISTMIN);       // TODO: Revisar si es DISTMIN (=5) o 2
        asteroide.set_vel_x(asteroide.get_vel_x() * -1);
    } else if (pos_y <= DISTMIN)
    {
        asteroide.set_pos_y(DISTMIN);       // TODO: Revisar si es DISTMIN (=5) o 2
        asteroide.set_vel_y(asteroide.get_vel_y() * -1);

    } else if (pos_x >= ANCHURA - DISTMIN)
    {
        asteroide.set_pos_x(ANCHURA - DISTMIN);       // TODO: Revisar si es DISTMIN (=5) o 2
        asteroide.set_vel_x(asteroide.get_vel_x() * -1);

    } else if (pos_y >= ALTURA - DISTMIN)
    {
        asteroide.set_pos_y(ALTURA - DISTMIN);       // TODO: Revisar si es DISTMIN (=5) o 2
        asteroide.set_vel_y(asteroide.get_vel_y() * -1);
    }
}


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


/* Cálculo de rebote entre asteroides.
    Recibe el vector de asteroides.
    No devuelve nada.
*/
void calc_rebote_asteroides(vector<Asteroide> asteroides){
    /* Copia temporal de vector asteroides para no perder las velocidades antes de los cambios */
    vector<Asteroide> asteroides_temp_copy;  

    for(size_t i = 0; i <= asteroides.size() - 1; ++i){
        Asteroide *asteroide_orig = &asteroides[i];
        Asteroide *asteroide_temp = clonar_asteroide(*asteroide_orig);        
        asteroides_temp_copy.push_back(*asteroide_temp);
    }

    /* Cálculo de intercambio de velocidades de los asteroides si estos rebotan (dist <= DISTMN */  // TODO: Revisar que es <= que DISTIM y no 0 o 2
    for(size_t i = 0; i <= asteroides.size() - 1; ++i){

        for(size_t j = 0; j <= asteroides_temp_copy.size() - 1; ++j){

            /* Comprueba que la distancia es menor a la mínima y que no está comparando el asteroide consigo mismo */
            if(asteroides_temp_copy[i].get_dist_asteroides().at(j) < DISTMIN && asteroides[i].get_pos_x() != asteroides_temp_copy[j].get_pos_x() && asteroides[i].get_pos_y() != asteroides_temp_copy[j].get_pos_y()){
                cout << "REBOTE ENTRE ASTEROIDES " << i << " con pos " << asteroides[i].get_pos_x() << ", " << asteroides[i].get_pos_y() << " y " << j << " con pos " << asteroides_temp_copy[j].get_pos_x() << ", " << asteroides_temp_copy[j].get_pos_y() << endl;
                /* Intercambias las velocidades*/
                // TODO: Revisar con la especificacion/enunciado que el cambio de velocidades está bien hecho (si es cruzado o no)
                asteroides[i].set_vel_x(-1 * asteroides_temp_copy[j].get_vel_x());
                asteroides[i].set_vel_y(-1 * asteroides_temp_copy[j].get_vel_y());
            }
        }
    }
}



