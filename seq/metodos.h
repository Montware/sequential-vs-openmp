#include <random>
#include <vector>
#include <math.h>
#include <iomanip>

#define GRAVITY 6.674E-5
#define PERIODO 0.1
#define DISTMIN 2.0
#define ANCHURA 200
#define ALTURA 200
#define MEDIADISTRIBUCIONMASAS 1000
#define DESVIACIONSDM 50

#define SPACE " "
#include "estructuras.h"

using namespace std;

//funciones generacion automatica de valores
uniform_real_distribution<double> xdist(0.0, std::nextafter(ANCHURA, std::numeric_limits<double>::max()));
uniform_real_distribution<double> ydist(0.0, std::nextafter(ALTURA, std::numeric_limits<double>::max()));
normal_distribution<double> mdist(MEDIADISTRIBUCIONMASAS, DESVIACIONSDM);

//funciones de la logica del programa, fuerza, distacias y mas cosas
vector<Asteriode> initAsteriodes(int asteriodes, int vsem){
    vector<Asteriode> asteriodes_fin;
    default_random_engine semilla{vsem};
    for(int i=0;i<=asteriodes-1;i++){
        double xd = xdist(semilla);
        double yd = ydist(semilla);
        double ma = mdist(semilla);
        Asteriode ast(xd, yd, ma);
        asteriodes_fin.push_back(ast); 
    }
    return asteriodes_fin;
}

vector<Planeta> initPlanetas(int planetas, int vsem){
    vector<Planeta> planetas_fin;
    default_random_engine semilla{vsem};
    double xd=0.0, yd=0.0, ma=0.0;
    for(int i=0;i<=planetas-1;i++){
        if(i%1==0){
            xd = xdist(semilla);
            yd = 0.0;
        }
        if(i%2==0){
            xd = ALTURA;
            yd = ydist(semilla);
        }
        if(i%3==0){
            xd = xdist(semilla);
            yd = ANCHURA;
        }
        if(i%4==0){
            xd = 0.0;
            yd = ydist(semilla);
        }
        ma = mdist(semilla);
        Planeta planeta(xd, yd, ma * 10);
        planetas_fin.push_back(planeta);
    }
    return planetas_fin;
}

void generateFile(string file, vector<Asteriode> gast, vector<Planeta> planetas, int num_asteroides, int num_iteraciones, int num_planetas, int semilla){
    ofstream initconf;
    initconf.open (file, ios::out | ios::app | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);
    initconf << num_asteroides << SPACE <<  num_iteraciones << SPACE <<  num_planetas << SPACE <<  semilla << "\n";
    for(unsigned int i=0;i<=gast.size()-1;i++){
        Asteriode a = gast.at(i);
        initconf << a.getCorx() << SPACE <<  a.getCory() << SPACE <<  a.getMasa() << "\n";
    }
    for(unsigned int i=0;i<=planetas.size()-1;i++){
        Planeta a = planetas.at(i);
        initconf << a.getCorx() << SPACE <<  a.getCory() << SPACE <<  a.getMasa() << "\n";
    }
    initconf.close();
}

void generateStepFile(string file, vector<Asteriode> gast, vector<Planeta> planetas, int ia){
    ofstream initconf;
    initconf.open (file, ios::out | ios::app | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);
    initconf << "*******ITERATION " << (ia + 1) << "*******\n";
    for(unsigned int i = 0; i <= gast.size() - 1; i++){
        Asteriode a = gast.at(i);
        initconf << "--- asteroid " << i << " vs asteroids ---\n";
        for(unsigned int j = 0; j <= gast.size() - 1; j++){
            initconf << i << SPACE << j << SPACE << a.getFuerzaX()[j] << SPACE << a.getAnguloInfluencia()[j] << "\n";
         }
        initconf << "--- asteroid " << i << " vs planets ---\n";
        for(unsigned int k = 0; k <= planetas.size() - 1; k++){
           initconf << i << SPACE << k << SPACE << a.getFuerzaX()[k + gast.size()] << SPACE << a.getAnguloInfluencia()[k + gast.size()] << "\n";
        }
    }
    initconf.close();
}

void generateOutFile(string file, vector<Asteriode> gast){
    ofstream initconf;
    initconf.open (file, ios::out | ios::app | ios::binary);
    initconf << std::fixed;
    initconf << std::setprecision(3);
    for(unsigned int i=0;i<=gast.size()-1;i++){
        Asteriode a = gast.at(i);
        initconf << a.getCorx() << SPACE <<  a.getCory() << SPACE << a.getVelx() << SPACE << a.getVely()  << SPACE <<  a.getMasa() << "\n";
    }
    initconf.close();
}

void calcularFuerzaX(Asteriode& ast, vector<Asteriode> asteroides, vector<Planeta> planetas){
    vector<double> distancia_ast = ast.getDistAsteroides();
    vector<double> distancia_pla = ast.getDistPlanetas();
    vector<double> angulo_inf = ast.getAnguloInfluencia();
    ast.removeFuerzasX();
    for(size_t i=0;i<=distancia_ast.size()-1;i++){
        if(ast.getDistAsteroides()[i] == 0){
            ast.addFuerzaX(0);
        }else{
            double fuerzax = ((GRAVITY * ast.getMasa() * asteroides[i].getMasa()) / pow(ast.getDistAsteroides()[i], 2)) * cos(angulo_inf[i]);
            ast.addFuerzaX(fuerzax);
        }
    }
    for(size_t i=0;i<=distancia_pla.size()-1;i++){
        double fuerzax1 = ((GRAVITY * ast.getMasa() * planetas[i].getMasa()) / pow(ast.getDistPlanetas()[i], 2)) * cos(angulo_inf[i + asteroides.size()]);
        ast.addFuerzaX(fuerzax1);
    }
}

void calcularFuerzaY(Asteriode &ast, vector<Asteriode> asteroides, vector<Planeta> planetas){
    vector<double> distancia_ast = ast.getDistAsteroides();
    vector<double> distancia_pla = ast.getDistPlanetas();
    vector<double> angulo_inf = ast.getAnguloInfluencia();
    ast.removeFuerzasY();
    for(size_t i=0;i<=distancia_ast.size()-1;i++){
        if(ast.getDistAsteroides()[i] == 0){
            ast.addFuerzaY(0);
        }else{
            double fuerzay = ((GRAVITY * ast.getMasa() * asteroides[i].getMasa()) / pow(ast.getDistAsteroides()[i], 2)) * sin(angulo_inf[i]);
            ast.addFuerzaY(fuerzay);
        }
    }
    for(size_t i=0;i<=distancia_pla.size()-1;i++){
        double fuerzay1 = ((GRAVITY * ast.getMasa() * planetas[i].getMasa()) / pow(ast.getDistPlanetas()[i], 2)) * sin(angulo_inf[i + asteroides.size()]);
        ast.addFuerzaY(fuerzay1);
    }
}

void calcularMovimientoNormal(Asteriode &ast, vector<Asteriode> asteroides, vector<Planeta> planetas){
    vector<double> pendiente = ast.getPendiente();
    vector<double> anguloInfluencia = ast.getAnguloInfluencia();
    vector<double> distancia_ast = ast.getDistAsteroides();
    vector<double> distancia_pla = ast.getDistAsteroides();
    ast.removeMovimientoNormal();
    for(size_t i=0;i<=asteroides.size()-1;i++){
        //calculo pendiente
        double pendiente_gen = (ast.getCory()-asteroides[i].getCory())/(ast.getCorx()-asteroides[i].getCorx());
        //corecion antes de guardar
        if (pendiente_gen > 1 || pendiente_gen < -1){
            pendiente_gen = pendiente_gen-((int)pendiente_gen/1);
        }
        if(isnan(pendiente_gen)>0){
            pendiente_gen=0.0;
        }
        ast.addPendiente(pendiente_gen);
        ast.addAnguloInfluencia(atan(pendiente_gen));
    }

    for(size_t i=0;i<=planetas.size()-1;i++){
        double pendiente_gen2 = (ast.getCory()-planetas[i].getCory())/(ast.getCorx()-planetas[i].getCorx());
        //corecion antes de guardar
        if (pendiente_gen2 > 1 || pendiente_gen2 < -1){
            pendiente_gen2 = pendiente_gen2-((int)pendiente_gen2/1);
        }
        if(isnan(pendiente_gen2)>0){
            pendiente_gen2=0.0;
        }
        ast.addPendiente(pendiente_gen2);
        ast.addAnguloInfluencia(atan(pendiente_gen2));
    }
}

void calcularDistancia(Asteriode &ast, vector<Asteriode> asteroides, vector<Planeta> planetas){
    ast.removeDistAsteroides();
    ast.removeDistPlanetas();
    for(size_t i=0;i<=asteroides.size()-1;i++){
        Asteriode tmp = asteroides[i];
        double d = sqrt(pow((ast.getCorx()-tmp.getCorx()),2) + pow((ast.getCory()-tmp.getCory()),2));
        ast.addDistAsteroides(d);
    }
    for(size_t j=0;j<=planetas.size()-1;j++){
        Planeta tmp = planetas[j];
        double d1 = sqrt(pow((ast.getCorx()-tmp.getCorx()),2) + pow((ast.getCory()-tmp.getCory()),2));
        ast.addDistPlanetas(d1);
    }
}

void calcularMovimientoAsteriode(Asteriode& ast){
    vector<double> fuerzasX = ast.getFuerzaX();
    vector<double> fuerzasY = ast.getFuerzaY();
    vector<double> angulo_inf = ast.getAnguloInfluencia();
    ast.setFuerzaTotx(0);
    ast.setFuerzaTotx(0);
    for(size_t i=0;i<=fuerzasX.size()-1;i++){
        ast.setFuerzaTotx(ast.getFuerzaTotx() + fuerzasX[i]);
    }
    for(size_t i=0;i<=fuerzasY.size()-1;i++){
        ast.setFuerzaToty(ast.getFuerzaToty() + fuerzasY[i]);
    }
    //calcular aceleracion
    ast.setAceleracionX(ast.getFuerzaTotx() / ast.getMasa());
    double aclx = ast.getAceleracionX();

    ast.setAceleracionY(ast.getFuerzaToty() / ast.getMasa());
    double acly = ast.getAceleracionY();

    //calcular velocidad
    double velx = ast.getVelx();
    velx = velx + aclx * PERIODO;
    double vely = ast.getVely();
    vely = vely + acly * PERIODO;
    ast.setVelx(velx);
    ast.setVely(vely);
    //calcular posicion
    double posx = ast.getCorx();
    posx = posx + velx * PERIODO;
    double posy = ast.getCory();
    posy = posy + vely * PERIODO;
    ast.setCorx(posx);
    ast.setCory(posy);
}

void calcularRebotePared(Asteriode& ast){
    double posx = ast.getCorx();
    double posy = ast.getCory();
    if(posx <= 0)
    {
        ast.setCorx(2);
        ast.setVelx(ast.getVelx() * -1);
    }else if (posy <= 0)
    {
        ast.setCory(2);
        ast.setVely(ast.getVely() * -1);

    }else if (posx >= ANCHURA)
    {
        ast.setCorx(ANCHURA - 2);
        ast.setVelx(ast.getVelx() * -1);

    }else if (posy >= ALTURA)
    {
        ast.setCory(ALTURA - 2);
        ast.setVely(ast.getVely() * -1);
    }
}

Asteriode* clonarAsteroide(const Asteriode& orig)
{
    Asteriode* temp = new Asteriode();
    *temp = orig;
    return temp;
}

void calcularReboteAsteroides(vector<Asteriode> asteroides){
    //Copia de vector asteroides para no perder las velocidades antes de los cambios 
    vector<Asteriode> asteroidesSaveAux;   
    for(size_t i = 0; i <= asteroides.size() - 1; i++){
        Asteriode *asteroideOrig = &asteroides[i];
        Asteriode *asteroideTemp = clonarAsteroide(*asteroideOrig);        
        asteroidesSaveAux.push_back(*asteroideTemp);
    }
    //Intercambia velocidades
    for(size_t i = 0; i <= asteroides.size() - 1; i++){
        for(size_t j = 0; j <= asteroidesSaveAux.size() - 1; j++){
            if(asteroidesSaveAux[i].getDistAsteroides().at(j) < 2 && asteroides[i].getCorx() != asteroidesSaveAux[j].getCorx() && asteroides[i].getCory() != asteroidesSaveAux[j].getCory()){
                cout << "REBOTE ENTRE ASTEROIDES " << i << " con pos " << asteroides[i].getCorx() << ", " << asteroides[i].getCory() << " y " << j << " con pos " << asteroidesSaveAux[j].getCorx() << ", " << asteroidesSaveAux[j].getCory() << endl;
                asteroides[i].setVelx(-1 * asteroidesSaveAux[j].getVelx());
                asteroides[i].setVely(-1 * asteroidesSaveAux[j].getVely());
            }
        }
    }
}



