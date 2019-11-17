#include <vector>
using namespace std;


class Planeta {
    private:
        double corx, cory, masa;

    public:
        Planeta(double corx1, double cory1, double masa1){
            corx = corx1;
            cory = cory1;
            masa = masa1;
        }//constructor por defecto
        virtual ~Planeta(){}; //destructor
        double getCorx(){
            return corx;
        }
        double getCory(){
            return cory;
        }
        double getMasa(){
            return masa;
        }
};

class Asteriode {
    private:
        double corx=0.0, cory=0.0, masa=0.0;
        vector <double> dist_asteriodes;
        vector <double> dist_planetas;
        vector <double> angulo_influencia, pendiente; 
        double angulo_movimiento=0.0; //angulo del moviento
        vector<double> fuerzax, fuerzay;//fuerzas de atrancion
        double fuerzax_tot=0.0, fuerzay_tot=0.0; //fuerza total
        double aceleracionx=0.0, aceleraciony=0.0; //acelacion en el instante pedido
        double velx=0.0, vely=0.0; //velocidad en el instante pedido

    public:
        Asteriode(){}
        Asteriode(double corx1, double cory1, double masa1){
            this->corx = corx1;
            this->cory = cory1;
            this->masa = masa1;
        };
        virtual ~Asteriode(){};
        double getCorx(){
            return this->corx;
        };
        double getCory(){
            return this->cory;
        };
        void setCorx(double corxq){
            this->corx=corxq;
        };
        void setCory(double coryq){
            this->cory=coryq;
        };
        double getMasa(){
            return this->masa;
        };
        void addDistAsteroides(double astw){
            this->dist_asteriodes.push_back(astw);
        };
        void addDistPlanetas(double plnw){
            this->dist_planetas.push_back(plnw);
        };
        void removeDistAsteroides(){
            this->dist_asteriodes.clear();
        };
        void removeDistPlanetas(){
            this->dist_planetas.clear();
        };
        void removeMovimientoNormal(){
            this->pendiente.clear();
            this->angulo_influencia.clear();
        };
        vector<double> getDistAsteroides(){
            return this->dist_asteriodes;
        };
        vector<double> getDistPlanetas(){
            return this->dist_planetas;
        };
        void addFuerzaX(double ast){
            this->fuerzax.push_back(ast);
        };
        void addFuerzaY(double pln){
            this->fuerzay.push_back(pln);
        };
        void removeFuerzasX(){
            this->fuerzax.clear();
        };
        void removeFuerzasY(){
            this->fuerzay.clear();
        };
        vector<double> getFuerzaX(){
            return this->fuerzax;
        };
        vector<double> getFuerzaY(){
            return this->fuerzay;
        };
        void addAnguloInfluencia(double asta){
            this->angulo_influencia.push_back(asta);
        };
        void addPendiente(double plnq){
            this->pendiente.push_back(plnq);
        };
        vector<double> getAnguloInfluencia(){
            return this->angulo_influencia;
        };
        vector<double> getPendiente(){
            return this->pendiente;
        }; 
        double getAceleracionX(){
            return this->aceleracionx;
        }; 
        double getAceleracionY(){
            return this->aceleraciony;
        };
        void setAceleracionX(double aclx){
            this->aceleracionx=aclx;
        };
        void setAceleracionY(double acly){
            this->aceleraciony=acly;
        }; 
        double getVelx(){
            return this->velx;
        };
        double getVely(){
            return this->vely;
        };
        void setVelx(double velx1){
            this->velx = velx1;
        };
        void setVely(double vely1){
            this->vely=vely1;
        };  
        double getFuerzaTotx(){
            return this->fuerzax_tot;
        };
        double getFuerzaToty(){
            return this->fuerzay_tot;
        };
        void setFuerzaTotx(double forx){
            this->fuerzax_tot=forx;
        };
        void setFuerzaToty(double fory){
            this->fuerzay_tot=fory;
        }; 
};