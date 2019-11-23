#include <vector>
using namespace std;



/* Clase para representar los objetos de tipo Planeta.
    Son elementos estáticos con masa que no se ven sujetos a la interacción gravitatoria.
*/ 
class Planeta {
    /* Variables */
    private:
        double pos_x, pos_y, masa;

    public:
        Planeta(){}
        /* Constructor por defecto */
        Planeta(double pos_x, double pos_y, double masa){
            this -> pos_x = pos_x;
            this -> pos_y = pos_y;
            this -> masa = masa;
        }
        
        /* Destructor */
        virtual ~Planeta(){};

        /* Getters y Setters */        
        double get_pos_x(){
            return pos_x;
        }
        
        double get_pos_y(){
            return pos_y;
        }
        
        double get_masa(){
            return masa;
        }
};


/* Clase para representar los objetos de tipo Asteroide.
    Son elementos dinámicos con masa que se ven sujetos a la interacción gravitatoria y los impactos.
*/ 
class Asteroide {
    private:
        double pos_x, pos_y, masa;
        vector <double> dist_asteroides;
        vector <double> dist_planetas;
        vector <double> angs_influ_asteroides;  // Ángulos de influencia con otros asteroides
        vector <double> angs_influ_planetas;  // Ángulos de influencia con los planetas
        double ang_movimiento;       // Ángulo del movimiento
        vector<double> fuerzas_x_asteroides, fuerzas_y_asteroides, fuerzas_x_planetas, fuerzas_y_planetas;        // Fuerzas de atracción
        double fuerza_x_tot, fuerza_y_tot;        // Fuerzas totales influyentes
        double acel_x, acel_y;      // Aceleraciones
        double vel_x, vel_y;      // Velocidades

    /* Constructor por defecto */
    public:
        Asteroide(){}
        Asteroide(double pos_x, double pos_y, double masa){
            this -> pos_x = pos_x;
            this -> pos_y = pos_y;
            this -> masa = masa;
        };
                
        /* Destructor */
        virtual ~Asteroide(){};

        /* Getters y Setters */        
        double get_pos_x(){
            return this -> pos_x;
        };

        void set_pos_x(double pos_x){
            this -> pos_x = pos_x;
        };

        double get_pos_y(){
            return this -> pos_y;
        };

        void set_pos_y(double pos_y){
            this -> pos_y = pos_y;
        };

        double get_masa(){
            return this -> masa;
        };

        vector<double> get_dist_asteroides(){
            return this -> dist_asteroides;
        };

        vector<double> get_dist_planetas(){
            return this -> dist_planetas;
        };

        vector<double> get_fuerzas_x_asteroides(){
            return this -> fuerzas_x_asteroides;
        };

        vector<double> get_fuerzas_y_asteroides(){
            return this -> fuerzas_y_asteroides;
        };

        vector<double> get_fuerzas_x_planetas(){
            return this -> fuerzas_x_planetas;
        };

        vector<double> get_fuerzas_y_planetas(){
            return this -> fuerzas_y_planetas;
        };
   
        vector<double> get_angs_influ_asteroides(){
            return this -> angs_influ_asteroides;
        };

        vector<double> get_angs_influ_planetas(){
            return this -> angs_influ_planetas;
        };

        double get_acel_x(){
            return this -> acel_x;
        }; 

        void set_acel_x(double acel_x){
            this -> acel_x = acel_x;
        };

        double get_acel_y(){
            return this -> acel_y;
        };

        void set_acel_y(double acel_y){
            this -> acel_y = acel_y;
        };

        double get_vel_x(){
            return this -> vel_x;
        };

        void set_vel_x(double vel_x){
            this -> vel_x = vel_x;
        };

        double get_vel_y(){
            return this -> vel_y;
        };

        void set_vel_y(double vel_y){
            this -> vel_y=vel_y;
        };  

        double get_fuerza_tot_x(){
            return this -> fuerza_x_tot;
        };
        double get_fuerza_tot_y(){
            return this -> fuerza_y_tot;
        };
        void set_fuerza_tot_x(double fuerza_x_tot){
            this -> fuerza_x_tot = fuerza_x_tot;
        };
        void set_fuerza_tot_y(double fuerza_y_tot){
            this -> fuerza_y_tot = fuerza_y_tot;
        }; 

        /* Helpers */
        /* Añade a un asteroide la distancia con otro asteroide */
        void add_dist_asteroides(double dist_asteroide){
            this -> dist_asteroides.push_back(dist_asteroide);
        };

        /* Vacía el vector de distancias de un asteroide con los demás */
        void clear_dists_asteroides(){
            this -> dist_asteroides.clear();
        };

        /* Añade a un asteroide la distancia con planeta */
        void add_dist_planetas(double dist_planeta){
            this -> dist_planetas.push_back(dist_planeta);
        };

        /* Vacía el vector de distancias de un asteroide con los planetas */
        void clear_dists_planetas(){
            this -> dist_planetas.clear();
        };
     
        /* Vacía el vector de movimientos normales incluyendo el ángulo de influencia */
        void clear_movs_norm_asteroides(){
            this -> angs_influ_asteroides.clear();
        };

        /* Vacía el vector de movimientos normales incluyendo el ángulo de influencia */
        void clear_movs_norm_planetas(){
            this -> angs_influ_planetas.clear();
        };

        /* Añade la componente X de una fuerza influyente sobre un asteroide por parte de otros asteroides */
        void add_fuerza_x_asteroides(double fuerzas_x){
            this -> fuerzas_x_asteroides.push_back(fuerzas_x);
        };

        /* Añade la componente Y de una fuerza influyente sobre un asteroide por parte de otros asteroides */
        void add_fuerza_y_asteroides(double fuerzas_y){
            this -> fuerzas_y_asteroides.push_back(fuerzas_y);
        };

        /* Añade la componente X de una fuerza influyente sobre un asteroide por parte de los planetas */
        void add_fuerza_x_planetas(double fuerzas_x){
            this -> fuerzas_x_planetas.push_back(fuerzas_x);
        };

        /* Añade la componente Y de una fuerza influyente sobre un asteroide por parte de los planetas */
        void add_fuerza_y_planetas(double fuerzas_y){
            this -> fuerzas_y_planetas.push_back(fuerzas_y);
        };

        /* Vacía el vector de componentes X de fuerzas influyentes en el asteroide por parte de otros asteroides */
        void clear_fuerzas_x_asteroides(){
            this -> fuerzas_x_asteroides.clear();
        };

        /* Vacía el vector de componentes Y de fuerzas influyentes en el asteroide por parte de otros asteroides */
        void clear_fuerzas_y_asteroides(){
            this -> fuerzas_y_asteroides.clear();
        };

        /* Vacía el vector de componentes X de fuerzas influyentes en el asteroide por parte de los planetas */
        void clear_fuerzas_x_planetas(){
            this -> fuerzas_x_planetas.clear();
        };

        /* Vacía el vector de componentes Y de fuerzas influyentes en el asteroide por parte de los planetas */
        void clear_fuerzas_y_planetas(){
            this -> fuerzas_y_planetas.clear();
        };

        /* Añade un ángulo de influencia sobre un asteroide respecto a otro asteroide */
        void add_ang_influ_asteroides(double ang_ingluencia){
            this -> angs_influ_asteroides.push_back(ang_ingluencia);
        };

        /* Añade un ángulo de influencia sobre un asteroide respecto a un planeta*/
        void add_ang_influ_planetas(double ang_ingluencia){
            this -> angs_influ_planetas.push_back(ang_ingluencia);
        };

};