#include <iostream>
#include <sstream>
#include <fstream> // pour gérer les fichiers
#include "Entite.h"
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <list>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// Déclaration des variables globales

vector<Entite> entites = {
    Entite ("A",500, 0, 0.5),
    Entite ("B", 500, 0, 0.5),
};


// Initialisation de la matrice stoechiométrique
double V=500.;
double p_renouvelé = 0.1;

// On définit la fonction pour la résolution des équations différentielles

void cycle(const vector<double>& y , vector<double> &dydt, double t){
    
    for (size_t k=0; k < entites.size(); k++){
        dydt[k]=0;
        dydt[k] += entites[k].concentration_ext*p_renouvelé - p_renouvelé*(y[k]);
    }
}


int main(int argc,char* argv[]) { // pour mettre des arguments
    
    double tmax = 1000; // par défaut c'est ce temps maximal
    if (argc > 1){
        tmax=stod(argv[1]);
    }
    
    vector<double> dydt(entites.size());
    
    // Definition du vecteur y
    
    vector<double> y(entites.size());
    for(size_t i=0;i<entites.size();i++){
        y[i] = entites[i].effectif/V; // concentration
    }
    
    
    ofstream out("resultats.csv");
    out << "Temps";
    for (size_t i = 0; i < entites.size(); i++) out << "," << entites[i].name;
    out << "\n";

    auto observer = [&](const vector<double>& y, double t) {
        out << t;
        for (double val : y) out << "," << val;
        out << "\n";
    };
    
    // Résolution avec la fonction integrate
    double dt = 0.1; // définition du pas de temps
    
    integrate_const(runge_kutta4<vector<double>>(),cycle,y,0.0,tmax,dt,observer);

    // Définition des vecteurs pour stocker les résultats de la dynamique


}
