#include <iostream>
#include <sstream>
#include <fstream>
#include "Reaction.h"
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
    Entite ("C", 500, 0, 0.5),
    Entite ("AB",50, 0, 0),
    Entite ("ABA",0, 0, 0),
    Entite ("ABAB",0, 0, 0),
    Entite ("ABC",0, 0, 0),
    Entite ("ABCB",0, 0, 0),
    Entite ("CB",0, 0, 0),
    Entite ("CBC",0, 0, 0),
    Entite ("CBCB",0, 0, 0),
};

//declarer des pointeurs vers les entites

Entite* A=&entites[0];
Entite* B=&entites[1];
Entite* C=&entites[2];
Entite* AB=&entites[3];
Entite* ABA=&entites[4];
Entite* ABAB=&entites[5];
Entite* ABC=&entites[6];
Entite* ABCB=&entites[7];
Entite* CB=&entites[8];
Entite* CBC=&entites[9];
Entite* CBCB=&entites[10];

// Création des objets de la classe Réaction
// Objets stockés dans un vecteur

// Reaction AB+A=ABA
// Reaction ABA+B=ABAB
// Reaction ABAB=2AB

vector<Reaction> reactions = {
    Reaction({A,AB}, {ABA}, 0.0),
    Reaction({ABA, B}, {ABAB}, 0.0),
    Reaction({ABAB}, {AB, AB}, 0.0),
    Reaction({AB, C}, {ABC}, 5.0), // réactions de mutation
    Reaction({ABC, B}, {ABCB}, 0.0),
    Reaction({ABCB}, {AB, CB}, 0.0),
    Reaction({CB, C}, {CBC}, 0.0),// réactions du cycle mutant
    Reaction({CBC,B}, {CBCB}, 0.0),
    Reaction({CBCB}, {CB,CB}, 0.0)
};

// Initialisation de la matrice stoechiométrique
vector<vector<int>> M (reactions.size(), vector<int> (entites.size(), 0));
double V=300;
double p_renouvelé=0.001;

// On définit la fonction pour la résolution des équations différentielles

void cycle(const vector<double>& y , vector<double> &dydt, double t){
    vector<double> v (reactions.size(), 0);
    size_t i;
    for (i=0; i< reactions.size(); i++){
        v[i]=reactions[i].vitesse(true, y, entites) - reactions[i].vitesse(false, y, entites) ;
    }
    
    for (size_t k=0; k < entites.size(); k++){
        dydt[k]=0;
        for (size_t j=0; j < reactions.size(); j++){
            dydt[k]+= M[j][k] * v[j];
        }
        dydt[k] += entites[k].concentration_ext*p_renouvelé - p_renouvelé*(y[k]);
    }
}


// Simulation cycle autocatalytique AB

int main(int argc,char* argv[]) { // pour mettre des arguments
    
    double tmax = 1000; // par défaut c'est ce temps maximal
    if (argc > 1){
        tmax=stod(argv[1]);
    }
    

    int n = 1; // pour la stoechiométrie
    Entite* lastreac = NULL;

    for (size_t i =0; i < reactions.size(); i++){
        for (const auto& r : reactions[i].reactifs){
            if (r==lastreac){
                n++;
            }
            //auto it = find(entites.begin(), entites.end(), r);
            auto it = find_if(entites.begin(), entites.end(),
                [&](const Entite &e){ return &e == r; });
            M[i][it - entites.begin()] += -n;
        }
        for (const auto& p : reactions[i].produits){
            if (p==lastreac){
                n++;
            }
            //auto it = find(entites.begin(), entites.end(), p);
            auto it = find_if(entites.begin(), entites.end(),
                [&](const Entite &e){ return &e == p; });
            M[i][it - entites.begin()] += n;
        }
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

}
