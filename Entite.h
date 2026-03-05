#ifndef ENTITE_H
#define ENTITE_H
#include <vector>
#include <string>
#include <map>

using namespace std;

class Entite {
    public :
    std::string name;
    double concentration;
    double energie_libre;
    double taux_creation;
    double taux_destruction;
    
    Entite(std::string n, double c, double e_l, double tc=0, double td=0);
};

#endif // ENTITE_H



