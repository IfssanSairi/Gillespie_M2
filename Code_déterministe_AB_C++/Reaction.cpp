#include "Reaction.h"
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <cmath>

// Définition des constructeurs

// On utilise la composition ici (nouvelle classe Reaction composée d'objets (les entités) de la classe existante Entite)

Reaction :: Reaction (vector<Entite*> r, vector<Entite*> p, double E){
    
    reactifs = r;
    produits = p;
    E_a = E;
    
    //kforward et kbackward sont bien des attributs, mais non explicités en arguments dans le constructeur
    
    if (DeltaG() <= 0){
        kforward = exp(-E_a);
        kbackward = exp(-E_a - abs(DeltaG())) ;
    }
    else {
        kforward = exp(-E_a - abs(DeltaG()));
        kbackward = exp(-E_a);
    }
}

double Reaction::vitesse(bool isForward, const vector<double>& y, const vector<Entite>& entites){
    double v = isForward ? kforward : kbackward;
    vector<Entite*> *vecEnt = isForward ? &reactifs : &produits;
    
    for (const auto& er : *vecEnt){
        auto it = find_if(entites.begin(), entites.end(),
                          [&](const Entite &e){ return &e == er; });
        
        size_t j = it - entites.begin();
        
        v *= y[j];
    }
    
    return v;
}

double Reaction :: DeltaG(){
    double G_initial=0;
    double G_final=0;
    for (const auto& er : reactifs){
        G_initial += er->energie_libre; // déférencer le pointeur
        
    }
    for (const auto& er : produits){
        G_final += er->energie_libre;
        
    }
    
    return (G_final - G_initial);
}










