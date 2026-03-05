#include "Reaction.h"
#include <iostream>
#include <vector>
#include <map>
#include <string>

// Définition des constructeurs

// On utilise la composition ici (nouvelle classe Reaction composée d'objets (les entités) de la classe existante Entite)

Reaction :: Reaction (vector<Entite*> r, vector<Entite*> p, double E){
    
    reactifs = r;
    produits = p;
    E_a = E;
    
    if (DeltaG() <= 0){
        kforward = exp(-E_a);
        kbackward = exp(-E_a - abs(DeltaG())) ;
    }
    else {
        kforward = exp(-E_a - abs(DeltaG()));
        kbackward = exp(-E_a);
    }
}


double Reaction::vitesse(bool isForward, double V){
    double v = isForward?kforward:kbackward;
    Entite* lastreac = NULL;
    int compteur =0;
    vector<Entite *> *vecEnt=isForward?&reactifs:&produits;
    for (const auto& er : *vecEnt){
        if (er==lastreac){ // we assume that identical reactants are consecutive
            compteur++;}
        else {
            compteur=0;
                
            }
            
        v*= (er->concentration - compteur)/V;
     
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










