#include <iostream>
#include <fstream> // pour gérer les fichiers
#include "Reaction.h"
#include "Entite.h"
#include <cmath>
#include <random>
#include <vector>
#include <map>
#include <string>
#include <list>


using namespace std;

// Simulation cycle autocatalytique AB

int main(int argc,char* argv[]) {
    // Par défaut, rand() fixe la graine (même résultat à chaque simulation)
    // Donc srand pour ne pas fixer la graine
    
    srand(time(NULL));
    
    // Déclaration et définition des variables propres à l'algo Gillespie
    
    double V=50.; // Volume réactionnel
    
    double tmax = 7000;
    if (argc >=1){
        tmax=stod(argv[1]);
        
    }
    
    // Liste des entités primaires
    // list<string> entites_p = {"A", "B"};
    
    // initialisation du vecteur des propensions de création des entités primaires
    //vector <double> c (entites_p.size());
    //fill(c.begin(), c.end(), 0.0);
    
    // initialisation du vecteur des propensions de destruction de toutes les entités
    //vector <double> d (entites.size());
    //fill(d.begin(), d.end(), 0.0);
    
    // Taux pour les propensions de création destruction
    //double taux_entree = 0.1;
    //double taux_sortie = 0.1;
    
    vector <double> r = {0.0, 0.0}; // initialisation pour les nombres aléatoires
    double atot, t, tau, somme;
    unsigned mu; // unsigned pour entiers non négatifs
    
    // Définition des vecteurs pour stocker les résultats de la dynamique
    
    vector <double> temps = {0.0};
    
    // Création des objets de la classe Entite
    
    vector<Entite> entites = {
        Entite ("A", 50, 0, 1, 0.1),
        Entite ("B", 50, 0, 0.1, 0.1),
        Entite ("AB", 1, 0),
        Entite ("ABA", 0, 0),
        Entite ("ABAB", 0, 0)
    };
    
    //declarer des pointeurs vers les entites
    Entite* A=&entites[0];
    Entite* B=&entites[1];
    Entite* AB=&entites[2];
    Entite* ABA=&entites[3];
    Entite* ABAB=&entites[4];
    
    vector<double> x ;
    
    for (const auto& e : entites){
        x.push_back(e.concentration);
    }
    
    vector <vector<double>> etats = {x};
    
    // Création des objets de la classe Réaction
    // Objets stockés dans un vecteur
    
    // On détaille les 2 sens des réactions
    // Reaction AB+A=ABA
    // Reaction ABA+B=ABAB
    // Reaction ABAB=2AB
    
    vector<Reaction> reactions = {
        //Reaction({ {entite.getEntityForName('AB'), 1}, {"A", 1}}, {{"ABA", 1}}, 0.0),
        Reaction({A,AB}, {ABA}, 0.0),
        Reaction({ABA, B}, {ABAB}, 0.0),
        Reaction({ABAB}, {AB, AB}, 0.0)
    };
    
    
    // initialisation du vecteur des propensions de réactions du cycle
    int nb_events = 2*reactions.size() + 2*entites.size();
    vector <double> a (nb_events);
    vector <vector<double>> propensions;
    
    // Algo Gillespie
    // Initialisation
    
    t = 0.0;
    
    while (t>=0 && t<tmax){
        
        // Calcul des propensions de création des entités primaires
        //for (int k=0; k < entites_p.size(); k++){
            //if (x[k]!=0){
                //c[k]=taux_entree*(1/x[k]); // propension à entrer inversement proportionnelle à la concentration
            //}
            //else {
                //c[k]=taux_entree
            //}
        //}
        
        // Calcul des propensions de destruction de toutes les entités
        //for (int k=0; k < entites.size(); k++){
            //d[k]=taux_sortie*x[k]; // propension à sortir proportionnelle à la concentration
        //}
        
        // Calcul des propensions de réaction (vitesses des réactions)
        // On distingue sens direct et sens inverse car on ne peut pas avoir de propensions négatives
        size_t i;
        for (i=0; i < reactions.size(); i++){ // car on cherche chaque objet dans la liste reactions
            //cout << i <<endl;
            a[2*i]=reactions[i].vitesse(true, V);
            a[2*i+1]=reactions[i].vitesse(false, V);
        }
        
        for (size_t j=0; j < entites.size(); j++){
            a[2*i+2*j]= entites[j].taux_creation; // on parcourt tout le tableau
            a[2*i+2*j+1]= entites[j].taux_destruction*entites[j].concentration/V;
            
        }
        
        atot=0.0;
        for(unsigned m=0; m < a.size(); ++m){
            atot += a[m]; // propension totale à réagir
        }
        
        // Génération de r1 et r2 dans une loi uniforme
        
        r.at(0)=(double)rand()/(double)RAND_MAX; // r1
        r.at(1)=(double)rand()/(double)RAND_MAX; // r2
        
        // Temps jusqu'à la prochaine simulation (distribution exponentielle) :
        
        tau = (-log(r[0])) / atot;
        t=t+tau;
        
        // Détermination de l'évènement :
        
        somme = 0.0;  // Initialisation de la somme
        
        for(unsigned j=0; j < a.size(); j++){
            somme = somme + a[j]; // on fait la somme mais ce qui compte c'est bien la longueur de l'intervalle que prend la réaction sur le segment total
            
            // puisque r tiré dans loi uniforme (0,1), en multipliant par atot, on "créé" l'intervalle (0,atot)
            if(somme >= atot*r[1] && somme - a[j] < atot*r[1]){ // Condition pour choisir la réaction
                mu=j; // mu représente l'indice de la réaction
                break;
            } // on sort de la boucle dès qu'on a trouvé la réaction à se produire
        }
        // Une fois la réaction choisie, on change le vecteur x du nombre d'entités
        
        if (mu < 2*reactions.size()){
            int dir=(mu%2==0)?1:-1;
            int idx = mu/2;
            
            for (const auto& reac : reactions[idx].reactifs){
                reac->concentration+=-dir;
                
            }
            
            for (const auto& prod : reactions[idx].produits){
                prod->concentration+=dir;
            }
            
        }
        else {
            int dir=(mu%2==0)?1:-1;// si pair, création sinon destruction
            int idx = mu - 2*reactions.size(); // on reprend à 0 après les réactions
            idx = idx/2; // on retrouve l'indice de l'entité
            entites[idx].concentration+=dir;
            cout << entites[idx].name <<","<<dir;
        }
        cout << endl;
        
        for (size_t i=0; i < entites.size(); i++){
            x[i]=entites[i].concentration;
            cout << x[i] <<",";
            
        }
        cout << endl;
        
        // Enregistrement de l'état :
        temps.push_back(t);
        etats.push_back(x);
        propensions.push_back(a);
        for (auto i : a){
            cout << i <<",";
            
        }
        cout << endl;
    }

    // Création d'un fichier csv
    
    ofstream file("gillespie.csv");
    
    file << "Temps,A,B,AB,ABA,ABAB\n"; // en-tête
    
    if (!file.is_open()) {
            cerr << "Failed to open file!" << endl;
            return 1;
        }
    
    vector<vector<double>> data ;
    
    for (size_t k=0; k< etats.size(); k++){
        
        vector <double> row;
        
        row.push_back(temps[k]);
        
        for (double x : etats[k]){
            row.push_back(x);
            
        }
        
        data.push_back(row);
        
    }
        
    for (const auto& row : data) {
            for (size_t i = 0; i < row.size(); ++i) {
                file << row[i];
                if (i != row.size() - 1) file << ",";
            }
                file << "\n";
        }
        
    
    file.close();
    cout << "CSV file created successfully." << endl;

    // Pour afficher des vecteurs, on doit faire des boucles (impossible de les afficher
    // directement)
    
    // Affichage Texte des résultats
    
    cout << "Temps : " << "Propensions : " << " Etats : \n ";
    
    
    for (size_t k=0; k < etats.size(); k++){
        cout << temps[k] << " ";
        
        cout << "{";
        for (double h : propensions[k]){
                    
                    cout << h << ",";
                }
                cout << "} ";
        
        cout << "{";

        for (double x : etats[k]){
            
            cout << x << ",";
        }
        cout << "} \n";
        
    }
    return 0;
}
