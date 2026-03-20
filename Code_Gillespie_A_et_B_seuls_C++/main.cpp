#include <iostream>
#include <sstream>
#include <fstream> // pour gérer les fichiers
#include "Entite.h"
#include <cmath>
#include <random>
#include <vector>
#include <map>
#include <string>
#include <list>


using namespace std;

int main(int argc,char* argv[]) { // pour mettre des arguments

    // File_entites to read
    vector<string> file_entites = {"entites - Feuille 2.csv"};

    // List to store all entites data
    vector<Entite> entites;


    for (const auto& filename : file_entites) {
    ifstream csv_file(filename);
    string line;

    // Skip the header line
    getline(csv_file, line);

    // Read rows with entites data
    while (getline(csv_file, line)) {
        stringstream ss(line);
        string name, number, free_energy, concentration_ext;

        // Read model and skip unwanted columns
        getline(ss, name, ',');
        getline(ss, number, ',');
        getline(ss, free_energy, ',');
        getline(ss, concentration_ext, ',');

        // Convert variables string to double
        double effectif = stod(number);
        double energie_libre = stod(free_energy);
        double con_ext = stod(concentration_ext);

        // Add entite to the vector
        entites.push_back(Entite (name, effectif, energie_libre, con_ext));
    }
}

    int nRuns = 1;
    vector<int> runs;
    vector <double> temps;
    vector <vector<double>> etats;
    vector<vector<vector<double>>> all_etats;
    vector<vector<double>> all_temps;
    vector<double> x ; // vecteur x pour stocker les concentrations des objets entites
    vector <vector<double>> propensions;
    
    vector<Entite> entites_init = entites;
    
    // Par défaut, rand() fixe la graine (même résultat à chaque simulation)
    // Donc srand pour ne pas fixer la graine
    
    srand(time(NULL));
    
    // Déclaration et définition des variables propres à l'algo Gillespie
    
    double V =500; // Volume total
    double p_renouvelé = 0.1; // part de volume renouvelé à l'entrée et à la sortie du système
    
    
    double tmax = 1000; // par défaut c'est ce temps maximal
    if (argc > 1){
        tmax=stod(argv[1]);
        
    }
    
    double atot, t, tau, somme;
    unsigned mu; // unsigned pour entiers non négatifs
        
        for (int i=0;i < nRuns; i++){
            
            // Initialisation des vecteurs pour les prochains runs
            
            entites = entites_init;
            
            x.clear();
            propensions.clear();
            etats.clear();
            
            vector <double> r = {0.0, 0.0}; // initialisation pour les nombres aléatoires
            
            temps = {0.0}; // initialisation du vecteur temps
            
            //  Vecteur x et vecteur etats exprimés en CONCENTRATION
            
            for (const auto& e : entites){
                x.push_back(e.effectif/V);
            }
            
            etats = {x}; // vecteur double états pour stocker tous les x
            
            
            // Initialisation du vecteur des propensions de réactions du cycle
            // Les propensions sont exprimées en NOMBRE D'ENTITES par unité de temps
            int nb_events =  2*entites.size();
            vector <double> a (nb_events);
            
            // Algo Gillespie
            // Initialisation
            
            t = 0.0;
            
            while (t>=0 && t<tmax){
                
                // Calcul des propensions de création(ou entrée) et de destruction(ou sortie)
                // exprimé en nombre d'entités par unité de temps
                
                for (size_t i=0; i < entites.size(); i++){
                    a[2*i]= entites[i].concentration_ext*V*p_renouvelé; // création
                    a[2*i+1]= V*p_renouvelé*entites[i].effectif/V; // destruction
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
                    
                    //  && somme - a[j] < atot*r[1]
                    
                    // puisque r tiré dans loi uniforme (0,1), en multipliant par atot, on "créé" l'intervalle (0,atot)
                    if(somme >= atot*r[1]){ // Condition pour choisir la réaction
                        mu=j; // mu représente l'indice de la réaction
                        break;
                    } // on sort de la boucle dès qu'on a trouvé la réaction à se produire
                }
                // Une fois la réaction choisie, on change le vecteur x du nombre d'entités
                
                int dir=(mu%2==0)?1:-1;// si pair, création sinon destruction
                int idx = mu/2; // on retrouve l'indice de l'entité
                entites[idx].effectif+=dir;
                
                for (size_t i=0; i < entites.size(); i++){
                    x[i]=entites[i].effectif/V;
                    
                }
                
                // Enregistrement de l'état :
                temps.push_back(t);
                etats.push_back(x);
                propensions.push_back(a);
            } // fin boucle while
            all_etats.push_back(etats);
            all_temps.push_back(temps);
            
        } // fin boucle Run
    
    //for (size_t run = 0; run < all_etats.size(); run++) {
        //const auto& etat_final = all_etats[run].back(); // on prend le dernier vecteur du vecteur de tous les états pour un run donné
        //if (etat_final[3] ==0 && etat_final[8] !=0) {
            //freq_fix_mut++;
            //temps_fix+= all_temps[run][];
        //}
    //}

        // calcul variance entre trajectoires

    // Création d'un fichier csv
    
    ofstream file("gillespie.csv");

    
    file << "Run,Temps,A,B \n"; // en-tête
    
    if (!file.is_open()) {
            cerr << "Failed to open file!" << endl;
            return 1;
        }
    
    //vector<vector<double>> data ;
    
    //for (size_t k=0; k< etats.size(); k++){
        //vector <double> row;
        //row.push_back(temps[k]);
        //for (double x : etats[k]){
            //row.push_back(x);
        //}
        //data.push_back(row);
    //}
        
    //for (const auto& row : data) {
            //for (size_t i = 0; i < row.size(); ++i) {
                //file << row[i];
                //if (i != row.size() - 1) file << ",";
            //}
                //file << "\n";
        //}
    
    for (size_t run = 0; run < all_etats.size(); run++){
        for (size_t k = 0; k < all_etats[run].size(); k++){
            file << run << "," << all_temps[run][k];
            for (double v : all_etats[run][k]){
                file << "," << v;
            }
            file << "\n";
        }
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
