#include <iostream>
#include <vector>
#include <array>
#include "cell.hpp"
#include "mesh.hpp"
#include "unsteady_state_solution.hpp"
#include "TransientConduction.hpp"



    
int main (){
    std::vector <double> vnorth = {0}; // στοιχεία συνοριακής συνθήκης στο άνω μέρος (αδιαβατικά τοιχώματα)
    std::vector <double> vsouth = {0}; // στοιχεία συνοριακής συνθήκης στο κάτω μέρος (αδιαβατικά τοιχώματα)
    std::vector <double> veast = {T_inf,h}; // στοιχεία συνοριακής συνθήκης στο δεξί μέρος (συναγωγή με ρευστό θερμοκρασίας T_inf και συναγωγιμότητας h)
    std::vector <double> vwest = {T_inf,h}; // στοιχεία συνοριακής συνθήκης στο αριστερό μέρος (συναγωγή με ρευστό θερμοκρασίας T_inf και συναγωγιμότητας h)

    // επεσήμανση είδους συνοριακής συνθήκης σε κάθε μια από τις 4 κατευθύνσεις
    boundary north(ADIABATIC, vnorth); 
    boundary south(ADIABATIC, vsouth);
    boundary east(CONVECTION, veast);
    boundary west(CONVECTION, vwest);
    // γέμισμα του πίνακα με τις συνοριακές συνθήκες
    boundary myboundaries[4] = {north, south, east, west};

    //μελέτη περίπτωσης, πρακτικά 1D προβλήματος, αφού έχω αδιαβατικά τοιχώματα πάνω και κάτω
    //Αριθμός διαίρεσης πλέγματος στην αξονική κατεύθυνση
    int M = 50;

    // Δημιουργία κλάσης του προβλήματος, όπου:
    //T_0: η αρχική θερμοκρασία κάθε κυψέλης για χρόνο t = 0
    //1e-4: επιθυμητή τιμή του residual
    UnsteadyState dummyProblem1(2, M, T_0, myboundaries,1e-4);

    //Υπολογιστική λύση στο σημέιο x = L/2
    dummyProblem1.RungeKutta1st_point("test.dat");
    //Αναλυτική λύση στο σημέιο x = L/2
    theoritical_transient_conduction(L/2,"theoritical_results.dat");

    


    return 0;
}