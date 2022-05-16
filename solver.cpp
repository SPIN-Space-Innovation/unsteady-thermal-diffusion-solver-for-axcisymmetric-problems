#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <Eigen/Dense>
#include "cell.hpp"
#include "mesh.hpp"
#include "coef.hpp"
#include "boundary.hpp"


//g++ -I C:\Users\VS0121\eigen-3.4.0 solver.cpp


int main() {

//εισαγωγή προβλήματος, με την επιθυμητή ακρίβεια για την λύση
boundary problem(100);
// εισαγωγή 4 συνοριακών συνθηκών σε κάθε κατεύθυνση
problem.convection(N,T_inf,K_steel,h);
problem.conduction(S,T_0);
problem.adiabatic(W);
problem.adiabatic(E);
// επίλυση του προβλήματος και αποθήκευση αποτελεσμάτων στο επιθυμητό αρχείο
problem.solver("results.txt");

return 0;
}

