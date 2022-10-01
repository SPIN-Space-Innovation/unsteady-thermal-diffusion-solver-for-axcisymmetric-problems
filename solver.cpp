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
#include "steady_state_solution.hpp"
#include "unsteady_state_solution.hpp"


//g++ -I C:\Users\VS0121\eigen-3.4.0 solver.cpp


int main() {

//εισαγωγή προβλήματος, με την επιθυμητή ακρίβεια για την λύση
SteadyState problem(20);
// εισαγωγή 4 συνοριακών συνθηκών σε κάθε κατεύθυνση
problem.convection(N,70,K_steel,h);
problem.convection(S,20,K_steel,h);
problem.adiabatic(W);
problem.adiabatic(E);
//problem.print_A();
// επίλυση του προβλήματος και αποθήκευση αποτελεσμάτων στο επιθυμητό αρχείο
problem.solver("results.dat");
//problem.solver_GMRES("results_gmres.txt",1e-2);

return 0;
}

