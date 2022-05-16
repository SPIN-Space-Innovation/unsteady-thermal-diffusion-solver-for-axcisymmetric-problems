#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <Eigen/Dense>
#include "cell.hpp"
#include "coef.hpp"


// κλάση που περιέχει την πληροφορία όλου του πλέγματος
class mesh{
    public:
    int N_; // αριθμός διαίρεσης του πλέγματος
    vector <cell> cells_; // διάνυσμα που περιέχει όλα τα κελία του πλέγματος (ουσιαστικά το σύνολο των πληροφοριών για όλο το πλέγμα)

// συνάρτηση που επιστρέφει το κελί με τα επιθυμητά χωρικά indexes
   const cell &at(int i, int j) const{
        return cells_.at(j + i*N_);
    }
// συνάρτηση που επιστρέφει το γραμμικό index (στο διάνυσμα cells_) της επιθυμητής κυψελής, δεδομένων των χωρικών indexes, i->r, j->x 
    int index(int i, int j){
        return j + i*N_;
    }

    mesh() = default;

// constructor της κλάσης mesh
    mesh(int N)
    :N_(N), cells_(N*N)
    {
        const array <double,DIM> dx{L/N, t/N}; // σταθερό διάνυσμα (δx,δr) ορισμένο από τις γεωμετρικές διαστάσεις του coef.hpp
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                cell a(array <double,DIM> {0.5*dx[X] + j*dx[X],   D/2 - t + 0.5*dx[R] + i*dx[R]},    dx);  // διάνυσμα θέσης κάθε κυψέλης ορισμένο από τις γεωμετρικές διαστάσεις του coef.hpp
                cells_.at(index(i,j)) = a;
            }
        }
    }
    
    mesh(const mesh &other) = default;
    ~mesh() = default;
};