#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include "cell.hpp"
#include "coef.hpp"


// κλάση που περιέχει την πληροφορία όλου του πλέγματος
class mesh{
    public:
    int N_; // αριθμός διαίρεσης του πλέγματος
    int M_; // αριθμός διαίρεσης του πλέγματος

    vector <cell> cells_; // διάνυσμα που περιέχει όλα τα κελία του πλέγματος (ουσιαστικά το σύνολο των πληροφοριών για όλο το πλέγμα)
// συνάρτηση που επιστρέφει το κελί με τα επιθυμητά χωρικά indexes
    // cell const &at(int i, int j) const{
    //     return cells_.at(j + i*N_);
    // }
    cell &at(int i, int j) {
        return cells_.at(j + i*M_);
    }
// συνάρτηση που επιστρέφει το γραμμικό index (στο διάνυσμα cells_) της επιθυμητής κυψελής, δεδομένων των χωρικών indexes, i->r, j->x 
    int index(int i, int j){
        return j + i*M_;
    }

    mesh() = default;

// constructor της κλάσης mesh για steady προβλημα
    mesh(int N, int M)
    :N_(N), M_(M), cells_(N*M)
    {
        const array <double,DIM> dx{L/M, t/N}; // σταθερό διάνυσμα (δx,δr) ορισμένο από τις γεωμετρικές διαστάσεις του coef.hpp
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                cell a(array <double,DIM> {0.5*dx[X] + j*dx[X],   D/2 - t + 0.5*dx[R] + i*dx[R]},    dx);  // διάνυσμα θέσης κάθε κυψέλης ορισμένο από τις γεωμετρικές διαστάσεις του coef.hpp
                cells_.at(index(i,j)) = a;
            }
        }
    }

    // constructor της κλάσης mesh για unsteady προβλημα
    mesh(int N, int M, double T_initial)
    :N_(N), M_(M), cells_(N*M)
    {
        const array <double,DIM> dx{L/M, t/N}; // σταθερό διάνυσμα (δx,δr) ορισμένο από τις γεωμετρικές διαστάσεις του coef.hpp
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                cell a(array <double,DIM> {0.5*dx[X] + j*dx[X],   D/2 - t + 0.5*dx[R] + i*dx[R]},    dx,    T_initial);  // διάνυσμα θέσης κάθε κυψέλης ορισμένο από τις γεωμετρικές διαστάσεις του coef.hpp
                at(i,j) = a;
            }
        }
    }
    
    mesh(const mesh &other) = default;
    ~mesh() = default;
};

double dt_min(int N){
    return density_steel*C_p_steel*L*L*t*t/(2*K_steel*(L*L+t*t)*N*N);
};

int maxInt(int a, int b){
    if(a>b){
        return a;
    }
    else{
        return b;
    }
};