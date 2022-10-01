#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <vector>

enum Coefs {N, S, E, W, P, COEFS}; // integers που δηλώνουν κατεύθυνση στον χώρο (βόρεια, νότια, ανατολικά, δυτικά), κεντρική θέση, πλήθος συντελεστών των τελικών γραμμικών εξισώσεων 
enum Dimension {X, R, DIM}; // integers που δηλώνουν τις διαστάσεις - κατευθύνσεις x, r του προβλήματος
//enum Solver_Category {steady, unsteady};

using namespace std;

// κλάση που περιέχει τις πληροφορίες των κυψελών του πλέγματος
class cell{
public:

    array <double,DIM> x_{0} ; // διάνυσμα θέσης (x,r) του κέντρου της κυψέλης
    array <double,DIM> dx_{0}; // διάνυσμα που περιέχει τις διαστάσεις (μήκος, πλατος) μιας κυψέλης (δx,δr)
    array <double,COEFS> a_{0};// διάνυσμα που περιέχει τους συντελεστές της εξίσωσης που προκύπτει από κάθε κυψέλη
    double T_;                 // Τοπική θερμοκρασία

    cell() = default;

    // constructor for steady problem
    cell(const array <double,2> &x, const array <double,2> &dx)
    :x_(x), dx_(dx) //γέμισμα διανυσμάτων (x,r) και (δx,δr)
    {
        // γέμισμα διανύσματος συντελεστών
        a_[N] = (x_[R]+dx_[R]/2)*dx_[X]/dx_[R]; // συντελεστής της θερμοκρασίας της βόρειας κυψέλης
        a_[S] = (x_[R]-dx_[R]/2)*dx_[X]/dx_[R]; //συντελεστής της θερμοκρασίας της νότιας κυψέλης
        a_[E] = x_[R]*dx_[R]/dx_[X]; //συντελεστής της θερμοκρασίας της ανατολικής κυψέλης
        a_[W] = x_[R]*dx_[R]/dx_[X]; //συντελεστής της θερμοκρασίας της δυτικής κυψέλης
        a_[P] = a_[N] + a_[S] + a_[E] + a_[W]; //συντελεστής της θερμοκρασίας της κυψέλης
    }
    // Construstor for unsteady problem
    cell(const array <double,2> &x, const array <double,2> &dx, double initial_temperature)
    :x_(x), dx_(dx), T_(initial_temperature) //γέμισμα διανυσμάτων (x,r) και (δx,δr)
    {
        // γέμισμα διανύσματος συντελεστών
        a_[N] = (x_[R]+dx_[R]/2)*dx_[X]/dx_[R]; // συντελεστής της θερμοκρασίας της βόρειας κυψέλης
        a_[S] = (x_[R]-dx_[R]/2)*dx_[X]/dx_[R]; //συντελεστής της θερμοκρασίας της νότιας κυψέλης
        a_[E] = x_[R]*dx_[R]/dx_[X]; //συντελεστής της θερμοκρασίας της ανατολικής κυψέλης
        a_[W] = x_[R]*dx_[R]/dx_[X]; //συντελεστής της θερμοκρασίας της δυτικής κυψέλης
        a_[P] = a_[N] + a_[S] + a_[E] + a_[W]; //συντελεστής της θερμοκρασίας της κυψέλης

    }
    
    cell(const cell &other) = default;

    ~cell() = default;
};
