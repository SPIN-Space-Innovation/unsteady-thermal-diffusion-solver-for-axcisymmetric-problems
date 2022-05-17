#pragma once
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
#include <Eigen/Dense> // Αναγκαία η εγκατάσταση της opensource βιβλιοθήλης Eigen (χρήσιμη για γραμμική άλγεβρα)


using namespace Eigen;


enum bound {CONDUCTION, CONVECTION, ADIABATIC};

// κλάση ορισμού του προβλήματος
class boundary{
    public:

    int N_; // αριθμός διαίρεσης του πλέγματος
    mesh mesh_; // πλέγμα
    MatrixXd A_; // πίνακας συντελεστών [Α](Ν*Ν,Ν*Ν) (συστήματος [Α]{T} = {Β})
    VectorXd B_; // διάνυσμα Β (Ν*Ν,1) (συστήματος [Α]{T} = {Β})

    //ΑΡΧΙΚΟ ΓΕΜΙΣΜΑ ΠΙΝΑΚΑ [Α]
    void fill(){
        
        int index;
        // Αρχικός μηδενισμός όλων των στοιχείων Αij, Βi
        for (int i = 0; i < N_*N_; i++)
        {
            B_(i) = 0;
            for (int j = 0; j < N_*N_; j++)
            {
                A_(i,j) = 0;
            }
            
        }
        
        for (int i = 0; i < N_; i++)
        {
            for (int j = 0; j < N_; j++)
            {
                index = j+i*N_;
                if (i == 0 && j == 0) //ΝΟΤΙΟΔΥΤΙΚΗ-ΑΡΙΣΤΕΡΗ ΚΑΤΩ ΓΩΝΙΑ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index + N_) = -mesh_.at(i,j).a_[N];
                    A_(index,index + 1) = -mesh_.at(i,j).a_[E];
                }
                if (i == 0 && j > 0 && j < N_-1) //ΝΟΤΙΟ-ΚΑΤΩ ΟΡΙΟ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index + N_) = -mesh_.at(i,j).a_[N];
                    A_(index,index + 1) = -mesh_.at(i,j).a_[E];
                    A_(index,index - 1) = -mesh_.at(i,j).a_[W];
                }
                if (i > 0 && i < N_-1 && j == 0) //ΔΥΤΙΚΟ-ΑΡΙΣΤΕΡΟ ΟΡΙΟ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index + N_) = -mesh_.at(i,j).a_[N];
                    A_(index,index - N_) = -mesh_.at(i,j).a_[S];
                    A_(index,index + 1) = -mesh_.at(i,j).a_[E];
                }
                if (i == N_-1 && j == 0) //ΒΟΡΕΙΟΑΝΑΤΟΛΙΚΗ-ΑΡΙΣΤΕΡΗ ΑΝΩ ΓΩΝΙΑ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index - N_) = -mesh_.at(i,j).a_[S];
                    A_(index,index + 1) = -mesh_.at(i,j).a_[E];
                }
                if (j > 0 && j < N_-1 && i == N_-1) //ΒΟΡΕΙΟ-ΑΝΩ ΟΡΙΟ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index - N_) = -mesh_.at(i,j).a_[S];
                    A_(index,index + 1) = -mesh_.at(i,j).a_[E];
                    A_(index,index - 1) = -mesh_.at(i,j).a_[W];
                }
                if (i == N_-1 && j == N_-1) //ΒΟΡΕΙΟΑΝΑΤΟΛΙΚΗ-ΔΕΞΙΑ ΑΝΩ ΓΩΝΙΑ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index - N_) = -mesh_.at(i,j).a_[S];
                    A_(index,index - 1) = -mesh_.at(i,j).a_[W];
                }
                if (i > 0 && i < N_-1 && j == N_-1) //ΑΝΑΤΟΛΙΚΟ-ΔΕΞΙ ΟΡΙΟ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index + N_) = -mesh_.at(i,j).a_[N];
                    A_(index,index - N_) = -mesh_.at(i,j).a_[S];
                    A_(index,index - 1) = -mesh_.at(i,j).a_[W];
                }
                if (i == 0 && j == N_-1) //ΝΟΤΙΟΑΝΑΤΟΛΙΚΗ-ΚΑΤΩ ΔΕΞΙΑ ΓΩΝΙΑ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index + N_) = -mesh_.at(i,j).a_[N];
                    A_(index,index - 1) = -mesh_.at(i,j).a_[W];
                }
                if (i > 0 && i < N_-1 && j > 0 && j < N_-1)//ΕΣΩΤΕΡΙΚΟ ΓΕΜΙΣΜΑ
                {
                    A_(index,index) = mesh_.at(i,j).a_[P];
                    A_(index,index + N_) = -mesh_.at(i,j).a_[N];
                    A_(index,index - N_) = -mesh_.at(i,j).a_[S];
                    A_(index,index + 1) = -mesh_.at(i,j).a_[E];
                    A_(index,index - 1) = -mesh_.at(i,j).a_[W];
                }
            }
        }
    }

    // constructor
    boundary(int N)
    :N_(N)
    {
        mesh mymesh(N);
        mesh_ = mymesh;

        A_ = MatrixXd::Zero(N*N,N*N);
        B_ = VectorXd::Zero(N*N);
        fill();
    }

    // συνάρτηση για την εισαγωγή συνοριακή συνθήκης αγωγής θερμότητας (σταθερή θερμοκρασία σε μία επιφάνεια όπως για παράδειγμα στις φλεγόμενες πλευρές των grains)
    void conduction(int direction_, double Temperature){

        int i,j,index;
        switch (direction_)
        {
            case S: // ΚΑΤΩ ΟΡΙΟ 
            for (i = 0,j = 0; j < N_; j++)
            {
                index = j+i*N_;
                A_(index,index) += mesh_.at(i,j).a_[S];

                B_(index) += 2*Temperature*mesh_.at(i,j).a_[S];
            }
            break;
            case N: //ΑΝΩ ΟΡΙΟ
            for (i = N_-1,j = 0; j < N_; j++)
            {
                index = j+i*N_;
                A_(index,index) += mesh_.at(i,j).a_[N];

                B_(index) += 2*Temperature*mesh_.at(i,j).a_[N];
            }
            break;
            case E: //ΔΕΞΙ ΟΡΙΟ
            for (j = N_-1,i = 0; i < N_; i++)
            {
                index = j+i*N_;
                A_(index,index) += mesh_.at(i,j).a_[E];

                B_(index) += 2*Temperature*mesh_.at(i,j).a_[E];
            }
            break;
            case W: //ΑΡΙΣΤΕΡΟ ΟΡΙΟ
            for (i = 0,j = 0; i < N_; i++)
            {
                index = j+i*N_;
                A_(index,index) += mesh_.at(i,j).a_[W];

                B_(index) = 2*Temperature*mesh_.at(i,j).a_[W];
            }
            break;
            default:
            cout<<"problem at conduction"<<"\n";
        }
    }
    // συνάρτηση για την εισαγωγή συνοριακής συνθήκης αδιαβατικού τοιχώματος (πρακτικά απουσία ροής ενέργειας, όπως για παρλάδειγμα στις περιπτώσεις τοιχωμάτων μεγάλου πάχους)
    void adiabatic(int direction_){

        int i,j,index;
        switch (direction_)
        {
            case S: //ΚΑΤΩ ΟΡΙΟ
            for (i = 0,j = 0; j < N_; j++)
            {
                index = j+i*N_;
                A_(index,index) -= mesh_.at(i,j).a_[S];
           }
            break;
            case N: //ΑΝΩ ΟΡΙΟ
            for (i = N_-1,j = 0; j < N_; j++)
            {
                index = j+i*N_;
                A_(index,index) -= mesh_.at(i,j).a_[N];
            }
            break;
            case E: //ΔΕΞΙ ΟΡΙΟ
            for (j = N_-1,i = 0; i < N_; i++)
            {
                index = j+i*N_;
                A_(index,index) -= mesh_.at(i,j).a_[E];
              }
            break;
            case W: //ΑΡΙΣΤΕΡΟ ΟΡΙΟ
            for (i = 0,j = 0; i < N_; i++)
            {
                index = j+i*N_;
                A_(index,index) -= mesh_.at(i,j).a_[W];
           }
            break;
            default:
            cout<<"problem at adiabatic"<<"   direction = "<<direction_<<"\n";
        }
    }
    // συνάρτηση για την εισαγωγή συνοριακή συνθήκης συναγωγής (μεταφορά θερμότητας λόγω της ροής κ΄ρου ή ζεστου ρευστού), Κ == συντελεστής αγωγιμότητας υλικού πλέγματος, Η == συντελεσής συναγωγιμότητας κινούμενου ρευστού
    void convection(int direction_, double Temperature, double K,  double H)
    {
        int i,j,index;
        double A;
        switch (direction_)
        {
            case S: // ΚΑΤΩ ΟΡΙΟ
            for (i = 0,j = 0; j < N_; j++)
            {
                index = j+i*N_;
                A = (mesh_.at(i,j).x_[R] - mesh_.at(i,j).dx_[R]/2)*mesh_.at(i,j).dx_[X];

                A_(index,index) -= (mesh_.at(i,j).a_[S] - H*A/K);

                B_(index) += Temperature*H*A/K;
            }
            break;
            case N: //ΑΝΩ ΟΡΙΟ
            for (i = N_-1,j = 0; j < N_; j++)
            {
                index = j+i*N_;
                A = (mesh_.at(i,j).x_[R] + mesh_.at(i,j).dx_[R]/2)*mesh_.at(i,j).dx_[X];

                A_(index,index) -= (mesh_.at(i,j).a_[N] - H*A/K);
                
                B_(index) += Temperature*H*A/K;
            }
            break;
            case E: //ΔΕΞΙ ΟΡΙΟ
            for (j = N_-1,i = 0; i < N_; i++)
            {
                index = j+i*N_;
                A = mesh_.at(i,j).x_[R]*mesh_.at(i,j).dx_[R];

                A_(index,index) -= (mesh_.at(i,j).a_[E] - H*A/K);
                
                B_(index) += Temperature*H*A/K;
            }
            break;
            case W: //ΑΡΙΣΤΕΡΟ ΟΡΙΟ
            for (i = 0,j = 0; i < N_; i++)
            {
                index = j+i*N_;
                A = mesh_.at(i,j).x_[R]*mesh_.at(i,j).dx_[R];

                A_(index,index) -= (mesh_.at(i,j).a_[W] - H*A/K);
                
                B_(index) += Temperature*H*A/K;
            }
            break;
            default:
            cout<<"problem at convection"<<"\n";
        }
    }
    // συνάρτηση που επιλύει το σύστημα [Α]Τ = Β και αποθηκεύει τα αποτελέσματα στο επιθυμητό αρχείο
    void solver(string textname){
        VectorXd T_ = A_.fullPivLu().solve(B_); // Mέθοδος επίλυσης του συστήματος [Α]{Τ} = {Β}
        ofstream Temp;
        Temp.open(textname);
        for (int i = 0; i < N_; ++i) {
            for (int j = 0; j < N_; ++j) {
                Temp<<mesh_.at(i,j).x_[X]<<"   "<<mesh_.at(i,j).x_[R]<<"   "<<T_(j+i*N_)<<"\n";
                if (j == N_-1){
                    Temp<<"\n";
                }
            }
        }
        Temp.close();
    }

    
};


