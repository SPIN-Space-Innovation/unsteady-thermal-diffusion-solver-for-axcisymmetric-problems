#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include "cell.hpp"
#include "mesh.hpp"
#include "coef.hpp"


enum bound {CONDUCTION, CONVECTION, ADIABATIC};
enum val {Temperature, convection_coef};




class boundary{
    public:

    int category_; // κατηγορία συνοριακής συνθήκης (CONDUCTION, CONVECTION, ADIABATIC)
    std::vector <double> values_; // Διάνυσμα με όλα τα απαραίτητα μεγέθη της συνοριακής συνθήκης (θερμοκρασια, συντελεστής συναγωγιμότητας)
    
    boundary(int category, std::vector <double> values) // constructor
    :category_(category), values_(values){}

    boundary() = default;
    boundary(const boundary &other) = default;

    ~boundary() = default;
};

class UnsteadyState{
    public:
    int N_; // αριθμός διαίρεσης του πλέγματος στην x κατεύθυνση
    int M_; // αριθμός διαίρεσης του πλέγματος στην r κατεύθυνση
    mesh mesh_; // πλέγμα
    double dt_; // χρονικό βήμα
    double residual_;
    double error_; // μικρότερη δυνατή τιμή του residual, την επιλέγει ο χρήστης
    boundary boundaries_[4]; //οι 4 συνοριακές συνθήκες


    UnsteadyState(int N, int M, double T_initial, boundary bound[4], double error) // constructor
    :N_(N), M_(M), dt_(0.05*dt_min(maxInt(N,M))), error_(error)
    {
        mesh mymesh(N,M,T_initial);
        mesh_ = mymesh;
        for (int i = 0; i < 4; i++)
        {
            boundaries_[i] = bound[i]; // γέμισμα των συνοριακών συνθηκών με array 'bound[]' το οποίο φτιάχνεται μέσα στην main
        }
        Residual();
    }

    UnsteadyState() = default;
    UnsteadyState(const UnsteadyState &other) = default;
    ~UnsteadyState() = default;

    double ap(int i, int j){
        return mesh_.at(i,j).a_[P];
    }

    double as(int i, int j){
        return mesh_.at(i,j).a_[S];
    }

    double an(int i, int j){
        return mesh_.at(i,j).a_[N];
    }

    double ae(int i, int j){
        return mesh_.at(i,j).a_[E];
    }

    double aw(int i, int j){
        return mesh_.at(i,j).a_[W];
    }

    double Tp(int i, int j){
        return mesh_.at(i,j).T_;
    }

    double Tn(int i, int j){
        return mesh_.at(i+1,j).T_;
    }

    double Ts(int i, int j){
        return mesh_.at(i-1,j).T_;
    }

    double Te(int i, int j){
        return mesh_.at(i,j+1).T_;
    }

    double Tw(int i, int j){
        return mesh_.at(i,j-1).T_;
    }
    
    double x(int i, int j){
        return mesh_.at(i,j).x_[X];
    }
    
    double r(int i, int j){
        return mesh_.at(i,j).x_[R];
    }
    
    double dx(int i, int j){
        return mesh_.at(i,j).dx_[X];
    }

    double dr(int i, int j){
        return mesh_.at(i,j).dx_[R];
    }

    double dV(int i, int j){
        return r(i,j)*dx(i,j)*dr(i,j);
    }
    
    double A(int i, int j, int Direction){
        switch (Direction)
        {
        case N:
            return dx(i,j)*(r(i,j) + dr(i,j)/2);
            break;
        case S:
            return dx(i,j)*(r(i,j) - dr(i,j)/2);
            break;
        case E:
            return r(i,j)*dr(i,j);
            break;
        case W:
            return r(i,j)*dr(i,j);  
            break;
        default:
            break;
        }
    }
    
    double F_dot(int Boundary,int i, int j, int Direction){

        switch (Boundary)
        {
        case CONDUCTION:
            return (-mesh_.at(i,j).a_[Direction]*Tp(i,j) + 2*boundaries_[Direction].values_[Temperature]*mesh_.at(i,j).a_[Direction])*(diffusivity_steel/dV(i,j));
            break;
        case ADIABATIC:
            return (mesh_.at(i,j).a_[Direction]*Tp(i,j))*(diffusivity_steel/dV(i,j));
            break;
        case CONVECTION:
            return (diffusivity_steel/dV(i,j))*(mesh_.at(i,j).a_[Direction]*Tp(i,j) - boundaries_[Direction].values_[convection_coef]*A(i,j,Direction)*Tp(i,j)/K_steel
            + boundaries_[Direction].values_[Temperature]*boundaries_[Direction].values_[convection_coef]*A(i,j,Direction)/K_steel);
            break;
        default:
            break;
        }
    }

    double F(int i, int j){

        double return_value = 0;

        // ΚΑΤΩ-ΑΡΙΣΤΕΡΗ ΓΩΝΙΑ
        if(i == 0 && j == 0){
            return_value += (diffusivity_steel/dV(i,j)) * (ae(i,j)*Te(i,j) + an(i,j)*Tn(i,j) - ap(i,j)*Tp(i,j));

            switch (boundaries_[S].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,S);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,S);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,S);
                break;
            default:
                break;
            }
            switch (boundaries_[W].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,W);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,W);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,W);
                break;
            default:
                break;
            }

        }
        // ΚΑΤΩ
        if (i == 0 && j>0 && j<M_-1)
        {
            return_value += (diffusivity_steel/dV(i,j)) * (ae(i,j)*Te(i,j) + an(i,j)*Tn(i,j) + aw(i,j)*Tw(i,j) - ap(i,j)*Tp(i,j));
            switch (boundaries_[S].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,S);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,S);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,S);
                break;
            default:
                break;
            }
            
        }
        // ΚΑΤΩ ΔΕΞΙΑ ΓΩΝΙΑ
        if(i == 0 && j == M_-1){
            return_value += (diffusivity_steel/dV(i,j)) * (an(i,j)*Tn(i,j) + aw(i,j)*Tw(i,j) - ap(i,j)*Tp(i,j));
            switch (boundaries_[S].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,S);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,S);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,S);
                break;
            default:
                break;
            }
            switch (boundaries_[E].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,E);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,E);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,E);
                break;
            default:
                break;
            }
            
        }
        // ΔΕΞΙΑ ΠΛΕΥΡΑ
        if (j == M_-1 && i>0 && i<N_-1)
        {
            return_value += (diffusivity_steel/dV(i,j)) * (aw(i,j)*Tw(i,j) + as(i,j)*Ts(i,j) + an(i,j)*Tn(i,j) - ap(i,j)*Tp(i,j));
            switch (boundaries_[E].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,E);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,E);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,E);
                break;
            default:
                break;
            }
            
        }
        // ΠΑΝΩ ΔΕΞΙΑ ΓΩΝΙΑ
        if(i == N_-1 && j == M_-1){
            return_value += (diffusivity_steel/dV(i,j)) * (aw(i,j)*Tw(i,j) + as(i,j)*Ts(i,j) - ap(i,j)*Tp(i,j));
            switch (boundaries_[N].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,N);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,N);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,N);
                break;
            default:
                break;
            }
            switch (boundaries_[E].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,E);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,E);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,E);
                break;
            default:
                break;
            }
            
        }
        // ΠΑΝΩ ΠΛΕΥΡΑ
        if (i == N_-1 && j>0 && j<M_-1)
        {
            return_value += (diffusivity_steel/dV(i,j)) * (aw(i,j)*Tw(i,j) + as(i,j)*Ts(i,j) + ae(i,j)*Te(i,j) - ap(i,j)*Tp(i,j));
            switch (boundaries_[N].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,N);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,N);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,N);
                break;
            default:
                break;
            }
            
        }
        // ΠΑΝΩ ΑΡΙΣΤΕΡΗ ΓΩΝΙΑ
        if(i == N_-1 && j == 0){
            return_value += (diffusivity_steel/dV(i,j)) * (as(i,j)*Ts(i,j) + ae(i,j)*Te(i,j) - ap(i,j)*Tp(i,j));
            switch (boundaries_[N].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,N);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,N);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,N);
                break;
            default:
                break;
            }
            switch (boundaries_[W].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,W);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,W);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,W);
                break;
            default:
                break;
            }
            
        }
        // ΑΡΙΣΤΕΡΗ ΠΛΕΥΡΑ
        if (j == 0 && i>0 && i<N_-1)
        {
            return_value += (diffusivity_steel/dV(i,j)) * (an(i,j)*Tn(i,j) + as(i,j)*Ts(i,j) + ae(i,j)*Te(i,j) - ap(i,j)*Tp(i,j));
            switch (boundaries_[W].category_)
            {
            case CONDUCTION:
                return_value += F_dot(CONDUCTION,i,j,W);
                break;
            case ADIABATIC:
                return_value += F_dot(ADIABATIC,i,j,W);
                break;
            case CONVECTION:
                return_value += F_dot(CONVECTION,i,j,W);
                break;
            default:
                break;
            }
            
        } 
        // ΕΣΩΤΕΡΙΚΟ
        if(i >0 && i<N_-1 && j>0 && j<M_-1){
            return_value += (diffusivity_steel/dV(i,j)) * (aw(i,j)*Tw(i,j) + as(i,j)*Ts(i,j) + ae(i,j)*Te(i,j) + an(i,j)*Tn(i,j) - ap(i,j)*Tp(i,j));
        }
    
        return return_value;
    }

    double Laplacian(int i, int j){
        return F(i,j)*dV(i,j)/diffusivity_steel;
    } 
    
    // συνάρτηση που τυπώνει όλη την πορεία των θερμοκρασιών όλων των κυψελών του πλέγματος για κάθε χρονική στιγμή, σε αρχείο .dat
    void RungeKutta1st(std::string datFile){
        std::ofstream results;
        results.open(datFile);
        double timeMoment = 0;
        if(results.is_open()){
            results<<" t = "<<timeMoment<<"   residual = "<<residual_<<"\n";
            //std::cout<<" t = "<<timeMoment<<"   residual = "<<residual_<<"\n";
            for (int  i = 0; i < N_; i++)
            {
                for (int j = 0; j < M_; j++)
                {
                results<<mesh_.at(i,j).T_ <<"   ";
                //std::cout<<mesh_.at(i,j).T_ <<"   ";
                if(j == M_-1){
                    results<<"\n";
                    //std::cout<<"\n";
                }
                else{
                    continue;
                }
                }
            }
            while (residual_>error_)
            {            
                results<<"\n\n t = "<<timeMoment<<"   residual = "<<residual_<<"\n";
                //std::cout<<"\n\n t = "<<timeMoment<<"   residual = "<<residual_<<"\n";
                for (int  i = 0; i < N_; i++)
                {
                    for (int j = 0; j < M_; j++)
                    {
                        mesh_.at(i,j).T_ += dt_*F(i,j);
                        results<<mesh_.at(i,j).T_<<"   ";
                        //std::cout<<mesh_.at(i,j).T_ <<"   ";
                        if(j == M_-1){
                            results<<"\n";
                            //std::cout<<"\n";
                        }
                        else{continue;}
                    }
                }
                timeMoment += dt_;
                Residual();
            }
            
         }
         results.close();    
    }

    void RungeKutta1st_point(std::string datFile){
        std::ofstream results;
        results.open(datFile);
        double timeMoment = 0;
        int k = 1;
        if(results.is_open()){
            results<<timeMoment<<"  "<<mesh_.at(0,int(floor(double(M_)/2))).T_<<"\n";
            while (residual_>error_ && timeMoment<1e4)
            {            
                for (int  i = 0; i < N_; i++)
                {
                    for (int j = 0; j < M_; j++)
                    {
                        mesh_.at(i,j).T_ += dt_*F(i,j);
                    }
                }
                if(k==10000)
                {
                    results<<timeMoment<<"  "<<mesh_.at(0,int(M_/2)).T_<<"\n";
                    k=1;
                }
                else{
                    k++;
                }
                timeMoment += dt_;
                Residual();
            }

         }
         results.close();    
    }



    void Residual(){
        double residual = 0;
        for (int i = 0; i < N_; i++)
        {
            for (int j = 0; j < M_; j++)
            {
                residual += pow(abs(Laplacian(i,j)),maxInt(N_,M_));
            } 
        }
        residual_ = pow(residual,1/double(maxInt(N_,M_)));
    }
    
    
    
};




    