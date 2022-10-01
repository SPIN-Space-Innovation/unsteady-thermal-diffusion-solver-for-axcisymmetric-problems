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




double z(int i){
    double z0 = 0.7 + 3*(i-1);
    double err = 1;
    double BI = h*L/(2*K_steel);
    while (abs(err)>1e-6)
    {
        err = (z0*tan(z0)-BI)/(tan(z0)+z0/pow(cos(z0),2));
        z0 = z0 -(z0*tan(z0)-BI)/(tan(z0)+z0/pow(cos(z0),2));
    }
    return z0;
    
};


void theoritical_transient_conduction(double x, std::string datFile){
    double Bi = h*L/(2*K_steel);
    std::cout<<"Bi = "<<Bi<<"\n";
    
    double Fo = 4*diffusivity_steel*t/pow(L,2);
    cout<<"Fo = "<<Fo<<"\n";
    
    double T[100];
    double dt = 100;
    double time = 0;
    std::vector<double> Zi(10+1,0);
    double Ci;
    for (int i = 1; i <= 10; i++)
    {
        Zi.at(i) = z(i);
    }

    std::ofstream results;
    results.open(datFile);
    if(results.is_open()){
        for (int i = 0; i < 1000; i++)
        {
            time = i*dt;
            Fo = 4*diffusivity_steel*time/pow(L,2);
            for (int m = 1; m <= 10; m++)  /////////////////                   //////////
            {
                Ci = 4*sin(Zi[m])/(2*Zi[m] + sin(2*Zi[m]));
                T[i] += Ci*exp(-pow(Zi[m],2)*Fo)*cos(Zi[m]*(2*x-L)/L);
            }
            T[i] = T[i]*(T_0-T_inf) + T_inf;
            
            if(abs(T[i]) < 1e3){
                results<<time<<"  "<<T[i]<<"\n";
            }
        }
    }
    
    
};
