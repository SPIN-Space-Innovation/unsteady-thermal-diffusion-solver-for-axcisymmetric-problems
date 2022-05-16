#pragma once
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <vector>
#include <Eigen/Dense>



double u_inf=150, // ταχύτητα ρευστού
        D = 100e-3, // εξωτερική διάμετρος σωλήνα
        T_inf = 0 , // θερμοκρασία κινούμενου ρευστού
        T_0 = 80 , // θερμοκρασία εσωτερικά του σωλήνα
        K_steel = 45, // συντελεστής αγωγιμότητας υλικού σωλήνα
        nu = 15.89e-6, // συνεκτικότητα ρευστού
        Pr = 0.707, // αριθμός Prandl του ρευστού
        k = 26.3e-3, // αγωγιμότητα ρευστού
        L = 0.4, // μήκος σωλήνα
        t = 45e-3; // πάχος σωλήνα

double Re = u_inf * D / nu; //αριθμός Reynolds του κινούμενου ρευστού
double Nu = 0.3 + 0.62 * pow(Re, 0.5) * pow(Pr, 1 / 3) * pow(1 + pow(Re / 282000, 5 / 8), 0.8) *
                  pow(1 + pow(0.4 / Pr, 0.25), -0.25); // αριθμός Nusselt του ρευστού
double h = Nu * k / D; // συντελεστής συναγωγιμότητας ρευστού

