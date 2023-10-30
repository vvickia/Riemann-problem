#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>

#include "headers/de_allocate.h"
#include "headers/iof.h"
#include "headers/configurations.h"
#include "headers/newton.h"
#include "headers/riemann.h"

using namespace std;

/*
    "The Riemann problem on the decay of an arbitrary discontinuity"

    Write a program that implements an exact solution for the decay of an arbitrary gas-dynamic 
    discontinuity in an ideal gas with adiabatic exponent 𝛾 = 5/3. 
    Initially, the discontinuity is located at the origin.
    Test the program.
*/

int main ()
{
// Add test selection for single compilation

    string test = "";

    while (true)
    {
        cout << "Choose a test!\nPossible input options: input1.txt, input2.txt, input3.txt\n";
        cin >> test;
        cout << '\n';

        if (!((test == "input1.txt") || (test == "input2.txt") || (test == "input3.txt")))
        {
            cerr << "Invalid input!\n";
        }
        else
        {
            break;
        }
    }

// Parameter initialization:
//    rho_L (rho_R) = gas density on the left (on the right),
//    v_L   (v_R)   = gas velocity on the left (on the right),
//    p_L   (p_R)   = gas pressure on the left (on the right),
//    time = a time point under study.

    double rho_L, v_L, p_L, rho_R, v_R, p_R, time;

    initialization(test, rho_L, v_L, p_L, rho_R, v_R, p_R, time);
    
/*
    Since there is a gamma() in C++, let’s replace the letter γ with the third letter of the
    Phoenician alphabet 𐤂 (gimel) that generates it. Bertrand Russell (although I don't respect
    him) posited that the original pictogram is a conventionalized image of a camel and means
    this animal in the Phoenician.
*/
    const double gimel = 5.0 / 3.0; // Ratio of specific heats (adiabatic exponent)

// Sound velocities (on the left & on the right)

    double c_L = sqrt(gimel * p_L / rho_L);
    double c_R = sqrt(gimel * p_R / rho_R);

// Configuration selection

    string selected_configuration = configuration_selection(gimel, rho_L, v_L, p_L, c_L, rho_R, v_R, p_R, c_R);

/* 
    Solving a nonlinear equation
        v*_L(p*) = v*_R(p*)
    using Newton's iterative method
*/

    double unknown_p = Newtons_iterative_method (rho_L, v_L, p_L, c_L, rho_R, v_R, p_R, c_R, gimel, selected_configuration);

// Allocation of memory to dynamic variables

    int N = 1000;
    double x_L = -0.5;
    double x_R = 0.5;

    double * x = create_vector(N);
    
    double * velocity = create_vector(N);
    double * density = create_vector(N);
    double * pressure = create_vector(N);

// Solving the Riemann problem
    
    riemann_solve (rho_L, v_L, p_L, c_L, rho_R, v_R, p_R, c_R, unknown_p, gimel, time, \
                   velocity, density, pressure, x, N, x_L, x_R, selected_configuration);

// Shutdown

    save_results(test, N, x, pressure, density, velocity);

// Freeing allocated memory

    free_vector(x);
    free_vector(velocity);
    free_vector(density);
    free_vector(pressure);

    return 0;
}
