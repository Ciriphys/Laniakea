#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>

#include "integrator.hpp"
#include "vec4.hpp"

int depth = 0;

double shoot(double c, double dx, int dim) {
    vec4 state = { 1.0, 0.0, c, 0.0 };
   
    bool isnan_flag = false;

    for(int i = 1; i < dim && !isnan_flag; i++) {
        state = runge_kutta_iv(state, i * dx, dx);
        if(std::isnan(state.f) || std::abs(state.f) > 1e6) 
            isnan_flag = true;
    }

    return state.f;
}

vec4 shoot_on_file(double c, double dx, int dim, std::ofstream& file, double* mass) {
    vec4 state = { 1.0, 0.0, c, 0.0 };

    file << "# r, f(r), phi(r)\n"; // I could precompute psi(x,0) with alpha, pos, vel.
    file << 0.0 << ", " << state.f << ", " << state.phi << "\n";

    for(int i = 1; i < dim; i++) {
        state = runge_kutta_iv(state, i * dx, dx);
        (*mass) += 4 * M_PI * SQR(state.f) * SQR(i * dx) * dx;
        file << i*dx << ", " << state.f << ", " << state.phi << "\n";
    }

    return state;
}

double shooting_estimate(double c_low, double c_high, double f_low, double f_high, double dx, int dim, double tol = 1e-6) {
    double c_mid = 0.5 * (c_low + c_high);
    double f_mid = shoot(c_mid, dx, dim);
    double abs_f_mid = std::abs(f_mid);

    std::cout << "|f_mid| = " << abs_f_mid << " (depth: " << depth << ") " << std::endl;
    if(++depth >= 100 || abs_f_mid < tol) { return c_mid; } 

    if(f_low * f_mid < 0) {
        return shooting_estimate(c_low, c_mid, f_low, f_mid, dx, dim, tol);
    } else {
        return shooting_estimate(c_mid, c_high, f_mid, f_high, dx, dim, tol);
    }   
}

int main(int argc, char** argv) {
    if(argc < 4) { 
        std::cout << "usage: <filename:string> <grid_size:int> <box_size:double>\n";
        return -1;
    }

    std::string filename = argv[1];
    int grid_size = std::atoi(argv[2]);
    double box_size = std::stod(argv[3]);

    double dx = box_size / (double)grid_size;

    std::ofstream file(filename);
    if(!file) { std::cout << "Error opening file!" << std::endl; return -1; }

    //  TODO : Add soliton data and save soliton data into a separate file: 
        // --- SOLITON 1 ---
        // mass : xxx 
        // pos : (xx, yy, zz)
        // vel : (vx, vy, vz)

    //  double alpha = 3.883;

    int dim = (int)(box_size / dx);
    double mass = 0.0;
    
    double c_low = -2.5, c_high = 0.0;
    double f_low = shoot(c_low, dx, dim);
    double f_high = shoot(c_high, dx, dim);

    if(f_low * f_high > 0) {
        std::cerr << "f_low * f_high: " << f_low * f_high << " No zero has been found. Try extending the resolution! (box_size/grid_size)\n";
        return -1;
    }

    double c = shooting_estimate(c_low, c_high, f_low, f_high, dx, dim, TOL);
    std::cout << "Shooting done, phi(0) = " << std::setprecision(12) << c << std::endl; 
    vec4 state = shoot_on_file(c, dx, dim, file, &mass);

    file.close();
    double beta = state.phi + box_size * 0.5 * state.dphi;
    std::cout << "Estimated beta: " << beta << "\n";
    std::cout << "Estimated mass: " << mass << "\n";

    return 0;
} 