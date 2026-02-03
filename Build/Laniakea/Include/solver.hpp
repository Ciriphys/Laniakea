#pragma once 

#include <fstream>
#include <vector>

#include "fftw3.h"
#include "vec3.hpp"

#define SQR(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))
#define DIST(x,y,z) std::sqrt((SQR(x)) + (SQR(y)) + (SQR(z)))
#define DIST2(x,y,z) ((SQR(x)) + (SQR(y)) + (SQR(z)))

struct poisson_schrodinger_state {
    fftw_complex psi;
    double phi;
};

struct soliton {
    soliton() = default;
    soliton(double mass, vec3 pos, vec3 vel, double phase);

    double mass;
    vec3 pos, vel;
    double phase;

    double alpha() const;
};

struct potential {
    potential(double mass = 0.0, std::string type = "", vec3 center = {});
    
    double compute(vec3 pos);

    double mass;
    vec3 center;
    std::string type;
};

struct settings {
    settings(int grid_size = -1.0, 
             double box_size = -1.0,   
             double max_time = -1.0,
             double vmin = 0.0,
             double vmax = 50.0,
             int threads = -1,
             int framerate = 60,
             int n_sol = 0,
             int n_pot = 0,
             bool ascii = false,
             bool infer = true,
             std::string soliton_file = "soliton.dat",
             std::string wisdom_file = "wisdom.fftw",
             std::string renderer = "ffmpeg",
             std::string video_filename = "output.mkv", 
             std::string ffmpeg_cmd = "-loglevel quiet",
             std::vector<soliton> solitons = {},
             std::vector<potential> potentials = {}
             );

    int grid_size; 
    double box_size;
    double max_time;
    double vmin;
    double vmax;
    int framerate;
    int threads;
    int n_sol;
    int n_pot;
    bool ascii;
    bool infer;
    
    std::string soliton_file;
    std::string wisdom_file;
    std::string renderer;
    std::string video_filename;
    std::string ffmpeg_cmd;

    std::vector<soliton> solitons;
    std::vector<potential> potentials;
};

//969

void expi(double value, fftw_complex result);
void mult(fftw_complex a, fftw_complex b, fftw_complex result);
double abs2(const fftw_complex& z);

std::string strip(const std::string& s, const std::string& chars=" ");

double init_grid(poisson_schrodinger_state* grid, double* fourier_grid, settings config, double& vmin, double& vmax, std::ofstream& en_file);

void pseudo_spectral_solver(poisson_schrodinger_state* grid, double* fourier_grid, settings config, double& vmin, double& vmax, std::ofstream& en_file, const double& E0);
void generate_fourier_grid(double* fourier_grid, int N, double L);

#ifdef FFMPEG 
    void print_ascii_slice(const poisson_schrodinger_state* grid, int frame, double vmin, double vmax);
#endif 