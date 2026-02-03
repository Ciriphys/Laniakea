#include <functional>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include "solver.hpp"
#include "fftw3.h"

#ifdef FFMPEG 
    #include "ffmpeg.hpp"
#endif 

// returns exp(i * value)
void expi(double value, fftw_complex result) {
    result[0] = std::cos(value);  
    result[1] = std::sin(value);  
}

void mult(fftw_complex a, fftw_complex b, fftw_complex result) {
    result[0] = a[0] * b[0] - a[1] * b[1];
    result[1] = a[0] * b[1] + a[1] * b[0];
}

double abs2(const fftw_complex& z) {
    return z[0] * z[0] + z[1] * z[1];
}

double dist(double x1, double y1, double z1, double x2 = 0.0, double y2 = 0.0, double z2 = 0.0) {
    return std::sqrt(SQR(x1 - x2) + SQR(y1 - y2) + SQR(z1 - z2));
}

potential::potential(double mass, std::string type, vec3 center)
     : mass(mass), type(type), center(center) {}

double potential::compute(vec3 pos) {
    double radius  = DIST (pos.x - center.x, pos.y - center.y, pos.z - center.z);
    double radius2 = DIST2(pos.x - center.x, pos.y - center.y, pos.z - center.z);

    double cent_fac = (type == "central") * -mass / radius;
    double quad_fac = (type == "quadratic") * 0.5 * radius2;

    return (radius > 1e-12 ? (cent_fac + quad_fac) : 0.0);
}

settings::settings(int grid_size, 
                   double box_size,   
                   double max_time,
                   double vmin,
                   double vmax,
                   int threads,
                   int framerate,
                   int n_sol,
                   int n_pot,
                   bool ascii,
                   bool infer,
                   std::string soliton_file,
                   std::string wisdom_file,
                   std::string renderer,
                   std::string video_filename,
                   std::string ffmpeg_cmd,
                   std::vector<soliton> solitons,
                   std::vector<potential> potentials) : 
                   grid_size(grid_size),
                   box_size(box_size),
                   max_time(max_time),
                   vmin(vmin),
                   vmax(vmax),
                   threads(threads),
                   framerate(framerate),
                   n_sol(n_sol),
                   n_pot(n_pot),
                   ascii(ascii),
                   infer(infer),
                   soliton_file(soliton_file),
                   wisdom_file(wisdom_file),
                   video_filename(video_filename),
                   ffmpeg_cmd(ffmpeg_cmd),
                   renderer(renderer),
                   solitons(solitons),
                   potentials(potentials) {}

double psi_norm2(const poisson_schrodinger_state* grid, int NNN) {
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < NNN; ++i) {
        sum += abs2(grid[i].psi);
    }
    return sum;
}


double psi_r2_avg(const poisson_schrodinger_state* grid, int L, int N) {
    double r2_sum = 0.0, norm = 0.0;
    const double dx = L / N;

    #pragma omp parallel for collapse(3) reduction(+:r2_sum, norm)
    for (int x = 0; x < N; x++) {
        double dx_pos = (x - (N - 1) / 2.0) * dx;
        for (int y = 0; y < N; y++) {
            double dy_pos = (y - (N - 1) / 2.0) * dx;
            for (int z = 0; z < N; z++) {
                double dz_pos = (z - (N - 1) / 2.0) * dx;
                int idx = z + y * N + x * N * N;

                double rho = abs2(grid[idx].psi);
                double r2 = dx_pos * dx_pos + dy_pos * dy_pos + dz_pos * dz_pos;

                r2_sum += r2 * rho;
                norm += rho;
            }
        }
    }

    return r2_sum / norm;
}

void debug_count(fftw_complex* arr, std::string msg, int NNN) {
    int zero_count = 0, nan_count = 0;

    for (int i = 0; i < NNN; ++i) {
        double re = arr[i][0];
        double im = arr[i][1];

        bool is_nan = std::isnan(re) || std::isnan(im);
        bool is_zero = (std::abs(re) < 1e-12 && std::abs(im) < 1e-12);

        if (is_nan) {
            // std::cout << "psi_K[" << i << "] = NaN!\n";
            ++nan_count;
        } else if (is_zero) {
            ++zero_count;
        }
    }

    std::cout << "\n -- " << msg << " -- \n";
    std::cout << "\tZeros: " << zero_count << " / " << NNN << "\n";
    std::cout << "\tNaNs : " << nan_count << " / " << NNN << "\n";
}

#ifdef FFMPEG 

void print_ascii_slice(const poisson_schrodinger_state* grid, int frame, double vmin, double vmax, int N) {
    constexpr const char* shades = " .:-=+*#%@";
    constexpr int n_shades = 10;
    int z = N / 2;

    std::cout << "Frame " << frame << " - [" << vmin << ":" << vmax << "]\n";
    std::cout << '+' << std::string(N, '-') << "+\n";

    for (int y = 0; y < N; ++y) {
        std::cout << '|';
        for (int x = 0; x < N; ++x) {
            int idx = z + y * N + x * N * N;
            double val = grid[idx].phi;
            double norm = (val - vmin) / (vmax - vmin);
            norm = std::clamp(norm, 0.0, 1.0);
            int shade_idx = static_cast<int>(norm * (n_shades - 1));
            std::cout << shades[shade_idx];
        }
        std::cout << "|\n";
    }

    std::cout << '+' << std::string(N, '-') << "+\n";
    std::cout << std::flush;
}
#endif 

// Thanks Peter Mortensen & Brian W.
std::string strip(const std::string& s, const std::string& chars) {
    size_t begin = 0;
    size_t end = s.size()-1;
    for(; begin < s.size(); begin++)
        if(chars.find_first_of(s[begin]) == std::string::npos)
            break;
    for(; end > begin; end--)
        if(chars.find_first_of(s[end]) == std::string::npos)
            break;
    return s.substr(begin, end-begin+1);
}

std::vector<double> read_csv_column(std::ifstream& file, int col) {
    std::vector<double> column;
    std::string line = "";

    while(std::getline(file, line)) {
        if(strip(line)[0] == '#') continue;
        std::istringstream sline(line);
        for(int i = 0; i <= col; i++)
            std::getline(sline, line, ','); 
        column.push_back(std::stod(line)); 
    }

    file.clear();
    file.seekg(0, std::ios::beg);
    return column;
}

double farthest_corner(double xx, double yy, double zz, int L) {
    double half_L = L / 2.0;
    double max_dist = 0.0;

    for (int x = 0; x <= 1; ++x) {
        for (int y = 0; y <= 1; ++y) {
            for (int z = 0; z <= 1; ++z) {
                double vx = (x ? half_L : -half_L);
                double vy = (y ? half_L : -half_L);
                double vz = (z ? half_L : -half_L);

                double dx = xx - vx;
                double dy = yy - vy;
                double dz = zz - vz;

                double distance = dist(dx, dy, dz);
                max_dist = std::max(max_dist, distance);
            }
        }
    }

    return max_dist;
} 

soliton::soliton(double mass, vec3 pos, vec3 vel, double phase) : mass(mass), pos(pos), vel(vel), phase(phase) {} 

double soliton::alpha() const {
    return SQR(mass / 3.88);
}

double init_grid(poisson_schrodinger_state* grid, double* fourier_grid, settings config, double& vmin, double& vmax, std::ofstream& en_file) {
    std::ifstream file(config.soliton_file);

    if(!file) { 
        std::cerr << "[ERROR] File handle (" << config.soliton_file << ") unresolved. Abort.";
        std::exit(-1);
    }

    std::vector<double> radii = read_csv_column(file, 0);
    std::vector<double> f_r = read_csv_column(file, 1);

    double L = config.box_size;
    int N = config.grid_size;
    int NNN = CUBE(N);

    // TODO -- Load Soliton Data, from soliton.info file.
    // Initialize soliton
    double dx = (double)(L / N);

    fftw_complex* psi_arr = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * NNN);
    fftw_complex* phi0_K = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N * (N/2 + 1));
    double* phi0_X = (double*)fftw_malloc(sizeof(double) * NNN);

    fftw_plan FT_psi = fftw_plan_dft_3d(N, N, N, psi_arr, psi_arr, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan IFT_psi = fftw_plan_dft_3d(N, N, N, psi_arr, psi_arr, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan FT_phi0  = fftw_plan_dft_r2c_3d(N, N, N, phi0_X, phi0_K, FFTW_ESTIMATE);
    fftw_plan IFT_phi0 = fftw_plan_dft_c2r_3d(N, N, N, phi0_K, phi0_X, FFTW_ESTIMATE);

    int max_rad = 0;

    #pragma omp parallel for collapse(3)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N;

                grid[idx].psi[0] = grid[idx].psi[1] = 0.0;

                double cx = (double)(x - (N - 1)/2.0) * dx;
                double cy = (double)(y - (N - 1)/2.0) * dx;
                double cz = (double)(z - (N - 1)/2.0) * dx;
                vec3 grid_loc = vec3(cx, cy, cz);

                for(soliton sol : config.solitons) {
                    // double farthest = farthest_corner(sol.pos.x, sol.pos.y, sol.pos.z);
                    double radius = std::sqrt(sol.alpha()) * dist(sol.pos.x, sol.pos.y, sol.pos.z, grid_loc.x, grid_loc.y, grid_loc.z);
                    int radius_idx = (int)(radius / 0.00001);
        
                    if(radius < 5.6) {
                        grid[idx].psi[0] += sol.alpha() * f_r[radius_idx] * std::cos(sol.phase + sol.vel * (grid_loc - sol.pos));
                        grid[idx].psi[1] += sol.alpha() * f_r[radius_idx] * std::sin(sol.phase + sol.vel * (grid_loc - sol.pos));
                    }
                }

                phi0_X[idx] = 4 * M_PI * abs2(grid[idx].psi);
            } 
        }      
    }

    fftw_execute(FT_phi0);

    #pragma omp parallel for collapse(3)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N/2 + 1; z++) {
                int idx = z + (y + x * N) * (N/2 + 1);
                int idx_f = z + y * N + x * N * N;

                if (fourier_grid[idx_f] > 1e-12) {
                    phi0_K[idx][0] *= -1.0 / fourier_grid[idx_f];
                    phi0_K[idx][1] *= -1.0 / fourier_grid[idx_f];
                } else {
                    phi0_K[idx][0] = 0.0;
                    phi0_K[idx][1] = 0.0;
                }
            } 
        }      
    }

    fftw_execute(IFT_phi0);

    #pragma omp parallel for collapse(3) 
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N;
                grid[idx].phi = phi0_X[idx] / NNN;

                double cx = (double)(x - (N - 1)/2.0) * dx;
                double cy = (double)(y - (N - 1)/2.0) * dx;
                double cz = (double)(z - (N - 1)/2.0) * dx;
                vec3 grid_loc = vec3(cx, cy, cz);

                for(potential pot : config.potentials) {
                    grid[idx].phi += pot.compute(grid_loc);
                }

                psi_arr[idx][0] = grid[idx].psi[0];
                psi_arr[idx][1] = grid[idx].psi[1];

                if(config.infer) {
                    vmin = std::min(vmin, abs2(grid[idx].psi));
                    vmax = std::max(vmax, abs2(grid[idx].psi));
                } else {
                    vmin = config.vmin;
                    vmax = config.vmax;
                }
            } 
        }
    }

    fftw_execute(FT_psi);

    double nrg = 0.0;
    double ek_eq = 0.0;
    double e_self = 0.0;
    double e_pot = 0.0;

    #pragma omp parallel for collapse(3) reduction (+:e_pot, e_self)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N; 

                double cx = (double)(x - (N - 1)/2.0) * dx;
                double cy = (double)(y - (N - 1)/2.0) * dx;
                double cz = (double)(z - (N - 1)/2.0) * dx;
                vec3 grid_loc = vec3(cx, cy, cz);

                double phi = grid[idx].phi;
                for(potential pot : config.potentials) {
                    double temp = pot.compute(grid_loc);
                    phi -= temp;
                    e_pot += temp * CUBE(dx) * abs2(grid[idx].psi);
                }

                e_self += 0.5 * CUBE(dx) * phi * abs2(grid[idx].psi);             

                psi_arr[idx][0] *= -fourier_grid[idx];
                psi_arr[idx][1] *= -fourier_grid[idx];
            }
        }
    }

    fftw_execute(IFT_psi);

    #pragma omp parallel for collapse(3) reduction (+:ek_eq)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N; 
                ek_eq += -0.5 * CUBE(dx) * (grid[idx].psi[0] * psi_arr[idx][0] + grid[idx].psi[1] * psi_arr[idx][1]) / NNN;             
            }
        }
    }

    nrg = ek_eq + e_pot + e_self;
    en_file << nrg << ", " << ek_eq << ", " << e_self << ", " << e_pot << "\n";

    fftw_destroy_plan(FT_psi);
    fftw_destroy_plan(IFT_psi);
    fftw_destroy_plan(FT_phi0);
    fftw_destroy_plan(IFT_phi0);

    fftw_free(phi0_X);
    fftw_free(phi0_K);
    fftw_free(psi_arr);

    std::cout << "[INFO] Suggested (vmin, vmax): (" << vmin << ", " << vmax << ")\n";

    return nrg;
}

void pseudo_spectral_solver(poisson_schrodinger_state* grid, double* fourier_grid, settings config, double& vmin, double& vmax, std::ofstream& en_file, const double& E0) {
    double L = config.box_size;
    int N = config.grid_size;
    int NNN = CUBE(N);
    double dx = (double)(L / N);
    double h = dx * dx / M_PI;

    #ifdef FFMPEG 
        double* slice = new double[N * N];
        unsigned char* rgb_frame = new unsigned char[N * N * 3];
    #endif 

    fftw_complex* psi_arr = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * NNN);
    fftw_complex* phi_K = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N * N * (N/2 + 1));
    double* phi_X = (double*)fftw_malloc(sizeof(double) * NNN);
   
    fftw_plan FT_psi  = fftw_plan_dft_3d(N, N, N, psi_arr, psi_arr, FFTW_FORWARD , FFTW_ESTIMATE);
    fftw_plan IFT_psi = fftw_plan_dft_3d(N, N, N, psi_arr, psi_arr, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_plan FT_phi  = fftw_plan_dft_r2c_3d(N, N, N, phi_X, phi_K, FFTW_ESTIMATE);
    fftw_plan IFT_phi = fftw_plan_dft_c2r_3d(N, N, N, phi_K, phi_X, FFTW_ESTIMATE);

    // Compute the half-step for Psi(x,t)
    #pragma omp parallel for collapse(3)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N;
                
                double _cos = std::cos(-0.5 * h * grid[idx].phi);
                double _sin = std::sin(-0.5 * h * grid[idx].phi);
                double _psi_re = grid[idx].psi[0];
                double _psi_im = grid[idx].psi[1];

                psi_arr[idx][0] = _psi_re * _cos - _psi_im * _sin;
                psi_arr[idx][1] = _psi_re * _sin + _psi_im * _cos;
            }
        }      
    }

    fftw_execute(FT_psi);
    // debug_count(psi_K, "fftw_execute(FT_psi);");

    #pragma omp parallel for collapse(3)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N;

                double _cos = std::cos(-0.5 * h * fourier_grid[idx]);
                double _sin = std::sin(-0.5 * h * fourier_grid[idx]);
                double _psi_re = psi_arr[idx][0];
                double _psi_im = psi_arr[idx][1];

                psi_arr[idx][0] = _psi_re * _cos - _psi_im * _sin;
                psi_arr[idx][1] = _psi_re * _sin + _psi_im * _cos; 
            }
        }      
    }

    fftw_execute(IFT_psi);
    // debug_count(psi_X, "fftw_execute(IFT_psi);");

    // Compute the step for Phi(x,t);
    #pragma omp parallel for collapse(3)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N;
                psi_arr[idx][0] /= NNN;
                psi_arr[idx][1] /= NNN;

                phi_X[idx] = 4 * M_PI * (psi_arr[idx][0] * psi_arr[idx][0] + psi_arr[idx][1] * psi_arr[idx][1]);
            }
        }      
    }

    fftw_execute(FT_phi); 
    // debug_count(phi_K, "fftw_execute(FT_phi);");

    #pragma omp parallel for collapse(3)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N/2 + 1; z++) {
                int idx = z + (y + x * N) * (N/2 + 1);
                int idx_f = z + y * N + x * N * N;

                if (fourier_grid[idx_f] > 1e-12) {
                    phi_K[idx][0] /= -fourier_grid[idx_f];
                    phi_K[idx][1] /= -fourier_grid[idx_f];
                } else {
                    phi_K[idx][0] = 0.0;
                    phi_K[idx][1] = 0.0;
                }
            }
        }      
    }

    fftw_execute(IFT_phi);
    // debug_count(phi_X, "fftw_execute(IFT_phi);");

    // Compute the final second-half step for Psi(x,t)
    #pragma omp parallel for collapse(3)
    for(int x = 0; x < N; x++) {
        for(int y = 0; y < N; y++) { 
            for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N;
                grid[idx].phi = phi_X[idx] / NNN;

                double cx = (double)(x - (N - 1)/2.0) * dx;
                double cy = (double)(y - (N - 1)/2.0) * dx;
                double cz = (double)(z - (N - 1)/2.0) * dx;
                vec3 grid_loc = {cx, cy, cz};

                for (potential pot : config.potentials) {
                    grid[idx].phi += pot.compute(grid_loc);
                }

                double _cos = std::cos(-0.5 * h * grid[idx].phi);
                double _sin = std::sin(-0.5 * h * grid[idx].phi);
                double _psi_re = psi_arr[idx][0];
                double _psi_im = psi_arr[idx][1];

                psi_arr[idx][0] = grid[idx].psi[0] = _psi_re * _cos - _psi_im * _sin;
                psi_arr[idx][1] = grid[idx].psi[1] = _psi_re * _sin + _psi_im * _cos; 

                #ifdef FFMPEG 
                    if(z == N/2) { 
                        double val = abs2(grid[idx].psi);
                        slice[x + y * N] = val; 
                       // vmax = std::max(vmax, val);
                    }
                #endif
            }
        }      
    }

    #ifdef FFMPEG
        apply_colormap_rgb_gamma(rgb_frame, slice, N * N, vmin, vmax, 1);
        // apply_colormap_rgb_log(rgb_frame, slice, N * N, vmin, vmax);
        stream_rgb_frame(rgb_frame, N, N);
        
        delete[] slice;
        delete[] rgb_frame;
    #endif

    // Energy code evaluation
    fftw_execute(FT_psi);

    double nrg = 0.0;
    double ek_eq = 0.0;
    double e_self = 0.0;
    double e_pot = 0.0;

    #pragma omp parallel for collapse(3) reduction (+:e_pot, e_self)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N; 

                double cx = (double)(x - (N - 1)/2.0) * dx;
                double cy = (double)(y - (N - 1)/2.0) * dx;
                double cz = (double)(z - (N - 1)/2.0) * dx;
                vec3 grid_loc = vec3(cx, cy, cz);

                double phi = grid[idx].phi;
                for(potential pot : config.potentials) {
                    double temp = pot.compute(grid_loc);
                    phi -= temp;
                    e_pot += temp * CUBE(dx) * abs2(grid[idx].psi);
                }

                e_self += 0.5 * CUBE(dx) * phi * abs2(grid[idx].psi);             

                psi_arr[idx][0] *= -fourier_grid[idx];
                psi_arr[idx][1] *= -fourier_grid[idx];
            }
        }
    }

    fftw_execute(IFT_psi);

    #pragma omp parallel for collapse(3) reduction (+:ek_eq)
    for(int x = 0; x < N; x++) {
         for(int y = 0; y < N; y++) {
             for(int z = 0; z < N; z++) {
                int idx = z + y * N + x * N * N; 
                ek_eq += -0.5 * CUBE(dx) * (grid[idx].psi[0] * psi_arr[idx][0] + grid[idx].psi[1] * psi_arr[idx][1]) / NNN;             
            }
        }
    }

    nrg = ek_eq + e_pot + e_self;
    en_file << std::setprecision(12) << (nrg - E0) / E0 << ", " << nrg << ", " << ek_eq << ", " << e_self << ", " << e_pot << "\n";

    fftw_destroy_plan(FT_psi);
    fftw_destroy_plan(IFT_psi);
    fftw_destroy_plan(FT_phi);
    fftw_destroy_plan(IFT_phi);

    fftw_free(psi_arr);
    fftw_free(phi_X);
    fftw_free(phi_K);
}

void generate_fourier_grid(double* fourier_grid, int N, double L) {
    for(int x = 0; x < N; x++) {
        double kx = 2.0 * M_PI * ((x < N / 2) ? x : x - N) / L;
        for(int y = 0; y < N; y++) {
            double ky = 2.0 * M_PI * ((y < N / 2) ? y : y - N) / L;
            for(int z = 0; z < N; z++) {
                double kz = 2.0 * M_PI * ((z < N / 2) ? z : z - N) / L;
                int idx = z + y * N + x * N * N;
                fourier_grid[idx] = kx*kx + ky*ky + kz*kz;
            }
        }      
    }
}