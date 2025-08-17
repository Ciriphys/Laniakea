#pragma once

#include <functional>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define SQR(x) ((x) * (x))
#define DIST(x,y,z) std::sqrt((SQR(x)) + (SQR(y)) + (SQR(z)))
#define DIST2(x,y,z) ((SQR(x)) + (SQR(y)) + (SQR(z)))

struct vec3 {
    vec3(double x = 0.0,
         double y = 0.0,
         double z = 0.0
        );

    union {
        double data[3];
        struct {
            double x, y, z;
        };
    };

    friend std::ostream& operator<<(std::ostream& os, const vec3& rhs);
};

struct soliton {
    soliton(double mass = 0.0, 
            vec3 pos = {0.0, 0.0, 0.0},
            vec3 vel = {0.0, 0.0, 0.0}, 
            double phase = 0.0
            );

    double mass; 
    vec3 pos, vel;
    double phase;
};

class potential {
    public: 
        potential(double mass = 0.0, std::string type = "", vec3 center = {});
        double compute(vec3 pos);

        double mass;
        vec3 center;
        std::string type;

    private:
        double inv_sqr(vec3 pos);
        double quadratic(vec3 pos);

        std::function<double(vec3)> cmp_ptr;  
};

struct settings {
    settings(int grid_size = -1.0, 
             double box_size = -1.0,   
             double max_time = -1.0,
             int threads = 0,
             int n_sol = 0,
             int n_pot = 0,
             bool ascii = false,
             std::string soliton_file = "soliton.dat",
             std::string wisdom_file = "wisdom.fftw",
             std::string renderer = "ffmpeg",
             std::string video_filename = "output.mkv", 
             std::vector<soliton> solitons = {},
             std::vector<potential> potentials = {}
             );

    int grid_size; 
    double box_size;
    double max_time;
    int threads;
    int n_sol;
    int n_pot;
    bool ascii;
    
    std::string soliton_file;
    std::string wisdom_file;
    std::string renderer;
    std::string video_ext;
    std::string video_filename;

    std::vector<soliton> solitons;
    std::vector<potential> potentials;
};

class parser {
    public: 
        parser(std::ifstream& handle);
        ~parser();

        settings parse();

    private:
        vec3 parse_vec3(std::string input);
        std::ifstream& handle;
};
