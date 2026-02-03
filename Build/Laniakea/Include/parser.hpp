#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "vec3.hpp"
#include "solver.hpp"

class parser {
    public: 
        parser(std::ifstream& handle);
        ~parser();

        settings parse();

    private:
        vec3 parse_vec3(std::string input);
        bool parse_bool(std::string input);
        std::ifstream& handle;
};
