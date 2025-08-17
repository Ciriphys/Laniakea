#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

std::string strip(const std::string& s, const std::string& chars=" ") {
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
        // std::cout << line << std::endl;
        std::istringstream sline(line);
        for(int i = 0; i <= col; i++) {
            std::getline(sline, line, ',');
            // std::cout << "\t" << line << std::endl;
        } 
        column.push_back(std::stod(line)); 
    }

    file.clear();
    file.seekg(0, std::ios::beg);
    return column;
}

int main() {

    std::ifstream file("../soliton.dat");
    std::vector<double> radii = read_csv_column(file, 0);
    std::vector<double> f_r = read_csv_column(file, 1);

    std::cout << radii.size() << "  " << f_r.size() << std::endl;

    double mean = 0.0;
    int N = 0;
    for(int i = 0; i < 10; i++) {
        mean *= N++;
        mean += (double)(i + 1);
        mean /= N;

        std::cout << mean << std::endl;
    }

    return 0;
} 
