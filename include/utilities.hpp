#ifndef UTILS_H
#define UTILS_H


#include "allIncludes.hpp"
#include <iostream>
#include <string>

#define DEFAULT_N 400
#define DEFAULT_ALPHA 10.0
#define DEFAULT_WIDTH 10.0
#define DEFAULT_LEVEL 2
#define DEFAULT_TEST 1

namespace Utils{

void Initialization_for_N(int argc, char** argv, size_t &N, double &alpha, double &width, unsigned char &level, unsigned char &functions_to_test);

template<class SpMat>
void saveMatrixOnFile(SpMat A, std::string fileName){

    std::ofstream file;
    file.open(fileName, std::ofstream::trunc);
    file<<A.rows()<<" "<<A.cols()<<" "<<A.nonZeros()<<std::endl;

    for(size_t i = 0; i < A.rows(); i++){
        std::vector<size_t> row = A.nonZerosInRow(i);
        for(const auto& j : row){
            file<<i<<" "<<j<<" "<<A.coeffRef(i,j)<<std::endl;
        }
    }
    file.close();
}

template<class Vector>
void saveVectorOnFile(Vector f, std::string fileName){
    
    std::ofstream file;
    file.open(fileName, std::ofstream::trunc);
    file<<f.size()<<std::endl;

    for(size_t i = 0; i < f.size(); i++){
        file<<f[i]<<std::endl;
    }
    file.close();
}

}

#endif