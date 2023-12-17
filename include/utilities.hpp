<<<<<<< HEAD
#ifndef UTILS_H
#define UTILS_H


#include "allIncludes.hpp"

namespace Utils{

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
=======
#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <numeric>
#include <memory>
#include <math.h>
#include <fstream>

// TODO: spiegazione cosa fa
#include "multigrid.hpp"

// TODO: spiegazione cosa fa

#include "solver.hpp"
// TODO: spiegazione cosa fa
#include "linear_system.hpp"
// TODO: spiegazione cosa fa
#include "domain.hpp"


namespace MultiGrid{
    // TODO: Spiegazione funzionalità metodo
    void Initialization_for_N(int argc, char** argv, size_t &N, double &alpha);

    // TODO: Spiegazione funzionalità metodo
    template<class SpMat>
    void saveMatrixOnFile(SpMat A, std::string fileName);

    // TODO: Spiegazione funzionalità metodo
    template<class Vector>
    void saveVectorOnFile(Vector f, std::string fileName);

    // TODO: Spiegazione funzionalità metodo
    std::vector<double> formatVector(std::vector<double> &in, Domain &domain);


    // TODO: Spiegazione funzionalità metodo
    inline double f(const double x, const double y){
        return -5.0 * exp(x) * exp(-2.0 * y);
    }
 // namespace UTL

};

#endif
>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb
