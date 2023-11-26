#include <iostream>
#include <tuple>
#include <fstream>
#include <cmath>
#include "classes.hpp"



template<class SpMat>
void saveMatrixOnFile(SpMat A, std::string fileName){

    std::ofstream file;
    file.open(fileName);
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
    file.open(fileName);
    file<<f.size()<<std::endl;

    for(size_t i = 0; i < f.size(); i++){
        file<<f[i]<<std::endl;
    }
    file.close();
}

double f(const double x, const double y){
    return 0.5;
}

double g(const double x, const double y){
    double r = sqrt(x*x + y*y);
    return sin(r);
}


int main(int argc, char** argv){
    size_t size = std::atoi(argv[1]);
    AMG::SquareDomain dominio(size,10.);
    AMG::PoissonMatrix<double> A(dominio);

    AMG::DataVector<double> fVec(dominio,f,g);
    
    //saveMatrixOnFile<AMG::PoissonMatrix<double>>(A,"out.mtx");
    //saveVectorOnFile<AMG::DataVector<double>>(fVec,"f.mtx");

    std::vector<double> x(fVec.size(),0.);      //our initial guess

    int niter = 1000;
    for(int i = 0; i < niter; i++){
        AMG::jacobiIteration(A,fVec,x);
    }

    saveVectorOnFile<std::vector<double>>(x,"x.mtx");

    
    return 0;
}