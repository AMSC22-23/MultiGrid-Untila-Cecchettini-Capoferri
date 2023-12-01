#include <iostream>
#include <tuple>
#include <fstream>
#include <cmath>
#include "classes.hpp"



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

std::vector<double> formatVector(std::vector<double> &in, AMG::Domain &domain){
    std::vector<double> temp;
    for(size_t i = 0; i < domain.N(); i++){
        temp.push_back(in[domain.mask(i)]);
    }
    return temp;
}

double f(const double x, const double y){
    //return 0.5;
    return -5.0 * exp(x) * exp(-2.0 * y);
}

double g(const double x, const double y){
    
    //return 0.;
    return exp(x) * exp(-2.0 * y);
}


int main(int argc, char** argv){
    size_t size = std::atoi(argv[1]);
    //size_t size = 13;
    double width = 1.0;
    AMG::SquareDomain dominio_h(size,width,0);
    AMG::SquareDomain dominio_2h(size,width,1);
    AMG::SquareDomain dominio_4h(size,width,2);


    AMG::PoissonMatrix<double> A_h(dominio_h,1.);
    AMG::PoissonMatrix<double> A_2h(dominio_2h,1.);
    AMG::PoissonMatrix<double> A_4h(dominio_4h,1.);
    
    AMG::DataVector<double> fvec(dominio_h, f, g);

    //initial guess
    std::vector<double> u(fvec.size(), 0.);
    
    //starting a V-cycle multigrid iteration

    int smoothIterations = 1000;

    //fine grid smoothing
    for(int i = 0; i < smoothIterations; i++){
        AMG::gaussSeidelIteration(A_2h,fvec,u,dominio_2h);
    }

    AMG::Interpolation(u,dominio_h,dominio_2h,fvec);

    //auto u4h = formatVector(u,dominio_4h);
    saveVectorOnFile(u,"x.mtx");

    return 0;
}
