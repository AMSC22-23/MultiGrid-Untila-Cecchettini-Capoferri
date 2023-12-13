#include <iostream>
#include <tuple>
#include <fstream>
#include <cmath>
#include <memory>
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
    //return -5.0 * exp(x) * exp(-2.0 * y);
    double r = sqrt(x*x + y*y);
    return -cos(r) / r + sin(r);
}

double g(const double x, const double y){
    //return 0.;
    //return exp(x) * exp(-2.0 * y);
    return sin(sqrt(x * x + y * y));
}


int main(int argc, char** argv){
    
    size_t size = std::atoi(argv[1]);
    double alpha = std::atof(argv[2]);
    double width = 10.0;

    //Create a vector to store residuals norm to plot them
    std::vector<double> hist;

    //Let's define the domains on the three levels
    AMG::SquareDomain dominio(size,width,0);
    AMG::SquareDomain dominio_2h(dominio);
    AMG::SquareDomain dominio_4h(dominio_2h);

    //Let's create the matrices
    AMG::PoissonMatrix<double> A(dominio,alpha);
    AMG::PoissonMatrix<double> A_2h(dominio_2h,alpha);
    AMG::PoissonMatrix<double> A_4h(dominio_4h,alpha);
    
    //Now we need to create the known vector
    AMG::DataVector<double> fvec(dominio, f, g);

    
    //Let's create the solution vector and ine for the residual (and for the error)
    std::vector<double> u(A.rows(),0.);
    std::vector<double> res(u.size(),0.);
    //std::vector<double> err(u.size(),0.);
    
   
    std::vector<AMG::PoissonMatrix<double>> matrici;
    matrici.push_back(A);
    matrici.push_back(A_2h);
    matrici.push_back(A_4h);

    AMG::SawtoothMGIteration<AMG::DataVector<double>,AMG::Gauss_Siedel_iteration<std::vector<double>>> MG(matrici,fvec);
    AMG::Residual<AMG::DataVector<double>> RES(matrici[0],fvec,res);

    u * RES;
    hist.push_back(RES.Norm());

    int mgiter = 20;
    for(int i = 0; i < mgiter; i++){
        u * MG * RES;
        hist.push_back(RES.Norm());
    }
    
    saveVectorOnFile(hist,"MGGS.txt");
    saveVectorOnFile(u,"x.mtx");

    return 0;
}
