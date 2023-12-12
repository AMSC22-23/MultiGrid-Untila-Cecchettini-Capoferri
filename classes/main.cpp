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
    return -5.0 * exp(x) * exp(-2.0 * y);
}

double g(const double x, const double y){
    //return 0.;
    return exp(x) * exp(-2.0 * y);
}


int main(int argc, char** argv){
    
    size_t size = std::atoi(argv[1]);
    double width = 1.0;

    //Create a vector to store residuals norm to plot them
    std::vector<double> hist;

    //Let's define the domains on the three levels
    AMG::SquareDomain dominio(size,width,0);
    AMG::SquareDomain dominio_2h(dominio);
    AMG::SquareDomain dominio_4h(dominio_2h);

    //Let's create the matrices
    AMG::PoissonMatrix<double> A(dominio,1.);
    AMG::PoissonMatrix<double> A_2h(dominio_2h,1.);
    AMG::PoissonMatrix<double> A_4h(dominio_4h,1.);
    
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

    /*
    AMG::Residual<AMG::DataVector<double>> RES(matrici[0],fvec,res);

    std::vector<double> coarse_res(u.size());
    AMG::Residual<std::vector<double>> COARSE_RES(matrici.back(),res,coarse_res);

    //create the smoothers
    std::vector<std::unique_ptr<AMG::Iteration<std::vector<double>>>> iterations(matrici.size());
    for(size_t j = 0; j < matrici.size(); j++){
            iterations[j] = std::make_unique<AMG::Gauss_Siedel_iteration<std::vector<double>>>(matrici[j],res);
    }
    
    //create solver on the coarsest grid
    AMG::Solver<std::vector<double>> COARSE_SOLVER((*iterations.back()),COARSE_RES,2000,1.e-6,1);

    //create the interpolators
    std::vector<std::unique_ptr<AMG::InterpolationClass>> interpolators(matrici.size()-1);
    for(size_t j = 0; j < matrici.size() - 1; j++){
        interpolators[j] = std::make_unique<AMG::InterpolationClass>(matrici[j+1],matrici[j]);
    }

    int mgIter1 = 0;
    int nu = 10;

    u * RES;
    hist.push_back(RES.Norm());

    for(int i = 0; i < mgIter1; i++){
        u * RES;
        COARSE_RES.refresh_normalization_constant();
        
        err * COARSE_SOLVER * COARSE_RES;
        std::cout<<COARSE_RES.Norm()<<std::endl;


        for(size_t j = matrici.size() - 1; j > 0; --j){
            err * (*interpolators[j-1]);
            for(int k = 0; k < nu; k++){
                err * (*iterations[j-1]);
            }
        }

        for(size_t j = 0; j < u.size(); j++){
            u[j] += err[j];
            err[j] = 0;
        }
        u * RES;
        hist.push_back(RES.Norm());
    }
    */

    AMG::SawtoothMGIteration<AMG::DataVector<double>,AMG::Gauss_Siedel_iteration<std::vector<double>>> MG(matrici,fvec);
    AMG::Residual<AMG::DataVector<double>> RES(matrici[0],fvec,res);

    u * RES;
    hist.push_back(RES.Norm());

    int mgiter = 15;
    for(int i = 0; i < mgiter; i++){
        u * MG * RES;
        hist.push_back(RES.Norm());
    }
    
    
    saveVectorOnFile(hist,"histMG3.txt");
    saveVectorOnFile(u,"x.mtx");

    return 0;
}
