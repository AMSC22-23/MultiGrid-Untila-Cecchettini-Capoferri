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
    return 0.5;
    //return -5.0 * exp(x) * exp(-2.0 * y);
}

double g(const double x, const double y){
    return 0.;
    //return exp(x) * exp(-2.0 * y);
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
    std::vector<double> err(u.size(),0.);

    //As a smoother we're gonna use Gauss Seidel
    AMG::Gauss_Siedel_iteration<AMG::DataVector<double>> SMOOTH(A,fvec);

    //We need a residual updater
    AMG::Residual<AMG::DataVector<double>> RES(A,fvec,res);
    u * RES;
    hist.push_back(RES.Norm());

    //Before starting the mg iterations we can define the interpolation operators
    AMG::InterpolationClass INTERPOLATE_4h(A_4h,A_2h);
    AMG::InterpolationClass INTERPOLATE_2h(A_2h,A);
    
    int mgIter = 20;
    int nu1 = 10;
    int nu2 = 10;
    int nu3 = 10;

    //***UNTIL WE SOLVE A PROBLEM WE NEED TO CREATE A NEW RESIDUAL VECTOR***
    std::vector<double> res_4h(u.size());

    for(int k = 0; k < mgIter; ++k){
        //V-cycle 3 levels multigrid iteration

        //Let's start by smoothing in the finest level
        for(int i = 0; i < nu1; i++){
            u = u * SMOOTH;
        }

        //Now we need to compute the residual
        u = u * RES;

        //Now we have to compute A * err = res on the coarsest grid
        AMG::Gauss_Siedel_iteration<std::vector<double>> GS_4h(A_4h,res);
        AMG::Residual<std::vector<double>> RES_4H(A_4h,res,res_4h);
        AMG::Solver<std::vector<double>> SOLVE(GS_4h,RES_4H,5000,1.e-6,10);

        err = err * SOLVE;
        std::cout<<"Coarse grid normalized residual "<<RES_4H.Norm()<<std::endl;

        //Once we solved the error on the coarsest level we can interpolate
        err = err * INTERPOLATE_4h;

        //After the interpolation we need to smooth the error
        AMG::Gauss_Siedel_iteration<std::vector<double>> GS_2H(A_2h,res);
        for(int i = 0; i < nu2; i++){
            err = err * GS_2H;
        }

        //We can now interpolate again
        err = err * INTERPOLATE_2h;

        //We can now sum the error to our solution
        for(size_t i = 0; i < u.size(); i++){
            u[i] += err[i];
            //set err to 0 to have a closer initial guess in the next iteration
            err[i] = 0.;
        }

        //now we can smooth again the solution
        for(int i = 0; i < nu3; i++){
            u = u * SMOOTH;
        }
        u = u * RES;
        hist.push_back(RES.Norm());
    }

    saveVectorOnFile(hist,"histMG3.txt");
    saveVectorOnFile(u,"x.mtx");

    return 0;
}
