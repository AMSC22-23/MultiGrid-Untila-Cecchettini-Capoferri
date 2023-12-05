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
    
    //exact solution
    std::vector<double> ue(fvec.size());
    for(size_t i = 0; i < fvec.size(); i++){
        auto [x, y] = dominio_h[i];
        ue[i] = g(x, y);
    }

    int mgIterations = 20;
    int nu1 = 10;
    int nu2 = 20;
    int nu3 = 30;
    int maxSolverIter = 300;
    double tol = 1.e-6;

    std::vector<double> res(u.size());
    std::vector<double> hist;
    hist.push_back(AMG::error(u,ue,dominio_h));

    for(int j = 0; j < mgIterations; j++){
        std::cout<<j<<std::endl;
        std::vector<double> err(u.size(),0.);
        // V-cycle multigrid iteration (3 levels)

        // Do nu1 iterations on A_h x_h = b_h
        for(int i = 0; i < nu1; i++){
            AMG::gaussSeidelIteration<AMG::DataVector<double>>(A_h,fvec,u);

        }

        //compute the fine grid residual
        AMG::residual<AMG::DataVector<double>>(u,fvec,A_h,res);

        //solve for A_4h e_4h = r_4h
        
        for(int i = 0; i < maxSolverIter; i++){
            AMG::gaussSeidelIteration<std::vector<double>>(A_4h,res,err);
        }

        
        //interpolate 4h -> 2h
        AMG::Interpolation(err,dominio_2h,dominio_4h);

        //smooth on level 2h
        for(int i = 0; i < nu2; i++){
            AMG::gaussSeidelIteration<std::vector<double>>(A_2h,res,err);
        }

        

        //interpolate 2h -> h
        AMG::Interpolation(err,dominio_h,dominio_2h);

        for(size_t i = 0; i < u.size(); i++){
            u[i] += err[i];
        }

        //smooth on the fine grid
        for(int i = 0; i < nu3; i++){
            AMG::gaussSeidelIteration<AMG::DataVector<double>>(A_h,fvec,u);
        }

        //compute the true error to plot it
        hist.push_back(AMG::error(u,ue,dominio_h));
    }


    
    saveVectorOnFile(hist,"provaMG.txt");

    saveVectorOnFile(u,"x.mtx");


    return 0;
}
