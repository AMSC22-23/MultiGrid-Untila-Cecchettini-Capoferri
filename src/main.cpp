#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <numeric>
#include <memory>

// TODO: spiegazione cosa fa
#include "multigrid.hpp"

// TODO: spiegazione cosa fa

#include "solver.hpp"
// TODO: spiegazione cosa fa
#include "utilities.hpp"
// TODO: spiegazione cosa fa
#include "linear_system.hpp"
// TODO: spiegazione cosa fa
#include "domain.hpp"

// to set, for the discretization of the domain, the N: -n number_of_elements


using namespace std;

int main(int argc, char** argv)
{
    size_t size;
    double alpha; 
    Initialization_for_N(argc, argv, size, alpha);
    double width = 1.0;

    //Create a vector to store residuals norm to plot them
    std::vector<double> hist;

    //Let's define the domains on the three levels
    MultiGrid::SquareDomain dominio(size,width,0);
    MultiGrid::SquareDomain dominio_2h(dominio);
    MultiGrid::SquareDomain dominio_4h(dominio_2h);
    MultiGrid::SquareDomain dominio_8h(dominio_4h);
    MultiGrid::SquareDomain dominio_16h(dominio_8h);


    //Let's create the matrices
    MultiGrid::PoissonMatrix<double> A(dominio,alpha);
    MultiGrid::PoissonMatrix<double> A_2h(dominio_2h,alpha);
    MultiGrid::PoissonMatrix<double> A_4h(dominio_4h,alpha);
    MultiGrid::PoissonMatrix<double> A_8h(dominio_8h,alpha);
    MultiGrid::PoissonMatrix<double> A_16h(dominio_16h,alpha);

    std::vector<MultiGrid::PoissonMatrix<double>> matrici;
    matrici.push_back(A);
    matrici.push_back(A_2h);
    matrici.push_back(A_4h);
    matrici.push_back(A_8h);
    matrici.push_back(A_16h);
    
    //Now we need to create the known vector
    MultiGrid::DataVector<double> fvec(dominio, f, g);

    
    //Let's create the solution vector and ine for the residual (and for the error)
    std::vector<double> u(A.rows(),0.);
    std::vector<double> res(u.size(),0.);
    

    MultiGrid::SawtoothMGIteration<MultiGrid::DataVector<double>,MultiGrid::Gauss_Siedel_iteration<std::vector<double>>> MG(matrici,fvec);
    MultiGrid::Residual<MultiGrid::DataVector<double>> RES(matrici[0],fvec,res);

    u * RES;
    hist.push_back(RES.Norm());

    

    
    int mgiter = 20;
    for(int i = 0; i < mgiter; i++){
        u * MG * RES;
        hist.push_back(RES.Norm());
    }
    
    
    saveVectorOnFile(hist,"MGGS4.txt");
    saveVectorOnFile(u,"x.mtx");

    return 0;
}