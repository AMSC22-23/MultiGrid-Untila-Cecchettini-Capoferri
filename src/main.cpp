#include "domain.hpp"
#include "linear_system.hpp"
#include "Main.hpp"

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
    AMG::SquareDomain dominio(size,width,0);
    AMG::SquareDomain dominio_2h(dominio);
    AMG::SquareDomain dominio_4h(dominio_2h);
    AMG::SquareDomain dominio_8h(dominio_4h);
    AMG::SquareDomain dominio_16h(dominio_8h);


    //Let's create the matrices
    AMG::PoissonMatrix<double> A(dominio,alpha);
    AMG::PoissonMatrix<double> A_2h(dominio_2h,alpha);
    AMG::PoissonMatrix<double> A_4h(dominio_4h,alpha);
    AMG::PoissonMatrix<double> A_8h(dominio_8h,alpha);
    AMG::PoissonMatrix<double> A_16h(dominio_16h,alpha);

    std::vector<AMG::PoissonMatrix<double>> matrici;
    matrici.push_back(A);
    matrici.push_back(A_2h);
    matrici.push_back(A_4h);
    matrici.push_back(A_8h);
    matrici.push_back(A_16h);
    
    //Now we need to create the known vector
    AMG::DataVector<double> fvec(dominio, f, g);

    
    //Let's create the solution vector and ine for the residual (and for the error)
    std::vector<double> u(A.rows(),0.);
    std::vector<double> res(u.size(),0.);
    

    AMG::SawtoothMGIteration<AMG::DataVector<double>,AMG::Gauss_Siedel_iteration<std::vector<double>>> MG(matrici,fvec);
    AMG::Residual<AMG::DataVector<double>> RES(matrici[0],fvec,res);

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