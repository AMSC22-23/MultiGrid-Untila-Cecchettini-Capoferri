<<<<<<< HEAD
#include "allIncludes.hpp"

// to set, for the discretization of the domain, the N: -n number_of_elements

/*
void Initialization_for_N(int argc, char** argv, unsigned int &N){
    string nn = "-n";

    if(argc < 2) cout<<"Inserted by default N = "<<DEFAULT_N<<endl;
    else{
        for(int i= 0; i < argc; i++)
        {
            if(nn.compare(argv[i]) == 0 && ( i+1 < argc ) )
            {
                try{
                    N = stoi(argv[i+1]);
                }catch(exception){
                    cout<<"Please, insert a number after -n"<<endl<<"Inserted by default N = 100"<<endl;
                    exit(1);
                }
            }
            else if(nn.compare(argv[i]) == 0 && ( i+1 == argc )) {
                cout<<"Please, insert a number after -n"<<endl<<"Inserted by default N = 100"<<endl;
                exit(1);
            }
        } 
    } 
}
=======
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <numeric>
#include <memory>

#include "main.hpp"

// to set, for the discretization of the domain, the N: -n number_of_elements

>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb

using namespace std*/

double f(const double x, const double y){
    //return 0.5;
    //return -5.0 * exp(x) * exp(-2.0 * y);
    double k = 2.;
    double r = sqrt(x*x + y*y);
    return -k*(cos(k * r) / r - k*sin(k * r));
}

double g(const double x, const double y){
    //return 0.;
    //return exp(x) * exp(-2.0 * y);
    return sin(2. * sqrt(x * x + y * y));
}

int main(int argc, char** argv)
{
<<<<<<< HEAD
    //unsigned int N;
    //Initialization_for_N(argc, argv, N);
=======
    size_t size;
    double alpha; 
    MultiGrid::Initialization_for_N(argc, argv, size, alpha);
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
    
    
    MultiGrid::saveVectorOnFile(hist,"MGGS4.txt");
    MultiGrid::saveVectorOnFile(u,"x.mtx");

    std::cout<<"Hello";
>>>>>>> 93daadda06d78812819edd8bb3abfea81b24b6fb

    return 0;
}