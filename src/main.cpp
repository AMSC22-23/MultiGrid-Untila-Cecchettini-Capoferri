#include "allIncludes.hpp"

// definition of solution function 
double f(const double x, const double y){
    return -5.0 * exp(x) * exp(-2.0 * y);
}

// definition of border solution function
double g(const double x, const double y){
    return exp(x) * exp(-2.0 * y);
}


int main(int argc, char** argv)
{
    
    // initialization of principal parameters from the user
    size_t size;
    double alpha;
    double width;
    Utils::Initialization_for_N(argc, argv, size, alpha, width);

    #ifdef _OPENMP
    std::cout<<"Openmp enabled"<<std::endl;
    #endif

    auto start = std::chrono::high_resolution_clock::now();
    //Create a vector to store residuals norm to plot them
    std::vector<double> hist;

    //Let's define the domains on all levels
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

    
    //Let's create the solution vector and one for the residual
    std::vector<double> u(A.rows(),0.);
    std::vector<double> res(u.size(),0.);
    

    MultiGrid::SawtoothMGIteration<MultiGrid::DataVector<double>,MultiGrid::Jacobi_iteration<std::vector<double>>> MG(matrici,fvec);
    MultiGrid::Residual<MultiGrid::DataVector<double>> RES(matrici[0],fvec,res);

    //We also need a smoother for the pre-smoothing
    MultiGrid::Gauss_Siedel_iteration<MultiGrid::DataVector<double>> GS(matrici.front(),fvec);

    
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> init_time = end - start;
    std::cout<<"Initialization time: "<<init_time.count()<<" seconds"<<std::endl;

    start = std::chrono::high_resolution_clock::now();
    
    u * RES;
    hist.push_back(RES.Norm());
    

    
    int mgiter = 20;
    for(int i = 0; i < mgiter; i++){
        u * GS * GS * MG;                   
        u * RES;                            
        hist.push_back(RES.Norm());
    }

    
    end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> solve_time = end - start;
    std::cout<<"Solving elapsed time: "<<solve_time.count()<<" seconds"<<std::endl;
    

    Utils::saveVectorOnFile(hist,"MGGS4.txt");
    Utils::saveVectorOnFile(u,"x.mtx");

    return 0;
}