#include "allIncludes.hpp"


int main(int argc, char** argv)
{
    
    // initialization of principal parameters from the user
    size_t size;
    double alpha;
    double width;
    int levels, test_functions;
    std::function<double(const double, const double)> f;
    std::function<double(const double, const double)> g;
    Utils::Initialization_for_N(argc, argv, size, alpha, width, levels, test_functions);
    Utils::init_test_functions(f, g, test_functions);


    #ifdef _OPENMP
    std::cout<<"Openmp enabled"<<std::endl;
    #endif

    auto start = std::chrono::high_resolution_clock::now();

    //Create a vector to store residuals norm to plot them
    std::vector<double> hist;

    
    std::vector<MultiGrid::SquareDomain> domains;

    for(int i = 0; i < levels; i++){
        domains.push_back(MultiGrid::SquareDomain(size,width,i));
    }
    
    std::vector<MultiGrid::PoissonMatrix<double>> matrici;

    for(auto &domain : domains){
        matrici.push_back(MultiGrid::PoissonMatrix<double>(domain,alpha));
    }
    
    //Now we need to create the known vector
    MultiGrid::DataVector<double> fvec(domains.front(), f, g);

    
    //Let's create the solution vector and one for the residual
    std::vector<double> u(matrici.front().rows(),0.);
    std::vector<double> res(u.size(),0.);
    

    MultiGrid::SawtoothMGIteration<MultiGrid::DataVector<double>,MultiGrid::Jacobi_iteration<std::vector<double>>> MG(matrici,fvec);
    MultiGrid::Residual<MultiGrid::DataVector<double>> RES(matrici.front(),fvec,res);

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