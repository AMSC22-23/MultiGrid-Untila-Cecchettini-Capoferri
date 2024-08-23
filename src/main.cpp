#include "allIncludes.hpp"


int main(int argc, char** argv)
{
    
    // initialization of principal parameters from the user
    size_t size;
    double alpha;
    double width;
    int levels, test_functions;
    SMOOTHERS smooter;
    std::function<double(const double, const double)> f;
    std::function<double(const double, const double)> g;
    Utils::Initialization_for_N(argc, argv, size, alpha, width, levels, test_functions, smooter);
    Utils::init_test_functions(f, g, test_functions);


    #ifdef _OPENMP
    std::cout<<"Openmp enabled"<<std::endl;
    #endif

    //Initialization of timer
    auto start = std::chrono::high_resolution_clock::now();

    //Create a vector to store residuals norm to plot them
    std::vector<double> hist;

    
    //Let's create the domains
    std::vector<MultiGrid::SquareDomain> domains;
    for(int i = 0; i < levels; i++){
        domains.push_back(MultiGrid::SquareDomain(size,width,i));
    }
    
    //Then we can create the matrices
    std::vector<MultiGrid::PoissonMatrix<double>> matrici;
    for(auto &domain : domains){
        matrici.push_back(MultiGrid::PoissonMatrix<double>(domain,alpha));
    }
    

    //Now we need to create the known vector (forzante)
    MultiGrid::DataVector<double> fvec(domains.front(), f, g);

    
    //Let's create the solution vector and one for the residual
    std::vector<double> u(matrici.front().rows(),0.);
    std::vector<double> res(u.size(),0.);
    

    //Now we can create a multigrid iteration and the residual calculator (just to report the history of convergence)
    MultiGrid::SawtoothMGIteration<MultiGrid::DataVector<double>,MultiGrid::Jacobi_iteration<std::vector<double>>> MG(matrici,fvec);
    MultiGrid::Residual<MultiGrid::DataVector<double>> RES(matrici.front(),fvec,res);

    //We also need a smoother for the pre-smoothing
    MultiGrid::Gauss_Seidel_iteration<MultiGrid::DataVector<double>> GS(matrici.front(),fvec);


    
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> init_time = end - start;
    std::cout<<"Initialization time: "<<init_time.count()<<" seconds"<<std::endl;

    start = std::chrono::high_resolution_clock::now();
    

    //At first let's compute the residual
    u * RES;
    hist.push_back(RES.Norm());
    
    //Then we can start our iterations
    int mgiter = 20;
    for(int i = 0; i < mgiter; i++){
        u * GS * GS * MG;                   
        u * RES;                            
        hist.push_back(RES.Norm());
    }

    
    end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> solve_time = end - start;
    std::cout<<"Solving elapsed time: "<<solve_time.count()<<" seconds"<<std::endl;
    

    //After solving the problem we can export the solution and the history
    Utils::saveVectorOnFile(hist,"MGGS4.txt");
    Utils::saveVectorOnFile(u,"x.mtx");

    return 0;
}