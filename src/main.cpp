#include "allIncludes.hpp"

template <typename T>
void printVector(std::vector<T> result){
    std::cout << "Vector elements: ";
    for (int i : result) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

int main(int argc, char** argv)
{
    size_t size = 25;
    double alpha = 1.;
    double width = 10.0;
    int levels = 1;

    std::vector<MultiGrid::SquareDomain> domains;
    for(int i = 0; i < levels; i++){
        domains.push_back(MultiGrid::SquareDomain(size,width,i));
    }
    
    //Then we can create the matrices
    std::vector<MultiGrid::PoissonMatrix<double>> matrici;
    for(auto &domain : domains){
        matrici.push_back(MultiGrid::PoissonMatrix<double>(domain,alpha));
    }

    MultiGrid::AMG M(matrici[0], 0.25);

    std::vector<double> a =  M.valueStrongConnection(80,1);
    printVector(a);
    /*
    
    // initialization of principal parameters from the user
    size_t size;
    double alpha;
    double width;
    int levels, test_functions;
    SMOOTHERS smoother;
    std::function<double(const double, const double)> f;
    std::function<double(const double, const double)> g;
    Utils::Initialization_for_N(argc, argv, size, alpha, width, levels, test_functions, smoother);
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
    MultiGrid::SawtoothMGIteration<MultiGrid::DataVector<double>,MultiGrid::Gauss_Seidel_iteration<std::vector<double>>> MG0(matrici,fvec);
    MultiGrid::SawtoothMGIteration<MultiGrid::DataVector<double>,MultiGrid::Jacobi_iteration<std::vector<double>>> MG1(matrici,fvec);
    MultiGrid::SawtoothMGIteration<MultiGrid::DataVector<double>,MultiGrid::BiCGSTAB<std::vector<double>>> MG2(matrici,fvec);

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
    int MaxIter = 20;


    switch (smoother)
    {

        case Gauss_Siedel:
            std::cout<<"GS iters"<<std::endl;   
            for(int i = 0; i < MaxIter; i++){
                u * GS * GS * MG0;                   
                u * RES;                            
                hist.push_back(RES.Norm());
                if(hist.back() <= TOL)
                    break;
            }
            break;
        case Jacobi:
            std::cout<<"Jacobi iters"<<std::endl;
            for(int i = 0; i < MaxIter; i++){
                u * GS * GS * MG1;                   
                u * RES;                            
                hist.push_back(RES.Norm());
                if(hist.back() <= TOL)
                    break;
            }
            break;
        
        case BiCGSTAB:
            std::cout<<"BiCGSTAB iters"<<std::endl;
            for(int i = 0; i < MaxIter; i++){
                u * GS * GS * MG1;                   
                u * RES;                            
                hist.push_back(RES.Norm());
                if(hist.back() <= TOL)
                    break;
            }
            break;
    
        default:
            break;
    }

    end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> solve_time = end - start;
    std::cout<<"||Solving elapsed time: "<<solve_time.count()<<" sec<br>"<<std::endl;
    std::cout<<"Tol: "<< TOL<<"<br>"<<std::endl;
    std::cout<<"Max iter: "<< MaxIter << "<br>"<< std::endl;
    

    //After solving the problem we can export the solution and the history
    Utils::saveVectorOnFile(hist,"MGGS4.txt");
    Utils::saveVectorOnFile(u,"x.mtx");
    */
    return 0;
}