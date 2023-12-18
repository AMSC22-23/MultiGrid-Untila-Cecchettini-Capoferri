#include "allClasses.hpp"


double f(const double x, const double y){
    //return 0.5;
    return -5.0 * exp(x) * exp(-2.0 * y);
    //double k = 30.;
    //double r = sqrt(x*x + y*y);
    //return -k*(cos(k * r) / r - k*sin(k * r));
}

double g(const double x, const double y){
    //return 0.;
    return exp(x) * exp(-2.0 * y);
    //return sin(30. * sqrt(x * x + y * y));
}


int main(int argc, char** argv){
    
    size_t size = std::atoi(argv[1]);
    double alpha = std::atof(argv[2]);
    double width = std::atof(argv[3]);

    //Create a vector to store residuals norm to plot them
    std::vector<double> hist;

    //Let's define the domains on all levels
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

    
    //Let's create the solution vector and one for the residual
    std::vector<double> u(A.rows(),0.);
    std::vector<double> res(u.size(),0.);
    

    AMG::SawtoothMGIteration<AMG::DataVector<double>,AMG::Jacobi_iteration<std::vector<double>>> MG(matrici,fvec);
    AMG::Residual<AMG::DataVector<double>> RES(matrici[0],fvec,res);

    u * RES;
    hist.push_back(RES.Norm());
    

    
    int mgiter = 20;
    for(int i = 0; i < mgiter; i++){
        u * MG * RES;
        hist.push_back(RES.Norm());
    }


    Utils::saveVectorOnFile(hist,"MGGS4.txt");
    Utils::saveVectorOnFile(u,"x.mtx");

    return 0;
}
