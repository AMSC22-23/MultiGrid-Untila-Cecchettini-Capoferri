#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "main.hpp"


namespace MultiGrid {
class InterpolationClass{
    private:
        PoissonMatrix<double> &m_A_inf, &m_A_sup;

    public:
        InterpolationClass(PoissonMatrix<double> &A_inf, PoissonMatrix<double> &A_sup): m_A_inf(A_inf), m_A_sup(A_sup){}

        void interpolate(std::vector<double> &vec);

        friend std::vector<double>& operator*(std::vector<double> &x_k, InterpolationClass &B)
        {
            B.interpolate(x_k);
            return x_k;
        }
    

};




template<class Vector, class Smoother>
class SawtoothMGIteration{
    private:
        std::vector<PoissonMatrix<double>> &A_level;
        Vector &b;

        std::vector<double> res;
        std::vector<double> err;

        std::vector<std::unique_ptr<SmootherClass<std::vector<double>>>> iterations;
        std::vector<std::unique_ptr<InterpolationClass>> interpolators;

        std::unique_ptr<Residual<Vector>> RES;
        std::unique_ptr<Residual<std::vector<double>>> COARSE_RES;

        std::unique_ptr<Solver<std::vector<double>>> COARSE_SOLVER;
        int nu = 5;

        //std::vector<int> &nu;

    public:
        SawtoothMGIteration(std::vector<PoissonMatrix<double>> &matrices, Vector &knownVec): A_level(matrices), b(knownVec) {
            res = std::vector<double>(b.size(),0.);
            err = std::vector<double>(b.size(),0.);

            
            for(size_t j = 0; j < A_level.size(); j++){
                iterations.push_back(std::make_unique<Smoother>(A_level[j],res));
            }
            
            
            for(size_t j = 0; j < A_level.size() - 1; j++){
                interpolators.push_back(std::make_unique<InterpolationClass>(A_level[j+1],A_level[j]));
            }
            

            RES = std::make_unique<Residual<Vector>>(A_level.at(0),b,res);
            COARSE_RES = std::make_unique<Residual<std::vector<double>>>(A_level.back(),res);

            COARSE_SOLVER = std::make_unique<AMG::Solver<std::vector<double>>>((*iterations.back()),(*COARSE_RES),2000,1.e-1,1);
        }


        void apply_iteration_to_vec(std::vector<double> &sol);

        friend std::vector<double>& operator*(std::vector<double> &x_k, SawtoothMGIteration &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        }
        ~SawtoothMGIteration(){
        }
};

};

#endif