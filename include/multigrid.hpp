#ifndef MG_H
#define MG_H

#include "allIncludes.hpp"

namespace MultiGrid{

class InterpolationClass{
    private:
        PoissonMatrix<double> &m_A_inf, &m_A_sup;

    public:
        InterpolationClass(PoissonMatrix<double> &A_inf, PoissonMatrix<double> &A_sup) : m_A_inf(A_inf), m_A_sup(A_sup){}
        
        void interpolate(std::vector<double> &vec);
        
        friend std::vector<double>& operator*(std::vector<double> &x_k, InterpolationClass &B)
        {
            B.interpolate(x_k);
            return x_k;
        }

};

template<typename T>
class AMG{
    private:
        PoissonMatrix<T> &m_A;
        std::vector<T> res;
        double epsilon;
    
    public:

        AMG(PoissonMatrix<T> &A, double epsilonI) : m_A(A), epsilon(epsilonI)  {}

        int getRandomInit(int max){
            std::random_device rd;

            // Initialize a random number generator with the random device
            std::mt19937 gen(rd());
            // Create a uniform distribution within the specified range
            std::uniform_int_distribution<> distr(0, max);

            // Generate a random number within the range
            return distr(gen);
        }

        std::vector<T> valueStrongConnection(int NodeIndex, unsigned char multi){
            T max = 0;
            std::vector<size_t> R = m_A.nonZerosInRow(NodeIndex);
            for(int i=0 ; i < R.size(); i++)
            {
                if(i != NodeIndex && std::abs(m_A.coeffRef(NodeIndex,R[i])) > max)
                    max = std::abs(m_A.coeffRef(NodeIndex,R[i]));
            }
            // compute strong connections
            std::vector<T> Ret; //m_a.cols
            for(int i=0; i < R.size(); i++)
            {
                if(i != NodeIndex && std::abs(m_A.coeffRef(NodeIndex,R[i])) >= epsilon*max)
                {
                    Ret.push_back(1*multi);
                }else if(i == NodeIndex || multi == 1){
                    Ret.push_back(0);
                }else{
                    Ret.push_back(1);
                }
            }
            return Ret;
        }

        T evaluateNode(int NodeIndex, std::vector<T> allNodes)
        {
            std::vector<T> R = valueStrongConnection(NodeIndex,1); // is a vector with position i-th equal to 1 if strongly connected, 0 otherwise

            if (allNodes.size() != R.size()) {
                std::cerr << "Vectors must be of the same length!" << std::endl;
                return -1;
            }

            return inner_product(R.begin(), R.end(), allNodes.begin(), 0);
        }

        void AMGIteration();

};

#ifndef CREATE_GIF
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

            COARSE_SOLVER = std::make_unique<MultiGrid::Solver<std::vector<double>>>((*iterations.back()),(*COARSE_RES),2000,1.e-1,1);
        }        

        void apply_iteration_to_vec(std::vector<double> &sol) {
            sol * (*RES);
            COARSE_RES->refresh_normalization_constant();

            err * (*COARSE_SOLVER) * (*COARSE_RES);
            std::cout<<"Achieved residual on coarse grid: "<<COARSE_RES->Norm()<<std::endl;


            for(size_t j = A_level.size() - 1; j > 0; --j){
                err * (*interpolators[j-1]);
                for(int i = 0; i < nu; i++){
                    err * (*iterations[j-1]);
                }
            }
            
            for(size_t j = 0; j < sol.size(); j++){
                sol[j] += err[j];
                err[j] = 0;
            }
        }

        friend std::vector<double>& operator*(std::vector<double> &x_k, SawtoothMGIteration &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        }



        ~SawtoothMGIteration(){
        }
};
#endif

#ifdef CREATE_GIF

using namespace Utils;

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
        int nu = 2;

        //I need this for the gif
        unsigned int counter = 0;
        std::vector<double> temp;
        std::string name;

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

            COARSE_SOLVER = std::make_unique<MultiGrid::Solver<std::vector<double>>>((*iterations.back()),(*COARSE_RES),2000,0.6,1);

            //for the gif
            temp = std::vector<double>(b.size());
        }


        
        std::vector<double> formatVector(std::vector<double> &in, PoissonMatrix<double> &A){
            std::vector<double> temp;
            for(size_t i = 0; i < A.rows(); i++){
                temp.push_back(in[A.mask(i)]);
            }
            return temp;
        }
        

        void apply_iteration_to_vec(std::vector<double> &sol) {
            //for the gif
            temp = sol;
            //********

            sol * (*RES);
            COARSE_RES->refresh_normalization_constant();

            //gif
            name = "./output/" + std::to_string(counter) + ".mtx";
            saveVectorOnFile(formatVector(temp,A_level.back()),name);
            counter++;
            //******

            err * (*COARSE_SOLVER) * (*COARSE_RES);
            std::cout<<"Achieved residual on coarse grid: "<<COARSE_RES->Norm()<<std::endl;

            //gif
            name = "./output/" + std::to_string(counter) + ".mtx";
            for(size_t j = 0; j < A_level.back().rows(); j++){
                auto id = A_level.back().mask(j);
                temp[id] += err[id];
            }
            saveVectorOnFile(formatVector(temp,A_level.back()),name);
            counter++;
            temp = sol;
            //******

            for(size_t j = A_level.size() - 1; j > 0; --j){
                err * (*interpolators[j-1]);
                //gif
                
                name = "./output/" + std::to_string(counter) + ".mtx";
                
                for(size_t k = 0; k < A_level.at(j-1).rows(); k++){
                    auto id = A_level.at(j-1).mask(k);
                    temp[id] += err[id];
                }
                

                saveVectorOnFile(formatVector(temp,A_level.at(j-1)),name);
                counter++;
                temp = sol;
                
                //******

                for(int i = 0; i < nu; i++){
                    err * (*iterations[j-1]);
                }

                //gif
                name = "./output/" + std::to_string(counter) + ".mtx";
                
                for(size_t k = 0; k < A_level.at(j-1).rows(); k++){
                    auto id = A_level.at(j-1).mask(k);
                    temp[id] += err[id];
                }
                

                saveVectorOnFile(formatVector(temp,A_level.at(j-1)),name);
                counter++;
                temp = sol;
                //******
            }
            
            for(size_t j = 0; j < sol.size(); j++){
                sol[j] += err[j];
                err[j] = 0;
            }
            
            //gif
            name = "./output/" + std::to_string(counter) + ".mtx";
            saveVectorOnFile(formatVector(sol,A_level.at(0)),name);
            counter++;
            //******

            
        }

        friend std::vector<double>& operator*(std::vector<double> &x_k, SawtoothMGIteration &B)
        {
            B.apply_iteration_to_vec(x_k);
            return x_k;
        }



        ~SawtoothMGIteration(){
        }
};

#endif


}


#endif