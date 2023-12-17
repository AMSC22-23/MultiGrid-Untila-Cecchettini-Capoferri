#ifndef CLASSES_H
#define CLASSES_H

#include <tuple>
#include <vector>
#include <functional>
#include <numeric>
#include <memory>
#include "utils.hpp"



namespace AMG{

/*
For generalization we decided to create an abstract class Domain containing the declaration of all the public methods
we're gonna use.
*/


/*
The goal of our project is to do some tests on a simple square domain, but the main idea was to generalize the abstract class in such
a way that we can easily make a more generic domain divided in triangles
*/


/*
We decided for the linear system to create a class that behaves as a matrix, without allocating memory for matrix entries
We wanted to make poissonMatrix fully compatible with Eigen::SparseMatrix, but we didn't found an equivalent of
nonZerosInRow for eigen matrices
*/




/*
// un risolutore temporaneo

template<class Vector>
void gaussSeidelIteration(PoissonMatrix<double> &A, Vector &f, std::vector<double> &x){
    for(size_t i = 0; i < A.rows(); i++){
        size_t index = A.mask(i);
        double sum = 0;
        for(const auto &id : A.nonZerosInRow(i)){
            if(id != i){
                sum += A.coeffRef(i,id) * x[A.mask(id)];
            }
        }
        x[index] = (f[index] - sum) / A.coeffRef(i,i);
    }
}


void Interpolation(std::vector<double> &sol, Domain &domain_sup, Domain &domain_inf){
    for(size_t i = 0; i < domain_inf.N() - domain_inf.getWidth(); i++){
        size_t index1 = domain_inf.mask(i);
        size_t index2 = domain_inf.mask(i + domain_inf.getWidth());
        size_t index3 = (index1 + index2) / 2;
        
        if(! domain_sup.isOnBoundary(index3)){
            sol[index3] = 0.5 * (sol[index1] + sol[index2]);
        }else{
            //sol[index3] = f[index3];
            continue;
        }
    }
    
    size_t width = domain_sup.getWidth();
    for(size_t i = 0; i < domain_sup.N() / width; i++){
        for(size_t j = i * width; j < (i + 1) * width - 1; j += 2){
            if(! domain_sup.isOnBoundary(domain_sup.mask(j+1)))
                sol[domain_sup.mask(j + 1)] = 0.5 * (sol[domain_sup.mask(j)] + sol[domain_sup.mask(j + 2)]);
            else
                //sol[domain_sup.mask(j + 1)] = f[domain_sup.mask(j + 1)];
                continue;
        }
    }
    
}

double error(std::vector<double> &vec, std::vector<double> &sol, Domain &domain){
    double squaredError = 0;
    double squaredNorm = 0;
    
    for(size_t i = 0; i < domain.N(); i++){
        size_t index = domain.mask(i);
        double solVal = sol[index];
        double errVal = vec[index] - solVal;
        squaredError += errVal * errVal;
        squaredNorm += solVal * solVal;
    }
    return sqrt(squaredError / squaredNorm); 
}


template<class Vector>
void residual(std::vector<double> &u, Vector &f, PoissonMatrix<double> &A, std::vector<double> &res){
    for(size_t i = 0; i < A.rows(); i++){
        size_t index = A.mask(i);
        double sum = 0;
        for(const auto &j : A.nonZerosInRow(i)){
            sum += A.coeffRef(i, j) * u[A.mask(j)];
        }
        res[index] = f[index] - sum;
    }
}

template<class Vector>
double residualNorm(std::vector<double> &u, Vector &f, PoissonMatrix<double> &A){
    double res = 0.;
    for(size_t i = 0; i < A.rows(); i++){
        size_t index = A.mask(i);
        double sum = 0.;
        for(const auto &j : A.nonZerosInRow(i)){
            sum += A.coeffRef(i, j) * u[A.mask(j)];
        }
        sum = f[index] - sum;
        res += sum * sum;
    }
    return sqrt(res);
}

template<class Vector>
double norm(Vector &u, PoissonMatrix<double> &A){
    double sum = 0;
    for(size_t i = 0; i < A.rows(); i++){
        double val = u[A.mask(i)];
        sum += val * val;
    }
    return sqrt(sum);
}

*/

// TODO: comment
















/*
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

            COARSE_SOLVER = std::make_unique<AMG::Solver<std::vector<double>>>((*iterations.back()),(*COARSE_RES),2000,0.5,1);

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
            Utils::saveVectorOnFile(formatVector(temp,A_level.back()),name);
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
            Utils::saveVectorOnFile(formatVector(temp,A_level.back()),name);
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
                

                Utils::saveVectorOnFile(formatVector(temp,A_level.at(j-1)),name);
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
                

                Utils::saveVectorOnFile(formatVector(temp,A_level.at(j-1)),name);
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
            Utils::saveVectorOnFile(formatVector(sol,A_level.at(0)),name);
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
*/





}


#endif

