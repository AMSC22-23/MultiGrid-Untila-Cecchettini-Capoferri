#include <iostream>
#include <vector>
#define DEFAULT_N 100

class DOMAIN {
    public:
        DOMAIN(int argc, char** argv) {
            using namespace std;
            // below... 
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

    void printN(){
        std::cout<<"And N is = "<< N<<std::endl;
    }

    private:
        std::size_t N = DEFAULT_N; // default value
        std::vector<int> Boundary = {0, 1000, 0, 1000};
};



template <typename T> 
class Jacobi1 : public System{

    
    public:

        Jacobi1(unsigned int max_N_iterations, double tol) : Max_N_iterations(max_N_iterations), Tol(tol) {}

        void Set_Tol(double tol){ Tol = tol; }
        double Get_Tol(){ return Tol; }

        void Set_Max_N_iterations(unsigned int max_N_iterations){ Max_N_iterations = max_N_iterations; }
        double Get_Max_N_iterations(){ return Max_N_iterations; }


        void iteration_method(std::vector<T>& sol) override{
            
            std::vector<T> x_sol;
            unsigned int k = 0;
            double sum;

            while(k <= Max_N_iterations){
                for(size_t i = 0; i < f.size(); i++){
                    sum = 0;
                    for(const auto &idx : A.nonZerosInRow(i)){
                        if(idx != i){
                            sum += A.coeffRef(i,idx) * x[idx];
                        }
                    }
                    x_sol[i] = (f[i] - sum) / A.coeffRef(i,i);
                }
                if(std::abs() < Tol)
                {
                    sol = x_sol;
                    break;
                }
                    
                k++;
            }
            
        }
    
    private:
        unsigned int Max_N_iterations = 1000; 
        double Tol; // tolerance
}